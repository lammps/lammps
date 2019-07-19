/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author (CAC) : Adrian Diaz (University of Florida)
------------------------------------------------------------------------- */

#include <cstring>
#include <cmath>
#include "comm_cac.h"
#include "comm_brick.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "output.h"
#include "dump.h"
#include "memory.h"
#include "error.h"
#include "nbin.h"
#include "nstencil.h"

using namespace LAMMPS_NS;

#define BUFFACTOR 1.5
#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000
#define EPSILON 1.0e-6
#define BOXEPSILON 1.0e-8
#define DELTA_PROCS 16

/* ---------------------------------------------------------------------- */

CommCAC::CommCAC(LAMMPS *lmp) : CommTiled(lmp)
{
  buf_send=NULL;
  buf_recv=NULL;
  pbc_flag = NULL;
  eboxes=NULL;
  foreign_eboxes=NULL;
  ebox_ref=NULL;
  maxebox=0;
  maxforeign_ebox=0;
  maxall=0;
  neboxes=0;
  local_neboxes=0;
  nforeign_eboxes=0;
  reset_array_flag=0;
  foreign_eprocs=NULL;
  foreign_image=NULL;
  foreign_swaps=NULL;
  work1=NULL;
  work2=NULL;
  overlap_repeat=NULL;
  proc2box=NULL;
  init_buffers();
}

/* ---------------------------------------------------------------------- */
//IMPORTANT: we *MUST* pass "*oldcomm" to the Comm initializer here, as
//           the code below *requires* that the (implicit) copy constructor
//           for Comm is run and thus creating a shallow copy of "oldcomm".
//           The call to Comm::copy_arrays() then converts the shallow copy
//           into a deep copy of the class with the new layout.

CommCAC::CommCAC(LAMMPS *lmp, Comm *oldcomm) : CommTiled(lmp, oldcomm)
{
  buf_send=NULL;
  buf_recv=NULL;
  eboxes=NULL;
  foreign_eboxes=NULL;
  ebox_ref=NULL;
  maxebox=0;
  maxforeign_ebox=0;
  maxall=0;
  neboxes=0;
  local_neboxes=0;
  nforeign_eboxes=0;
  reset_array_flag=0;
  foreign_eprocs=NULL;
  foreign_image=NULL;
  foreign_swaps=NULL;
  work1=NULL;
  work2=NULL;
  overlap_repeat=NULL;
  proc2box=NULL;
  comm_style = (const char *) "cac";
  Comm::copy_arrays(oldcomm);
  init_buffers();
}

/* ---------------------------------------------------------------------- */

CommCAC::~CommCAC()
{
  
  memory->destroy(buf_send);
  memory->destroy(buf_recv);
  memory->destroy(overlap);
  deallocate_swap(nswap);
  memory->sfree(rcbinfo);
   if (mode == Comm::MULTI) {
    memory->destroy(cutghostmulti);
  }
  memory->destroy(work1);
  memory->destroy(work2);
  atom->CAC_comm_flag=0;
}

/* ----------------------------------------------------------------------
   initialize comm buffers and other data structs local to CommTiled
------------------------------------------------------------------------- */

void CommCAC::init_buffers()
{
  sendbox_multi = NULL;
  cutghostmulti = NULL;
  // bufextra = max size of one exchanged atom
  //          = allowed overflow of sendbuf in exchange()
  // atomvec, fix reset these 2 maxexchange values if needed
  // only necessary if their size > BUFEXTRA
  
  maxexchange = 0;
  bufextra = maxexchange + BUFEXTRA;

  maxsend = BUFMIN;
  memory->create(buf_send,maxsend+bufextra,"comm:buf_send");
  maxrecv = BUFMIN;
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
  maxall = BUFMIN;
  memory->create(recv_flag,maxall,"commCAC:recv_flag");
  memory->create(lo2_set,8,3,"commCAC:lo2_set");
  memory->create(hi2_set,8,3,"commCAC:hi2_set");
  memory->create(current_pbc_set,8,3,"commCAC:current_pbc_set");
  max_image=DELTA_PROCS;
  maxoverlap = 0;
  maxoverlap_box = 0;
  overlap = NULL;
  overlap_pbc = NULL;
  memory->grow(work1,nprocs,"commCAC:work1");
  memory->grow(work2,nprocs,"commCAC:work2");
  memory->grow(overlap_repeat,nprocs,"commCAC:overlap_repeat");
  if(domain->dimension==2) nswap_border=8;
  if(domain->dimension==3) nswap_border=26;
  nswap = domain->dimension*2;
  
  allocate_swap(nswap_border);

  rcbinfo = NULL;
}

/* ---------------------------------------------------------------------- */

void CommCAC::init()
{
  Comm::init();
  
    // temporary restrictions

  if (triclinic)
    error->all(FLERR,"Cannot yet use comm_style cac with triclinic box");
  if (mode == Comm::MULTI)
    error->all(FLERR,"Cannot yet use comm_style cac with multi comm option");
  // memory for multi-style communication
 int maxswap = 26;
  if (mode == Comm::MULTI) {
    memory->create(cutghostmulti,atom->ntypes+1,3,"comm:cutghostmulti");
  }

  if (mode == Comm::SINGLE && cutghostmulti) {
    memory->destroy(cutghostmulti);
  }




  //CAC atom style check
  if (!atom->CAC_flag)
  error->all(FLERR,"Cannot use comm_style cac with non CAC atom style");

  atom->CAC_comm_flag=1;
  
}

/* ----------------------------------------------------------------------
   setup spatial-decomposition communication patterns
   function of neighbor cutoff(s) & cutghostuser & current box size and tiling
------------------------------------------------------------------------- */

void CommCAC::setup()
{
  int i,j,n;
  int ntypes = atom->ntypes;
  double lamda_temp[3];
  double nodal_temp[3];
	double ****nodal_positions=atom->nodal_positions;
	double current_distancesq;
  double search_radius;
  double dx,dy,dz;
	int *element_type = atom->element_type;
	int *poly_count = atom->poly_count;
	double ebounding_boxlo[3];
	double ebounding_boxhi[3];
	int *nodes_per_element_list = atom->nodes_per_element_list;
  double max_distancesq;
  int nodetotal;
  max_search_range= 0;
  // domain properties used in setup method and methods it calls
  size_forward = atom->avec->size_velocity;
  dimension = domain->dimension;
  prd = domain->prd;
  boxlo = domain->boxlo;
  boxhi = domain->boxhi;
  sublo = domain->sublo;
  subhi = domain->subhi;

    if (LAYOUT_NONUNIFORM==layout)
    error->all(FLERR,"Can only use comm style CAC with brick and rcb decompositions");
    
  int *periodicity = domain->periodicity;
  
  //check if CAC_pair style was invoked
  if(!atom->CAC_pair_flag)
  error->one(FLERR,"Cannot use the CAC comm style without a CAC pair style");

  // set function pointers

  if (layout != Comm::LAYOUT_TILED) {
    box_drop = &CommCAC::box_drop_brick;
    box_drop_full = &CommCAC::box_drop_brick_full;
    box_other = &CommCAC::box_other_brick;
    box_other_full = &CommCAC::box_other_brick_full;
    box_touch = &CommCAC::box_touch_brick;
    point_drop = &CommCAC::point_drop_brick;
  } else {
    box_drop = &CommCAC::box_drop_tiled;
    box_drop_full = &CommCAC::box_drop_tiled_full;
    box_other = &CommCAC::box_other_tiled;
    box_other_full = &CommCAC::box_other_tiled;
    box_touch = &CommCAC::box_touch_tiled;
    point_drop = &CommCAC::point_drop_tiled;
  }

  // if RCB decomp exists and just changed, gather needed global RCB info

  if (layout == Comm::LAYOUT_TILED) coord2proc_setup();

  // set cutoff for comm forward and comm reverse
  // check that cutoff < any periodic box length

  double cut;

    if (mode == Comm::MULTI) {
      double *cuttype = neighbor->cuttype;
      for (i = 1; i <= ntypes; i++) {
        cut = 0.0;
        if (cutusermulti) cut = cutusermulti[i];
        cutghostmulti[i][0] = MAX(cut,cuttype[i]);
        cutghostmulti[i][1] = MAX(cut,cuttype[i]);
        cutghostmulti[i][2] = MAX(cut,cuttype[i]);
      }
    }
  //cut = MAX(neighbor->cutneighmax,cutghostuser);
  /*cutoff for ghost information is the sum of the force cutoff
  and the difference between the limiting element bounding boxes
  and the local simulation box limits*/ 
  cut=neighbor->cutneighmax;

  //define the maximum search range
  int error_scale=1.10; //essentially a fudge factor for the search range 

	for(int element_index=0; element_index < atom->nlocal; element_index++)
	{
		 max_distancesq=0;
	   nodetotal=nodes_per_element_list[element_type[element_index]];
		//compute search radius using maximum distance between nodes of an element;
		//not the most rigorous approach; feel free to improve :).
	
		for(int ipoly=0; ipoly < poly_count[element_index]; ipoly++){
		for(i=0; i<nodetotal; i++){
			for (j=i+1; j<nodetotal; j++){
				dx=nodal_positions[element_index][i][ipoly][0]-nodal_positions[element_index][j][ipoly][0];
				dy=nodal_positions[element_index][i][ipoly][1]-nodal_positions[element_index][j][ipoly][1];
				dz=nodal_positions[element_index][i][ipoly][2]-nodal_positions[element_index][j][ipoly][2];
				current_distancesq=dx*dx+dy*dy+dz*dz;
				if(current_distancesq>max_distancesq) max_distancesq=current_distancesq;
			}
		}
    }
		search_radius=sqrt(max_distancesq);
		search_radius*=error_scale;
		if(search_radius>atom->max_search_range) atom->max_search_range=search_radius;
		
	}
  
   
	
	MPI_Allreduce(&atom->max_search_range,&max_search_range,1,MPI_DOUBLE,MPI_MAX,world);
  atom->max_search_range=max_search_range;
  if(cut>max_search_range)
  atom->max_search_range=max_search_range=cut;
  
  //double element_overlap_range[6];
  element_overlap_range[0]=element_overlap_range[1]=element_overlap_range[2]=
  element_overlap_range[3]=element_overlap_range[4]=element_overlap_range[5]=max_search_range+BOXEPSILON;
  cutghost[0] = cutghost[1] = cutghost[2] = cut+BOXEPSILON;
  
  if ((periodicity[0] && cutghost[0] > prd[0]) ||
      (periodicity[1] && cutghost[1] > prd[1]) ||
      (dimension == 3 && periodicity[2] && cutghost[2] > prd[2]))
    error->all(FLERR,"Communication cutoff for comm_style cac "
               "cannot exceed periodic box length");

  // if cut = 0.0, set to epsilon to induce nearest neighbor comm
  // this is b/c sendproc is used below to infer touching exchange procs
  // exchange procs will be empty (leading to lost atoms) if sendproc = 0
  // will reset sendproc/etc to 0 after exchange is setup, down below

  int cutzero = 0;
  if (cut == 0.0) {
    cutzero = 1;
    max_search_range = MIN(prd[0],prd[1]);
    if (dimension == 3) max_search_range = MIN(max_search_range,prd[2]);
    max_search_range *= EPSILON*EPSILON;
  }
  
  if(cutzero){
    element_overlap_range[0]=element_overlap_range[1]=element_overlap_range[2]=
  element_overlap_range[3]=element_overlap_range[4]=element_overlap_range[5]=0;
  }
  // setup forward/reverse communication
  // loop over 6 swap directions
  // determine which procs I will send to and receive from in each swap
  // done by intersecting ghost box with all proc sub-boxes it overlaps
  // sets nsendproc, nrecvproc, sendproc, recvproc
  // sets sendother, recvother, sendself, pbc_flag, pbc, sendbox
  // resets nprocmax

  int noverlap1,indexme;
  double lo1[3],hi1[3],lo2[3],hi2[3];
  int one,two;

  int iswap = 0;
  for (int idim = 0; idim < dimension; idim++) {
    for (int idir = 0; idir < 2; idir++) {

      // one = first ghost box in same periodic image
      // two = second ghost box wrapped across periodic boundary
      // either may not exist

      one = 1;
      lo1[0] = sublo[0]; lo1[1] = sublo[1]; lo1[2] = sublo[2];
      hi1[0] = subhi[0]; hi1[1] = subhi[1]; hi1[2] = subhi[2];
      if (idir == 0) {
        lo1[idim] = sublo[idim] - cutghost[idim]-element_overlap_range[idim];
        hi1[idim] = sublo[idim];
      } else {
        lo1[idim] = subhi[idim];
        hi1[idim] = subhi[idim] + cutghost[idim]+element_overlap_range[3+idim];
       
      }

      two = 0;
      if (idir == 0 && periodicity[idim] && lo1[idim] < boxlo[idim]) two = 1;
      if (idir == 1 && periodicity[idim] && hi1[idim] > boxhi[idim]) two = 1;

      if (two) {
        lo2[0] = sublo[0]; lo2[1] = sublo[1]; lo2[2] = sublo[2];
        hi2[0] = subhi[0]; hi2[1] = subhi[1]; hi2[2] = subhi[2];
        if (idir == 0) {
          lo2[idim] = lo1[idim] + prd[idim];
          hi2[idim] = boxhi[idim];
          if (sublo[idim] == boxlo[idim]) one = 0;
        } else {
          lo2[idim] = boxlo[idim];
          hi2[idim] = hi1[idim] - prd[idim];
          if (subhi[idim] == boxhi[idim]) one = 0;
        }
      }

      if (one) {
        if (idir == 0) lo1[idim] = MAX(lo1[idim],boxlo[idim]);
        else hi1[idim] = MIN(hi1[idim],boxhi[idim]);
        if (lo1[idim] == hi1[idim]) one = 0;
      }

      // noverlap = # of overlaps of box1/2 with procs via box_drop()
      // overlap = list of overlapping procs
      // if overlap with self, indexme = index of me in list

      indexme = -1;
      noverlap = 0;
      if (one) (this->*box_drop)(idim,lo1,hi1,indexme);
      noverlap1 = noverlap;
      if (two) (this->*box_drop)(idim,lo2,hi2,indexme);

      // if self is in overlap list, move it to end of list

      if (indexme >= 0) {
        int tmp = overlap[noverlap-1];
        overlap[noverlap-1] = overlap[indexme];
        overlap[indexme] = tmp;
      }

      // reallocate 2nd dimensions of all send/recv arrays, based on noverlap
      // # of sends of this swap = # of recvs of iswap +/- 1

      if (noverlap > nprocmax[iswap]) {
        int oldmax = nprocmax[iswap];
        while (nprocmax[iswap] < noverlap) nprocmax[iswap] += DELTA_PROCS;
        grow_swap_send(iswap,nprocmax[iswap],oldmax);
        if (idir == 0) {
          nrecv_procmax[iswap+1]=nprocmax[iswap];
          grow_swap_recv(iswap+1,nprocmax[iswap],oldmax);}
        else {
          nrecv_procmax[iswap-1]=nprocmax[iswap];
          grow_swap_recv(iswap-1,nprocmax[iswap],oldmax);}
        
      }

      // overlap how has list of noverlap procs
      // includes PBC effects

      if (noverlap && overlap[noverlap-1] == me) sendself[iswap] = 1;
      else sendself[iswap] = 0;
      if (noverlap && noverlap-sendself[iswap]) sendother[iswap] = 1;
      else sendother[iswap] = 0;

      nsendproc[iswap] = noverlap;
      for (i = 0; i < noverlap; i++) sendproc[iswap][i] = overlap[i];

      if (idir == 0) {
        recvother[iswap+1] = sendother[iswap];
        nrecvproc[iswap+1] = noverlap;
        for (i = 0; i < noverlap; i++) recvproc[iswap+1][i] = overlap[i];
      } else {
        recvother[iswap-1] = sendother[iswap];
        nrecvproc[iswap-1] = noverlap;
        for (i = 0; i < noverlap; i++) recvproc[iswap-1][i] = overlap[i];
      }
      
      iswap++;
    }
  }

  // setup exchange communication = subset of forward/reverse comm procs
  // loop over dimensions
  // determine which procs I will exchange with in each dimension
  // subset of procs that touch my proc in forward/reverse comm
  // sets nexchproc & exchproc, resets nexchprocmax

  int proc;

  for (int idim = 0; idim < dimension; idim++) {

    // overlap = list of procs that touch my sub-box in idim
    // proc can appear twice in list if touches in both directions
    // 2nd add-to-list checks to insure each proc appears exactly once

    noverlap = 0;
    iswap = 2*idim;
    n = nsendproc[iswap];
    for (i = 0; i < n; i++) {
      proc = sendproc[iswap][i];
      if (proc == me) continue;
      if ((this->*box_touch)(proc,idim,0)) {
        if (noverlap == maxoverlap) {
          maxoverlap += DELTA_PROCS;
          memory->grow(overlap,maxoverlap,"comm:overlap");
        }
        overlap[noverlap++] = proc;
      }
    }
    noverlap1 = noverlap;
    iswap = 2*idim+1;
    n = nsendproc[iswap];

    //MPI_Barrier(world);

    for (i = 0; i < n; i++) {
      proc = sendproc[iswap][i];
      if (proc == me) continue;
      if ((this->*box_touch)(proc,idim,1)) {
        for (j = 0; j < noverlap1; j++)
          if (overlap[j] == proc) break;
        if (j < noverlap1) continue;
        if (noverlap == maxoverlap) {
          maxoverlap += DELTA_PROCS;
          memory->grow(overlap,maxoverlap,"comm:overlap");
        }
        overlap[noverlap++] = proc;
      }
    }

    //MPI_Barrier(world);

    // reallocate exchproc and exchnum if needed based on noverlap

    if (noverlap > nexchprocmax[idim]) {
      while (nexchprocmax[idim] < noverlap) nexchprocmax[idim] += DELTA_PROCS;
      delete [] exchproc[idim];
      exchproc[idim] = new int[nexchprocmax[idim]];
      delete [] exchnum[idim];
      exchnum[idim] = new int[nexchprocmax[idim]];
    }

    nexchproc[idim] = noverlap;
    for (i = 0; i < noverlap; i++) exchproc[idim][i] = overlap[i];
  }

  // reset sendproc/etc to 0 if cut is really 0.0

  if (cutzero) {
    for (i = 0; i < nswap; i++) {
      nsendproc[i] = nrecvproc[i] =
        sendother[i] = recvother[i] = sendself[i] = 0;
    }
  }
  
  //setup borders for full corner and edge swaps now that exchange has been defined

      // one = first ghost box in same periodic image
      // two = second ghost box wrapped across periodic boundary
      // either may not exist

      one = 1;
      iswap=0;
      indexme = -1;
      noverlap = 0;
      //initialize repeat check array
      for(int init=0; init<nprocs; init++){
        overlap_repeat[init]=0;
      }
      for(int box_count=0; box_count<26; box_count++){
      lo1[0] = sublo[0]; lo1[1] = sublo[1]; lo1[2] = sublo[2];
      hi1[0] = subhi[0]; hi1[1] = subhi[1]; hi1[2] = subhi[2];
         one=1;
          //faces
          if(box_count==0){
          lo1[0] = sublo[0] - cutghost[0]-element_overlap_range[0];
          hi1[0] = sublo[0];
          }
          if(box_count==1){
          lo1[1] = sublo[1] - cutghost[1]-element_overlap_range[1];
          hi1[1] = sublo[1];
          }
          if(box_count==2){
          lo1[2] = sublo[2] - cutghost[2]-element_overlap_range[2];
          hi1[2] = sublo[2];
          }
          if(box_count==3){
          lo1[0] = subhi[0];
          hi1[0] = subhi[0]+ cutghost[0]+element_overlap_range[0];
          }
          if(box_count==4){
          lo1[1] = subhi[1];
          hi1[1] = subhi[1]+ cutghost[1]+element_overlap_range[1];
          }
          if(box_count==5){
          lo1[2] = subhi[2];
          hi1[2] = subhi[2]+ cutghost[2]+element_overlap_range[2];
          }
          //z axis edges
          if(box_count==6){
          hi1[0] = sublo[0];
          hi1[1] = sublo[1];
          lo1[0] = sublo[0] - cutghost[0]-element_overlap_range[0];
          lo1[1] = sublo[1] - cutghost[1]-element_overlap_range[1];
          }
          if(box_count==7){
          hi1[0] = sublo[0];
          lo1[1] = subhi[1];
          lo1[0] = sublo[0] - cutghost[0]-element_overlap_range[0];
          hi1[1] = subhi[1] + cutghost[1]+element_overlap_range[1];
          }
          if(box_count==8){
          lo1[0] = subhi[0];
          lo1[1] = subhi[1];
          hi1[0] = subhi[0] + cutghost[0]+element_overlap_range[0];
          hi1[1] = subhi[1] + cutghost[1]+element_overlap_range[1];
          }
          if(box_count==9){
          lo1[0] = subhi[0];
          hi1[1] = sublo[1];
          hi1[0] = subhi[0] + cutghost[0]+element_overlap_range[0];
          lo1[1] = sublo[1] - cutghost[1]-element_overlap_range[1];
          }
          //y axis edges
          if(box_count==10){
          hi1[0] = sublo[0];
          hi1[2] = sublo[2];
          lo1[0] = sublo[0] - cutghost[0]-element_overlap_range[0];
          lo1[2] = sublo[2] - cutghost[2]-element_overlap_range[2];
          }
          if(box_count==11){
          hi1[0] = sublo[0];
          lo1[2] = subhi[2];
          lo1[0] = sublo[0] - cutghost[0]-element_overlap_range[0];
          hi1[2] = subhi[2] + cutghost[2]+element_overlap_range[2];
          }
          if(box_count==12){
          lo1[0] = subhi[0];
          lo1[2] = subhi[2];
          hi1[0] = subhi[0] + cutghost[0]+element_overlap_range[0];
          hi1[2] = subhi[2] + cutghost[2]+element_overlap_range[2];
          }
          if(box_count==13){
          lo1[0] = subhi[0];
          hi1[2] = sublo[2];
          hi1[0] = subhi[0] + cutghost[0]+element_overlap_range[0];
          lo1[2] = sublo[2] - cutghost[2]-element_overlap_range[2];
          }
          //x axis edges
          if(box_count==14){
          hi1[1] = sublo[1];
          hi1[2] = sublo[2];
          lo1[1] = sublo[1] - cutghost[1]-element_overlap_range[1];
          lo1[2] = sublo[2] - cutghost[2]-element_overlap_range[2];
          }
          if(box_count==15){
          hi1[1] = sublo[1];
          lo1[2] = subhi[2];
          lo1[1] = sublo[1] - cutghost[1]-element_overlap_range[1];
          hi1[2] = subhi[2] + cutghost[2]+element_overlap_range[2];
          }
          if(box_count==16){
          lo1[1] = subhi[1];
          lo1[2] = subhi[2];
          hi1[1] = subhi[1] + cutghost[1]+element_overlap_range[1];
          hi1[2] = subhi[2] + cutghost[2]+element_overlap_range[2];
          }
          if(box_count==17){
          lo1[1] = subhi[1];
          hi1[2] = sublo[2];
          hi1[1] = subhi[1] + cutghost[1]+element_overlap_range[1];
          lo1[2] = sublo[2] - cutghost[2]-element_overlap_range[2];
          }
          //corners
          if(box_count==18){
          hi1[0]=sublo[0];
          hi1[1]=sublo[1];
          hi1[2]=sublo[2];
          lo1[0] = sublo[0] - cutghost[0]-element_overlap_range[0];
          lo1[1] = sublo[1] - cutghost[1]-element_overlap_range[1];
          lo1[2] = sublo[2] - cutghost[2]-element_overlap_range[2];
          }
          if(box_count==19){
          hi1[0]=sublo[0];
          lo1[1]=subhi[1];
          hi1[2]=sublo[2];
          lo1[0] = sublo[0] - cutghost[0]-element_overlap_range[0];
          hi1[1] = subhi[1] + cutghost[1]+element_overlap_range[1];
          lo1[2] = sublo[2] - cutghost[2]-element_overlap_range[2];
          }
          if(box_count==20){
          lo1[0]=subhi[0];
          lo1[1]=subhi[1];
          hi1[2]=sublo[2];
          hi1[0] = subhi[0] + cutghost[0]+element_overlap_range[0];
          hi1[1] = subhi[1] + cutghost[1]+element_overlap_range[1];
          lo1[2] = sublo[2] - cutghost[2]-element_overlap_range[2];
          }
          if(box_count==21){
          lo1[0]=subhi[0];
          hi1[1]=sublo[1];
          hi1[2]=sublo[2];
          hi1[0] = subhi[0] + cutghost[0]+element_overlap_range[0];
          lo1[1] = sublo[1] - cutghost[1]-element_overlap_range[1];
          lo1[2] = sublo[2] - cutghost[2]-element_overlap_range[2];
          }
          if(box_count==22){
          hi1[0]=sublo[0];
          hi1[1]=sublo[1];
          lo1[2]=subhi[2];
          lo1[0] = sublo[0] - cutghost[0]-element_overlap_range[0];
          lo1[1] = sublo[1] - cutghost[1]-element_overlap_range[1];
          hi1[2] = subhi[2] + cutghost[2]+element_overlap_range[2];
          }
          if(box_count==23){
          hi1[0]=sublo[0];
          lo1[1]=subhi[1];
          lo1[2]=subhi[2];
          lo1[0] = sublo[0] - cutghost[0]-element_overlap_range[0];
          hi1[1] = subhi[1] + cutghost[1]+element_overlap_range[1];
          hi1[2] = subhi[2] + cutghost[2]+element_overlap_range[2];
          }
          if(box_count==24){
          lo1[0]=subhi[0];
          lo1[1]=subhi[1];
          lo1[2]=subhi[2];
          hi1[0] = subhi[0] + cutghost[0]+element_overlap_range[0];
          hi1[1] = subhi[1] + cutghost[1]+element_overlap_range[1];
          hi1[2] = subhi[2] + cutghost[2]+element_overlap_range[2];
          }
          if(box_count==25){
          lo1[0]=subhi[0];
          hi1[1]=sublo[1];
          lo1[2]=subhi[2];
          hi1[0] = subhi[0] + cutghost[0]+element_overlap_range[0];
          lo1[1] = sublo[1] - cutghost[1]-element_overlap_range[1];
          hi1[2] = subhi[2] + cutghost[2]+element_overlap_range[2];
          }
      
      
      two = 0;
      current_pbc[0]=current_pbc[1]=current_pbc[2]=0;
      int pbc_loop[3];
      pbc_loop[2]=pbc_loop[1]=pbc_loop[0]=0;
      for(int idim=0; idim<domain->dimension; idim++){
      if (periodicity[idim] && lo1[idim] < boxlo[idim]) current_pbc[idim] = 1;
      if (periodicity[idim] && hi1[idim] > boxhi[idim]) current_pbc[idim] = -1;
      if(current_pbc[idim]) {
        two=1;
        pbc_loop[idim]=1;
        }
      }
      int image_count=0;
      int reduced_image_count=0;
      //compute how many and which images of lo1 to test for overlap
  
      if (two) {
        
        for(int pbx=0; pbx<=pbc_loop[0]; pbx++)
          for(int pby=0; pby<=pbc_loop[1]; pby++)
            for(int pbz=0; pbz<=pbc_loop[2]; pbz++){
               if(pbx==0&&pby==0&&pbz==0) continue;
               lo2_set[image_count][0] = lo1[0]; 
               lo2_set[image_count][1] = lo1[1]; 
               lo2_set[image_count][2] = lo1[2];
               hi2_set[image_count][0] = hi1[0];
               hi2_set[image_count][1] = hi1[1]; 
               hi2_set[image_count][2] = hi1[2];
               current_pbc_set[image_count][0]=0;
               current_pbc_set[image_count][1]=0;
               current_pbc_set[image_count][2]=0;
               if(pbx) {
                 lo2_set[image_count][0]+=prd[0]*current_pbc[0];
                 hi2_set[image_count][0]+=prd[0]*current_pbc[0];
                 current_pbc_set[image_count][0]=current_pbc[0];
               }
               if(pby) {
                 lo2_set[image_count][1]+=prd[1]*current_pbc[1];
                 hi2_set[image_count][1]+=prd[1]*current_pbc[1];
                 current_pbc_set[image_count][1]=current_pbc[1];
               }
               if(pbz) {
                 lo2_set[image_count][2]+=prd[2]*current_pbc[2];
                 hi2_set[image_count][2]+=prd[2]*current_pbc[2];
                 current_pbc_set[image_count][2]=current_pbc[2];
               }
               image_count+=1;
            }
          
      }
       
      if (nprocs==1) one = 0;
      if (one) {
        for(int idim=0; idim<domain->dimension; idim++){ 
         
         lo1[idim] = MAX(lo1[idim],boxlo[idim]);
         hi1[idim] = MIN(hi1[idim],boxhi[idim]);
         
        }
      }
      
      if (two) {
        for(int image_loop=0; image_loop<image_count; image_loop++){

        
        for(int idim=0; idim<domain->dimension; idim++){ 
         
         lo2_set[image_loop][idim] = MAX(lo2_set[image_loop][idim],boxlo[idim]);
         hi2_set[image_loop][idim] = MIN(hi2_set[image_loop][idim],boxhi[idim]);
         
        }
        }
      }
      
        
      
      // noverlap = # of overlaps of box1/2 with procs via box_drop()
      // overlap = list of overlapping procs
      // if overlap with self, indexme = index of me in list
      for(int idim=0; idim<domain->dimension; idim++){ 

         if(hi1[idim]==lo1[idim]) one=0; 
         
        }
      pbc_overlap=0;
      int overlap_find=1;
      if (one) (this->*box_drop_full)(0,lo1,hi1,indexme);
      noverlap1 = noverlap;
      pbc_overlap=1;
      if (two){
        for(int image_loop=0; image_loop<image_count; image_loop++){
        overlap_find=1;
        current_pbc[0]=current_pbc_set[image_loop][0];
        current_pbc[1]=current_pbc_set[image_loop][1];
        current_pbc[2]=current_pbc_set[image_loop][2];
        for(int idim=0; idim<domain->dimension; idim++){ 
         
         if(hi2_set[image_loop][idim]==lo2_set[image_loop][idim]) overlap_find=0; 
         
        }
        if(overlap_find)
        (this->*box_drop_full)(0,lo2_set[image_loop],hi2_set[image_loop],indexme); 
        } 
        
        }
      }

      sendself[iswap]=0;
      int tmp[10];
      int swapindex=noverlap-1;
      //move all overlaps with me to the end of the list with their respective pbc offsets and other overlap arrays
      for(int overlap_scan=0; overlap_scan<noverlap; overlap_scan++){
        if(overlap[overlap_scan]==me){
        while(overlap[swapindex]==me) swapindex-=1;
        if(overlap_scan>=swapindex) break;
        tmp[0]=overlap[swapindex];
        tmp[1]=overlap_pbc[swapindex][0];
        tmp[2]=overlap_pbc[swapindex][1];
        tmp[3]=overlap_pbc[swapindex][2];
        tmp[4]=proc2box[swapindex][0];
        tmp[5]=proc2box[swapindex][1];
        tmp[6]=proc2box[swapindex][2];
        tmp[7]=proc2box[swapindex][3];
        tmp[8]=proc2box[swapindex][4];
        tmp[9]=proc2box[swapindex][5];
        overlap[swapindex]=overlap[overlap_scan];
        overlap_pbc[swapindex][0]=overlap_pbc[overlap_scan][0];
        overlap_pbc[swapindex][1]=overlap_pbc[overlap_scan][1];
        overlap_pbc[swapindex][2]=overlap_pbc[overlap_scan][2];
        proc2box[swapindex][0]=proc2box[overlap_scan][0];
        proc2box[swapindex][1]=proc2box[overlap_scan][1];
        proc2box[swapindex][2]=proc2box[overlap_scan][2];
        proc2box[swapindex][3]=proc2box[overlap_scan][3];
        proc2box[swapindex][4]=proc2box[overlap_scan][4];
        proc2box[swapindex][5]=proc2box[overlap_scan][5];
        overlap[overlap_scan]=tmp[0];
        overlap_pbc[overlap_scan][0]=tmp[1];
        overlap_pbc[overlap_scan][1]=tmp[2];
        overlap_pbc[overlap_scan][2]=tmp[3];
        proc2box[overlap_scan][0]=tmp[4];
        proc2box[overlap_scan][1]=tmp[5];
        proc2box[overlap_scan][2]=tmp[6];
        proc2box[overlap_scan][3]=tmp[7];
        proc2box[overlap_scan][4]=tmp[8];
        proc2box[overlap_scan][5]=tmp[9];
        sendself[iswap]+=1;
        }
        
      }

      // reallocate 2nd dimensions of all send/recv arrays, based on noverlap
      // # of sends of this swap = # of recvs of iswap +/- 1
      
      if (noverlap > nprocmax[iswap]) {
        int oldmax = nprocmax[iswap];
        while (nprocmax[iswap] < noverlap) nprocmax[iswap] += DELTA_PROCS;
        grow_swap_send(iswap,nprocmax[iswap],oldmax);
        nrecv_procmax[iswap]=nprocmax[iswap];
        grow_swap_recv(iswap,nprocmax[iswap],oldmax);
     
        
      }

      // overlap how has list of noverlap procs
      // includes PBC effects

      if (noverlap && noverlap-sendself[iswap]) sendother[iswap] = 1;
      else sendother[iswap] = 0;
      
     
      nsendproc[iswap] = noverlap;
      for (i = 0; i < noverlap; i++) sendproc[iswap][i] = overlap[i];
      
        recvother[iswap] = sendother[iswap];
        nrecvproc[iswap] = noverlap;
        for (i = 0; i < noverlap; i++) recvproc[iswap][i] = overlap[i];
       
       
      
       

      // compute sendbox for each of my sends
      // obox = intersection of ghostbox with other proc's sub-domain
      // sbox = what I need to send to other proc
      //      = sublo to MIN(sublo+cut,subhi) in idim, for idir = 0
      //      = MIN(subhi-cut,sublo) to subhi in idim, for idir = 1
      //      = obox in other 2 dims
      // if sbox touches other proc's sub-box boundaries in lower dims,
      //   extend sbox in those lower dims to include ghost atoms

      double oboxlo[3],oboxhi[3],sbox[6],sbox_multi[6];
      double eoboxlo[3],eoboxhi[3];
      
      for(int init=0; init<noverlap; init++){
        sendbox_flag[iswap][init]=1;
        repeatsend_flag[iswap][init]=0;
      }
      if (mode == Comm::SINGLE) {
      for (i = 0; i < noverlap; i++) {
        pbc_flag[iswap][i] = 0;
        pbc[iswap][i][0] = pbc[iswap][i][1] = pbc[iswap][i][2] =
          pbc[iswap][i][3] = pbc[iswap][i][4] = pbc[iswap][i][5] = 0;
        if(overlap_pbc[i][0]||overlap_pbc[i][1]||overlap_pbc[i][2])
        {
          pbc_flag[iswap][i]=1;
          pbc[iswap][i][0]=overlap_pbc[i][0];
          pbc[iswap][i][1]=overlap_pbc[i][1];
          pbc[iswap][i][2]=overlap_pbc[i][2];
        }
        overlap_counter=i;
        (this->*box_other_full)(0,0,overlap[i],oboxlo,oboxhi);
        oboxlo[0]-=pbc[iswap][i][0]*prd[0];
        oboxlo[1]-=pbc[iswap][i][1]*prd[1];
        oboxlo[2]-=pbc[iswap][i][2]*prd[2];
        oboxhi[0]-=pbc[iswap][i][0]*prd[0];
        oboxhi[1]-=pbc[iswap][i][1]*prd[1];
        oboxhi[2]-=pbc[iswap][i][2]*prd[2];
        eoboxlo[0]=oboxlo[0]-cutghost[0];
        eoboxlo[1]=oboxlo[1]-cutghost[1];
        eoboxlo[2]=oboxlo[2]-cutghost[2];
        eoboxhi[0]=oboxhi[0]+cutghost[0];
        eoboxhi[1]=oboxhi[1]+cutghost[1];
        eoboxhi[2]=oboxhi[2]+cutghost[2];

         for(int idim=0; idim<dimension; idim++){
          if(eoboxlo[idim]<=sublo[idim]){
            sbox[3+idim]=MIN(eoboxhi[idim],subhi[idim]);
            sbox[idim]=sublo[idim];
            if(sbox[3+idim]<=sbox[idim])sbox[3+idim]=sbox[idim];
          }
          else{
            sbox[3+idim]=MIN(eoboxhi[idim],subhi[idim]);
            sbox[idim]=eoboxlo[idim];
            if(sbox[idim]>=sbox[3+idim])sbox[idim]=sbox[3+idim];
          }
          
        }
         
      //determine if sendbox has zero thickness due to a lack of overlap with force cutoff radius
      for(int boxcheck=0; boxcheck<dimension; boxcheck++){
        if(sbox[boxcheck]==sbox[boxcheck+3]) sendbox_flag[iswap][i]=0;
      }
        
        memcpy(overlap_sendbox[iswap][i],sbox,6*sizeof(double));
        memcpy(sendbox[iswap][i],sbox,6*sizeof(double));
      }
      }
      else{
        for (i = 0; i < noverlap; i++) {
        overlap_counter=i;
          
        pbc_flag[iswap][i] = 0;
        pbc[iswap][i][0] = pbc[iswap][i][1] = pbc[iswap][i][2] =
          pbc[iswap][i][3] = pbc[iswap][i][4] = pbc[iswap][i][5] = 0;

        (this->*box_other_full)(0,0,overlap[i],oboxlo,oboxhi);
        }
      }
      
  // reallocate MPI Requests and Statuses as needed

  int nmax = 0;
  for (i = 0; i < nswap; i++) nmax = MAX(nmax,nprocmax[i]);
  for (i = 0; i < dimension; i++) nmax = MAX(nmax,nexchprocmax[i]);
  if (nmax > maxreqstat) {
    maxreqstat = nmax;
    delete [] requests;
    requests = new MPI_Request[maxreqstat];
  }
}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommCAC::forward_comm(int /*dummy*/)
{
  int i,irecv,n,nsend,nrecv;
  AtomVec *avec = atom->avec;
  double **x = atom->x;
  
  // exchange data with another set of procs in each swap
  // post recvs from all procs except self
  // send data to all procs except self
  // copy data to self if sendself is set
  // wait on all procs except self and unpack received data
  // if comm_x_only set, exchange or copy directly to x, don't unpack
 //send first set of ghosts
    int iswap = 0; 
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];

     
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++){
          MPI_Irecv(&buf_recv[recvoffset[iswap][i]],
                    recvsize[iswap][i],
                    MPI_DOUBLE,recvproc[iswap][i],9,world,&requests[i]);
                  
        }
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          n = avec->pack_comm_vel(sendnum[iswap][i],sendlist[iswap][i],
                              buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],9,world);
        }
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
          avec->unpack_comm_vel(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                            &buf_recv[recvoffset[iswap][irecv]]);
        }
      }
      if (sendself[iswap]) {
        for(int selfcount=nsendproc[iswap]-sendself[iswap]; selfcount<nsendproc[iswap]; selfcount++){
        avec->pack_comm_vel(sendnum[iswap][selfcount],sendlist[iswap][selfcount],
                        buf_send,pbc_flag[iswap][selfcount],pbc[iswap][selfcount]);
        avec->unpack_comm_vel(recvnum[iswap][selfcount],firstrecv[iswap][selfcount],
                          buf_send);
        }
      }
    
  
}


/* ----------------------------------------------------------------------
   exchange: move atoms to correct processors
   atoms exchanged with procs that touch sub-box in each of 3 dims
   send out atoms that have left my box, receive ones entering my box
   atoms will be lost if not inside a touching proc's box
     can happen if atom moves outside of non-periodic bounary
     or if atom moves more than one proc away
   this routine called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before exchange is called
------------------------------------------------------------------------- */

void CommCAC::exchange()
{
  int i,m,nexch,nsend,nrecv,nlocal,proc,offset;
  double lo,hi,value,lo_ep,hi_ep;
  double **x = atom->x;
  double ****nodal_positions = atom->nodal_positions;
  double ****initial_nodal_positions = atom->initial_nodal_positions;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int *nodes_count_list = atom->nodes_per_element_list;	
  int nodes_per_element;
  double xcom[3];
  double dx[3];
  nlocal = atom->nlocal;
  int pbc_sign;
  AtomVec *avec = atom->avec;
  // domain properties used in exchange method and methods it calls
  // subbox bounds for orthogonal or triclinic

  prd = domain->prd;
  boxlo = domain->boxlo;
  boxhi = domain->boxhi;

  //resize avec arrays if they're much larger than nlocal
  if(atom->nmax > 1.4*nlocal)
  avec->shrink_array(1.4*nlocal);

  //check for pbc remaps and set nodal positions
  for(i=0; i<nlocal; i++){
  //compute finite element centroid
  nodes_per_element = nodes_count_list[element_type[i]];
  	xcom[0] = 0;
		xcom[1] = 0;
		xcom[2] = 0;
    for(int k=0; k<nodes_per_element; k++){
		  for (int poly_counter = 0; poly_counter < poly_count[i];poly_counter++) {
			
				xcom[0] += nodal_positions[i][k][poly_counter][0];
				xcom[1] += nodal_positions[i][k][poly_counter][1];
				xcom[2] += nodal_positions[i][k][poly_counter][2];

			}
		}
	xcom[0] = xcom[0] / nodes_per_element / poly_count[i];
	xcom[1] = xcom[1] / nodes_per_element / poly_count[i];
	xcom[2] = xcom[2] / nodes_per_element / poly_count[i];

  //test the difference 
  for(int dim=0; dim < dimension; dim++){
  dx[dim] = x[i][dim]-xcom[dim];
  if(dx[dim]>0) pbc_sign = 1;
  else pbc_sign = -1;
  //if the difference exceeds the skin it was almost certainly remapped
  if(dx[dim]>neighbor->skin||dx[dim]<-neighbor->skin){
    for(int k=0; k<nodes_per_element; k++){
		  for (int poly_counter = 0; poly_counter < poly_count[i];poly_counter++) {
		   nodal_positions[i][k][poly_counter][dim] += pbc_sign*prd[dim];
       initial_nodal_positions[i][k][poly_counter][dim] += pbc_sign*prd[dim];
			}
		}
  }
  
  }
  }

  //end of pbc corrections

  // clear global->local map for owned and ghost atoms
  // b/c atoms migrate to new procs in exchange() and
  //   new ghosts are created in borders()
  // map_set() is done at end of borders()
  // clear ghost count and any ghost bonus data internal to AtomVec

  if (map_style) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();

  // insure send buf is large enough for single atom
  // bufextra = max size of one atom = allowed overflow of sendbuf
  // fixes can change per-atom size requirement on-the-fly

  int bufextra_old = bufextra;
  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;
  if (bufextra > bufextra_old)
    memory->grow(buf_send,maxsend+bufextra,"comm:buf_send");

  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // loop over dimensions

  dimension = domain->dimension;

  for (int dim = 0; dim < dimension; dim++) {

    // fill buffer with atoms leaving my box, using < and >=
    // when atom is deleted, fill it in with last atom
    x = atom->x;
    lo = sublo[dim];
    hi = subhi[dim];
    lo_ep = sublo[dim] - BOXEPSILON;
    hi_ep = subhi[dim] + BOXEPSILON;
    nlocal = atom->nlocal;
    i = nsend = 0;

    while (i < nlocal) {
      if (x[i][dim] < lo_ep || x[i][dim] >= hi_ep) {
        if (nsend > maxsend) grow_send(nsend,1);
        proc = (this->*point_drop)(dim,x[i]);
        if (proc != me) {
          buf_send[nsend++] = proc;
          nsend += avec->pack_exchange(i,&buf_send[nsend]);
        }
        avec->copy(nlocal-1,i,1);
        nlocal--;
      } else i++;
    }
    atom->nlocal = nlocal;

    // send and recv atoms from neighbor procs that touch my sub-box in dim
    // no send/recv with self
    // send size of message first
    // receiver may receive multiple messages, realloc buf_recv if needed

    nexch = nexchproc[dim];
    if (!nexch) continue;

    for (m = 0; m < nexch; m++)
      MPI_Irecv(&exchnum[dim][m],1,MPI_INT,
                exchproc[dim][m],7,world,&requests[m]);
    for (m = 0; m < nexch; m++)
      MPI_Send(&nsend,1,MPI_INT,exchproc[dim][m],7,world);
    MPI_Waitall(nexch,requests,MPI_STATUS_IGNORE);

    nrecv = 0;
    for (m = 0; m < nexch; m++) nrecv += exchnum[dim][m];
    if (nrecv > maxrecv) grow_recv(nrecv);

    offset = 0;
    for (m = 0; m < nexch; m++) {
      MPI_Irecv(&buf_recv[offset],exchnum[dim][m],
                MPI_DOUBLE,exchproc[dim][m],8,world,&requests[m]);
      offset += exchnum[dim][m];
    }
    for (m = 0; m < nexch; m++)
      MPI_Send(buf_send,nsend,MPI_DOUBLE,exchproc[dim][m],8,world);
    MPI_Waitall(nexch,requests,MPI_STATUS_IGNORE);

    // check incoming atoms to see if I own it and they are in my box
    // if so, add to my list
    // box check is only for this dimension,
    //   atom may be passed to another proc in later dims

    m = 0;
    while (m < nrecv) {
      proc = static_cast<int> (buf_recv[m++]);
      if (proc == me) {
        value = buf_recv[m+dim+1];
        if (value >= lo && value < hi) {
          m += avec->unpack_exchange(&buf_recv[m]);
          continue;
        }
      }
      m += static_cast<int> (buf_recv[m]);
    }
  }

  if (atom->firstgroupname) atom->first_reorder();
}

/* ----------------------------------------------------------------------
   compute_eboxes: called by borders to compute the current set of eboxes 
------------------------------------------------------------------------- */

void CommCAC::compute_eboxes(int iswap){
int i,j,n,m;
  int ntypes = atom->ntypes;
  double lamda_temp[3];
  double nodal_temp[3];
	double ***nodal_positions;
	int idir, idim;
	int *element_type = atom->element_type;
	int *poly_count = atom->poly_count;
	int *nodes_per_element_list = atom->nodes_per_element_list;
  //double max_search_range= atom->max_search_range;
  nebox_considered=0;
  neboxes=0;
  local_neboxes=0;
  // domain properties used in setup method and methods it calls
  
  dimension = domain->dimension;
  prd = domain->prd;
  boxlo = domain->boxlo;
  boxhi = domain->boxhi;
  sublo = domain->sublo;
  subhi = domain->subhi;
  elimit=atom->nlocal+atom->nghost;
  int *periodicity = domain->periodicity;
 

  //grow bounding boxes if needed
   ebox_ref=memory->grow(atom->ebox_ref,atom->nlocal+atom->nghost,"commCAC: ebox_ref");
  //find maximum element overlap length in each swap direction
  //NOTE: CONVERT TO LAMBDA COORDS FOR TRICLINIC COMPATIBILITY
  for(int element_index=0; element_index<elimit; element_index++){
  	if(element_type[element_index]){
   nodal_positions = atom->nodal_positions[element_index];
	//int current_poly_count = poly_count[element_index];
    int current_poly_count = poly_count[element_index];
	int nodes_per_element = nodes_per_element_list[element_type[element_index]];
  if(neboxes == maxebox){
  maxebox+=BUFEXTRA;
  eboxes=memory->grow(atom->eboxes,maxebox,6,"commCAC: eboxes");
  }



	//initialize bounding box values
	eboxes[neboxes][0] = nodal_positions[0][0][0];
	eboxes[neboxes][1] = nodal_positions[0][0][1];
	eboxes[neboxes][2] = nodal_positions[0][0][2];
	eboxes[neboxes][3] = nodal_positions[0][0][0];
	eboxes[neboxes][4] = nodal_positions[0][0][1];
	eboxes[neboxes][5] = nodal_positions[0][0][2];
  ebox_ref[element_index]=neboxes;
     //define the bounding box for the element being considered as a neighbor
	
	for (int poly_counter = 0; poly_counter < current_poly_count; poly_counter++) {
		for (int kkk = 0; kkk < nodes_per_element; kkk++) {
			for (int dim = 0; dim < 3; dim++) {
				if (nodal_positions[kkk][poly_counter][dim] < eboxes[neboxes][dim]) {
					eboxes[neboxes][dim] = nodal_positions[kkk][poly_counter][dim];
				}
				if (nodal_positions[kkk][poly_counter][dim] > eboxes[neboxes][3+dim]) {
					eboxes[neboxes][3+dim] = nodal_positions[kkk][poly_counter][dim];
				}
			}
		}
	}
  eboxes[neboxes][0] -= cutghost[0];
	eboxes[neboxes][1] -= cutghost[1];
	eboxes[neboxes][2] -= cutghost[2];
	eboxes[neboxes][3] += cutghost[0];
	eboxes[neboxes][4] += cutghost[1];
	eboxes[neboxes][5] += cutghost[2];
 
  neboxes++;
  if(element_index<atom->nlocal) local_neboxes++;
	}
  }
  


  // setup forward/reverse communication
  // loop over 6 swap directions
  // determine which procs I will send to and receive from in each swap
  // done by intersecting ghost box with all proc sub-boxes it overlaps
  // sets nsendproc, nrecvproc, sendproc, recvproc
  // sets sendother, recvother, sendself, pbc_flag, pbc, sendbox
  // resets nprocmax
  atom->neboxes=neboxes;
  atom->local_neboxes=local_neboxes;
}



/* ----------------------------------------------------------------------
   overlap_element_comm: called by borders to obtain information about 
   the bounding box for other procs in recvproc that have been expanded to 
   include element bounding boxes
------------------------------------------------------------------------- */

void CommCAC::get_aug_oboxes(int iswap){
//compute expanded subbox for me
int i,j,n,m;
  int ntypes = atom->ntypes;
  double lamda_temp[3];
  double nodal_temp[3];
	double ***nodal_positions;
	int idir, idim;
	int *element_type = atom->element_type;
	int *poly_count = atom->poly_count;
	int *nodes_per_element_list = atom->nodes_per_element_list;
  //double max_search_range= atom->max_search_range;
  double* current_ebox;
  //double  expanded_subbox[6];
  int nsend,nrecv;
  // domain properties used in setup method and methods it calls
  aug_box[0]=sublo[0];
  aug_box[1]=sublo[1];
  aug_box[2]=sublo[2];
  aug_box[3]=subhi[0];
  aug_box[4]=subhi[1];
  aug_box[5]=subhi[2];
  dimension = domain->dimension;
  prd = domain->prd;
  boxlo = domain->boxlo;
  boxhi = domain->boxhi;
  sublo = domain->sublo;
  subhi = domain->subhi;
  int *periodicity = domain->periodicity;

  // send sendnum counts to procs who recv from me except self
 // copy data to self if sendself is set

    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];
  //initialize expanded subbox

  //find maximum element overlap length in each swap direction
  //NOTE: CONVERT TO LAMBDA COORDS FOR TRICLINIC COMPATIBILITY
  for(int iebox=0; iebox<neboxes; iebox++){
   current_ebox=eboxes[iebox];
   

        for(int idim=0; idim<dimension; idim++){
          if(current_ebox[idim]<aug_box[idim])
          aug_box[idim]=current_ebox[idim];
          if(current_ebox[3+idim]>aug_box[3+idim])
          aug_box[3+idim]=current_ebox[3+idim];
        }
    
    
  }
  
  //send my expanded subbox to all my sendprocs (except self send cases)
    if (recvother[iswap])
    for (m = 0; m < nrecv; m++)
        MPI_Irecv(&aug_oboxes[iswap][m][0],6,MPI_DOUBLE,
                  recvproc[iswap][m],0,world,&requests[m]);
    if (sendother[iswap])
      for (m = 0; m < nsend; m++)
        MPI_Send(&aug_box[0],6,MPI_DOUBLE,sendproc[iswap][m],0,world);
    
    if (recvother[iswap]) MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);
    
    if(sendself[iswap])
    for(int selfcount=nsendproc[iswap]-sendself[iswap]; selfcount<nsendproc[iswap]; selfcount++){
      aug_oboxes[iswap][selfcount][0]=aug_box[0];
      aug_oboxes[iswap][selfcount][1]=aug_box[1];
      aug_oboxes[iswap][selfcount][2]=aug_box[2];
      aug_oboxes[iswap][selfcount][3]=aug_box[3];
      aug_oboxes[iswap][selfcount][4]=aug_box[4];
      aug_oboxes[iswap][selfcount][5]=aug_box[5];
    }
 
}

/* ----------------------------------------------------------------------
   overlap_element_comm: called by borders to obtain information about 
   element bounding boxes that overlap with me but belong to other procs; 
   this information is then used to shape the send box into a more complicated
   structure to reduce the amount of extra ghosts in the send list
------------------------------------------------------------------------- */

void CommCAC::overlap_element_comm(int iswap){
  double lamda_temp[3];
  double nodal_temp[3];
	double ***nodal_positions;
	
	int *element_type = atom->element_type;
	int *poly_count = atom->poly_count;
	double ebounding_boxlo[3];
	double ebounding_boxhi[3];
	int *nodes_per_element_list = atom->nodes_per_element_list;
  //double max_search_range= atom->max_search_range;
  double box_limit;
  int i,m,n,nlast,nsend,nrecv,ngroup,ncount;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double *bbox;
  
  double **x;
  int box_overlap_flag[3];
  
  int iswap_old=iswap;
  double* current_ebox;
  double oboxlo[3],oboxhi[3],sbox[6],sbox_multi[6];

  // domain properties used in setup method and methods it calls
  
  dimension = domain->dimension;
  
  //initialize sendnums to zero
  for (m = 0; m < nsendproc[iswap]; m++) {
    overlap_sendnum[iswap][m]=0;
  }
  
  for (m = 0; m < nrecvproc[iswap]; m++) {
    overlap_recvnum[iswap][m]=0;
  }

  sublo = domain->sublo;
  subhi = domain->subhi;
//send elements that exceed the subbox limit for this swap

 ebox_limit=neboxes;
 for (m = 0; m < nsendproc[iswap]; m++) {
   overlap_counter=m;
         

        //(this->*box_other_full)(0,0,sendproc[iswap][m],oboxlo,oboxhi);
   oboxlo[0]=aug_oboxes[iswap][m][0];
   oboxlo[1]=aug_oboxes[iswap][m][1];
   oboxlo[2]=aug_oboxes[iswap][m][2];
   oboxhi[0]=aug_oboxes[iswap][m][3];
   oboxhi[1]=aug_oboxes[iswap][m][4];
   oboxhi[2]=aug_oboxes[iswap][m][5];
   oboxlo[0]-=pbc[iswap][m][0]*prd[0];
   oboxlo[1]-=pbc[iswap][m][1]*prd[1];
   oboxlo[2]-=pbc[iswap][m][2]*prd[2];
   oboxhi[0]-=pbc[iswap][m][0]*prd[0];
   oboxhi[1]-=pbc[iswap][m][1]*prd[1];
   oboxhi[2]-=pbc[iswap][m][2]*prd[2];

   //determine if aug_obox has any overlap with my aug_box
   box_overlap_flag[2]=box_overlap_flag[1]=box_overlap_flag[0]=0;
     for(int idim=0; idim<dimension; idim++){
          
          if(aug_box[idim]>=oboxlo[idim]&&aug_box[idim]<=oboxhi[idim])
          box_overlap_flag[idim]=1;
          if(aug_box[3+idim]>=oboxlo[idim]&&aug_box[3+idim]<=oboxhi[idim])
          box_overlap_flag[idim]=1;
          //case where sendbox is contained along this dimension within this dimension of ebox
          if(aug_box[idim]<=oboxlo[idim]&&aug_box[3+idim]>=oboxhi[idim])
          box_overlap_flag[idim]=1;
        }
     if(!(box_overlap_flag[0]==1&&box_overlap_flag[1]==1&&box_overlap_flag[2]==1)) continue; 


   for(int iebox=0; iebox<ebox_limit; iebox++){
   current_ebox=eboxes[iebox];
   

	//test if this bounding box exceeds local sub box
  
    //test if element ebox overlaps the sendbox for this comm pair
         box_overlap_flag[2]=box_overlap_flag[1]=box_overlap_flag[0]=0;
        

        for(int idim=0; idim<dimension; idim++){
          
          if(current_ebox[idim]>=oboxlo[idim]&&current_ebox[idim]<=oboxhi[idim])
          box_overlap_flag[idim]=1;
          if(current_ebox[3+idim]>=oboxlo[idim]&&current_ebox[3+idim]<=oboxhi[idim])
          box_overlap_flag[idim]=1;
          //case where sendbox is contained along this dimension within this dimension of ebox
          if(current_ebox[idim]<=oboxlo[idim]&&current_ebox[3+idim]>=oboxhi[idim])
          box_overlap_flag[idim]=1;
        }
        if(box_overlap_flag[0]==1&&box_overlap_flag[1]==1&&box_overlap_flag[2]==1){
         
        if (overlap_sendnum[iswap][m] == overlap_maxsendlist[iswap][m]) overlap_grow_list(iswap,m,overlap_sendnum[iswap][m]);
            overlap_sendlist[iswap][m][overlap_sendnum[iswap][m]++] = iebox;
        
        }
    }
    
  }


    // send sendnum counts to procs who recv from me except self
    // copy data to self if sendself is set

    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];

    if (recvother[iswap])
      for (m = 0; m < nrecv; m++)
        MPI_Irecv(&overlap_recvnum[iswap][m],1,MPI_INT,
                  recvproc[iswap][m],1,world,&requests[m]);
    if (sendother[iswap])
      for (m = 0; m < nsend; m++)
        MPI_Send(&overlap_sendnum[iswap][m],1,MPI_INT,sendproc[iswap][m],1,world);
    //if (sendself[iswap]){ overlap_recvnum[iswap][nrecv] = overlap_sendnum[iswap][nsend];}
    if (sendself[iswap]) 
        for(int selfcount=nsendproc[iswap]-sendself[iswap]; selfcount<nsendproc[iswap]; selfcount++)
        overlap_recvnum[iswap][selfcount]=overlap_sendnum[iswap][selfcount];
    if (recvother[iswap]) MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);

    // setup other per swap/proc values from sendnum and recvnum

    for (m = 0; m < nsendproc[iswap]; m++) {
      size_reverse_recv[iswap][m] = overlap_sendnum[iswap][m]*size_reverse;
      if (m == 0) reverse_recv_offset[iswap][0] = 0;
      else reverse_recv_offset[iswap][m] =
             reverse_recv_offset[iswap][m-1] + overlap_sendnum[iswap][m-1];
    }

   
    for (m = 0; m < nrecvproc[iswap]; m++) {
     
      if (m == 0) {
        overlap_firstrecv[iswap][0] = nforeign_eboxes;
      } else {
        overlap_firstrecv[iswap][m] = overlap_firstrecv[iswap][m-1] + overlap_recvnum[iswap][m-1];
      }
    }
    
    // insure send/recv buffers are large enough for this border comm swap

    //if (smaxone*size_border > maxsend) grow_send(smaxone*size_border,0);
    //if (rmaxall*size_border > maxrecv) grow_recv(rmaxall*size_border);

    // swap atoms with other procs using pack_border(), unpack_border()
    // use Waitall() instead of Waitany() because calls to unpack_border()
    //   must increment per-atom arrays in ascending order

    //compute buffer sizes for each of the nsend communications for this task
     int send_accumulation=0;
    
     if (sendother[iswap]) {
        for (m = 0; m < nsend; m++) {
          overlap_sendoffset[iswap][m]=send_accumulation;
          for (int sendcounter = 0; sendcounter < overlap_sendnum[iswap][m]; sendcounter++) {
          
          
          if (send_accumulation+overlap_sendsize[iswap][m] > maxsend) grow_send(send_accumulation+overlap_sendsize[iswap][m],1);
          overlap_sendsize[iswap][m] += pack_eboxes(1,&overlap_sendlist[iswap][m][sendcounter],
                                &buf_send[overlap_sendoffset[iswap][m]+overlap_sendsize[iswap][m]],pbc_flag[iswap][m],pbc[iswap][m],iswap);
          }
          //compute offsets in buffer index for each proc to send to; i.e. there is one buffer for all sends
          send_accumulation+=overlap_sendsize[iswap][m];
        }
      }
      //receive buffer sizes
      
      if (recvother[iswap]) {
        for (m = 0; m < nrecv; m++)
          MPI_Irecv(&overlap_recvsize[iswap][m],
                    1, MPI_INT,recvproc[iswap][m],2,world,&requests[m]);
      }
      //send buffer sizes to recv procs
     if (sendother[iswap]) {
        for (m = 0; m < nsend; m++) {
          MPI_Send(&overlap_sendsize[iswap][m],1,
          MPI_INT,sendproc[iswap][m],2,world);
        }
      }
      //wait for sizes to be received and resize recv buffer accordingly
      int total_recvsize=0;
      if (recvother[iswap]) {
        MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);
        for (m = 0; m < nrecv; m++){
          overlap_recvoffset[iswap][m]=total_recvsize;
          total_recvsize+=overlap_recvsize[iswap][m];
          
        }
         if (total_recvsize > maxrecv) grow_recv(total_recvsize); 
      }
      if (recvother[iswap]) {
        for (m = 0; m < nrecv; m++)
          MPI_Irecv(&buf_recv[overlap_recvoffset[iswap][m]],
                    overlap_recvsize[iswap][m],
                    MPI_DOUBLE,recvproc[iswap][m],3,world,&requests[m]);
      }
      if (sendother[iswap]) {
        for (m = 0; m < nsend; m++) {
          MPI_Send(&buf_send[overlap_sendoffset[iswap][m]],overlap_sendsize[iswap][m],
          MPI_DOUBLE,sendproc[iswap][m],3,world);
        }
      }
      if (recvother[iswap]) {
        MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);
        for (m = 0; m < nrecv; m++)
          unpack_eboxes(overlap_recvnum[iswap][m],overlap_firstrecv[iswap][m],
                              &buf_recv[overlap_recvoffset[iswap][m]]);
      }
     
      if (sendself[iswap]) {
        
        for(int selfcount=nsendproc[iswap]-sendself[iswap]; selfcount<nsendproc[iswap]; selfcount++){
          int self_accumulation=0;
          for (int sendcounter = 0; sendcounter < overlap_sendnum[iswap][selfcount]; sendcounter++) {
          
          
          if (self_accumulation+overlap_sendsize[iswap][selfcount] > maxsend) grow_send(self_accumulation+overlap_sendsize[iswap][selfcount],1);
          overlap_sendsize[iswap][selfcount] += pack_eboxes(1,&overlap_sendlist[iswap][selfcount][sendcounter],
                                &buf_send[sendsize[iswap][selfcount]],pbc_flag[iswap][selfcount],pbc[iswap][selfcount],iswap);
          
          
          }
          //compute offsets in buffer index for each proc to send to; i.e. there is one buffer for all sends
          self_accumulation+=overlap_sendsize[iswap][m];
        unpack_eboxes(overlap_recvnum[iswap][selfcount],overlap_firstrecv[iswap][selfcount],
                            buf_send);
        }
      }
}

/* ----------------------------------------------------------------------
   borders: list nearby atoms to send to neighboring procs at every timestep
   one list is created per swap/proc that will be made
   as list is made, actually do communication
   this does equivalent of a forward_comm(), so don't need to explicitly
     call forward_comm() on reneighboring timestep
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
------------------------------------------------------------------------- */

void CommCAC::borders()
{
  int i,m,n,nlast,nsend,nrecv,ngroup,ncount,ncountall;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double *bbox;
  double **x;
  
  int max_sendaccumulation=0;
  int max_recvaccumulation=0;
  AtomVec *avec = atom->avec;

  //flag indicating whether to zero out arrays for other comms routines
  reset_array_flag=1;
  
  // send/recv max one = max # of atoms in single send/recv for any swap
  // send/recv max all = max # of atoms in all sends/recvs within any swap

  smaxone = smaxall = 0;
  rmaxone = rmaxall = 0;
  //zero out arrays
   int iswap = 0; 
   for (m = 0; m < nsendproc[iswap]; m++) {
       sendsize[iswap][m]=0;
       
       sendnum[iswap][m]=0;
      
       sendoffset[iswap][m]=0;
       
       overlap_sendsize[iswap][m]=0;
      
       overlap_sendoffset[iswap][m]=0;
     }
   
   
    //zero out arrays
   
     for (m = 0; m < nrecvproc[iswap]; m++) {
       
       recvsize[iswap][m]=0;
       
       recvnum[iswap][m]=0;
       
       recvoffset[iswap][m]=0;
       
       overlap_recvsize[iswap][m]=0;
       
       overlap_recvoffset[iswap][m]=0;
       
     }
  //check if buffer is sized large enough 
  int bufextra_old = bufextra;
  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;
  if (bufextra > bufextra_old)
   memory->grow(buf_send,maxsend+bufextra,"comm:buf_send");
  
  nforeign_eboxes=0;
  //communicate all overlapping elements in all 6 swaps first to modify all sendboxes as needed
  compute_eboxes(iswap);

  //communicate overlapping elements
  get_aug_oboxes(iswap);
  overlap_element_comm(iswap);
  atom->nforeign_eboxes=nforeign_eboxes;
  atom->foreign_eboxes=foreign_eboxes;
  atom->bin_foreign=1;
  bin_pointer->bin_atoms(); 
  atom->bin_foreign=0;
  stencil_pointer->post_create_setup();
  stencil_pointer->post_create();
  bin_ncontent=bin_pointer->bin_ncontent;
  bin_content=bin_pointer->bin_content;
  nstencil=stencil_pointer->nstencil;
  stencil=stencil_pointer->stencil;
  nbin_element_overlap=bin_pointer->nbin_element_overlap;  //array storing the number of bins this element overlaps
  bin_element_overlap=bin_pointer->bin_element_overlap;  //set of bins this element overlaps
  //received ghosts may change necessary proc overlaps in subsequent swaps
  
    // find atoms within sendboxes using >= and <
    // hi test with ">" is important b/c don't want to send an atom
    //   in lower dim (on boundary) that a proc will recv again in higher dim
    // for x-dim swaps, check owned atoms
    // for yz-dim swaps, check owned and ghost atoms
    // store sent atom indices in sendlist for use in future timesteps
    // NOTE: assume SINGLE mode, add logic for MULTI mode later
    
    x = atom->x;
    nlast = atom->nlocal;
     //check is recv flag array is sized properly
    
    ncountall = 0;
    for (m = 0; m < nsendproc[iswap]; m++) {
        
      if (mode == Comm::SINGLE) {
      ncount = 0; 
      if((!sendbox_flag[iswap][m]&&overlap_recvnum[iswap][m]==0)&&overlap_sendnum==0) continue;
      bbox = sendbox[iswap][m];
      xlo = bbox[0]; ylo = bbox[1]; zlo = bbox[2];
      xhi = bbox[3]; yhi = bbox[4]; zhi = bbox[5];

      

      if (!bordergroup) {
        for (i = 0; i < nlast; i++) {
          if (sendbox_include(iswap, m, i)) {
            if (ncount == maxsendlist[iswap][m]) grow_list(iswap,m,ncount);
            sendlist[iswap][m][ncount++] = i;
          }
        }
      } else {
        ngroup = atom->nfirst;
        for (i = 0; i < ngroup; i++) {
          if (sendbox_include(iswap, m, i)) {
            if (ncount == maxsendlist[iswap][m]) grow_list(iswap,m,ncount);
            sendlist[iswap][m][ncount++] = i;
          }
        }
        for (i = atom->nlocal; i < nlast; i++) {
          if (sendbox_include(iswap, m, i)) {
            if (ncount == maxsendlist[iswap][m]) grow_list(iswap,m,ncount);
            sendlist[iswap][m][ncount++] = i;
          }
        }
      }

      sendnum[iswap][m] = ncount;
      smaxone = MAX(smaxone,ncount);
      ncountall += ncount;
    }
    else{
      int* type=atom->type;
      int itype;
      ncount = 0;

      if (!bordergroup) {
        for (i = 0; i < nlast; i++) {
          itype=type[i];    
          bbox = sendbox_multi[iswap][m][itype];
          xlo = bbox[0]; ylo = bbox[1]; zlo = bbox[2];
          xhi = bbox[3]; yhi = bbox[4]; zhi = bbox[5];
          if (sendbox_include(iswap, m, i)) {
            if (ncount == maxsendlist[iswap][m]) grow_list(iswap,m,ncount);
            sendlist[iswap][m][ncount++] = i;
          }
          
        }
      } else {
        ngroup = atom->nfirst;
        for (i = 0; i < ngroup; i++) {
          itype=type[i];    
          bbox = sendbox_multi[iswap][m][itype];
          xlo = bbox[0]; ylo = bbox[1]; zlo = bbox[2];
          xhi = bbox[3]; yhi = bbox[4]; zhi = bbox[5];
          if (sendbox_include(iswap, m, i)) {
           if (ncount == maxsendlist[iswap][m]) grow_list(iswap,m,ncount);
            sendlist[iswap][m][ncount++] = i;
          }
        }
        for (i = atom->nlocal; i < nlast; i++) {
           itype=type[i];    
           bbox = sendbox_multi[iswap][m][itype];
           xlo = bbox[0]; ylo = bbox[1]; zlo = bbox[2];
           xhi = bbox[3]; yhi = bbox[4]; zhi = bbox[5];
          if (sendbox_include(iswap, m, i)) {
            if (ncount == maxsendlist[iswap][m]) grow_list(iswap,m,ncount);
            sendlist[iswap][m][ncount++] = i;
          }
        }
      }

      sendnum[iswap][m] = ncount;
      smaxone = MAX(smaxone,ncount);
      ncountall += ncount;
    }
    }
    smaxall = MAX(smaxall,ncountall);

    // send sendnum counts to procs who recv from me except self
    // copy data to self if sendself is set

    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];

    if (recvother[iswap])
      for (m = 0; m < nrecv; m++)
        MPI_Irecv(&recvnum[iswap][m],1,MPI_INT,
                  recvproc[iswap][m],4,world,&requests[m]);
    if (sendother[iswap])
      for (m = 0; m < nsend; m++)
        MPI_Send(&sendnum[iswap][m],1,MPI_INT,sendproc[iswap][m],4,world);
    //if (sendself[iswap]) recvnum[iswap][nrecv] = sendnum[iswap][nsend];//look at pbc later
    if (sendself[iswap]) 
        for(int selfcount=nsendproc[iswap]-sendself[iswap]; selfcount<nsendproc[iswap]; selfcount++)
        recvnum[iswap][selfcount]=sendnum[iswap][selfcount];
    if (recvother[iswap]) MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);

    // setup other per swap/proc values from sendnum and recvnum
    
    //check on this later
    for (m = 0; m < nsendproc[iswap]; m++) {
      size_reverse_recv[iswap][m] = sendnum[iswap][m]*size_reverse;
      if (m == 0) reverse_recv_offset[iswap][0] = 0;
      else reverse_recv_offset[iswap][m] =
             reverse_recv_offset[iswap][m-1] + sendnum[iswap][m-1];
    }

    ncountall = 0;
    
    for (m = 0; m < nrecvproc[iswap]; m++) {
      ncount = recvnum[iswap][m];
      rmaxone = MAX(rmaxone,ncount);
      ncountall += ncount;
      
      if (m == 0) {
        firstrecv[iswap][0] = atom->nlocal + atom->nghost;
        forward_recv_offset[iswap][0] = 0;
      } else {
        firstrecv[iswap][m] = firstrecv[iswap][m-1] + recvnum[iswap][m-1];
        forward_recv_offset[iswap][m] =
          forward_recv_offset[iswap][m-1] + recvnum[iswap][m-1];
      }
    }
    rmaxall = MAX(rmaxall,ncountall);

    // insure send/recv buffers are large enough for this border comm swap
    // swap atoms with other procs using pack_border(), unpack_border()
    // use Waitall() instead of Waitany() because calls to unpack_border()
    //   must increment per-atom arrays in ascending order

    //compute buffer sizes for each of the nsend communications for this task
     int send_accumulation=0;
    
     if (sendother[iswap]) {
        for (m = 0; m < nsend; m++) {
          sendoffset[iswap][m]=send_accumulation;
          for (int sendcounter = 0; sendcounter < sendnum[iswap][m]; sendcounter++) {
          
          
          if (send_accumulation+sendsize[iswap][m] > maxsend) grow_send(send_accumulation+sendsize[iswap][m],1);
          sendsize[iswap][m] += avec->pack_border(1,&sendlist[iswap][m][sendcounter],
                                &buf_send[sendoffset[iswap][m]+sendsize[iswap][m]],pbc_flag[iswap][m],pbc[iswap][m]);
          
          
          }
          //compute offsets in buffer index for each proc to send to; i.e. there is one buffer for all sends
          send_accumulation+=sendsize[iswap][m];
          
        }
      }
      //receive buffer sizes
      
      if (recvother[iswap]) {
        for (m = 0; m < nrecv; m++)
          MPI_Irecv(&recvsize[iswap][m],
                    1, MPI_INT,recvproc[iswap][m],5,world,&requests[m]);
      }
      //send buffer sizes to recv procs
     if (sendother[iswap]) {
        for (m = 0; m < nsend; m++) {
          MPI_Send(&sendsize[iswap][m],1,
          MPI_INT,sendproc[iswap][m],5,world);
        }
      }
      //wait for sizes to be received and resize recv buffer accordingly
      int total_recvsize=0;
      if (recvother[iswap]) {
        MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);
        for (m = 0; m < nrecv; m++){
          recvoffset[iswap][m]=total_recvsize;
          total_recvsize+=recvsize[iswap][m];
          
        }
         if (total_recvsize > maxrecv) grow_recv(total_recvsize); 
      }
      if (recvother[iswap]) {
        for (m = 0; m < nrecv; m++)
          MPI_Irecv(&buf_recv[recvoffset[iswap][m]],
                    recvsize[iswap][m],
                    MPI_DOUBLE,recvproc[iswap][m],6,world,&requests[m]);
      }
      if (sendother[iswap]) {
        for (m = 0; m < nsend; m++) {
          MPI_Send(&buf_send[sendoffset[iswap][m]],sendsize[iswap][m],
          MPI_DOUBLE,sendproc[iswap][m],6,world);
        }
      }
      if (recvother[iswap]) {
        MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);
        for (m = 0; m < nrecv; m++)
          avec->unpack_border(recvnum[iswap][m],firstrecv[iswap][m],
                              &buf_recv[recvoffset[iswap][m]]);
      }
      if (sendself[iswap]) {
        
        for(int selfcount=nsendproc[iswap]-sendself[iswap]; selfcount<nsendproc[iswap]; selfcount++){
          int self_accumulation=0;
          for (int sendcounter = 0; sendcounter < sendnum[iswap][selfcount]; sendcounter++) {
          
          
          if (self_accumulation+sendsize[iswap][selfcount] > maxsend) grow_send(self_accumulation+sendsize[iswap][selfcount],1);
          sendsize[iswap][selfcount] += avec->pack_border(1,&sendlist[iswap][selfcount][sendcounter],
                                &buf_send[sendsize[iswap][selfcount]],pbc_flag[iswap][selfcount],pbc[iswap][selfcount]);
          
          
          }
          //compute offsets in buffer index for each proc to send to; i.e. there is one buffer for all sends
          self_accumulation+=sendsize[iswap][m];
        avec->unpack_border(recvnum[iswap][selfcount],firstrecv[iswap][selfcount],
                            buf_send);
        }
      }
    
    //compute maximum buffer size so far
    max_sendaccumulation = MAX(max_sendaccumulation,send_accumulation);
    max_recvaccumulation = MAX(max_recvaccumulation,total_recvsize);
    // increment ghost atoms

    n = nrecvproc[iswap];
    if (n){
      int old_nghost=atom->nghost;
      atom->nghost += forward_recv_offset[iswap][n-1] + recvnum[iswap][n-1];
      if(atom->nlocal + atom->nghost>maxall){ 
       maxall = atom->nlocal + atom->nghost+BUFEXTRA;
       memory->grow(recv_flag,maxall,"commCAC:recv_flag");
      }
    for(int m=0; m<nrecvproc[iswap]; m++){
    for(int k=atom->nlocal+forward_recv_offset[iswap][m]; k<atom->nlocal+forward_recv_offset[iswap][m]+recvnum[iswap][m]; k++){
      recv_flag[k]=recvproc[iswap][m];
    }
    }
    }
    

  // insure send/recv buffers are long enough for all forward & reverse comm
  // send buf is for one forward or reverse sends to one proc
  // recv buf is for all forward or reverse recvs in one swap

  int max = MAX(maxsend,max_sendaccumulation);
  if (max > maxsend) grow_send(max,0);
  max = MAX(maxrecv,max_recvaccumulation);
  if (max > maxrecv) grow_recv(max);

  
  double lamda_temp[3];
  double nodal_temp[3];
	double ***nodal_positions;

	int *element_type = atom->element_type;
	int *poly_count = atom->poly_count;
	int *nodes_per_element_list = atom->nodes_per_element_list;
  
  nebox_considered=0;
  neboxes=0;
  local_neboxes=0;
  ebox_ref=memory->grow(atom->ebox_ref,atom->nlocal+atom->nghost,"commCAC: ebox_ref");
  //compute final set of eboxes in me
  //find maximum element overlap length in each swap direction
  //NOTE: CONVERT TO LAMBDA COORDS FOR TRICLINIC COMPATIBILITY
  for(int element_index=0; element_index<atom->nlocal+atom->nghost; element_index++){
  	if(element_type[element_index]){
   nodal_positions = atom->nodal_positions[element_index];
    int current_poly_count = poly_count[element_index];
	int nodes_per_element = nodes_per_element_list[element_type[element_index]];
  if(neboxes == maxebox){
  maxebox+=BUFEXTRA;
  eboxes=memory->grow(atom->eboxes,maxebox,6,"commCAC: eboxes");
  }



	//initialize bounding box values
	eboxes[neboxes][0] = nodal_positions[0][0][0];
	eboxes[neboxes][1] = nodal_positions[0][0][1];
	eboxes[neboxes][2] = nodal_positions[0][0][2];
	eboxes[neboxes][3] = nodal_positions[0][0][0];
	eboxes[neboxes][4] = nodal_positions[0][0][1];
	eboxes[neboxes][5] = nodal_positions[0][0][2];
  ebox_ref[element_index]=neboxes;
     //define the bounding box for the element being considered as a neighbor
	
	for (int poly_counter = 0; poly_counter < current_poly_count; poly_counter++) {
		for (int kkk = 0; kkk < nodes_per_element; kkk++) {
			for (int dim = 0; dim < 3; dim++) {
				if (nodal_positions[kkk][poly_counter][dim] < eboxes[neboxes][dim]) {
					eboxes[neboxes][dim] = nodal_positions[kkk][poly_counter][dim];
				}
				if (nodal_positions[kkk][poly_counter][dim] > eboxes[neboxes][3+dim]) {
					eboxes[neboxes][3+dim] = nodal_positions[kkk][poly_counter][dim];
				}
			}
		}
	}
  
  eboxes[neboxes][0] -= cutghost[0];
	eboxes[neboxes][1] -= cutghost[1];
	eboxes[neboxes][2] -= cutghost[2];
	eboxes[neboxes][3] += cutghost[0];
	eboxes[neboxes][4] += cutghost[1];
	eboxes[neboxes][5] += cutghost[2];

  neboxes++;
  if(element_index<atom->nlocal) local_neboxes++;
	}
  }
  atom->neboxes=neboxes;
  atom->local_neboxes=local_neboxes;
 
  // reset global->local map

  if (map_style) atom->map_set();
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Pair
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommCAC::forward_comm_pair(Pair *pair)
{
  int i,irecv,n,nsend,nrecv;

  int nsize = pair->comm_forward;

  int iswap = 0;
  nsend = nsendproc[iswap] - sendself[iswap];
  nrecv = nrecvproc[iswap] - sendself[iswap];

  //zero out size and offset arrays
  if(reset_array_flag) pair_comm_setup(pair);
  
  if(!reset_array_flag){
    if (recvother[iswap]) {
      for (i = 0; i < nrecv; i++)
        MPI_Irecv(&buf_recv[pair_recvoffset[iswap][i]],
                  pair_recvsize[iswap][i],
                  MPI_DOUBLE,recvproc[iswap][i],12,world,&requests[i]);
    }

    if (sendother[iswap]) {
      for (i = 0; i < nsend; i++) {
        n = pair->pack_forward_comm(sendnum[iswap][i],sendlist[iswap][i],
                                    buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
        MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],12,world);
      }
    }

    if (sendself[iswap]) {
        for(int selfcount=nsendproc[iswap]-sendself[iswap]; selfcount<nsendproc[iswap]; selfcount++){
        pair->pack_forward_comm(sendnum[iswap][selfcount],sendlist[iswap][selfcount],
                        buf_send,pbc_flag[iswap][selfcount],pbc[iswap][selfcount]);
        pair->unpack_forward_comm(recvnum[iswap][selfcount],firstrecv[iswap][selfcount],
                          buf_send);
        }
    }
  }
    if (recvother[iswap]) {
      for (i = 0; i < nrecv; i++) {
        MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
        pair->unpack_forward_comm(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                                  &buf_recv[pair_recvoffset[iswap][irecv]]);
      }
    }
   reset_array_flag=0;
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Pair
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommCAC::pair_comm_setup(Pair *pair)
{
  int i,irecv,n, m, nsend,nrecv;
  AtomVec *avec = atom->avec;
  double **x = atom->x;
  
  // exchange data with another set of procs in each swap
  // post recvs from all procs except self
  // send data to all procs except self
  // copy data to self if sendself is set
  // wait on all procs except self and unpack received data
  // if comm_x_only set, exchange or copy directly to x, don't unpack
  //send first set of ghosts
    int iswap = 0; 
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];
  for (m = 0; m < nsendproc[iswap]; m++) {
    pair_sendsize[iswap][m]=0;
    pair_sendoffset[iswap][m]=0;
  }
   
  for (m = 0; m < nrecvproc[iswap]; m++) {
    pair_recvsize[iswap][m]=0;
    pair_recvoffset[iswap][m]=0; 
  } 

  //compute buffer sizes for each of the nsend communications for this task
     int send_accumulation=0;
    
     if (sendother[iswap]) {
       for (m = 0; m < nsend; m++) {
         pair_sendoffset[iswap][m]=send_accumulation;
           for (int sendcounter = 0; sendcounter < sendnum[iswap][m]; sendcounter++) {
           if (send_accumulation+pair_sendsize[iswap][m] > maxsend) grow_send(send_accumulation+pair_sendsize[iswap][m],1);
           pair_sendsize[iswap][m] += pair->pack_forward_comm(1,&sendlist[iswap][m][sendcounter],
             &buf_send[pair_sendoffset[iswap][m]+pair_sendsize[iswap][m]],pbc_flag[iswap][m],pbc[iswap][m]);
           }
          //compute offsets in buffer index for each proc to send to; i.e. there is one buffer for all sends
          send_accumulation+=pair_sendsize[iswap][m];
        }
      }
      //receive buffer sizes
      
      if (recvother[iswap]) {
        for (m = 0; m < nrecv; m++)
          MPI_Irecv(&pair_recvsize[iswap][m],
                    1, MPI_INT,recvproc[iswap][m],10,world,&requests[m]);
      }
      //send buffer sizes to recv procs
     if (sendother[iswap]) {
        for (m = 0; m < nsend; m++) {
          MPI_Send(&pair_sendsize[iswap][m],1,
          MPI_INT,sendproc[iswap][m],10,world);
        }
      }
      //wait for sizes to be received and resize recv buffer accordingly
      int total_recvsize=0;
      if (recvother[iswap]) {
        MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);
        for (m = 0; m < nrecv; m++){
          pair_recvoffset[iswap][m]=total_recvsize;
          total_recvsize+=pair_recvsize[iswap][m];
          
        }
         if (total_recvsize > maxrecv) grow_recv(total_recvsize); 
      }
      if (recvother[iswap]) {
        for (m = 0; m < nrecv; m++)
          MPI_Irecv(&buf_recv[pair_recvoffset[iswap][m]],
                    pair_recvsize[iswap][m],
                    MPI_DOUBLE,recvproc[iswap][m],11,world,&requests[m]);
      }
      if (sendother[iswap]) {
        for (m = 0; m < nsend; m++) {
          MPI_Send(&buf_send[pair_sendoffset[iswap][m]],pair_sendsize[iswap][m],
          MPI_DOUBLE,sendproc[iswap][m],11,world);
        }
      } 

      if (sendself[iswap]) {
        
        for(int selfcount=nsendproc[iswap]-sendself[iswap]; selfcount<nsendproc[iswap]; selfcount++){
          int self_accumulation=0;
          for (int sendcounter = 0; sendcounter < sendnum[iswap][selfcount]; sendcounter++) {
          if (self_accumulation+pair_sendsize[iswap][selfcount] > maxsend) grow_send(self_accumulation+pair_sendsize[iswap][selfcount],1);
          pair_sendsize[iswap][selfcount] += pair->pack_forward_comm(1,&sendlist[iswap][selfcount][sendcounter],
                                &buf_send[pair_sendsize[iswap][selfcount]],pbc_flag[iswap][selfcount],pbc[iswap][selfcount]);
          }
          //compute offsets in buffer index for each proc to send to; i.e. there is one buffer for all sends
          self_accumulation+=pair_sendsize[iswap][m];
          pair->unpack_forward_comm(recvnum[iswap][selfcount],firstrecv[iswap][selfcount],
                            buf_send);
        }
      }

}

/* ----------------------------------------------------------------------
   decides if element or atom should be sent to a neighboring proc in swap
------------------------------------------------------------------------- */

int CommCAC::sendbox_include(int iswap, int m, int current_element)
{
	int flag=0;
	int *element_type = atom->element_type;
	int *poly_count = atom->poly_count;
	double ebounding_boxlo[3];
	double ebounding_boxhi[3];
	int *nodes_per_element_list = atom->nodes_per_element_list;
  //double max_search_range= atom->max_search_range;
  double box_limit;
  int i,n,nlast,nsend,nrecv,ngroup,ncount,ncountall;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double *bbox;
  double **x=atom->x;
  int box_overlap_flag[3];
  double exp_ebox[6];
  double reduced_ebox[6];
  double contained_ebox[6];
  double oboxlo[3],oboxhi[3];
  double eoboxlo[3],eoboxhi[3];
  double x_send[3];
  int lo[3],hi[3];
  int ebox_id;
  int iebox;
  int swapdim=iswap/2;
  overlap_counter=m;
  
  if(element_type[current_element]!=0){
  ebox_id=ebox_ref[current_element];
  reduced_ebox[0]=eboxes[ebox_id][0]+cutghost[0];
  reduced_ebox[1]=eboxes[ebox_id][1]+cutghost[1];
  reduced_ebox[2]=eboxes[ebox_id][2]+cutghost[2];
  reduced_ebox[3]=eboxes[ebox_id][3]-cutghost[0];
  reduced_ebox[4]=eboxes[ebox_id][4]-cutghost[1];
  reduced_ebox[5]=eboxes[ebox_id][5]-cutghost[2];
  }
   bbox = sendbox[iswap][m];
      xlo = bbox[0]-BOXEPSILON; ylo = bbox[1]-BOXEPSILON; zlo = bbox[2]-BOXEPSILON;
      xhi = bbox[3]+BOXEPSILON; yhi = bbox[4]+BOXEPSILON; zhi = bbox[5]+BOXEPSILON;
  
  (this->*box_other_full)(0,0,sendproc[iswap][m],oboxlo,oboxhi);
  oboxlo[0]-=pbc[iswap][m][0]*prd[0];
  oboxlo[1]-=pbc[iswap][m][1]*prd[1];
  oboxlo[2]-=pbc[iswap][m][2]*prd[2];
  oboxhi[0]-=pbc[iswap][m][0]*prd[0];
  oboxhi[1]-=pbc[iswap][m][1]*prd[1];
  oboxhi[2]-=pbc[iswap][m][2]*prd[2];
  eoboxlo[0]=oboxlo[0]-cutghost[0];
  eoboxlo[1]=oboxlo[1]-cutghost[1];
  eoboxlo[2]=oboxlo[2]-cutghost[2];
  eoboxhi[0]=oboxhi[0]+cutghost[0];
  eoboxhi[1]=oboxhi[1]+cutghost[1];
  eoboxhi[2]=oboxhi[2]+cutghost[2];
      
      //partial test to help make sure eboxes are placed correctly; use max ebox overlap in me

      //determine if current element ebox or point particle position overlaps or lies cutghost sendbox
      if(!element_type[current_element]){
      if (x[current_element][0] >= xlo && x[current_element][0] < xhi &&
          x[current_element][1] >= ylo && x[current_element][1] < yhi &&
          x[current_element][2] >= zlo && x[current_element][2] < zhi) return 1;
      }
      else{
        
        box_overlap_flag[2]=box_overlap_flag[1]=box_overlap_flag[0]=0;
        for(int idim=0; idim<dimension; idim++){
          if(reduced_ebox[idim]>=eoboxlo[idim]&&reduced_ebox[idim]<=eoboxhi[idim])
          box_overlap_flag[idim]=1;
          if(reduced_ebox[3+idim]>=eoboxlo[idim]&&reduced_ebox[3+idim]<=eoboxhi[idim])
          box_overlap_flag[idim]=1;
          if(reduced_ebox[idim]<=eoboxlo[idim]&&reduced_ebox[3+idim]>=eoboxhi[idim])
          box_overlap_flag[idim]=1;
        }
        if(box_overlap_flag[0]==1&&box_overlap_flag[1]==1&&box_overlap_flag[2]==1){
         return 1;  
        }
      }
    
     //loop through communicated eboxes to determine if this particle is needed as a ghost for that ebox
     //belonging to another task
     
     for (int ibin_counter=0; ibin_counter<nbin_element_overlap[current_element]; ibin_counter++){
       int ibin=bin_element_overlap[current_element][ibin_counter];
     for (int k = 0; k < nstencil; k++) {
     for (int jj = 0; jj < bin_ncontent[ibin + stencil[k]]; jj++) {
			if(ibin + stencil[k]<0) error->one(FLERR," negative bin index");
			if(ibin + stencil[k]>=bin_pointer->mbins) error->one(FLERR," excessive bin index");
		  iebox = bin_content[ibin + stencil[k]][jj];
      
			
       //check if this ebox is relevant to this sendbox first
       if(foreign_eprocs[iebox]!=sendproc[iswap][m]) continue; 
       if(foreign_image[iebox][0]!=-pbc[iswap][m][0]||foreign_image[iebox][1]!=-pbc[iswap][m][1]
       ||foreign_image[iebox][2]!=-pbc[iswap][m][2]) continue; 
       exp_ebox[0]=foreign_eboxes[iebox][0];
       exp_ebox[1]=foreign_eboxes[iebox][1];
       exp_ebox[2]=foreign_eboxes[iebox][2];
       exp_ebox[3]=foreign_eboxes[iebox][3];
       exp_ebox[4]=foreign_eboxes[iebox][4];
       exp_ebox[5]=foreign_eboxes[iebox][5];
      
       if(!element_type[current_element]){
         box_overlap_flag[2]=box_overlap_flag[1]=box_overlap_flag[0]=0;
         for(int idim=0; idim<dimension; idim++){
           if(x[current_element][idim]>=exp_ebox[idim]&&x[current_element][idim]<=exp_ebox[3+idim])
           box_overlap_flag[idim]=1;;
         }
         if(box_overlap_flag[0]==1&&box_overlap_flag[1]==1&&box_overlap_flag[2]==1)
         return 1;
       }
       else{
        box_overlap_flag[2]=box_overlap_flag[1]=box_overlap_flag[0]=0;
        for(int idim=0; idim<dimension; idim++){
          if(reduced_ebox[idim]>=exp_ebox[idim]&&reduced_ebox[idim]<=exp_ebox[3+idim])
          box_overlap_flag[idim]=1;
          if(reduced_ebox[3+idim]>=exp_ebox[idim]&&reduced_ebox[3+idim]<=exp_ebox[3+idim])
          box_overlap_flag[idim]=1;
          if(reduced_ebox[idim]<=exp_ebox[idim]&&reduced_ebox[3+idim]>=exp_ebox[3+idim])
          box_overlap_flag[idim]=1;
        }
        if(box_overlap_flag[0]==1&&box_overlap_flag[1]==1&&box_overlap_flag[2]==1)
        return 1; 
       }

     }
     }
     }
     
     
 return flag;
}

/* ----------------------------------------------------------------------
   determine overlap list of Noverlap procs the lo/hi box overlaps
   overlap = non-zero area in common between box and proc sub-domain
   box is owned by me and extends in dim
------------------------------------------------------------------------- */

void CommCAC::box_drop_brick(int idim, double *lo, double *hi, int &indexme)
{
  // NOTE: this is not triclinic compatible
  // NOTE: these error messages are internal sanity checks
  //       should not occur, can be removed at some point

  int index=-1,dir;
  if (hi[idim] == sublo[idim]) {
    index = myloc[idim] - 1;
    dir = -1;
  } else if (lo[idim] == subhi[idim]) {
    index = myloc[idim] + 1;
    dir = 1;
  } else if (hi[idim] == boxhi[idim]) {
    index = procgrid[idim] - 1;
    dir = -1;
  } else if (lo[idim] == boxlo[idim]) {
    index = 0;
    dir = 1;
  } else error->one(FLERR,"Comm cac mis-match in box drop brick");

  int other1,other2,proc;
  double lower,upper;
  double *split;

  if (idim == 0) {
    other1 = myloc[1]; other2 = myloc[2];
    split = xsplit;
  } else if (idim == 1) {
    other1 = myloc[0]; other2 = myloc[2];
    split = ysplit;
  } else {
    other1 = myloc[0]; other2 = myloc[1];
    split = zsplit;
  }

  if (index < 0 || index > procgrid[idim])
    error->one(FLERR,"Comm cac invalid index in box drop brick");

  while (1) {
    lower = boxlo[idim] + prd[idim]*split[index];
    if (index < procgrid[idim]-1)
      upper = boxlo[idim] + prd[idim]*split[index+1];
    else upper = boxhi[idim];
    if (lower >= hi[idim] || upper <= lo[idim]) break;

    if (idim == 0) proc = grid2proc[index][other1][other2];
    else if (idim == 1) proc = grid2proc[other1][index][other2];
    else proc = grid2proc[other1][other2][index];

    if (noverlap == maxoverlap) {
      maxoverlap += DELTA_PROCS;
      memory->grow(overlap,maxoverlap,"comm:overlap");
    }

    if (proc == me) indexme = noverlap;
    overlap[noverlap++] = proc;
    index += dir;
    if (index < 0 || index >= procgrid[idim]) break;
  }
}

/* ----------------------------------------------------------------------
   determine overlap list of Noverlap procs the lo/hi box overlaps
   overlap = non-zero area in common between box and proc sub-domain
   box is owned by me and extends in dim
------------------------------------------------------------------------- */

void CommCAC::box_drop_brick_full(int idim, double *lo, double *hi, int &indexme)
{
  // NOTE: this is not triclinic compatible
  // NOTE: these error messages are internal sanity checks
  //       should not occur, can be removed at some point

  int index=-1,dir;
  int pbc_bit;
  int dim_range[6];
  double subbox_size[3];
  double procbox_lo[3];
  double procbox_hi[3];
  int other1,other2,proc;
  double lower,upper;
  double *split;
  double *split_array[3];
  int index_init[3];
  subbox_size[0]=subhi[0]-sublo[0];
  subbox_size[1]=subhi[1]-sublo[1];
  subbox_size[2]=subhi[2]-sublo[2];
  dim_range[0]=0;
  dim_range[1]=0;
  dim_range[2]=0;
  dim_range[3]=0;
  dim_range[4]=0;
  dim_range[5]=0;
  int xi[3],xf[3];
  split_array[0]=xsplit;
  split_array[1]=ysplit;
  split_array[2]=zsplit;
  for(int i=0; i<dimension; i++){
    xi[i]=(unsigned int) (lo[i]-boxlo[i])/subbox_size[i];
    xf[i]=(unsigned int) (hi[i]-boxlo[i])/subbox_size[i];
    
     if(hi[i]==boxhi[i]) xf[i]=procgrid[i]-1;
     if(lo[i]==boxlo[i]) xi[i]=0;
  }
 
    for(int x=xi[0]; x<=xf[0]; x++){
      for(int y=xi[1]; y<=xf[1]; y++){
        for(int z=xi[2]; z<=xf[2]; z++){
         if (noverlap == maxoverlap) {
         maxoverlap += DELTA_PROCS;
         memory->grow(overlap,maxoverlap,"comm_CAC:overlap");
         }
         if (noverlap >= maxoverlap_box) {
         maxoverlap_box += DELTA_PROCS;
         memory->grow(proc2box,maxoverlap_box,6,"comm_CAC:proc2box");
         memory->grow(overlap_pbc,maxoverlap_box,3,"comm_CAC:proc2box");
         }
         for(int boxdim=0; boxdim < domain->dimension; boxdim++){
           int split_index;
           if(boxdim==0) split_index=x;
           if(boxdim==1) split_index=y;
           if(boxdim==2) split_index=z;
           procbox_lo[boxdim] = boxlo[boxdim] + prd[boxdim]*split_array[boxdim][split_index];
           if (split_index < procgrid[boxdim]-1)
           procbox_hi[boxdim] = boxlo[boxdim] + prd[boxdim]*split_array[boxdim][split_index+1];
           else procbox_hi[boxdim] = boxhi[boxdim];
         }
         proc = grid2proc[x][y][z];
         //remove repeats from overlap
        pbc_bit=((current_pbc[0]+1)+3*(current_pbc[1]+1)+9*(current_pbc[2]+1));
        if(!pbc_overlap&&proc==me) continue;
        else if(!pbc_overlap) pbc_bit=13;
         if(!(overlap_repeat[proc]&(1<<pbc_bit))) overlap_repeat[proc]+=1<<pbc_bit;
         else continue;
         
         //if (proc == me) indexme = noverlap;
         proc2box[noverlap][0]=procbox_lo[0];
         proc2box[noverlap][1]=procbox_lo[1];
         proc2box[noverlap][2]=procbox_lo[2];
         proc2box[noverlap][3]=procbox_hi[0];
         proc2box[noverlap][4]=procbox_hi[1];
         proc2box[noverlap][5]=procbox_hi[2];
         if(pbc_overlap){
         overlap_pbc[noverlap][0]=current_pbc[0];
         overlap_pbc[noverlap][1]=current_pbc[1];
         overlap_pbc[noverlap][2]=current_pbc[2];
         }
         else{
         overlap_pbc[noverlap][0]=0;
         overlap_pbc[noverlap][1]=0;
         overlap_pbc[noverlap][2]=0;  
         }
         overlap[noverlap++] = proc;
         
         }
      }

    }

}

/* ----------------------------------------------------------------------
   determine overlap list of Noverlap procs the lo/hi box overlaps
   overlap = non-zero area in common between box and proc sub-domain
   recursive method for traversing an RCB tree of cuts
   no need to split lo/hi box as recurse b/c OK if box extends outside RCB box
------------------------------------------------------------------------- */

void CommCAC::box_drop_tiled(int /*idim*/, double *lo, double *hi, int &indexme)
{
  box_drop_tiled_recurse(lo,hi,0,nprocs-1,indexme);
}

void CommCAC::box_drop_tiled_recurse(double *lo, double *hi,
                                       int proclower, int procupper,
                                       int &indexme)
{
  // end recursion when partition is a single proc
  // add proc to overlap list

  if (proclower == procupper) {
    if (noverlap == maxoverlap) {
      maxoverlap += DELTA_PROCS;
      memory->grow(overlap,maxoverlap,"comm:overlap");
    }
    if (proclower == me) indexme = noverlap;
    overlap[noverlap++] = proclower;
    return;
  }

  // drop box on each side of cut it extends beyond
  // use > and < criteria so does not include a box it only touches
  // procmid = 1st processor in upper half of partition
  //         = location in tree that stores this cut
  // dim = 0,1,2 dimension of cut
  // cut = position of cut

  int procmid = proclower + (procupper - proclower) / 2 + 1;
  int idim = rcbinfo[procmid].dim;
  double cut = boxlo[idim] + prd[idim]*rcbinfo[procmid].cutfrac;

  if (lo[idim] < cut)
    box_drop_tiled_recurse(lo,hi,proclower,procmid-1,indexme);
  if (hi[idim] > cut)
    box_drop_tiled_recurse(lo,hi,procmid,procupper,indexme);
}

/* ----------------------------------------------------------------------
   determine overlap list of Noverlap procs the lo/hi box overlaps
   overlap = non-zero area in common between box and proc sub-domain
   recursive method for traversing an RCB tree of cuts
   no need to split lo/hi box as recurse b/c OK if box extends outside RCB box
------------------------------------------------------------------------- */

void CommCAC::box_drop_tiled_full(int /*idim*/, double *lo, double *hi, int &indexme)
{
  box_drop_tiled_recurse_full(lo,hi,0,nprocs-1,indexme);
}

void CommCAC::box_drop_tiled_recurse_full(double *lo, double *hi,
                                       int proclower, int procupper,
                                       int &indexme)
{
  // end recursion when partition is a single proc
  // add proc to overlap list
  int pbc_bit;
  if (proclower == procupper) {
    if (noverlap == maxoverlap) {
      maxoverlap += DELTA_PROCS;
      memory->grow(overlap,maxoverlap,"comm:overlap");
    }
    if (noverlap >= maxoverlap_box) {
         maxoverlap_box += DELTA_PROCS;
         memory->grow(proc2box,maxoverlap_box,6,"comm_CAC:proc2box");
         memory->grow(overlap_pbc,maxoverlap_box,3,"comm_CAC:proc2box");
    }
        pbc_bit=((current_pbc[0]+1)+3*(current_pbc[1]+1)+9*(current_pbc[2]+1));
        if(!pbc_overlap) pbc_bit=13;
        if(!pbc_overlap&&proclower==me) return;
         if(!(overlap_repeat[proclower]&(1<<pbc_bit))) overlap_repeat[proclower]+=1<<pbc_bit;
         else return;
          
    //if (proclower == me) indexme = noverlap;
    if(pbc_overlap){
         overlap_pbc[noverlap][0]=current_pbc[0];
         overlap_pbc[noverlap][1]=current_pbc[1];
         overlap_pbc[noverlap][2]=current_pbc[2];
         }
         else{
         overlap_pbc[noverlap][0]=0;
         overlap_pbc[noverlap][1]=0;
         overlap_pbc[noverlap][2]=0;  
         }
    overlap[noverlap++] = proclower;
    return;
  }

  // drop box on each side of cut it extends beyond
  // use > and < criteria so does not include a box it only touches
  // procmid = 1st processor in upper half of partition
  //         = location in tree that stores this cut
  // dim = 0,1,2 dimension of cut
  // cut = position of cut

  int procmid = proclower + (procupper - proclower) / 2 + 1;
  int idim = rcbinfo[procmid].dim;
  double cut = boxlo[idim] + prd[idim]*rcbinfo[procmid].cutfrac;

  if (lo[idim] < cut)
    box_drop_tiled_recurse_full(lo,hi,proclower,procmid-1,indexme);
  if (hi[idim] > cut)
    box_drop_tiled_recurse_full(lo,hi,procmid,procupper,indexme);
}

/* ----------------------------------------------------------------------
   return other box owned by proc as lo/hi corner pts
------------------------------------------------------------------------- */

void CommCAC::box_other_brick(int idim, int idir,
                                int proc, double *lo, double *hi)
{
  lo[0] = sublo[0]; lo[1] = sublo[1]; lo[2] = sublo[2];
  hi[0] = subhi[0]; hi[1] = subhi[1]; hi[2] = subhi[2];

  int other1,other2,oproc;
  double *split;

  if (idim == 0) {
    other1 = myloc[1]; other2 = myloc[2];
    split = xsplit;
  } else if (idim == 1) {
    other1 = myloc[0]; other2 = myloc[2];
    split = ysplit;
  } else {
    other1 = myloc[0]; other2 = myloc[1];
    split = zsplit;
  }

  int dir = -1;
  if (idir) dir = 1;
  int index = myloc[idim];
  int n = procgrid[idim];

  for (int i = 0; i < n; i++) {
    index += dir;
    if (index < 0) index = n-1;
    else if (index >= n) index = 0;

    if (idim == 0) oproc = grid2proc[index][other1][other2];
    else if (idim == 1) oproc = grid2proc[other1][index][other2];
    else oproc = grid2proc[other1][other2][index];

    if (proc == oproc) {
      lo[idim] = boxlo[idim] + prd[idim]*split[index];
      if (split[index+1] < 1.0)
        hi[idim] = boxlo[idim] + prd[idim]*split[index+1];
      else hi[idim] = boxhi[idim];
      return;
    }
  }
}

/* ----------------------------------------------------------------------
   return other box owned by proc as lo/hi corner pts
------------------------------------------------------------------------- */

void CommCAC::box_other_brick_full(int idim, int idir,
                                int proc, double *lo, double *hi)
{
  lo[0]=proc2box[overlap_counter][0];
  lo[1]=proc2box[overlap_counter][1];
  lo[2]=proc2box[overlap_counter][2];
  hi[0]=proc2box[overlap_counter][3];
  hi[1]=proc2box[overlap_counter][4];
  hi[2]=proc2box[overlap_counter][5];

}

/* ----------------------------------------------------------------------
   return other box owned by proc as lo/hi corner pts
------------------------------------------------------------------------- */

void CommCAC::box_other_tiled(int /*idim*/, int /*idir*/,
                                int proc, double *lo, double *hi)
{
  double (*split)[2] = rcbinfo[proc].mysplit;

  lo[0] = boxlo[0] + prd[0]*split[0][0];
  if (split[0][1] < 1.0) hi[0] = boxlo[0] + prd[0]*split[0][1];
  else hi[0] = boxhi[0];

  lo[1] = boxlo[1] + prd[1]*split[1][0];
  if (split[1][1] < 1.0) hi[1] = boxlo[1] + prd[1]*split[1][1];
  else hi[1] = boxhi[1];

  lo[2] = boxlo[2] + prd[2]*split[2][0];
  if (split[2][1] < 1.0) hi[2] = boxlo[2] + prd[2]*split[2][1];
  else hi[2] = boxhi[2];
}

/* ----------------------------------------------------------------------
   return 1 if proc's box touches me, else 0
   procneigh stores 6 procs that touch me
------------------------------------------------------------------------- */

int CommCAC::box_touch_brick(int proc, int idim, int idir)
{
  if (procneigh[idim][idir] == proc) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return 1 if proc's box touches me, else 0
------------------------------------------------------------------------- */

int CommCAC::box_touch_tiled(int proc, int idim, int idir)
{
  // sending to left
  // only touches if proc hi = my lo, or if proc hi = boxhi and my lo = boxlo

  if (idir == 0) {
    if (rcbinfo[proc].mysplit[idim][1] == rcbinfo[me].mysplit[idim][0])
      return 1;
    else if (rcbinfo[proc].mysplit[idim][1] == 1.0 &&
             rcbinfo[me].mysplit[idim][0] == 0.0)
      return 1;

  // sending to right
  // only touches if proc lo = my hi, or if proc lo = boxlo and my hi = boxhi

  } else {
    if (rcbinfo[proc].mysplit[idim][0] == rcbinfo[me].mysplit[idim][1])
      return 1;
    else if (rcbinfo[proc].mysplit[idim][0] == 0.0 &&
             rcbinfo[me].mysplit[idim][1] == 1.0)
      return 1;
  }

  return 0;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int CommCAC::point_drop_brick(int idim, double *x)
{
  if (closer_subbox_edge(idim,x)) return procneigh[idim][1];
  return procneigh[idim][0];
}

/* ----------------------------------------------------------------------
   determine which proc owns point x via recursion thru RCB tree
------------------------------------------------------------------------- */

int CommCAC::point_drop_tiled(int idim, double *x)
{
  double xnew[3];
  int dim1_other, dim2_other;
  xnew[0] = x[0]; xnew[1] = x[1]; xnew[2] = x[2];

  if (idim==0) {
    dim1_other=1;
    dim2_other=2;
  }
  if (idim==1) {
    dim1_other=0;
    dim2_other=2;
  }
  if (idim==2) {
    dim1_other=0;
    dim2_other=1;
  }

    if (xnew[dim1_other] < sublo[dim1_other] || xnew[dim1_other] > subhi[dim1_other]) {
      if (closer_subbox_edge(dim1_other,x)) xnew[dim1_other] = subhi[dim1_other];
      else xnew[dim1_other] = sublo[dim1_other];
    }
    if (xnew[dim2_other] < sublo[dim2_other] || xnew[dim2_other] > subhi[dim2_other]) {
      if (closer_subbox_edge(dim2_other,x)) xnew[dim2_other] = subhi[dim2_other];
      else xnew[dim2_other] = sublo[dim2_other];
    }
  


  int proc = point_drop_tiled_recurse(xnew,0,nprocs-1);//set to 6 here
  if (proc == me) return me;

    /*since point tiled recurse includes the closure of the split interval (i.e an == case)
    //perturb the particle in that dimension lower to bring it out of the closure that normally returns an adjacent proc
    to that which is intended */
    int done = 1;
    if (rcbinfo[proc].mysplit[dim1_other][0] == rcbinfo[me].mysplit[dim1_other][1]) {
      xnew[dim1_other] -= EPSILON * (subhi[dim1_other]-sublo[dim1_other]);
      done = 0;
    }
    if (rcbinfo[proc].mysplit[dim2_other][0] == rcbinfo[me].mysplit[dim2_other][1]) {
      xnew[dim2_other] -= EPSILON * (subhi[dim2_other]-sublo[dim2_other]);
      done = 0;
    }
    
    if (!done) {
      proc = point_drop_tiled_recurse(xnew,0,nprocs-1);
      done = 1;
      //sanity catch in case epsilon wasn't enough and floating point precision still catches an adjacent proc
      if (rcbinfo[proc].mysplit[dim1_other][0] == rcbinfo[me].mysplit[dim1_other][1]) {
        xnew[dim1_other] -= EPSILON * (subhi[dim1_other]-sublo[dim1_other]);
        done = 0;
      }
      if (rcbinfo[proc].mysplit[dim2_other][0] == rcbinfo[me].mysplit[dim2_other][1]) {
        xnew[dim2_other] -= EPSILON * (subhi[dim2_other]-sublo[dim2_other]);
        done = 0;
      }
      if (!done) proc = point_drop_tiled_recurse(xnew,0,nprocs-1);
    }
  

  return proc;
}

/* ----------------------------------------------------------------------
   recursive point drop thru RCB tree
------------------------------------------------------------------------- */

int CommCAC::point_drop_tiled_recurse(double *x,
                                        int proclower, int procupper)
{
  // end recursion when partition is a single proc
  // return proc

  if (proclower == procupper) return proclower;

  // drop point on side of cut it is on
  // use < criterion so point is not on high edge of proc sub-domain
  // procmid = 1st processor in upper half of partition
  //         = location in tree that stores this cut
  // dim = 0,1,2 dimension of cut
  // cut = position of cut

  int procmid = proclower + (procupper - proclower) / 2 + 1;
  int idim = rcbinfo[procmid].dim;
  double cut = boxlo[idim] + prd[idim]*rcbinfo[procmid].cutfrac;

  if (x[idim] < cut) return point_drop_tiled_recurse(x,proclower,procmid-1);
  else return point_drop_tiled_recurse(x,procmid,procupper);
}

/* ----------------------------------------------------------------------
   assume x[idim] is outside subbox bounds in same dim
------------------------------------------------------------------------- */

int CommCAC::closer_subbox_edge(int idim, double *x)
{
  double deltalo,deltahi;

  if (sublo[idim] == boxlo[idim])
    deltalo = fabs(x[idim]-prd[idim] - sublo[idim]);
  else deltalo = fabs(x[idim] - sublo[idim]);

  if (subhi[idim] == boxhi[idim])
    deltahi = fabs(x[idim]+prd[idim] - subhi[idim]);
  else deltahi = fabs(x[idim] - subhi[idim]);

  if (deltalo < deltahi) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   if RCB decomp exists and just changed, gather needed global RCB info
------------------------------------------------------------------------- */

void CommCAC::coord2proc_setup()
{
  if (!rcbnew) return;

  if (!rcbinfo)
    rcbinfo = (RCBinfo *)
      memory->smalloc(nprocs*sizeof(RCBinfo),"comm:rcbinfo");
  rcbnew = 0;
  RCBinfo rcbone;
  memcpy(&rcbone.mysplit[0][0],&mysplit[0][0],6*sizeof(double));
  rcbone.cutfrac = rcbcutfrac;
  rcbone.dim = rcbcutdim;
  MPI_Allgather(&rcbone,sizeof(RCBinfo),MPI_CHAR,
                rcbinfo,sizeof(RCBinfo),MPI_CHAR,world);
}

/* ----------------------------------------------------------------------
   determine which proc owns atom with coord x[3] based on current decomp
   x will be in box (orthogonal) or lamda coords (triclinic)
   if layout = UNIFORM or NONUNIFORM, invoke parent method
   if layout = TILED, use point_drop_recurse()
   return owning proc ID, ignore igx,igy,igz
------------------------------------------------------------------------- */

int CommCAC::coord2proc(double *x, int &igx, int &igy, int &igz)
{
  if (layout != Comm::LAYOUT_TILED) return Comm::coord2proc(x,igx,igy,igz);
  return point_drop_tiled_recurse(x,0,nprocs-1);
}

/* ---------------------------------------------------------------------- */

int CommCAC::pack_eboxes(int n, int *list, double *buf,
                               int pbc_flag, int *pbc,int iswap)
{
  int i,j,m;
  double dx,dy,dz;
  double *current_ebox;
  m = 0;
  
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      if(j<neboxes) current_ebox=eboxes[j];
      else current_ebox=foreign_eboxes[j-neboxes];
      //buf[m++] = x[j][0];
      buf[m++]=me;
      buf[m++]=iswap;
      buf[m++] = 0;
			buf[m++] = 0;
			buf[m++] = 0;
	 	  buf[m++] = current_ebox[0];
			buf[m++] = current_ebox[1];
			buf[m++] = current_ebox[2];
			buf[m++] = current_ebox[3];
			buf[m++] = current_ebox[4];
			buf[m++] = current_ebox[5];
		  
	  
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      if(j<neboxes) current_ebox=eboxes[j];
      else current_ebox=foreign_eboxes[j-neboxes];
      buf[m++]=me;
      buf[m++]=iswap;
      buf[m++] = pbc[0];
			buf[m++] = pbc[1];
			buf[m++] = pbc[2];
			buf[m++] = current_ebox[0]+dx;
			buf[m++] = current_ebox[1]+dy;
			buf[m++] = current_ebox[2]+dz;
			buf[m++] = current_ebox[3]+dx;
			buf[m++] = current_ebox[4]+dy;
			buf[m++] = current_ebox[5]+dz;
			
	  
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void CommCAC::unpack_eboxes(int n, int first, double *buf)
{
  int i,m,last;
  int *nodes_count_list = atom->nodes_per_element_list;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    
    if(neboxes>=maxebox) {
      maxebox+=BUFEXTRA;
      eboxes=memory->grow(atom->eboxes,maxebox,6,"commCAC: eboxes");
      }
    if(nforeign_eboxes>=maxforeign_ebox){
    maxforeign_ebox+=BUFEXTRA;
    foreign_eprocs=memory->grow(foreign_eprocs,maxforeign_ebox,"commCAC: foreign_eprocs");
    foreign_swaps=memory->grow(foreign_swaps,maxforeign_ebox,"commCAC: foreign_swaps");
    foreign_eboxes=memory->grow(foreign_eboxes,maxforeign_ebox,6, "commCAC: foreign_eboxes");
    foreign_image=memory->grow(foreign_image,maxforeign_ebox,3, "commCAC: foreign_image");
    }
    foreign_eprocs[i]=buf[m++];
    foreign_swaps[i]=buf[m++];
    foreign_image[i][0] = buf[m++];
    foreign_image[i][1] = buf[m++];
    foreign_image[i][2] = buf[m++];
    foreign_eboxes[i][0] = buf[m++];
    foreign_eboxes[i][1] = buf[m++];
    foreign_eboxes[i][2] = buf[m++];
    foreign_eboxes[i][3] = buf[m++];
    foreign_eboxes[i][4] = buf[m++];
    foreign_eboxes[i][5] = buf[m++];
    nforeign_eboxes++;
	
  }

}
/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR and bufextra
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void CommCAC::grow_send(int n, int flag)
{
  maxsend = static_cast<int> (BUFFACTOR * n);
  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;
  if (flag)
    memory->grow(buf_send,maxsend+bufextra,"comm:buf_send");
  else {
    memory->destroy(buf_send);
    memory->create(buf_send,maxsend+bufextra,"comm:buf_send");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommCAC::grow_recv(int n)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  memory->destroy(buf_recv);
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
}


/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommCAC::grow_list(int iswap, int iwhich, int n)
{
  maxsendlist[iswap][iwhich] = static_cast<int> (BUFFACTOR * n);
  memory->grow(sendlist[iswap][iwhich],maxsendlist[iswap][iwhich],
               "comm:sendlist[i][j]");
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap element overlap sendlist as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommCAC::overlap_grow_list(int iswap, int iwhich, int n)
{
  overlap_maxsendlist[iswap][iwhich] = static_cast<int> (BUFFACTOR * n);
  memory->grow(overlap_sendlist[iswap][iwhich],overlap_maxsendlist[iswap][iwhich],
               "comm:overlap_sendlist[i][j]");
}


/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommCAC::grow_sent_list(int iswap, int iwhich, int n)
{
  int oldsize=maxsent[iswap][iwhich];
  maxsent[iswap][iwhich] = static_cast<int> (BUFFACTOR * n);
  memory->grow(sent_flag[iswap][iwhich],maxsent[iswap][iwhich],
               "comm:sent_flag[i][j]");
  //initialize new portion of array to zero
  for(int init=oldsize;init<maxsent[iswap][iwhich];init++)
  sent_flag[iswap][iwhich][init]=-1;
}

/* ----------------------------------------------------------------------
   allocation of swap info
------------------------------------------------------------------------- */

void CommCAC::allocate_swap(int n)
{
  nsendproc = new int[n];
  nrecvproc = new int[n];
  sendother = new int[n];
  recvother = new int[n];
  sendself = new int[n];
  nprocmax = new int[n];
  nrecv_procmax = new int[n];
  sendproc = new int*[n];
  sendbox_flag = new int*[n];
  repeatsend_flag = new int*[n];
  recvproc = new int*[n];
  sendnum = new int*[n];
  overlap_sendnum = new int*[n];
  sendsize = new int*[n];
  recvsize = new int*[n];
  sendoffset = new int*[n];
  recvoffset = new int*[n];
  overlap_sendsize = new int*[n];
  overlap_recvsize = new int*[n];
  overlap_sendoffset = new int*[n];
  overlap_recvoffset = new int*[n];
  pair_sendsize = new int*[n];
  pair_recvsize = new int*[n];
  pair_recvoffset = new int*[n];
  pair_sendoffset = new int*[n];
  recvnum = new int*[n];
  overlap_recvnum = new int*[n];
  size_forward_recv = new int*[n];
  firstrecv = new int*[n];
  overlap_firstrecv = new int*[n];
  size_reverse_send = new int*[n];
  size_reverse_recv = new int*[n];
  forward_recv_offset = new int*[n];
  reverse_recv_offset = new int*[n];

  pbc_flag = new int*[n];
  pbc = new int**[n];
  sendbox = new double**[n];
  overlap_sendbox = new double**[n];
  sendbox_multi = new double***[n];
  maxsendlist = new int*[n];
  overlap_maxsendlist = new int*[n];
  maxsent = new int*[n];
  
  sendlist = new int**[n];
  overlap_sendlist = new int**[n];
  aug_oboxes = new double**[n];

 

  for (int i = 0; i < n; i++) {
    sendproc[i] = recvproc[i] = NULL;
    sendbox_flag[i] = NULL;
    repeatsend_flag[i] = NULL;
    sendnum[i] = recvnum[i] = overlap_sendnum[i] = overlap_recvnum[i] = NULL;
    sendsize[i] = NULL;
    recvsize[i] = NULL;
    sendoffset[i] = NULL;
    recvoffset[i] = NULL;
    overlap_sendsize[i] = NULL;
    overlap_recvsize[i] = NULL;
    overlap_sendoffset[i] = NULL;
    overlap_recvoffset[i] = NULL;
    pair_sendsize[i] = NULL;
    pair_recvsize[i] = NULL;
    pair_sendoffset[i] = NULL;
    pair_recvoffset[i] = NULL;
    size_forward_recv[i] = firstrecv[i] = overlap_firstrecv[i] = NULL;
    size_reverse_send[i] = size_reverse_recv[i] = NULL;
    forward_recv_offset[i] = reverse_recv_offset[i] = NULL;

    pbc_flag[i] = NULL;
    pbc[i] = NULL;
    sendbox[i] = NULL;
    overlap_sendbox[i] = NULL;
    sendbox_multi[i] = NULL;
  
    maxsendlist[i] = NULL;
    overlap_maxsendlist[i] = NULL;
    maxsent[i] = NULL;
    sendlist[i] = NULL;
    aug_oboxes[i] = NULL;
    overlap_sendlist[i] = NULL;

  }

  maxreqstat = 0;
  requests = NULL;

  for (int i = 0; i < n; i++) {
    nprocmax[i] = DELTA_PROCS;
    nrecv_procmax[i] = DELTA_PROCS;
    grow_swap_send(i,DELTA_PROCS,0);
    grow_swap_recv(i,DELTA_PROCS,0);
  }

  nexchproc = new int[n/2];
  nexchprocmax = new int[n/2];
  exchproc = new int*[n/2];
  exchnum = new int*[n/2];

  for (int i = 0; i < n/2; i++) {
    nexchprocmax[i] = DELTA_PROCS;
    exchproc[i] = new int[DELTA_PROCS];
    exchnum[i] = new int[DELTA_PROCS];
  }
}

/* ----------------------------------------------------------------------
   grow info for swap I, to allow for N procs to communicate with
   ditto for complementary recv for swap I+1 or I-1, as invoked by caller
------------------------------------------------------------------------- */

void CommCAC::grow_swap_send(int i, int n, int nold)
{
  delete [] sendproc[i];
  sendproc[i] = new int[n];
  delete [] sendbox_flag[i];
  sendbox_flag[i] = new int[n];
  delete [] repeatsend_flag[i];
  repeatsend_flag[i] = new int[n];
  delete [] sendnum[i];
  sendnum[i] = new int[n];
  delete [] overlap_sendnum[i];
  overlap_sendnum[i] = new int[n];
  delete [] sendsize[i];
  sendsize[i] = new int[n];
  delete [] sendoffset[i];
  sendoffset[i] = new int[n];
  delete [] overlap_sendsize[i];
  overlap_sendsize[i] = new int[n];
  delete [] overlap_sendoffset[i];
  overlap_sendoffset[i] = new int[n];
  delete [] pair_sendsize[i];
  pair_sendsize[i] = new int[n];
  delete [] pair_sendoffset[i];
  pair_sendoffset[i] = new int[n];

  delete [] size_reverse_recv[i];
  size_reverse_recv[i] = new int[n];
  delete [] reverse_recv_offset[i];
  reverse_recv_offset[i] = new int[n];

  delete [] pbc_flag[i];
  pbc_flag[i] = new int[n];
  memory->destroy(pbc[i]);
  memory->create(pbc[i],n,6,"comm:pbc_flag");
  memory->destroy(sendbox[i]);
  memory->create(sendbox[i],n,6,"comm:sendbox");
  memory->destroy(overlap_sendbox[i]);
  memory->create(overlap_sendbox[i],n,6,"comm:sendbox");
 
  memory->destroy(sendbox_multi[i]);
  memory->create(sendbox_multi[i],n,atom->ntypes+1,6,"comm:sendbox_multi");
  

  delete [] maxsendlist[i];
  maxsendlist[i] = new int[n];
  delete [] overlap_maxsendlist[i];
  overlap_maxsendlist[i] = new int[n];
  delete [] maxsent[i];
  maxsent[i] = new int[n];



  for (int j = 0; j < nold; j++){ 
    memory->destroy(sendlist[i][j]);
    memory->destroy(overlap_sendlist[i][j]);
    memory->destroy(aug_oboxes[i][j]);
    
    }
  delete [] sendlist[i];
  sendlist[i] = new int*[n];
  delete [] overlap_sendlist[i];
  overlap_sendlist[i] = new int*[n];
  delete [] aug_oboxes[i];
  aug_oboxes[i] = new double*[n];
  
  for (int j = 0; j < n; j++) {
    maxsendlist[i][j] = BUFMIN;
    overlap_maxsendlist[i][j] = BUFMIN;
    maxsent[i][j] = BUFMIN;
    
    memory->create(sendlist[i][j],BUFMIN,"comm_CAC:sendlist[i][j]");
    memory->create(overlap_sendlist[i][j],BUFMIN,"comm_CAC:sendlist[i][j]");
    memory->create(aug_oboxes[i][j],6,"comm_CAC:aug_oboxes[i][j]");
    
  }
}

void CommCAC::grow_swap_recv(int i, int n, int nold)
{
  delete [] recvproc[i];
  recvproc[i] = new int[n];
  delete [] recvnum[i];
  recvnum[i] = new int[n];
  delete [] overlap_recvnum[i];
  overlap_recvnum[i] = new int[n];
  delete [] size_forward_recv[i];
  size_forward_recv[i] = new int[n];
  delete [] firstrecv[i];
  firstrecv[i] = new int[n];
  delete [] overlap_firstrecv[i];
  overlap_firstrecv[i] = new int[n];
  delete [] forward_recv_offset[i];
  forward_recv_offset[i] = new int[n];
  delete [] recvoffset[i];
  recvoffset[i] = new int[n];
  delete [] recvsize[i];
  recvsize[i] = new int[n];
  delete [] overlap_recvsize[i];
  overlap_recvsize[i] = new int[n];
  delete [] overlap_recvoffset[i];
  overlap_recvoffset[i] = new int[n];
  delete [] pair_recvsize[i];
  pair_recvsize[i] = new int[n];
  delete [] pair_recvoffset[i];
  pair_recvoffset[i] = new int[n];
  delete [] size_reverse_send[i];
  size_reverse_send[i] = new int[n];
  
}

/* ----------------------------------------------------------------------
   deallocate swap info
------------------------------------------------------------------------- */

void CommCAC::deallocate_swap(int n)
{
  delete [] nsendproc;
  delete [] nrecvproc;
  delete [] sendother;
  delete [] recvother;
  delete [] sendself;

  for (int i = 0; i < n; i++) {
    delete [] sendproc[i];
    delete [] sendbox_flag[i];
    delete [] repeatsend_flag[i];
    delete [] recvproc[i];
    delete [] sendnum[i];
    delete [] overlap_sendnum[i];
    delete [] sendsize[i];
    delete [] recvsize[i];
    delete [] sendoffset[i];
    delete [] recvoffset[i];
    delete [] overlap_sendsize[i];
    delete [] overlap_recvsize[i];
    delete [] overlap_sendoffset[i];
    delete [] overlap_recvoffset[i];
    delete [] pair_sendsize[i];
    delete [] pair_recvsize[i];
    delete [] pair_sendoffset[i];
    delete [] pair_recvoffset[i];
    delete [] recvnum[i];
    delete [] overlap_recvnum[i];
    delete [] size_forward_recv[i];
    delete [] firstrecv[i];
    delete [] overlap_firstrecv[i];
    delete [] size_reverse_send[i];
    delete [] size_reverse_recv[i];
    delete [] forward_recv_offset[i];
    delete [] reverse_recv_offset[i];

    delete [] pbc_flag[i];
    memory->destroy(pbc[i]);
    memory->destroy(sendbox[i]);
    memory->destroy(overlap_sendbox[i]);
    memory->destroy(sendbox_multi[i]);
 
    delete [] maxsendlist[i];
    delete [] overlap_maxsendlist[i];
    delete [] maxsent[i];
  

    for (int j = 0; j < nprocmax[i]; j++){ 
      memory->destroy(sendlist[i][j]);
      memory->destroy(overlap_sendlist[i][j]);
      memory->destroy(aug_oboxes[i][j]);
  
      }
    delete [] sendlist[i];
    delete [] overlap_sendlist[i];
    delete [] aug_oboxes[i];
  
  }

  delete [] sendproc;
  delete [] sendbox_flag;
  delete [] repeatsend_flag;
  delete [] recvproc;
  delete [] sendnum;
  delete [] overlap_sendnum;
  delete [] sendsize;
  delete [] recvsize;
  delete [] sendoffset;
  delete [] recvoffset;
  delete [] overlap_sendsize;
  delete [] overlap_recvsize;
  delete [] overlap_sendoffset;
  delete [] overlap_recvoffset;
  delete [] pair_sendsize;
  delete [] pair_recvsize;
  delete [] pair_sendoffset;
  delete [] pair_recvoffset;
  delete [] recvnum;
  delete [] overlap_recvnum;
  delete [] size_forward_recv;
  delete [] firstrecv;
  delete [] overlap_firstrecv;
  delete [] size_reverse_send;
  delete [] size_reverse_recv;
  delete [] forward_recv_offset;
  delete [] reverse_recv_offset;

  delete [] pbc_flag;
  delete [] pbc;
  delete [] sendbox;
  delete [] overlap_sendbox;
  delete [] sendbox_multi;
  
  delete [] maxsendlist;
  delete [] overlap_maxsendlist;
  delete [] maxsent;
  delete [] sendlist;
  delete [] overlap_sendlist;
  delete [] aug_oboxes;
  


  delete [] requests;

  delete [] nprocmax;
  delete [] nrecv_procmax;
  delete [] nexchproc;
  delete [] nexchprocmax;

  for (int i = 0; i < n/2; i++) {
    delete [] exchproc[i];
    delete [] exchnum[i];
  }

  delete [] exchproc;
  delete [] exchnum;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint CommCAC::memory_usage()
{
  bigint bytes = 0;
  return bytes;
}
