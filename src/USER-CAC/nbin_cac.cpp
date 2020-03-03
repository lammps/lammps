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

#include "nbin_cac.h"
#include "neighbor.h"
#include "atom.h"
#include "group.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

enum{NSQ,BIN,MULTI};       // also in Neighbor

#define SMALL 1.0e-6
#define CUT2BIN_RATIO 100
#define MAXBINCONTENT 100  //used to bound local memory
#define EXPAND 100
#define MAXBINOVERLAP 10

/* ---------------------------------------------------------------------- */

NBinCAC::NBinCAC(LAMMPS *lmp) : NBin(lmp) {
  bin_ncontent=NULL;
  bin_content=NULL;
  atom2bin=NULL;
  bin_expansion_counts=NULL;
  nbin_element_overlap=NULL;
  bin_element_overlap=NULL;
  max_nbin_overlap=NULL;
  first_alloc=0;
  max_bin_expansion_count=0;
  nmax=0;
  max_nall=0;
}

/* ---------------------------------------------------------------------- */

NBinCAC::~NBinCAC() {
  memory->destroy(bin_ncontent);
  memory->destroy(bin_expansion_counts);

  for(int i=0 ; i<maxbin; i++)
  memory->destroy(bin_content[i]);
  memory->sfree(bin_content);
  memory->destroy(nbin_element_overlap);
  for(int i=0 ; i<max_nall; i++)
  memory->destroy(bin_element_overlap[i]);
  memory->sfree(bin_element_overlap);
}

/* ----------------------------------------------------------------------
   setup for bin_atoms()
------------------------------------------------------------------------- */
void NBinCAC::bin_atoms_setup(int nall)
{

}

/* ----------------------------------------------------------------------
   setup required for CAc binning strategy
------------------------------------------------------------------------- */

void NBinCAC::CAC_bin_atoms_setup(int nall)
{

  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **element_scale = atom->element_scale;
  double ****nodal_positions = atom->nodal_positions;

  // binhead = per-bin vector, mbins in length
  // add 1 bin for USER-INTEL package

  if (mbins > maxbin) {

    if(!first_alloc){
    first_alloc=1;
    memory->grow(bin_expansion_counts,mbins,"bin_CAC:bin_ncontent");
      for (int init = 0; init < mbins; init++){
      bin_expansion_counts[init]=max_bin_expansion_count;
    }
    bin_content= (int **) memory->smalloc(mbins*sizeof(int *), "bin_CAC:bin_content");
    for (int init = 0; init < mbins; init++)
    memory->create(bin_content[init],MAXBINCONTENT,"bin_CAC:bin_content");
    }
    else{
    //compute maximum expansion count used in previous allocations
    for (int init = 0; init < maxbin; init++){
      if(bin_expansion_counts[init]>max_bin_expansion_count) max_bin_expansion_count=bin_expansion_counts[init];
    }
    bin_content= (int **) memory->srealloc(bin_content,mbins*sizeof(int *), "bin_CAC:bin_content");
    for (int init = 0; init < mbins; init++){
    if(init>=maxbin) bin_content[init]=NULL;
    memory->grow(bin_content[init],MAXBINCONTENT+max_bin_expansion_count*EXPAND,"neigh:bin_content");
    }
    memory->grow(bin_expansion_counts,mbins,"neigh:bin_ncontent");
    for (int init = 0; init < mbins; init++)
      bin_expansion_counts[init]=max_bin_expansion_count;

    }
    maxbin = mbins;
    memory->grow(bin_ncontent,maxbin,"neigh:bin_ncontent");
  }


  // bins and atom2bin = per-atom vectors
  // for both local and ghost atoms

  if (nall > maxatom) {
    maxatom = nall;
    memory->destroy(bins);
    memory->create(bins,maxatom,"neigh:bins");
    memory->destroy(atom2bin);
    memory->create(atom2bin,maxatom,"neigh:atom2bin");
  }
}

/* ----------------------------------------------------------------------
   setup neighbor binning geometry
   bin numbering in each dimension is global:
     0 = 0.0 to binsize, 1 = binsize to 2*binsize, etc
     nbin-1,nbin,etc = bbox-binsize to bbox, bbox to bbox+binsize, etc
     -1,-2,etc = -binsize to 0.0, -2*binsize to -binsize, etc
   code will work for any binsize
     since next(xyz) and stencil extend as far as necessary
     binsize = 1/2 of cutoff is roughly optimal
   for orthogonal boxes:
     a dim must be filled exactly by integer # of bins
     in periodic, procs on both sides of PBC must see same bin boundary
     in non-periodic, coord2bin() still assumes this by use of nbin xyz
   for triclinic boxes:
     tilted simulation box cannot contain integer # of bins
     stencil & neigh list built differently to account for this
   mbinlo = lowest global bin any of my ghost atoms could fall into
   mbinhi = highest global bin any of my ghost atoms could fall into
   mbin = number of bins I need in a dimension
------------------------------------------------------------------------- */
void NBinCAC::setup_bins(int style)
{
comm->bin_pointer=this;
}

void NBinCAC::CAC_setup_bins(int style)
{
  // bbox = size of bbox of entire domain
  // bsubbox lo/hi = bounding box of my subdomain extended by comm->cutghost
  // for triclinic:
  //   bbox bounds all 8 corners of tilted box
  //   subdomain is in lamda coords
  //   include dimension-dependent extension via comm->cutghost
  //   domain->bbox() converts lamda extent to box coords and computes bbox

  double bbox[3];
  double *cutghost = comm->cutghost;
  double lamda_temp[3];
  double nodal_temp[3];
  double ***nodal_positions;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  double **x = atom->x;
  double ebounding_boxlo[3];
  double ebounding_boxhi[3];

  int *nodes_per_element_list = atom->nodes_per_element_list;
  double max_search_range=atom->max_search_range;
  double cut_max=neighbor->cutneighmax;
  double **eboxes=atom->eboxes;
  int *ebox_ref=atom->ebox_ref;
  //expand sub-box sizes to contain all element bounding boxes
  //in general initial box wont since only x, element centroid, must lie
  //inside box for lammps comms

  if (triclinic == 0) {
    bsubboxlo[0] = domain->sublo[0];
    bsubboxlo[1] = domain->sublo[1];
    bsubboxlo[2] = domain->sublo[2];
    bsubboxhi[0] = domain->subhi[0];
    bsubboxhi[1] = domain->subhi[1];
    bsubboxhi[2] = domain->subhi[2];
    bsubboxlo[0] -= cut_max;
    bsubboxlo[1] -= cut_max;
    bsubboxlo[2] -= cut_max;
    bsubboxhi[0] += cut_max;
    bsubboxhi[1] += cut_max;
    bsubboxhi[2] += cut_max;
  //loop through elements to compute bounding boxes and test
  //whether they should stretch the local bounding box
  for(int element_index=0; element_index < atom->nlocal; element_index++){
  if(element_type[element_index]){
  nodal_positions = atom->nodal_positions[element_index];

  double *current_ebox;

  current_ebox = eboxes[ebox_ref[element_index]];
  //define this elements bounding box
  ebounding_boxlo[0] = current_ebox[0];
  ebounding_boxlo[1] = current_ebox[1];
  ebounding_boxlo[2] = current_ebox[2];
  ebounding_boxhi[0] = current_ebox[3];
  ebounding_boxhi[1] = current_ebox[4];
  ebounding_boxhi[2] = current_ebox[5];

  //test if this bounding box exceeds local sub box
  for(int dim=0; dim < dimension; dim++){
  if(ebounding_boxhi[dim]>bsubboxhi[dim])
  bsubboxhi[dim]=ebounding_boxhi[dim];
  if(ebounding_boxlo[dim]<bsubboxlo[dim])
  bsubboxlo[dim]=ebounding_boxlo[dim];
  }
  }
  else{
  for(int dim=0; dim < dimension; dim++){
  if(x[element_index][dim]>bsubboxhi[dim])
  bsubboxhi[dim]=x[element_index][dim];
  if(x[element_index][dim]<bsubboxlo[dim])
  bsubboxlo[dim]=x[element_index][dim];
  }
  }
  }
  }
  else{
  double lo[3],hi[3];
  lo[0] = domain->sublo_lamda[0];
  lo[1] = domain->sublo_lamda[1];
  lo[2] = domain->sublo_lamda[2];
  hi[0] = domain->subhi_lamda[0];
  hi[1] = domain->subhi_lamda[1];
  hi[2] = domain->subhi_lamda[2];
  //loop through elements to compute bounding boxes and test
  //whether they should stretch the local bounding box
  for(int element_index=0; element_index < atom->nlocal+atom->nghost; element_index++){
  if(element_type[element_index]){
    nodal_positions = atom->nodal_positions[element_index];
  //int current_poly_count = poly_count[element_index];
    int current_poly_count = poly_count[element_index];
    int nodes_per_element = nodes_per_element_list[element_type[element_index]];

    //initialize bounding box values
    ebounding_boxlo[0] = nodal_positions[0][0][0];
    ebounding_boxlo[1] = nodal_positions[0][0][1];
    ebounding_boxlo[2] = nodal_positions[0][0][2];
    ebounding_boxhi[0] = nodal_positions[0][0][0];
    ebounding_boxhi[1] = nodal_positions[0][0][1];
    ebounding_boxhi[2] = nodal_positions[0][0][2];
     //define the bounding box for the element being considered as a neighbor

    for (int poly_counter = 0; poly_counter < current_poly_count; poly_counter++) {
    for (int kkk = 0; kkk < nodes_per_element; kkk++) {

    nodal_temp[0]=nodal_positions[poly_counter][kkk][0];
    nodal_temp[1]=nodal_positions[poly_counter][kkk][1];
    nodal_temp[2]=nodal_positions[poly_counter][kkk][2];
    domain->x2lamda(nodal_temp, lamda_temp);
      //test if this node lies outside local box and stretch box
    for(int dim=0; dim < dimension; dim++){
      if(lamda_temp[dim]>hi[dim])
      hi[dim]=lamda_temp[dim];
      if(lamda_temp[dim]<lo[dim])
      lo[dim]=lamda_temp[dim];
    }
    }
  }

  }
  else if(element_index>=atom->nlocal){
    for(int dim=0; dim < dimension; dim++){
     if(x[element_index][dim]>hi[dim])
     hi[dim]=x[element_index][dim];
     if(x[element_index][dim]<lo[dim])
     lo[dim]=x[element_index][dim];
    }
  }
  }

    lo[0] -= cutghost[0];
    lo[1] -= cutghost[1];
    lo[2] -= cutghost[2];
    hi[0] += cutghost[0];
    hi[1] += cutghost[1];
    hi[2] += cutghost[2];
    domain->bbox(lo,hi,bsubboxlo,bsubboxhi);
  }



  bbox[0] = bboxhi[0] - bboxlo[0];
  bbox[1] = bboxhi[1] - bboxlo[1];
  bbox[2] = bboxhi[2] - bboxlo[2];
  bsubbox[0]=bsubboxhi[0]-bsubboxlo[0];
  bsubbox[1]=bsubboxhi[1]-bsubboxlo[1];
  bsubbox[2]=bsubboxhi[2]-bsubboxlo[2];
  // optimal bin size is roughly 1/2 the cutoff
  // for BIN style, binsize = 1/2 of max neighbor cutoff
  // for MULTI style, binsize = 1/2 of min neighbor cutoff
  // special case of all cutoffs = 0.0, binsize = box size

  //compute optimal binsize using the average element ebox size for me
  //only use if greater than half the force cutoff radius
  double average_size=0;
  for(int element_index=0; element_index < atom->nlocal+atom->nghost; element_index++){
    if(element_type[element_index]){
    double *current_ebox;

  current_ebox = eboxes[ebox_ref[element_index]];
  //define this elements bounding box
  ebounding_boxlo[0] = current_ebox[0];
  ebounding_boxlo[1] = current_ebox[1];
  ebounding_boxlo[2] = current_ebox[2];
  ebounding_boxhi[0] = current_ebox[3];
  ebounding_boxhi[1] = current_ebox[4];
  ebounding_boxhi[2] = current_ebox[5];
  for(int idim=0; idim < dimension; idim++){
  if(ebounding_boxhi[idim]-ebounding_boxlo[idim]>cut_max)
  average_size+=ebounding_boxhi[idim]-ebounding_boxlo[idim];
  else
  average_size+=cut_max;
  }
  }
  else{
  average_size+=dimension*cut_max;
  }
  }
  if(atom->nlocal+atom->nghost!=0)
  average_size=average_size/(dimension*(atom->nlocal+atom->nghost));
  if(average_size<=cut_max) average_size=cut_max;
  double binsize_optimal;
  if (binsizeflag) binsize_optimal = binsize_user;
  else if (style == BIN) binsize_optimal = 0.5*cutneighmax;
  else binsize_optimal = 0.5*cutneighmin;
  if (binsize_optimal == 0.0) binsize_optimal = bbox[0];
  binsize_optimal=0.5*average_size;
  double binsizeinv = 1.0/binsize_optimal;

  // test for too many global bins in any dimension due to huge global domain

  if (bbox[0]*binsizeinv > MAXSMALLINT || bbox[1]*binsizeinv > MAXSMALLINT ||
      bbox[2]*binsizeinv > MAXSMALLINT)
    error->all(FLERR,"Domain too large for neighbor bins");

  // create actual bins
  // always have one bin even if cutoff > bbox
  // for 2d, nbinz = 1

  nbinx = static_cast<int> (bsubbox[0]*binsizeinv);
  nbiny = static_cast<int> (bsubbox[1]*binsizeinv);
  if (dimension == 3) nbinz = static_cast<int> (bsubbox[2]*binsizeinv);
  else nbinz = 1;

  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;

  // compute actual bin size for nbins to fit into box exactly
  // error if actual bin size << cutoff, since will create a zillion bins
  // this happens when nbin = 1 and box size << cutoff
  // typically due to non-periodic, flat system in a particular dim
  // in that extreme case, should use NSQ not BIN neighbor style

  binsizex = bsubbox[0]/nbinx;
  binsizey = bsubbox[1]/nbiny;
  binsizez = bsubbox[2]/nbinz;

  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;

  if (binsize_optimal*bininvx > CUT2BIN_RATIO ||
      binsize_optimal*bininvy > CUT2BIN_RATIO ||
      binsize_optimal*bininvz > CUT2BIN_RATIO)
    error->all(FLERR,"Cannot use neighbor bins - box size << cutoff");

  // mbinlo/hi = lowest and highest global bins my ghost atoms could be in
  // coord = lowest and highest values of coords for my ghost atoms
  // static_cast(-1.5) = -1, so subract additional -1
  // add in SMALL for round-off safety

  mbinzlo=mbinylo=mbinxlo=-1;

  // extend bins by 1 to insure stencil extent is included
  // for 2d, only 1 bin in z

  mbinx = nbinx+2;

  mbiny = nbiny+2;

  mbinz = nbinz+2;

  bigint bbin = ((bigint) mbinx) * ((bigint) mbiny) * ((bigint) mbinz) + 1;
  if (bbin > MAXSMALLINT) error->one(FLERR,"Too many neighbor bins");
  if (bbin<=0) error->one(FLERR,"Too few neighbor bins");
  mbins = bbin;
}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms
------------------------------------------------------------------------- */

void NBinCAC::bin_atoms()
{
  int i,ibin;
  if(atom->CAC_comm_flag==0)
  error->all(FLERR,"Cannot use the CAC method without the CAC comm style");

  //grow element bin overlap arrays if needed
  if(atom->nlocal+atom->nghost>max_nall)
  expand_overlap_arrays(atom->nlocal+atom->nghost);
  //initialize bin overlap arrays
  for(int i=0; i<atom->nlocal+atom->nghost; i++)
  nbin_element_overlap[i]=0;
  
  //for (i = 0; i < mbins; i++) binhead[i] = -1;

  // bin in reverse order so linked list will be in forward order
  // also puts ghost atoms at end of list, which is necessary

  double **x = atom->x;
  int **element_scale = atom->element_scale;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int current_bin_ncount;
  int *element_type= atom->element_type;
  int current_bin;
  CAC_setup_bins(0);
  CAC_bin_atoms_setup(atom->nlocal+atom->nghost);
  //initialize bin counts
  for(int init=0; init<maxbin; init++)
  bin_ncontent[init]=0;

  if (includegroup) {
    int bitmask = group->bitmask[includegroup];
    foreign_boxes=0;
    for (i = nlocal; i < nall; i++) {
      if (mask[i] & bitmask) {
      //returns the set of bins an elements bounding box overlaps
      //atoms are binned typically
      ibin = element2bins(i);
      if(element_type[i]!=0) rboundingbox2bins(i);
      atom2bin[i] = ibin;
      if(element_type[i]==0){
            if(!atom->bin_foreign){
            current_bin_ncount=bin_ncontent[ibin];
            if(bin_ncontent[ibin]==MAXBINCONTENT+bin_expansion_counts[ibin]*EXPAND){
            bin_expansion_counts[ibin]++;
            memory->grow(bin_content[ibin],MAXBINCONTENT+bin_expansion_counts[ibin]*EXPAND,"neigh:bin_content grow");
            }
            bin_content[ibin][current_bin_ncount] = i;
            bin_ncontent[ibin]++;
            }

            bin_element_overlap[i][nbin_element_overlap[i]]=ibin;
            nbin_element_overlap[i]++;
      }
      else{
      for(int overlapx=bin_overlap_limits[0]; overlapx<=bin_overlap_limits[3]; overlapx++){
        for(int overlapy=bin_overlap_limits[1]; overlapy<=bin_overlap_limits[4]; overlapy++){
          for(int overlapz=bin_overlap_limits[2]; overlapz<=bin_overlap_limits[5]; overlapz++){
            current_bin=(overlapz-mbinzlo)*mbiny*mbinx + (overlapy-mbinylo)*mbinx + (overlapx-mbinxlo);
            if(!atom->bin_foreign){
            current_bin_ncount=bin_ncontent[current_bin];
            //check bin memory allocation is large enough and grow if needed
            if(bin_ncontent[current_bin]==MAXBINCONTENT+bin_expansion_counts[current_bin]*EXPAND){
            bin_expansion_counts[current_bin]++;
            memory->grow(bin_content[current_bin],MAXBINCONTENT+bin_expansion_counts[current_bin]*EXPAND,"neigh:bin_content grow");
            }
            bin_content[current_bin][current_bin_ncount] = i;
            bin_ncontent[current_bin]++;
            }
            if(overlapx>rbin_overlap_limits[3]||overlapx<rbin_overlap_limits[0]) continue;
            if(overlapy>rbin_overlap_limits[4]||overlapy<rbin_overlap_limits[1]) continue;
            if(overlapz>rbin_overlap_limits[5]||overlapz<rbin_overlap_limits[2]) continue;
            if(nbin_element_overlap[i]>=max_nbin_overlap[i]){
              max_nbin_overlap[i]+=MAXBINOVERLAP;
              memory->grow(bin_element_overlap[i],max_nbin_overlap[i],"nbin_CAC:bin_element_overlap[i] grow");
            }
            bin_element_overlap[i][nbin_element_overlap[i]]=current_bin;
            nbin_element_overlap[i]++;
          }
        }
      }
      }
      }
    }
    for (i = 0; i < atom->nfirst; i++) {
      //returns the set of bins an elements bounding box overlaps
      //atoms are binned typically
      ibin = element2bins(i);
      if(element_type[i]!=0) rboundingbox2bins(i);
      atom2bin[i] = ibin;
      if(element_type[i]==0){
            if(!atom->bin_foreign){
            current_bin_ncount=bin_ncontent[ibin];
            if(bin_ncontent[ibin]==MAXBINCONTENT+bin_expansion_counts[ibin]*EXPAND){
            bin_expansion_counts[ibin]++;
            memory->grow(bin_content[ibin],MAXBINCONTENT+bin_expansion_counts[ibin]*EXPAND,"neigh:bin_content grow");
            }
            bin_content[ibin][current_bin_ncount] = i;
            bin_ncontent[ibin]++;
            }

            bin_element_overlap[i][nbin_element_overlap[i]]=ibin;
            nbin_element_overlap[i]++;
      }
      else{
      for(int overlapx=bin_overlap_limits[0]; overlapx<=bin_overlap_limits[3]; overlapx++){
        for(int overlapy=bin_overlap_limits[1]; overlapy<=bin_overlap_limits[4]; overlapy++){
          for(int overlapz=bin_overlap_limits[2]; overlapz<=bin_overlap_limits[5]; overlapz++){
            current_bin=(overlapz-mbinzlo)*mbiny*mbinx + (overlapy-mbinylo)*mbinx + (overlapx-mbinxlo);
            if(!atom->bin_foreign){
            current_bin_ncount=bin_ncontent[current_bin];
            //check bin memory allocation is large enough and grow if needed
            if(bin_ncontent[current_bin]==MAXBINCONTENT+bin_expansion_counts[current_bin]*EXPAND){
            bin_expansion_counts[current_bin]++;
            memory->grow(bin_content[current_bin],MAXBINCONTENT+bin_expansion_counts[current_bin]*EXPAND,"neigh:bin_content grow");
            }
            bin_content[current_bin][current_bin_ncount] = i;
            bin_ncontent[current_bin]++;
            }
            if(overlapx>rbin_overlap_limits[3]||overlapx<rbin_overlap_limits[0]) continue;
            if(overlapy>rbin_overlap_limits[4]||overlapy<rbin_overlap_limits[1]) continue;
            if(overlapz>rbin_overlap_limits[5]||overlapz<rbin_overlap_limits[2]) continue;
            if(nbin_element_overlap[i]>=max_nbin_overlap[i]){
              max_nbin_overlap[i]+=MAXBINOVERLAP;
              memory->grow(bin_element_overlap[i],max_nbin_overlap[i],"nbin_CAC:bin_element_overlap[i] grow");
            }
            bin_element_overlap[i][nbin_element_overlap[i]]=current_bin;
            nbin_element_overlap[i]++;
          }
        }
      }
      }
    }

  }
  else {
    foreign_boxes=0;
    for (i = 0; i < nall; i++) {

      //returns the set of bins an elements bounding box overlaps
      //atoms are binned typically
      ibin = element2bins(i);
      if(element_type[i]!=0) rboundingbox2bins(i);
      if(ibin<0) error->one(FLERR," negative bin index");
      if(ibin>=mbins) error->one(FLERR," excessive bin index");
      if(element_type[i]==0){
            if(!atom->bin_foreign){
            current_bin_ncount=bin_ncontent[ibin];
            if(bin_ncontent[ibin]==MAXBINCONTENT+bin_expansion_counts[ibin]*EXPAND){
            bin_expansion_counts[ibin]++;
            memory->grow(bin_content[ibin],MAXBINCONTENT+bin_expansion_counts[ibin]*EXPAND,"neigh:bin_content grow");
            }
            bin_content[ibin][current_bin_ncount] = i;
            bin_ncontent[ibin]++;
            }

            bin_element_overlap[i][nbin_element_overlap[i]]=ibin;
            nbin_element_overlap[i]++;
      }
      if(element_type[i]!=0){
      for(int overlapx=bin_overlap_limits[0]; overlapx<=bin_overlap_limits[3]; overlapx++){
        for(int overlapy=bin_overlap_limits[1]; overlapy<=bin_overlap_limits[4]; overlapy++){
          for(int overlapz=bin_overlap_limits[2]; overlapz<=bin_overlap_limits[5]; overlapz++){
                        current_bin=(overlapz-mbinzlo)*mbiny*mbinx + (overlapy-mbinylo)*mbinx + (overlapx-mbinxlo);
            if(current_bin<0) error->one(FLERR," negative bin index");
            if(current_bin>=mbins) error->one(FLERR," excessive bin index");
            if(!atom->bin_foreign){
                        current_bin_ncount=bin_ncontent[current_bin];
            //check bin memory allocation is large enough and grow if needed
            if(bin_ncontent[current_bin]==MAXBINCONTENT+bin_expansion_counts[current_bin]*EXPAND){
            bin_expansion_counts[current_bin]++;
            memory->grow(bin_content[current_bin],MAXBINCONTENT+bin_expansion_counts[current_bin]*EXPAND,"neigh:bin_content grow");
            }
            bin_content[current_bin][current_bin_ncount] = i;
            bin_ncontent[current_bin]++;
            }
            if(overlapx>rbin_overlap_limits[3]||overlapx<rbin_overlap_limits[0]) continue;
            if(overlapy>rbin_overlap_limits[4]||overlapy<rbin_overlap_limits[1]) continue;
            if(overlapz>rbin_overlap_limits[5]||overlapz<rbin_overlap_limits[2]) continue;
            if(nbin_element_overlap[i]>=max_nbin_overlap[i]){
              max_nbin_overlap[i]+=MAXBINOVERLAP;
              memory->grow(bin_element_overlap[i],max_nbin_overlap[i],"nbin_CAC:bin_element_overlap[i] grow");
            }
            bin_element_overlap[i][nbin_element_overlap[i]]=current_bin;
            nbin_element_overlap[i]++;
          }
        }
      }
      }
    }
  }

  //bin foreign eboxes for Comm CAC if needed
  if(atom->bin_foreign){
    foreign_boxes=1;
    for(int i=0; i<atom->nforeign_eboxes; i++){
    ibin = element2bins(i);
     for(int overlapx=bin_overlap_limits[0]; overlapx<=bin_overlap_limits[3]; overlapx++)
        for(int overlapy=bin_overlap_limits[1]; overlapy<=bin_overlap_limits[4]; overlapy++)
          for(int overlapz=bin_overlap_limits[2]; overlapz<=bin_overlap_limits[5]; overlapz++){
            current_bin=(overlapz-mbinzlo)*mbiny*mbinx + (overlapy-mbinylo)*mbinx + (overlapx-mbinxlo);
            if(current_bin<0) error->one(FLERR," negative bin index");
            if(current_bin>=mbins) error->one(FLERR," excessive bin index");
            current_bin_ncount=bin_ncontent[current_bin];
            //check bin memory allocation is large enough and grow if needed
            if(bin_ncontent[current_bin]==MAXBINCONTENT+bin_expansion_counts[current_bin]*EXPAND){
            bin_expansion_counts[current_bin]++;
            memory->grow(bin_content[current_bin],MAXBINCONTENT+bin_expansion_counts[current_bin]*EXPAND,"neigh:bin_content grow");
            }
            bin_content[current_bin][current_bin_ncount] = i;
            bin_ncontent[current_bin]++;

          }
    }
  }
}

/* ----------------------------------------------------------------------
   convert atom and CAC element coords into local bin #
   for orthogonal, only ghost atoms will have coord >= bboxhi or coord < bboxlo
     take special care to insure ghosts are in correct bins even w/ roundoff
     hi ghost atoms = nbin,nbin+1,etc
     owned atoms = 0 to nbin-1
     lo ghost atoms = -1,-2,etc
     this is necessary so that both procs on either side of PBC
       treat a pair of atoms straddling the PBC in a consistent way
   for triclinic, doesn't matter since stencil & neigh list built differently
------------------------------------------------------------------------- */

 int NBinCAC::element2bins(int element_index)
{
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  double bounding_boxlo[3];
  double bounding_boxhi[3];
  double *cutghost = comm->cutghost;

  double **eboxes=atom->eboxes;
  double **foreign_eboxes=atom->foreign_eboxes;
  int *ebox_ref=atom->ebox_ref;
  double *x;

  int *nodes_per_element_list = atom->nodes_per_element_list;
  int ix,iy,iz;
  int ixl,iyl,izl,ixh,iyh,izh;
  if(!foreign_boxes){
   x = atom->x[element_index];
  }
  //typical binning for atoms
  if(!foreign_boxes){
  if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2]))
    error->one(FLERR,"Non-numeric positions - simulation unstable");

  if (x[0] > bsubboxhi[0])
    ix = mbinx-2;
  else if (x[0] >= bsubboxlo[0]) {
    ix = static_cast<int> ((x[0]-bsubboxlo[0])*bininvx);
    ix = MIN(ix,mbinx-3);
  } else
    ix = -1;

if (x[1] > bsubboxhi[1])
    iy = mbiny-2;
  else if (x[1] >= bsubboxlo[1]) {
    iy = static_cast<int> ((x[1]-bsubboxlo[1])*bininvy);
    iy = MIN(iy,mbiny-3);
  } else
    iy = -1;

if (x[2] > bsubboxhi[2])
    iz = mbinz-2;
  else if (x[2] >= bsubboxlo[2]) {
    iz = static_cast<int> ((x[2]-bsubboxlo[2])*bininvz);
    iz = MIN(iz,mbinz-3);
  } else
    iz = -1;
  }

  int decision_flag;
  if (foreign_boxes) decision_flag=1;
  else decision_flag=element_type[element_index];

  //calculate the set of bins this element's bounding box overlaps
  if(decision_flag){

  double *current_ebox;
  if(!foreign_boxes)
  current_ebox = eboxes[ebox_ref[element_index]];
  else
  current_ebox = foreign_eboxes[element_index];
  bounding_boxlo[0] = current_ebox[0];
  bounding_boxlo[1] = current_ebox[1];
  bounding_boxlo[2] = current_ebox[2];
  bounding_boxhi[0] = current_ebox[3];
  bounding_boxhi[1] = current_ebox[4];
  bounding_boxhi[2] = current_ebox[5];

  if(!foreign_boxes&&atom->bin_foreign){
  bounding_boxlo[0] += cutghost[0];
  bounding_boxlo[1] += cutghost[1];
  bounding_boxlo[2] += cutghost[2];
  bounding_boxhi[0] -= cutghost[0];
  bounding_boxhi[1] -= cutghost[1];
  bounding_boxhi[2] -= cutghost[2];
  }

  //compute lowest overlap id for each dimension
  if (bounding_boxlo[0] > bsubboxhi[0])
    ixl = mbinx - 2;
  else if (bounding_boxlo[0] >= bsubboxlo[0]) {
    ixl = static_cast<int> ((bounding_boxlo[0]-bsubboxlo[0])*bininvx);
    ixl = MIN(ixl,mbinx-3);
  }
  else
    ixl = -1;

  if (bounding_boxlo[1] > bsubboxhi[1])
    iyl = mbiny-2;
  else if (bounding_boxlo[1] >= bsubboxlo[1]) {
    iyl = static_cast<int> ((bounding_boxlo[1]-bsubboxlo[1])*bininvy);
    iyl = MIN(iyl,mbiny-3);
  }
  else
    iyl = -1;

  if (bounding_boxlo[2] > bsubboxhi[2])
    izl = mbinz-2;
  else if (bounding_boxlo[2] >= bsubboxlo[2]) {
    izl = static_cast<int> ((bounding_boxlo[2]-bsubboxlo[2])*bininvz);
    izl = MIN(izl,mbinz-3);
  }
  else
    izl = -1;

  //compute highest overlap id for each dimension
  if (bounding_boxhi[0] > bsubboxhi[0])
    ixh = mbinx - 2;
  else if (bounding_boxhi[0] >= bsubboxlo[0]) {
    ixh = static_cast<int> ((bounding_boxhi[0]-bsubboxlo[0])*bininvx);
    ixh = MIN(ixh,mbinx-3);
  }
  else
    ixh = -1;

  if (bounding_boxhi[1] > bsubboxhi[1])
    iyh = mbiny-2;
  else if (bounding_boxhi[1] >= bsubboxlo[1]) {
    iyh = static_cast<int> ((bounding_boxhi[1]-bsubboxlo[1])*bininvy);
    iyh = MIN(iyh,mbiny-3);
  }
  else
    iyh = -1;

  if (bounding_boxhi[2] > bsubboxhi[2])
    izh = mbinz-2;
  else if (bounding_boxhi[2] >= bsubboxlo[2]) {
    izh = static_cast<int> ((bounding_boxhi[2]-bsubboxlo[2])*bininvz);
    izh = MIN(izh,mbinz-3);
  }
  else
    izh = -1;

  bin_overlap_limits[0]=ixl;
  bin_overlap_limits[1]=iyl;
  bin_overlap_limits[2]=izl;
  bin_overlap_limits[3]=ixh;
  bin_overlap_limits[4]=iyh;
  bin_overlap_limits[5]=izh;
  }
  return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
}

/* ----------------------------------------------------------------------
   convert reduced element bounding boxes to overlapped bins
------------------------------------------------------------------------- */

 void NBinCAC::rboundingbox2bins(int element_index)
{
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  double bounding_boxlo[3];
  double bounding_boxhi[3];
  double *cutghost = comm->cutghost;

  double **eboxes=atom->eboxes;
  double **foreign_eboxes=atom->foreign_eboxes;
  int *ebox_ref=atom->ebox_ref;
  double *x;

  int *nodes_per_element_list = atom->nodes_per_element_list;
  int ix,iy,iz;
  int ixl,iyl,izl,ixh,iyh,izh;


  //calculate the set of bins this element's bounding box overlaps

  double *current_ebox;
  current_ebox = eboxes[ebox_ref[element_index]];
  bounding_boxlo[0] = current_ebox[0];
  bounding_boxlo[1] = current_ebox[1];
  bounding_boxlo[2] = current_ebox[2];
  bounding_boxhi[0] = current_ebox[3];
  bounding_boxhi[1] = current_ebox[4];
  bounding_boxhi[2] = current_ebox[5];
  bounding_boxlo[0] += cutghost[0];
  bounding_boxlo[1] += cutghost[1];
  bounding_boxlo[2] += cutghost[2];
  bounding_boxhi[0] -= cutghost[0];
  bounding_boxhi[1] -= cutghost[1];
  bounding_boxhi[2] -= cutghost[2];
  

  //compute lowest overlap id for each dimension
  if (bounding_boxlo[0] > bsubboxhi[0])
    ixl = mbinx - 2;
  else if (bounding_boxlo[0] >= bsubboxlo[0]) {
    ixl = static_cast<int> ((bounding_boxlo[0]-bsubboxlo[0])*bininvx);
    ixl = MIN(ixl,mbinx-3);
  }
  else
    ixl = -1;

  if (bounding_boxlo[1] > bsubboxhi[1])
    iyl = mbiny-2;
  else if (bounding_boxlo[1] >= bsubboxlo[1]) {
    iyl = static_cast<int> ((bounding_boxlo[1]-bsubboxlo[1])*bininvy);
    iyl = MIN(iyl,mbiny-3);
  }
  else
    iyl = -1;

  if (bounding_boxlo[2] > bsubboxhi[2])
    izl = mbinz-2;
  else if (bounding_boxlo[2] >= bsubboxlo[2]) {
    izl = static_cast<int> ((bounding_boxlo[2]-bsubboxlo[2])*bininvz);
    izl = MIN(izl,mbinz-3);
  }
  else
    izl = -1;

  //compute highest overlap id for each dimension
  if (bounding_boxhi[0] > bsubboxhi[0])
    ixh = mbinx - 2;
  else if (bounding_boxhi[0] >= bsubboxlo[0]) {
    ixh = static_cast<int> ((bounding_boxhi[0]-bsubboxlo[0])*bininvx);
    ixh = MIN(ixh,mbinx-3);
  }
  else
    ixh = -1;

  if (bounding_boxhi[1] > bsubboxhi[1])
    iyh = mbiny-2;
  else if (bounding_boxhi[1] >= bsubboxlo[1]) {
    iyh = static_cast<int> ((bounding_boxhi[1]-bsubboxlo[1])*bininvy);
    iyh = MIN(iyh,mbiny-3);
  }
  else
    iyh = -1;

  if (bounding_boxhi[2] > bsubboxhi[2])
    izh = mbinz-2;
  else if (bounding_boxhi[2] >= bsubboxlo[2]) {
    izh = static_cast<int> ((bounding_boxhi[2]-bsubboxlo[2])*bininvz);
    izh = MIN(izh,mbinz-3);
  }
  else
    izh = -1;

  rbin_overlap_limits[0]=ixl;
  rbin_overlap_limits[1]=iyl;
  rbin_overlap_limits[2]=izl;
  rbin_overlap_limits[3]=ixh;
  rbin_overlap_limits[4]=iyh;
  rbin_overlap_limits[5]=izh;
  
}

//resize bin element overlap arrays
void NBinCAC::expand_overlap_arrays(int size){
if(max_nall==0){
max_nall=size+100;
bin_element_overlap= (int **) memory->smalloc(max_nall*sizeof(int *), "bin_CAC:bin_element_overlap");
for(int i=0 ; i<max_nall; i++)
memory->create(bin_element_overlap[i],MAXBINOVERLAP,"bin_CAC:bin_element_overlap[i]");
memory->grow(nbin_element_overlap,max_nall,"bin_CAC:nbin_element_overlap");
memory->grow(max_nbin_overlap,max_nall,"bin_CAC:max_nbin_overlap");
for(int i=0; i<max_nall; i++){
  max_nbin_overlap[i]=MAXBINOVERLAP;
}
}
else{
size+=100;
bin_element_overlap= (int **) memory->srealloc(bin_element_overlap,size*sizeof(int *), "bin_CAC:bin_element_overlap");
for(int i=max_nall ; i<size; i++)
memory->create(bin_element_overlap[i],MAXBINOVERLAP,"bin_CAC:bin_element_overlap[i]");
memory->grow(nbin_element_overlap,size,"bin_CAC:nbin_element_overlap");
memory->grow(max_nbin_overlap,size,"bin_CAC:max_nbin_overlap");
//initialize bin overlap arrays
for(int i=max_nall; i<size; i++){
  max_nbin_overlap[i]=MAXBINOVERLAP;
}
max_nall=size;
}
}

/* ----------------------------------------------------------------------
   convert quadrature point coordinate into local bin #
   for orthogonal, only ghost atoms will have coord >= bboxhi or coord < bboxlo
     take special care to insure ghosts are in correct bins even w/ roundoff
     hi ghost atoms = nbin,nbin+1,etc
     owned atoms = 0 to nbin-1
     lo ghost atoms = -1,-2,etc
     this is necessary so that both procs on either side of PBC
       treat a pair of atoms straddling the PBC in a consistent way
   for triclinic, doesn't matter since stencil & neigh list built differently
------------------------------------------------------------------------- */

 int NBinCAC::coord2bin(double *x)
{
  int ix,iy,iz;
  
  //binning for a quadrature point
  
  if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2]))
    error->one(FLERR,"Non-numeric positions - simulation unstable");

  if (x[0] >= bsubboxhi[0])
    error->one(FLERR,"quadrature_point outside of bin domain");
  else if (x[0] >= bsubboxlo[0]) {
    ix = static_cast<int> ((x[0]-bsubboxlo[0])*bininvx);
    ix = MIN(ix,mbinx-3);
  } else
    error->one(FLERR,"quadrature_point outside of bin domain");

  if (x[1] >= bsubboxhi[1])
    error->one(FLERR,"quadrature_point outside of bin domain");
  else if (x[1] >= bsubboxlo[1]) {
    iy = static_cast<int> ((x[1]-bsubboxlo[1])*bininvy);
    iy = MIN(iy,mbiny-3);
  } else
    error->one(FLERR,"quadrature_point outside of bin domain");

  if (x[2] >= bsubboxhi[2])
    error->one(FLERR,"quadrature_point outside of bin domain");
  else if (x[2] >= bsubboxlo[2]) {
    iz = static_cast<int> ((x[2]-bsubboxlo[2])*bininvz);
    iz = MIN(iz,mbinz-3);
  } else
    error->one(FLERR,"quadrature_point outside of bin domain");


  //return the bin id of this quadrature point
  
  return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
}