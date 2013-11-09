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
   Contributing author (triclinic) : Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "comm.h"
#include "universe.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "neighbor.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "output.h"
#include "dump.h"
#include "procmap.h"
#include "math_extra.h"
#include "error.h"
#include "memory.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace LAMMPS_NS;

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000
#define BIG 1.0e20

enum{SINGLE,MULTI};
enum{MULTIPLE};                   // same as in ProcMap
enum{ONELEVEL,TWOLEVEL,NUMA,CUSTOM};
enum{CART,CARTREORDER,XYZ};

/* ----------------------------------------------------------------------
   setup MPI and allocate buffer space
------------------------------------------------------------------------- */

Comm::Comm(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  user_procgrid[0] = user_procgrid[1] = user_procgrid[2] = 0;
  coregrid[0] = coregrid[1] = coregrid[2] = 1;
  gridflag = ONELEVEL;
  mapflag = CART;
  customfile = NULL;
  recv_from_partition = send_to_partition = -1;
  otherflag = 0;
  outfile = NULL;

  grid2proc = NULL;

  bordergroup = 0;
  style = SINGLE;
  uniform = 1;
  xsplit = ysplit = zsplit = NULL;
  multilo = multihi = NULL;
  cutghostmulti = NULL;
  cutghostuser = 0.0;
  ghost_velocity = 0;

  // use of OpenMP threads
  // query OpenMP for number of threads/process set by user at run-time
  // if the OMP_NUM_THREADS environment variable is not set, we default
  // to using 1 thread. This follows the principle of the least surprise,
  // while practically all OpenMP implementations violate it by using
  // as many threads as there are (virtual) CPU cores by default.

  nthreads = 1;
#ifdef _OPENMP
  if (getenv("OMP_NUM_THREADS") == NULL) {
    nthreads = 1;
    if (me == 0)
      error->warning(FLERR,"OMP_NUM_THREADS environment is not set.");
  } else {
    nthreads = omp_get_max_threads();
  }

  // enforce consistent number of threads across all MPI tasks

  MPI_Bcast(&nthreads,1,MPI_INT,0,world);
  omp_set_num_threads(nthreads);

  if (me == 0) {
    if (screen)
      fprintf(screen,"  using %d OpenMP thread(s) per MPI task\n",nthreads);
    if (logfile)
      fprintf(logfile,"  using %d OpenMP thread(s) per MPI task\n",nthreads);
  }
#endif

  // initialize comm buffers & exchange memory
  // NOTE: allow for AtomVec to set maxexchange_atom, e.g. for atom_style body

  maxexchange_atom = maxexchange_fix = 0;
  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;

  maxsend = BUFMIN;
  memory->create(buf_send,maxsend+bufextra,"comm:buf_send");
  maxrecv = BUFMIN;
  memory->create(buf_recv,maxrecv,"comm:buf_recv");

  maxswap = 6;
  allocate_swap(maxswap);

  sendlist = (int **) memory->smalloc(maxswap*sizeof(int *),"comm:sendlist");
  memory->create(maxsendlist,maxswap,"comm:maxsendlist");
  for (int i = 0; i < maxswap; i++) {
    maxsendlist[i] = BUFMIN;
    memory->create(sendlist[i],BUFMIN,"comm:sendlist[i]");
  }
}

/* ---------------------------------------------------------------------- */

Comm::~Comm()
{
  memory->destroy(xsplit);
  memory->destroy(ysplit);
  memory->destroy(zsplit);

  delete [] customfile;
  delete [] outfile;

  memory->destroy(grid2proc);

  free_swap();
  if (style == MULTI) {
    free_multi();
    memory->destroy(cutghostmulti);
  }

  if (sendlist) for (int i = 0; i < maxswap; i++) memory->destroy(sendlist[i]);
  memory->sfree(sendlist);
  memory->destroy(maxsendlist);

  memory->destroy(buf_send);
  memory->destroy(buf_recv);
}

/* ----------------------------------------------------------------------
   create 3d grid of procs based on Nprocs and box size & shape
   map processors to grid, setup xyz split for a uniform grid
------------------------------------------------------------------------- */

void Comm::set_proc_grid(int outflag)
{
  // recv 3d proc grid of another partition if my 3d grid depends on it

  if (recv_from_partition >= 0) {
    MPI_Status status;
    if (me == 0) {
      MPI_Recv(other_procgrid,3,MPI_INT,
               universe->root_proc[recv_from_partition],0,
               universe->uworld,&status);
      MPI_Recv(other_coregrid,3,MPI_INT,
               universe->root_proc[recv_from_partition],0,
               universe->uworld,&status);
    }
    MPI_Bcast(other_procgrid,3,MPI_INT,0,world);
    MPI_Bcast(other_coregrid,3,MPI_INT,0,world);
  }

  // create ProcMap class to create 3d grid and map procs to it

  ProcMap *pmap = new ProcMap(lmp);

  // create 3d grid of processors
  // produces procgrid and coregrid (if relevant)

  if (gridflag == ONELEVEL) {
    pmap->onelevel_grid(nprocs,user_procgrid,procgrid,
                        otherflag,other_style,other_procgrid,other_coregrid);

  } else if (gridflag == TWOLEVEL) {
    pmap->twolevel_grid(nprocs,user_procgrid,procgrid,
                        ncores,user_coregrid,coregrid,
                        otherflag,other_style,other_procgrid,other_coregrid);

  } else if (gridflag == NUMA) {
    pmap->numa_grid(nprocs,user_procgrid,procgrid,coregrid);

  } else if (gridflag == CUSTOM) {
    pmap->custom_grid(customfile,nprocs,user_procgrid,procgrid);
  }

  // error check on procgrid
  // should not be necessary due to ProcMap

  if (procgrid[0]*procgrid[1]*procgrid[2] != nprocs)
    error->all(FLERR,"Bad grid of processors");
  if (domain->dimension == 2 && procgrid[2] != 1)
    error->all(FLERR,"Processor count in z must be 1 for 2d simulation");

  // grid2proc[i][j][k] = proc that owns i,j,k location in 3d grid

  if (grid2proc) memory->destroy(grid2proc);
  memory->create(grid2proc,procgrid[0],procgrid[1],procgrid[2],
                 "comm:grid2proc");

  // map processor IDs to 3d processor grid
  // produces myloc, procneigh, grid2proc

  if (gridflag == ONELEVEL) {
    if (mapflag == CART)
      pmap->cart_map(0,procgrid,myloc,procneigh,grid2proc);
    else if (mapflag == CARTREORDER)
      pmap->cart_map(1,procgrid,myloc,procneigh,grid2proc);
    else if (mapflag == XYZ)
      pmap->xyz_map(xyz,procgrid,myloc,procneigh,grid2proc);

  } else if (gridflag == TWOLEVEL) {
    if (mapflag == CART)
      pmap->cart_map(0,procgrid,ncores,coregrid,myloc,procneigh,grid2proc);
    else if (mapflag == CARTREORDER)
      pmap->cart_map(1,procgrid,ncores,coregrid,myloc,procneigh,grid2proc);
    else if (mapflag == XYZ)
      pmap->xyz_map(xyz,procgrid,ncores,coregrid,myloc,procneigh,grid2proc);

  } else if (gridflag == NUMA) {
    pmap->numa_map(0,coregrid,myloc,procneigh,grid2proc);

  } else if (gridflag == CUSTOM) {
    pmap->custom_map(procgrid,myloc,procneigh,grid2proc);
  }

  // print 3d grid info to screen and logfile

  if (outflag && me == 0) {
    if (screen) {
      fprintf(screen,"  %d by %d by %d MPI processor grid\n",
              procgrid[0],procgrid[1],procgrid[2]);
      if (gridflag == NUMA || gridflag == TWOLEVEL)
        fprintf(screen,"  %d by %d by %d core grid within node\n",
                coregrid[0],coregrid[1],coregrid[2]);
    }
    if (logfile) {
      fprintf(logfile,"  %d by %d by %d MPI processor grid\n",
              procgrid[0],procgrid[1],procgrid[2]);
      if (gridflag == NUMA || gridflag == TWOLEVEL)
        fprintf(logfile,"  %d by %d by %d core grid within node\n",
                coregrid[0],coregrid[1],coregrid[2]);
    }
  }

  // print 3d grid details to outfile

  if (outfile) pmap->output(outfile,procgrid,grid2proc);

  // free ProcMap class

  delete pmap;

  // set xsplit,ysplit,zsplit for uniform spacings

  memory->destroy(xsplit);
  memory->destroy(ysplit);
  memory->destroy(zsplit);

  memory->create(xsplit,procgrid[0]+1,"comm:xsplit");
  memory->create(ysplit,procgrid[1]+1,"comm:ysplit");
  memory->create(zsplit,procgrid[2]+1,"comm:zsplit");

  for (int i = 0; i < procgrid[0]; i++) xsplit[i] = i * 1.0/procgrid[0];
  for (int i = 0; i < procgrid[1]; i++) ysplit[i] = i * 1.0/procgrid[1];
  for (int i = 0; i < procgrid[2]; i++) zsplit[i] = i * 1.0/procgrid[2];

  xsplit[procgrid[0]] = ysplit[procgrid[1]] = zsplit[procgrid[2]] = 1.0;

  // set lamda box params after procs are assigned
  // only set once unless load-balancing occurs

  if (domain->triclinic) domain->set_lamda_box();

  // send my 3d proc grid to another partition if requested

  if (send_to_partition >= 0) {
    if (me == 0) {
      MPI_Send(procgrid,3,MPI_INT,
               universe->root_proc[send_to_partition],0,
               universe->uworld);
      MPI_Send(coregrid,3,MPI_INT,
               universe->root_proc[send_to_partition],0,
               universe->uworld);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Comm::init()
{
  triclinic = domain->triclinic;
  map_style = atom->map_style;

  // comm_only = 1 if only x,f are exchanged in forward/reverse comm
  // comm_x_only = 0 if ghost_velocity since velocities are added

  comm_x_only = atom->avec->comm_x_only;
  comm_f_only = atom->avec->comm_f_only;
  if (ghost_velocity) comm_x_only = 0;

  // set per-atom sizes for forward/reverse/border comm
  // augment by velocity and fix quantities if needed

  size_forward = atom->avec->size_forward;
  size_reverse = atom->avec->size_reverse;
  size_border = atom->avec->size_border;

  if (ghost_velocity) size_forward += atom->avec->size_velocity;
  if (ghost_velocity) size_border += atom->avec->size_velocity;

  for (int i = 0; i < modify->nfix; i++)
    size_border += modify->fix[i]->comm_border;
  
  // maxexchange = max # of datums/atom in exchange communication
  // maxforward = # of datums in largest forward communication
  // maxreverse = # of datums in largest reverse communication
  // query pair,fix,compute,dump for their requirements
  // pair style can force reverse comm even if newton off

  maxexchange = BUFMIN + maxexchange_fix;
  maxforward = MAX(size_forward,size_border);
  maxreverse = size_reverse;

  if (force->pair) maxforward = MAX(maxforward,force->pair->comm_forward);
  if (force->pair) maxreverse = MAX(maxreverse,force->pair->comm_reverse);

  for (int i = 0; i < modify->nfix; i++) {
    maxforward = MAX(maxforward,modify->fix[i]->comm_forward);
    maxreverse = MAX(maxreverse,modify->fix[i]->comm_reverse);
  }

  for (int i = 0; i < modify->ncompute; i++) {
    maxforward = MAX(maxforward,modify->compute[i]->comm_forward);
    maxreverse = MAX(maxreverse,modify->compute[i]->comm_reverse);
  }

  for (int i = 0; i < output->ndump; i++) {
    maxforward = MAX(maxforward,output->dump[i]->comm_forward);
    maxreverse = MAX(maxreverse,output->dump[i]->comm_reverse);
  }

  if (force->newton == 0) maxreverse = 0;
  if (force->pair) maxreverse = MAX(maxreverse,force->pair->comm_reverse_off);

  // memory for multi-style communication

  if (style == MULTI && multilo == NULL) {
    allocate_multi(maxswap);
    memory->create(cutghostmulti,atom->ntypes+1,3,"comm:cutghostmulti");
  }
  if (style == SINGLE && multilo) {
    free_multi();
    memory->destroy(cutghostmulti);
  }
}

/* ----------------------------------------------------------------------
   setup spatial-decomposition communication patterns
   function of neighbor cutoff(s) & cutghostuser & current box size
   single style sets slab boundaries (slablo,slabhi) based on max cutoff
   multi style sets type-dependent slab boundaries (multilo,multihi)
------------------------------------------------------------------------- */

void Comm::setup()
{
  // cutghost[] = max distance at which ghost atoms need to be acquired
  // for orthogonal:
  //   cutghost is in box coords = neigh->cutghost in all 3 dims
  // for triclinic:
  //   neigh->cutghost = distance between tilted planes in box coords
  //   cutghost is in lamda coords = distance between those planes
  // for multi:
  //   cutghostmulti = same as cutghost, only for each atom type

  int i;
  int ntypes = atom->ntypes;
  double *prd,*sublo,*subhi;

  double cut = MAX(neighbor->cutneighmax,cutghostuser);

  if (triclinic == 0) {
    prd = domain->prd;
    sublo = domain->sublo;
    subhi = domain->subhi;
    cutghost[0] = cutghost[1] = cutghost[2] = cut;

    if (style == MULTI) {
      double *cuttype = neighbor->cuttype;
      for (i = 1; i <= ntypes; i++)
        cutghostmulti[i][0] = cutghostmulti[i][1] = cutghostmulti[i][2] =
          cuttype[i];
    }

  } else {
    prd = domain->prd_lamda;
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
    double *h_inv = domain->h_inv;
    double length0,length1,length2;
    length0 = sqrt(h_inv[0]*h_inv[0] + h_inv[5]*h_inv[5] + h_inv[4]*h_inv[4]);
    cutghost[0] = cut * length0;
    length1 = sqrt(h_inv[1]*h_inv[1] + h_inv[3]*h_inv[3]);
    cutghost[1] = cut * length1;
    length2 = h_inv[2];
    cutghost[2] = cut * length2;

    if (style == MULTI) {
      double *cuttype = neighbor->cuttype;
      for (i = 1; i <= ntypes; i++) {
        cutghostmulti[i][0] = cuttype[i] * length0;
        cutghostmulti[i][1] = cuttype[i] * length1;
        cutghostmulti[i][2] = cuttype[i] * length2;
      }
    }
  }

  // recvneed[idim][0/1] = # of procs away I recv atoms from, within cutghost
  //   0 = from left, 1 = from right
  //   do not cross non-periodic boundaries, need[2] = 0 for 2d
  // sendneed[idim][0/1] = # of procs away I send atoms to
  //   0 = to left, 1 = to right
  //   set equal to recvneed[idim][1/0] of neighbor proc
  // maxneed[idim] = max procs away any proc recvs atoms in either direction
  // uniform = 1 = uniform sized sub-domains:
  //   maxneed is directly computable from sub-domain size
  //     limit to procgrid-1 for non-PBC
  //   recvneed = maxneed except for procs near non-PBC
  //   sendneed = recvneed of neighbor on each side
  // uniform = 0 = non-uniform sized sub-domains:
  //   compute recvneed via updown() which accounts for non-PBC
  //   sendneed = recvneed of neighbor on each side
  //   maxneed via Allreduce() of recvneed

  int *periodicity = domain->periodicity;
  int left,right;

  if (uniform) {
    maxneed[0] = static_cast<int> (cutghost[0] * procgrid[0] / prd[0]) + 1;
    maxneed[1] = static_cast<int> (cutghost[1] * procgrid[1] / prd[1]) + 1;
    maxneed[2] = static_cast<int> (cutghost[2] * procgrid[2] / prd[2]) + 1;
    if (domain->dimension == 2) maxneed[2] = 0;
    if (!periodicity[0]) maxneed[0] = MIN(maxneed[0],procgrid[0]-1);
    if (!periodicity[1]) maxneed[1] = MIN(maxneed[1],procgrid[1]-1);
    if (!periodicity[2]) maxneed[2] = MIN(maxneed[2],procgrid[2]-1);

    if (!periodicity[0]) {
      recvneed[0][0] = MIN(maxneed[0],myloc[0]);
      recvneed[0][1] = MIN(maxneed[0],procgrid[0]-myloc[0]-1);
      left = myloc[0] - 1;
      if (left < 0) left = procgrid[0] - 1;
      sendneed[0][0] = MIN(maxneed[0],procgrid[0]-left-1);
      right = myloc[0] + 1;
      if (right == procgrid[0]) right = 0;
      sendneed[0][1] = MIN(maxneed[0],right);
    } else recvneed[0][0] = recvneed[0][1] =
             sendneed[0][0] = sendneed[0][1] = maxneed[0];

    if (!periodicity[1]) {
      recvneed[1][0] = MIN(maxneed[1],myloc[1]);
      recvneed[1][1] = MIN(maxneed[1],procgrid[1]-myloc[1]-1);
      left = myloc[1] - 1;
      if (left < 0) left = procgrid[1] - 1;
      sendneed[1][0] = MIN(maxneed[1],procgrid[1]-left-1);
      right = myloc[1] + 1;
      if (right == procgrid[1]) right = 0;
      sendneed[1][1] = MIN(maxneed[1],right);
    } else recvneed[1][0] = recvneed[1][1] =
             sendneed[1][0] = sendneed[1][1] = maxneed[1];

    if (!periodicity[2]) {
      recvneed[2][0] = MIN(maxneed[2],myloc[2]);
      recvneed[2][1] = MIN(maxneed[2],procgrid[2]-myloc[2]-1);
      left = myloc[2] - 1;
      if (left < 0) left = procgrid[2] - 1;
      sendneed[2][0] = MIN(maxneed[2],procgrid[2]-left-1);
      right = myloc[2] + 1;
      if (right == procgrid[2]) right = 0;
      sendneed[2][1] = MIN(maxneed[2],right);
    } else recvneed[2][0] = recvneed[2][1] =
             sendneed[2][0] = sendneed[2][1] = maxneed[2];

  } else {
    recvneed[0][0] = updown(0,0,myloc[0],prd[0],periodicity[0],xsplit);
    recvneed[0][1] = updown(0,1,myloc[0],prd[0],periodicity[0],xsplit);
    left = myloc[0] - 1;
    if (left < 0) left = procgrid[0] - 1;
    sendneed[0][0] = updown(0,1,left,prd[0],periodicity[0],xsplit);
    right = myloc[0] + 1;
    if (right == procgrid[0]) right = 0;
    sendneed[0][1] = updown(0,0,right,prd[0],periodicity[0],xsplit);

    recvneed[1][0] = updown(1,0,myloc[1],prd[1],periodicity[1],ysplit);
    recvneed[1][1] = updown(1,1,myloc[1],prd[1],periodicity[1],ysplit);
    left = myloc[1] - 1;
    if (left < 0) left = procgrid[1] - 1;
    sendneed[1][0] = updown(1,1,left,prd[1],periodicity[1],ysplit);
    right = myloc[1] + 1;
    if (right == procgrid[1]) right = 0;
    sendneed[1][1] = updown(1,0,right,prd[1],periodicity[1],ysplit);

    if (domain->dimension == 3) {
      recvneed[2][0] = updown(2,0,myloc[2],prd[2],periodicity[2],zsplit);
      recvneed[2][1] = updown(2,1,myloc[2],prd[2],periodicity[2],zsplit);
      left = myloc[2] - 1;
      if (left < 0) left = procgrid[2] - 1;
      sendneed[2][0] = updown(2,1,left,prd[2],periodicity[2],zsplit);
      right = myloc[2] + 1;
      if (right == procgrid[2]) right = 0;
      sendneed[2][1] = updown(2,0,right,prd[2],periodicity[2],zsplit);
    } else recvneed[2][0] = recvneed[2][1] =
             sendneed[2][0] = sendneed[2][1] = 0;

    int all[6];
    MPI_Allreduce(&recvneed[0][0],all,6,MPI_INT,MPI_MAX,world);
    maxneed[0] = MAX(all[0],all[1]);
    maxneed[1] = MAX(all[2],all[3]);
    maxneed[2] = MAX(all[4],all[5]);
  }

  // allocate comm memory

  nswap = 2 * (maxneed[0]+maxneed[1]+maxneed[2]);
  if (nswap > maxswap) grow_swap(nswap);

  // setup parameters for each exchange:
  // sendproc = proc to send to at each swap
  // recvproc = proc to recv from at each swap
  // for style SINGLE:
  //   slablo/slabhi = boundaries for slab of atoms to send at each swap
  //   use -BIG/midpt/BIG to insure all atoms included even if round-off occurs
  //   if round-off, atoms recvd across PBC can be < or > than subbox boundary
  //   note that borders() only loops over subset of atoms during each swap
  //   treat all as PBC here, non-PBC is handled in borders() via r/s need[][]
  // for style MULTI:
  //   multilo/multihi is same, with slablo/slabhi for each atom type
  // pbc_flag: 0 = nothing across a boundary, 1 = something across a boundary
  // pbc = -1/0/1 for PBC factor in each of 3/6 orthogonal/triclinic dirs
  // for triclinic, slablo/hi and pbc_border will be used in lamda (0-1) coords
  // 1st part of if statement is sending to the west/south/down
  // 2nd part of if statement is sending to the east/north/up

  int dim,ineed;

  int iswap = 0;
  for (dim = 0; dim < 3; dim++) {
    for (ineed = 0; ineed < 2*maxneed[dim]; ineed++) {
      pbc_flag[iswap] = 0;
      pbc[iswap][0] = pbc[iswap][1] = pbc[iswap][2] =
        pbc[iswap][3] = pbc[iswap][4] = pbc[iswap][5] = 0;

      if (ineed % 2 == 0) {
        sendproc[iswap] = procneigh[dim][0];
        recvproc[iswap] = procneigh[dim][1];
        if (style == SINGLE) {
          if (ineed < 2) slablo[iswap] = -BIG;
          else slablo[iswap] = 0.5 * (sublo[dim] + subhi[dim]);
          slabhi[iswap] = sublo[dim] + cutghost[dim];
        } else {
          for (i = 1; i <= ntypes; i++) {
            if (ineed < 2) multilo[iswap][i] = -BIG;
            else multilo[iswap][i] = 0.5 * (sublo[dim] + subhi[dim]);
            multihi[iswap][i] = sublo[dim] + cutghostmulti[i][dim];
          }
        }
        if (myloc[dim] == 0) {
          pbc_flag[iswap] = 1;
          pbc[iswap][dim] = 1;
          if (triclinic) {
            if (dim == 1) pbc[iswap][5] = 1;
            else if (dim == 2) pbc[iswap][4] = pbc[iswap][3] = 1;
          }
        }

      } else {
        sendproc[iswap] = procneigh[dim][1];
        recvproc[iswap] = procneigh[dim][0];
        if (style == SINGLE) {
          slablo[iswap] = subhi[dim] - cutghost[dim];
          if (ineed < 2) slabhi[iswap] = BIG;
          else slabhi[iswap] = 0.5 * (sublo[dim] + subhi[dim]);
        } else {
          for (i = 1; i <= ntypes; i++) {
            multilo[iswap][i] = subhi[dim] - cutghostmulti[i][dim];
            if (ineed < 2) multihi[iswap][i] = BIG;
            else multihi[iswap][i] = 0.5 * (sublo[dim] + subhi[dim]);
          }
        }
        if (myloc[dim] == procgrid[dim]-1) {
          pbc_flag[iswap] = 1;
          pbc[iswap][dim] = -1;
          if (triclinic) {
            if (dim == 1) pbc[iswap][5] = -1;
            else if (dim == 2) pbc[iswap][4] = pbc[iswap][3] = -1;
          }
        }
      }

      iswap++;
    }
  }
}

/* ----------------------------------------------------------------------
   walk up/down the extent of nearby processors in dim and dir
   loc = myloc of proc to start at
   dir = 0/1 = walk to left/right
   do not cross non-periodic boundaries
   is not called for z dim in 2d
   return how many procs away are needed to encompass cutghost away from loc
------------------------------------------------------------------------- */

int Comm::updown(int dim, int dir, int loc,
                 double prd, int periodicity, double *split)
{
  int index,count;
  double frac,delta;

  if (dir == 0) {
    frac = cutghost[dim]/prd;
    index = loc - 1;
    delta = 0.0;
    count = 0;
    while (delta < frac) {
      if (index < 0) {
        if (!periodicity) break;
        index = procgrid[dim] - 1;
      }
      count++;
      delta += split[index+1] - split[index];
      index--;
    }

  } else {
    frac = cutghost[dim]/prd;
    index = loc + 1;
    delta = 0.0;
    count = 0;
    while (delta < frac) {
      if (index >= procgrid[dim]) {
        if (!periodicity) break;
        index = 0;
      }
      count++;
      delta += split[index+1] - split[index];
      index++;
    }
  }

  return count;
}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void Comm::forward_comm(int dummy)
{
  int n;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;
  double **x = atom->x;
  double *buf;

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap] != me) {
      if (comm_x_only) {
        if (size_forward_recv[iswap]) buf = x[firstrecv[iswap]];
        else buf = NULL;
        if (size_forward_recv[iswap])
          MPI_Irecv(buf,size_forward_recv[iswap],MPI_DOUBLE,
                    recvproc[iswap],0,world,&request);
        n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
                            buf_send,pbc_flag[iswap],pbc[iswap]);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        if (size_forward_recv[iswap]) MPI_Wait(&request,&status);
      } else if (ghost_velocity) {
        if (size_forward_recv[iswap])
          MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,
                    recvproc[iswap],0,world,&request);
        n = avec->pack_comm_vel(sendnum[iswap],sendlist[iswap],
                                buf_send,pbc_flag[iswap],pbc[iswap]);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        if (size_forward_recv[iswap]) MPI_Wait(&request,&status);
        avec->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],buf_recv);
      } else {
        if (size_forward_recv[iswap])
          MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,
                    recvproc[iswap],0,world,&request);
        n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
                            buf_send,pbc_flag[iswap],pbc[iswap]);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        if (size_forward_recv[iswap]) MPI_Wait(&request,&status);
        avec->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_recv);
      }

    } else {
      if (comm_x_only) {
        if (sendnum[iswap])
          n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
                              x[firstrecv[iswap]],pbc_flag[iswap],
                              pbc[iswap]);
      } else if (ghost_velocity) {
        n = avec->pack_comm_vel(sendnum[iswap],sendlist[iswap],
                                buf_send,pbc_flag[iswap],pbc[iswap]);
        avec->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],buf_send);
      } else {
        n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
                            buf_send,pbc_flag[iswap],pbc[iswap]);
        avec->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_send);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication of forces on atoms every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void Comm::reverse_comm()
{
  int n;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;
  double **f = atom->f;
  double *buf;

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_f_only set, exchange or copy directly from f, don't pack

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (sendproc[iswap] != me) {
      if (comm_f_only) {
        if (size_reverse_recv[iswap])
          MPI_Irecv(buf_recv,size_reverse_recv[iswap],MPI_DOUBLE,
                    sendproc[iswap],0,world,&request);
        if (size_reverse_send[iswap]) buf = f[firstrecv[iswap]];
        else buf = NULL;
        if (size_reverse_send[iswap])
          MPI_Send(buf,size_reverse_send[iswap],MPI_DOUBLE,
                   recvproc[iswap],0,world);
        if (size_reverse_recv[iswap]) MPI_Wait(&request,&status);
      } else {
        if (size_reverse_recv[iswap])
          MPI_Irecv(buf_recv,size_reverse_recv[iswap],MPI_DOUBLE,
                    sendproc[iswap],0,world,&request);
        n = avec->pack_reverse(recvnum[iswap],firstrecv[iswap],buf_send);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,world);
        if (size_reverse_recv[iswap]) MPI_Wait(&request,&status);
      }
      avec->unpack_reverse(sendnum[iswap],sendlist[iswap],buf_recv);

    } else {
      if (comm_f_only) {
        if (sendnum[iswap])
            avec->unpack_reverse(sendnum[iswap],sendlist[iswap],
                                f[firstrecv[iswap]]);
      } else {
        n = avec->pack_reverse(recvnum[iswap],firstrecv[iswap],buf_send);
        avec->unpack_reverse(sendnum[iswap],sendlist[iswap],buf_send);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   exchange: move atoms to correct processors
   atoms exchanged with all 6 stencil neighbors
   send out atoms that have left my box, receive ones entering my box
   atoms will be lost if not inside some proc's box
     can happen if atom moves outside of non-periodic bounary
     or if atom moves more than one proc away
   this routine called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before exchange is called
------------------------------------------------------------------------- */

void Comm::exchange()
{
  int i,m,nsend,nrecv,nrecv1,nrecv2,nlocal;
  double lo,hi,value;
  double **x;
  double *sublo,*subhi,*buf;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;

  // clear global->local map for owned and ghost atoms
  // b/c atoms migrate to new procs in exchange() and
  //   new ghosts are created in borders()
  // map_set() is done at end of borders()
  // clear ghost count and any ghost bonus data internal to AtomVec

  if (map_style) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();

  // insure send buf is large enough for single atom
  // fixes can change per-atom size requirement on-the-fly

  int bufextra_old = bufextra;
  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;
  if (bufextra > bufextra_old)
    memory->grow(buf_send,maxsend+bufextra,"comm:buf_send");

  // subbox bounds for orthogonal or triclinic

  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // loop over dimensions

  for (int dim = 0; dim < 3; dim++) {

    // fill buffer with atoms leaving my box, using < and >=
    // when atom is deleted, fill it in with last atom

    x = atom->x;
    lo = sublo[dim];
    hi = subhi[dim];
    nlocal = atom->nlocal;
    i = nsend = 0;

    while (i < nlocal) {
      if (x[i][dim] < lo || x[i][dim] >= hi) {
        if (nsend > maxsend) grow_send(nsend,1);
        nsend += avec->pack_exchange(i,&buf_send[nsend]);
        avec->copy(nlocal-1,i,1);
        nlocal--;
      } else i++;
    }
    atom->nlocal = nlocal;

    // send/recv atoms in both directions
    // if 1 proc in dimension, no send/recv, set recv buf to send buf
    // if 2 procs in dimension, single send/recv
    // if more than 2 procs in dimension, send/recv to both neighbors

    if (procgrid[dim] == 1) {
      nrecv = nsend;
      buf = buf_send;

    } else {
      MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,
                   &nrecv1,1,MPI_INT,procneigh[dim][1],0,world,&status);
      nrecv = nrecv1;
      if (procgrid[dim] > 2) {
        MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,
                     &nrecv2,1,MPI_INT,procneigh[dim][0],0,world,&status);
        nrecv += nrecv2;
      }
      if (nrecv > maxrecv) grow_recv(nrecv);

      MPI_Irecv(buf_recv,nrecv1,MPI_DOUBLE,procneigh[dim][1],0,
                world,&request);
      MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][0],0,world);
      MPI_Wait(&request,&status);

      if (procgrid[dim] > 2) {
        MPI_Irecv(&buf_recv[nrecv1],nrecv2,MPI_DOUBLE,procneigh[dim][0],0,
                  world,&request);
        MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][1],0,world);
        MPI_Wait(&request,&status);
      }

      buf = buf_recv;
    }

    // check incoming atoms to see if they are in my box
    // if so, add to my list

    m = 0;
    while (m < nrecv) {
      value = buf[m+dim+1];
      if (value >= lo && value < hi) m += avec->unpack_exchange(&buf[m]);
      else m += static_cast<int> (buf[m]);
    }
  }

  if (atom->firstgroupname) atom->first_reorder();
}

/* ----------------------------------------------------------------------
   borders: list nearby atoms to send to neighboring procs at every timestep
   one list is created for every swap that will be made
   as list is made, actually do swaps
   this does equivalent of a communicate, so don't need to explicitly
     call communicate routine on reneighboring timestep
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
------------------------------------------------------------------------- */

void Comm::borders()
{
  int i,n,itype,iswap,dim,ineed,twoneed,smax,rmax;
  int nsend,nrecv,sendflag,nfirst,nlast,ngroup;
  double lo,hi;
  int *type;
  double **x;
  double *buf,*mlo,*mhi;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;

  // do swaps over all 3 dimensions

  iswap = 0;
  smax = rmax = 0;

  for (dim = 0; dim < 3; dim++) {
    nlast = 0;
    twoneed = 2*maxneed[dim];
    for (ineed = 0; ineed < twoneed; ineed++) {

      // find atoms within slab boundaries lo/hi using <= and >=
      // check atoms between nfirst and nlast
      //   for first swaps in a dim, check owned and ghost
      //   for later swaps in a dim, only check newly arrived ghosts
      // store sent atom indices in list for use in future timesteps

      x = atom->x;
      if (style == SINGLE) {
        lo = slablo[iswap];
        hi = slabhi[iswap];
      } else {
        type = atom->type;
        mlo = multilo[iswap];
        mhi = multihi[iswap];
      }
      if (ineed % 2 == 0) {
        nfirst = nlast;
        nlast = atom->nlocal + atom->nghost;
      }

      nsend = 0;

      // sendflag = 0 if I do not send on this swap
      // sendneed test indicates receiver no longer requires data
      // e.g. due to non-PBC or non-uniform sub-domains

      if (ineed/2 >= sendneed[dim][ineed % 2]) sendflag = 0;
      else sendflag = 1;

      // find send atoms according to SINGLE vs MULTI
      // all atoms eligible versus atoms in bordergroup
      // only need to limit loop to bordergroup for first sends (ineed < 2)
      // on these sends, break loop in two: owned (in group) and ghost

      if (sendflag) {
        if (!bordergroup || ineed >= 2) {
          if (style == SINGLE) {
            for (i = nfirst; i < nlast; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
          } else {
            for (i = nfirst; i < nlast; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          }

        } else {
          if (style == SINGLE) {
            ngroup = atom->nfirst;
            for (i = 0; i < ngroup; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            for (i = atom->nlocal; i < nlast; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
          } else {
            ngroup = atom->nfirst;
            for (i = 0; i < ngroup; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
            for (i = atom->nlocal; i < nlast; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          }
        }
      }

      // pack up list of border atoms

      if (nsend*size_border > maxsend) grow_send(nsend*size_border,0);
      if (ghost_velocity)
        n = avec->pack_border_vel(nsend,sendlist[iswap],buf_send,
                                  pbc_flag[iswap],pbc[iswap]);
      else
        n = avec->pack_border(nsend,sendlist[iswap],buf_send,
                              pbc_flag[iswap],pbc[iswap]);

      // swap atoms with other proc
      // no MPI calls except SendRecv if nsend/nrecv = 0
      // put incoming ghosts at end of my atom arrays
      // if swapping with self, simply copy, no messages

      if (sendproc[iswap] != me) {
        MPI_Sendrecv(&nsend,1,MPI_INT,sendproc[iswap],0,
                     &nrecv,1,MPI_INT,recvproc[iswap],0,world,&status);
        if (nrecv*size_border > maxrecv) grow_recv(nrecv*size_border);
        if (nrecv) MPI_Irecv(buf_recv,nrecv*size_border,MPI_DOUBLE,
                             recvproc[iswap],0,world,&request);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        if (nrecv) MPI_Wait(&request,&status);
        buf = buf_recv;
      } else {
        nrecv = nsend;
        buf = buf_send;
      }

      // unpack buffer

      if (ghost_velocity)
        avec->unpack_border_vel(nrecv,atom->nlocal+atom->nghost,buf);
      else
        avec->unpack_border(nrecv,atom->nlocal+atom->nghost,buf);

      // set all pointers & counters

      smax = MAX(smax,nsend);
      rmax = MAX(rmax,nrecv);
      sendnum[iswap] = nsend;
      recvnum[iswap] = nrecv;
      size_forward_recv[iswap] = nrecv*size_forward;
      size_reverse_send[iswap] = nrecv*size_reverse;
      size_reverse_recv[iswap] = nsend*size_reverse;
      firstrecv[iswap] = atom->nlocal + atom->nghost;
      atom->nghost += nrecv;
      iswap++;
    }
  }

  // insure send/recv buffers are long enough for all forward & reverse comm

  int max = MAX(maxforward*smax,maxreverse*rmax);
  if (max > maxsend) grow_send(max,0);
  max = MAX(maxforward*rmax,maxreverse*smax);
  if (max > maxrecv) grow_recv(max);

  // reset global->local map

  if (map_style) atom->map_set();
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Pair
------------------------------------------------------------------------- */

void Comm::forward_comm_pair(Pair *pair)
{
  int iswap,n;
  double *buf;
  MPI_Request request;
  MPI_Status status;

  for (iswap = 0; iswap < nswap; iswap++) {

    // pack buffer

    n = pair->pack_comm(sendnum[iswap],sendlist[iswap],
                        buf_send,pbc_flag[iswap],pbc[iswap]);

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      if (recvnum[iswap])
        MPI_Irecv(buf_recv,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
                  world,&request);
      if (sendnum[iswap])
        MPI_Send(buf_send,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,world);
      if (recvnum[iswap]) MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    pair->unpack_comm(recvnum[iswap],firstrecv[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Pair
------------------------------------------------------------------------- */

void Comm::reverse_comm_pair(Pair *pair)
{
  int iswap,n;
  double *buf;
  MPI_Request request;
  MPI_Status status;

  for (iswap = nswap-1; iswap >= 0; iswap--) {

    // pack buffer

    n = pair->pack_reverse_comm(recvnum[iswap],firstrecv[iswap],buf_send);

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      if (sendnum[iswap])
        MPI_Irecv(buf_recv,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
                  world,&request);
      if (recvnum[iswap])
        MPI_Send(buf_send,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,world);
      if (sendnum[iswap]) MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    pair->unpack_reverse_comm(sendnum[iswap],sendlist[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void Comm::forward_comm_fix(Fix *fix)
{
  int iswap,n;
  double *buf;
  MPI_Request request;
  MPI_Status status;

  for (iswap = 0; iswap < nswap; iswap++) {

    // pack buffer

    n = fix->pack_comm(sendnum[iswap],sendlist[iswap],
                       buf_send,pbc_flag[iswap],pbc[iswap]);

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      if (recvnum[iswap])
        MPI_Irecv(buf_recv,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
                  world,&request);
      if (sendnum[iswap])
        MPI_Send(buf_send,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,world);
      if (recvnum[iswap]) MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    fix->unpack_comm(recvnum[iswap],firstrecv[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void Comm::reverse_comm_fix(Fix *fix)
{
  int iswap,n;
  double *buf;
  MPI_Request request;
  MPI_Status status;

  for (iswap = nswap-1; iswap >= 0; iswap--) {

    // pack buffer

    n = fix->pack_reverse_comm(recvnum[iswap],firstrecv[iswap],buf_send);

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      if (sendnum[iswap])
        MPI_Irecv(buf_recv,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
                  world,&request);
      if (recvnum[iswap])
        MPI_Send(buf_send,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,world);
      if (sendnum[iswap]) MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    fix->unpack_reverse_comm(sendnum[iswap],sendlist[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   n = total datums for all atoms, allows for variable number/atom
------------------------------------------------------------------------- */

void Comm::forward_comm_variable_fix(Fix *fix)
{
  int iswap,n;
  double *buf;
  MPI_Request request;
  MPI_Status status;

  for (iswap = 0; iswap < nswap; iswap++) {

    // pack buffer

    n = fix->pack_comm(sendnum[iswap],sendlist[iswap],
                       buf_send,pbc_flag[iswap],pbc[iswap]);

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      if (recvnum[iswap])
        MPI_Irecv(buf_recv,maxrecv,MPI_DOUBLE,recvproc[iswap],0,
                  world,&request);
      if (sendnum[iswap])
        MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
      if (recvnum[iswap]) MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    fix->unpack_comm(recvnum[iswap],firstrecv[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   n = total datums for all atoms, allows for variable number/atom
------------------------------------------------------------------------- */

void Comm::reverse_comm_variable_fix(Fix *fix)
{
  int iswap,n;
  double *buf;
  MPI_Request request;
  MPI_Status status;

  for (iswap = nswap-1; iswap >= 0; iswap--) {

    // pack buffer

    n = fix->pack_reverse_comm(recvnum[iswap],firstrecv[iswap],buf_send);

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      if (sendnum[iswap])
        MPI_Irecv(buf_recv,maxrecv,MPI_DOUBLE,sendproc[iswap],0,
                  world,&request);
      if (recvnum[iswap])
        MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,world);
      if (sendnum[iswap]) MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    fix->unpack_reverse_comm(sendnum[iswap],sendlist[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Compute
------------------------------------------------------------------------- */

void Comm::forward_comm_compute(Compute *compute)
{
  int iswap,n;
  double *buf;
  MPI_Request request;
  MPI_Status status;

  for (iswap = 0; iswap < nswap; iswap++) {

    // pack buffer

    n = compute->pack_comm(sendnum[iswap],sendlist[iswap],
                           buf_send,pbc_flag[iswap],pbc[iswap]);

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      if (recvnum[iswap])
        MPI_Irecv(buf_recv,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
                  world,&request);
      if (sendnum[iswap])
        MPI_Send(buf_send,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,world);
      if (recvnum[iswap]) MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    compute->unpack_comm(recvnum[iswap],firstrecv[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Compute
------------------------------------------------------------------------- */

void Comm::reverse_comm_compute(Compute *compute)
{
  int iswap,n;
  double *buf;
  MPI_Request request;
  MPI_Status status;

  for (iswap = nswap-1; iswap >= 0; iswap--) {

    // pack buffer

    n = compute->pack_reverse_comm(recvnum[iswap],firstrecv[iswap],buf_send);

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      if (sendnum[iswap])
        MPI_Irecv(buf_recv,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
                  world,&request);
      if (recvnum[iswap])
        MPI_Send(buf_send,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,world);
      if (sendnum[iswap]) MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    compute->unpack_reverse_comm(sendnum[iswap],sendlist[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Dump
------------------------------------------------------------------------- */

void Comm::forward_comm_dump(Dump *dump)
{
  int iswap,n;
  double *buf;
  MPI_Request request;
  MPI_Status status;

  for (iswap = 0; iswap < nswap; iswap++) {

    // pack buffer

    n = dump->pack_comm(sendnum[iswap],sendlist[iswap],
                        buf_send,pbc_flag[iswap],pbc[iswap]);

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      if (recvnum[iswap])
        MPI_Irecv(buf_recv,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
                  world,&request);
      if (sendnum[iswap])
        MPI_Send(buf_send,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,world);
      if (recvnum[iswap]) MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    dump->unpack_comm(recvnum[iswap],firstrecv[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Dump
------------------------------------------------------------------------- */

void Comm::reverse_comm_dump(Dump *dump)
{
  int iswap,n;
  double *buf;
  MPI_Request request;
  MPI_Status status;

  for (iswap = nswap-1; iswap >= 0; iswap--) {

    // pack buffer

    n = dump->pack_reverse_comm(recvnum[iswap],firstrecv[iswap],buf_send);

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      if (sendnum[iswap])
        MPI_Irecv(buf_recv,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
                  world,&request);
      if (recvnum[iswap])
        MPI_Send(buf_send,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,world);
      if (sendnum[iswap]) MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    dump->unpack_reverse_comm(sendnum[iswap],sendlist[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   forward communication of N values in array
------------------------------------------------------------------------- */

void Comm::forward_comm_array(int n, double **array)
{
  int i,j,k,m,iswap,last;
  double *buf;
  MPI_Request request;
  MPI_Status status;

  // NOTE: check that buf_send and buf_recv are big enough

  for (iswap = 0; iswap < nswap; iswap++) {

    // pack buffer

    m = 0;
    for (i = 0; i < sendnum[iswap]; i++) {
      j = sendlist[iswap][i];
      for (k = 0; k < n; k++)
        buf_send[m++] = array[j][k];
    }

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      if (recvnum[iswap])
        MPI_Irecv(buf_recv,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
                  world,&request);
      if (sendnum[iswap])
        MPI_Send(buf_send,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,world);
      if (recvnum[iswap]) MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    m = 0;
    last = firstrecv[iswap] + recvnum[iswap];
    for (i = firstrecv[iswap]; i < last; i++)
      for (k = 0; k < n; k++)
        array[i][k] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   communicate inbuf around full ring of processors with messtag
   nbytes = size of inbuf = n datums * nper bytes
   callback() is invoked to allow caller to process/update each proc's inbuf
   if self=1 (default), then callback() is invoked on final iteration
     using original inbuf, which may have been updated
   for non-NULL outbuf, final updated inbuf is copied to it
     outbuf = inbuf is OK
------------------------------------------------------------------------- */

void Comm::ring(int n, int nper, void *inbuf, int messtag,
                void (*callback)(int, char *), void *outbuf, int self)
{
  MPI_Request request;
  MPI_Status status;

  int nbytes = n*nper;
  int maxbytes;
  MPI_Allreduce(&nbytes,&maxbytes,1,MPI_INT,MPI_MAX,world);

  char *buf,*bufcopy;
  memory->create(buf,maxbytes,"comm:buf");
  memory->create(bufcopy,maxbytes,"comm:bufcopy");
  memcpy(buf,inbuf,nbytes);

  int next = me + 1;
  int prev = me - 1;
  if (next == nprocs) next = 0;
  if (prev < 0) prev = nprocs - 1;

  for (int loop = 0; loop < nprocs; loop++) {
    if (me != next) {
      MPI_Irecv(bufcopy,maxbytes,MPI_CHAR,prev,messtag,world,&request);
      MPI_Send(buf,nbytes,MPI_CHAR,next,messtag,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_CHAR,&nbytes);
      memcpy(buf,bufcopy,nbytes);
    }
    if (self || loop < nprocs-1) callback(nbytes/nper,buf);
  }

  if (outbuf) memcpy(outbuf,buf,nbytes);

  memory->destroy(buf);
  memory->destroy(bufcopy);
}

/* ----------------------------------------------------------------------
   proc 0 reads Nlines from file into buf and bcasts buf to all procs
   caller allocates buf to max size needed
   each line is terminated by newline, even if last line in file is not
   return 0 if successful, 1 if get EOF error before read is complete
------------------------------------------------------------------------- */

int Comm::read_lines_from_file(FILE *fp, int nlines, int maxline, char *buf)
{
  int m;

  if (me == 0) {
    m = 0;
    for (int i = 0; i < nlines; i++) {
      if (!fgets(&buf[m],maxline,fp)) {
	m = 0;
	break;
      }
      m += strlen(&buf[m]);
    }
    if (m) {
      if (buf[m-1] != '\n') strcpy(&buf[m++],"\n");
      m++;
    }
  }

  MPI_Bcast(&m,1,MPI_INT,0,world);
  if (m == 0) return 1;
  MPI_Bcast(buf,m,MPI_CHAR,0,world);
  return 0;
}

/* ----------------------------------------------------------------------
   proc 0 reads Nlines from file into buf and bcasts buf to all procs
   caller allocates buf to max size needed
   each line is terminated by newline, even if last line in file is not
   return 0 if successful, 1 if get EOF error before read is complete
------------------------------------------------------------------------- */

int Comm::read_lines_from_file_universe(FILE *fp, int nlines, int maxline,
                                        char *buf)
{
  int m;

  int me_universe = universe->me;
  MPI_Comm uworld = universe->uworld;

  if (me_universe == 0) {
    m = 0;
    for (int i = 0; i < nlines; i++) {
      if (!fgets(&buf[m],maxline,fp)) {
	m = 0;
	break;
      }
      m += strlen(&buf[m]);
    }
    if (m) {
      if (buf[m-1] != '\n') strcpy(&buf[m++],"\n");
      m++;
    }
  }

  MPI_Bcast(&m,1,MPI_INT,0,uworld);
  if (m == 0) return 1;
  MPI_Bcast(buf,m,MPI_CHAR,0,uworld);
  return 0;
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR and bufextra
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void Comm::grow_send(int n, int flag)
{
  maxsend = static_cast<int> (BUFFACTOR * n);
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

void Comm::grow_recv(int n)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  memory->destroy(buf_recv);
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
------------------------------------------------------------------------- */

void Comm::grow_list(int iswap, int n)
{
  maxsendlist[iswap] = static_cast<int> (BUFFACTOR * n);
  memory->grow(sendlist[iswap],maxsendlist[iswap],"comm:sendlist[iswap]");
}

/* ----------------------------------------------------------------------
   realloc the buffers needed for swaps
------------------------------------------------------------------------- */

void Comm::grow_swap(int n)
{
  free_swap();
  allocate_swap(n);
  if (style == MULTI) {
    free_multi();
    allocate_multi(n);
  }

  sendlist = (int **)
    memory->srealloc(sendlist,n*sizeof(int *),"comm:sendlist");
  memory->grow(maxsendlist,n,"comm:maxsendlist");
  for (int i = maxswap; i < n; i++) {
    maxsendlist[i] = BUFMIN;
    memory->create(sendlist[i],BUFMIN,"comm:sendlist[i]");
  }
  maxswap = n;
}

/* ----------------------------------------------------------------------
   allocation of swap info
------------------------------------------------------------------------- */

void Comm::allocate_swap(int n)
{
  memory->create(sendnum,n,"comm:sendnum");
  memory->create(recvnum,n,"comm:recvnum");
  memory->create(sendproc,n,"comm:sendproc");
  memory->create(recvproc,n,"comm:recvproc");
  memory->create(size_forward_recv,n,"comm:size");
  memory->create(size_reverse_send,n,"comm:size");
  memory->create(size_reverse_recv,n,"comm:size");
  memory->create(slablo,n,"comm:slablo");
  memory->create(slabhi,n,"comm:slabhi");
  memory->create(firstrecv,n,"comm:firstrecv");
  memory->create(pbc_flag,n,"comm:pbc_flag");
  memory->create(pbc,n,6,"comm:pbc");
}

/* ----------------------------------------------------------------------
   allocation of multi-type swap info
------------------------------------------------------------------------- */

void Comm::allocate_multi(int n)
{
  multilo = memory->create(multilo,n,atom->ntypes+1,"comm:multilo");
  multihi = memory->create(multihi,n,atom->ntypes+1,"comm:multihi");
}

/* ----------------------------------------------------------------------
   free memory for swaps
------------------------------------------------------------------------- */

void Comm::free_swap()
{
  memory->destroy(sendnum);
  memory->destroy(recvnum);
  memory->destroy(sendproc);
  memory->destroy(recvproc);
  memory->destroy(size_forward_recv);
  memory->destroy(size_reverse_send);
  memory->destroy(size_reverse_recv);
  memory->destroy(slablo);
  memory->destroy(slabhi);
  memory->destroy(firstrecv);
  memory->destroy(pbc_flag);
  memory->destroy(pbc);
}

/* ----------------------------------------------------------------------
   free memory for multi-type swaps
------------------------------------------------------------------------- */

void Comm::free_multi()
{
  memory->destroy(multilo);
  memory->destroy(multihi);
}

/* ----------------------------------------------------------------------
   set communication style
   invoked from input script by communicate command
------------------------------------------------------------------------- */

void Comm::set(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal communicate command");

  if (strcmp(arg[0],"single") == 0) style = SINGLE;
  else if (strcmp(arg[0],"multi") == 0) style = MULTI;
  else error->all(FLERR,"Illegal communicate command");

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal communicate command");
      bordergroup = group->find(arg[iarg+1]);
      if (bordergroup < 0)
        error->all(FLERR,"Invalid group in communicate command");
      if (bordergroup && (atom->firstgroupname == NULL ||
                          strcmp(arg[iarg+1],atom->firstgroupname) != 0))
        error->all(FLERR,"Communicate group != atom_modify first group");
      iarg += 2;
    } else if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal communicate command");
      cutghostuser = force->numeric(FLERR,arg[iarg+1]);
      if (cutghostuser < 0.0)
        error->all(FLERR,"Invalid cutoff in communicate command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"vel") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal communicate command");
      if (strcmp(arg[iarg+1],"yes") == 0) ghost_velocity = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) ghost_velocity = 0;
      else error->all(FLERR,"Illegal communicate command");
      iarg += 2;
    } else error->all(FLERR,"Illegal communicate command");
  }
}

/* ----------------------------------------------------------------------
   set dimensions for 3d grid of processors, and associated flags
   invoked from input script by processors command
------------------------------------------------------------------------- */

void Comm::set_processors(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal processors command");

  if (strcmp(arg[0],"*") == 0) user_procgrid[0] = 0;
  else user_procgrid[0] = force->inumeric(FLERR,arg[0]);
  if (strcmp(arg[1],"*") == 0) user_procgrid[1] = 0;
  else user_procgrid[1] = force->inumeric(FLERR,arg[1]);
  if (strcmp(arg[2],"*") == 0) user_procgrid[2] = 0;
  else user_procgrid[2] = force->inumeric(FLERR,arg[2]);

  if (user_procgrid[0] < 0 || user_procgrid[1] < 0 || user_procgrid[2] < 0)
    error->all(FLERR,"Illegal processors command");

  int p = user_procgrid[0]*user_procgrid[1]*user_procgrid[2];
  if (p && p != nprocs)
    error->all(FLERR,"Specified processors != physical processors");

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"grid") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal processors command");

      if (strcmp(arg[iarg+1],"onelevel") == 0) {
        gridflag = ONELEVEL;

      } else if (strcmp(arg[iarg+1],"twolevel") == 0) {
        if (iarg+6 > narg) error->all(FLERR,"Illegal processors command");
        gridflag = TWOLEVEL;

        ncores = force->inumeric(FLERR,arg[iarg+2]);
        if (strcmp(arg[iarg+3],"*") == 0) user_coregrid[0] = 0;
        else user_coregrid[0] = force->inumeric(FLERR,arg[iarg+3]);
        if (strcmp(arg[iarg+4],"*") == 0) user_coregrid[1] = 0;
        else user_coregrid[1] = force->inumeric(FLERR,arg[iarg+4]);
        if (strcmp(arg[iarg+5],"*") == 0) user_coregrid[2] = 0;
        else user_coregrid[2] = force->inumeric(FLERR,arg[iarg+5]);

        if (ncores <= 0 || user_coregrid[0] < 0 ||
            user_coregrid[1] < 0 || user_coregrid[2] < 0)
          error->all(FLERR,"Illegal processors command");
        iarg += 4;

      } else if (strcmp(arg[iarg+1],"numa") == 0) {
        gridflag = NUMA;

      } else if (strcmp(arg[iarg],"custom") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal processors command");
        gridflag = CUSTOM;
        delete [] customfile;
        int n = strlen(arg[iarg+2]) + 1;
        customfile = new char[n];
        strcpy(customfile,arg[iarg+2]);
        iarg += 1;

      } else error->all(FLERR,"Illegal processors command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"map") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal processors command");
      if (strcmp(arg[iarg+1],"cart") == 0) mapflag = CART;
      else if (strcmp(arg[iarg+1],"cart/reorder") == 0) mapflag = CARTREORDER;
      else if (strcmp(arg[iarg+1],"xyz") == 0 ||
               strcmp(arg[iarg+1],"xzy") == 0 ||
               strcmp(arg[iarg+1],"yxz") == 0 ||
               strcmp(arg[iarg+1],"yzx") == 0 ||
               strcmp(arg[iarg+1],"zxy") == 0 ||
               strcmp(arg[iarg+1],"zyx") == 0) {
        mapflag = XYZ;
        strcpy(xyz,arg[iarg+1]);
      } else error->all(FLERR,"Illegal processors command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"part") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal processors command");
      if (universe->nworlds == 1)
        error->all(FLERR,
                   "Cannot use processors part command "
                   "without using partitions");
      int isend = force->inumeric(FLERR,arg[iarg+1]);
      int irecv = force->inumeric(FLERR,arg[iarg+2]);
      if (isend < 1 || isend > universe->nworlds ||
          irecv < 1 || irecv > universe->nworlds || isend == irecv)
        error->all(FLERR,"Invalid partitions in processors part command");
      if (isend-1 == universe->iworld) {
        if (send_to_partition >= 0)
          error->all(FLERR,
                     "Sending partition in processors part command "
                     "is already a sender");
        send_to_partition = irecv-1;
      }
      if (irecv-1 == universe->iworld) {
        if (recv_from_partition >= 0)
          error->all(FLERR,
                     "Receiving partition in processors part command "
                     "is already a receiver");
        recv_from_partition = isend-1;
      }

      // only receiver has otherflag dependency

      if (strcmp(arg[iarg+3],"multiple") == 0) {
        if (universe->iworld == irecv-1) {
          otherflag = 1;
          other_style = MULTIPLE;
        }
      } else error->all(FLERR,"Illegal processors command");
      iarg += 4;

    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal processors command");
      delete [] outfile;
      int n = strlen(arg[iarg+1]) + 1;
      outfile = new char[n];
      strcpy(outfile,arg[iarg+1]);
      iarg += 2;

    } else error->all(FLERR,"Illegal processors command");
  }

  // error checks

  if (gridflag == NUMA && mapflag != CART)
    error->all(FLERR,"Processors grid numa and map style are incompatible");
  if (otherflag && (gridflag == NUMA || gridflag == CUSTOM))
    error->all(FLERR,
               "Processors part option and grid style are incompatible");
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint Comm::memory_usage()
{
  bigint bytes = 0;
  bytes += nprocs * sizeof(int);    // grid2proc
  for (int i = 0; i < nswap; i++)
    bytes += memory->usage(sendlist[i],maxsendlist[i]);
  bytes += memory->usage(buf_send,maxsend+bufextra);
  bytes += memory->usage(buf_recv,maxrecv);
  return bytes;
}
