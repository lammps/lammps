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
#include "math_extra.h"
#include "error.h"
#include "memory.h"

#include <map>
#include <string>

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace LAMMPS_NS;

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000
#define BIG 1.0e20

enum{SINGLE,MULTI};
enum{NONE,MULTIPLE};
enum{CART,CARTREORDER,XYZ,XZY,YXZ,YZX,ZXY,ZYX};

/* ----------------------------------------------------------------------
   setup MPI and allocate buffer space 
------------------------------------------------------------------------- */

Comm::Comm(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  user_procgrid[0] = user_procgrid[1] = user_procgrid[2] = 0;
  grid2proc = NULL;
  layoutflag = CART;
  numa_nodes = 0;

  bordergroup = 0;
  style = SINGLE;
  multilo = multihi = NULL;
  cutghostmulti = NULL;
  cutghostuser = 0.0;
  ghost_velocity = 0;
  recv_from_partition = send_to_partition = -1;
  other_partition_style = NONE;

  // use of OpenMP threads
  // query OpenMP for number of threads/process set by user at run-time
  // need to be in a parallel area for this operation

  nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel default(shared)
  {
#pragma omp master
    { nthreads = omp_get_num_threads(); }
  }
  if (me == 0) {
    if (screen)
      fprintf(screen,"  using %d OpenMP thread(s) per MPI task\n",nthreads);
    if (logfile)
      fprintf(logfile,"  using %d OpenMP thread(s) per MPI task\n",nthreads);
  }
#endif

  // initialize comm buffers & exchange memory

  maxsend = BUFMIN;
  memory->create(buf_send,maxsend+BUFEXTRA,"comm:buf_send");
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
  if (grid2proc) memory->destroy(grid2proc);

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
   create 3d grid of procs based on box size
------------------------------------------------------------------------- */

void Comm::set_proc_grid()
{
  // recv proc layout of another partition if my layout depends on it

  if (recv_from_partition >= 0) {
    MPI_Status status;
    if (me == 0) MPI_Recv(other_procgrid,3,MPI_INT,
			  universe->root_proc[recv_from_partition],0,
			  universe->uworld,&status);
    MPI_Bcast(other_procgrid,3,MPI_INT,0,world);
  }

  // use NUMA routines if numa_nodes > 0
  // if NUMA routines fail, just continue

  if (numa_nodes)
    if (numa_set_proc_grid()) return;

  // create layout of procs mapped to simulation box
  // can fail (on one partition) if constrained by other_partition_style

  int flag = 
    procs2box(nprocs,user_procgrid,procgrid,1,1,1,other_partition_style);
  if (!flag) error->all(FLERR,"Could not layout grid of processors",1);

  // error check

  if (procgrid[0]*procgrid[1]*procgrid[2] != nprocs)
    error->all(FLERR,"Bad grid of processors");
  if (domain->dimension == 2 && procgrid[2] != 1)
    error->all(FLERR,"Processor count in z must be 1 for 2d simulation");

  // grid2proc[i][j][k] = proc that owns i,j,k location in grid

  if (grid2proc) memory->destroy(grid2proc);
  memory->create(grid2proc,procgrid[0],procgrid[1],procgrid[2],
		 "comm:grid2proc");

  // use MPI Cartesian routines to layout 3d grid of procs
  // MPI may do layout in machine-optimized fashion

  if (layoutflag == CART || layoutflag == CARTREORDER) {
    int reorder = 0;
    if (layoutflag == CARTREORDER) reorder = 1;
    int periods[3];
    periods[0] = periods[1] = periods[2] = 1;
    MPI_Comm cartesian;
      
    MPI_Cart_create(world,3,procgrid,periods,reorder,&cartesian);
    MPI_Cart_get(cartesian,3,procgrid,periods,myloc);
    MPI_Cart_shift(cartesian,0,1,&procneigh[0][0],&procneigh[0][1]);
    MPI_Cart_shift(cartesian,1,1,&procneigh[1][0],&procneigh[1][1]);
    MPI_Cart_shift(cartesian,2,1,&procneigh[2][0],&procneigh[2][1]);

    int coords[3];
    int i,j,k;
    for (i = 0; i < procgrid[0]; i++)
      for (j = 0; j < procgrid[1]; j++)
	for (k = 0; k < procgrid[2]; k++) {
	  coords[0] = i; coords[1] = j; coords[2] = k;
	  MPI_Cart_rank(cartesian,coords,&grid2proc[i][j][k]);
	}
    
    MPI_Comm_free(&cartesian);

  // layout 3d grid of procs explicitly, via user-requested XYZ ordering

  } else {
    int i,j,k,ime,jme,kme;
    for (i = 0; i < procgrid[0]; i++)
      for (j = 0; j < procgrid[1]; j++)
	for (k = 0; k < procgrid[2]; k++) {
	  if (layoutflag == XYZ)
	    grid2proc[i][j][k] = k*procgrid[1]*procgrid[0] + j*procgrid[0] + i;
	  else if (layoutflag == XZY)
	    grid2proc[i][j][k] = j*procgrid[2]*procgrid[0] + k*procgrid[0] + i;
	  else if (layoutflag == YXZ)
	    grid2proc[i][j][k] = k*procgrid[0]*procgrid[1] + i*procgrid[1] + j;
	  else if (layoutflag == YZX)
	    grid2proc[i][j][k] = i*procgrid[2]*procgrid[1] + k*procgrid[1] + j;
	  else if (layoutflag == ZXY)
	    grid2proc[i][j][k] = j*procgrid[0]*procgrid[2] + i*procgrid[2] + k;
	  else if (layoutflag == ZYX)
	    grid2proc[i][j][k] = i*procgrid[1]*procgrid[2] + j*procgrid[2] + k;

	  if (grid2proc[i][j][k] == me) {
	    ime = i; jme = j, kme = k;
	  }
	}

    myloc[0] = ime;
    myloc[1] = jme;
    myloc[2] = kme;

    i = ime-1;
    if (i < 0) i = procgrid[0]-1;
    procneigh[0][0] = grid2proc[i][jme][kme];
    i = ime+1;
    if (i == procgrid[0]) i = 0;
    procneigh[0][1] = grid2proc[i][jme][kme];

    j = jme-1;
    if (j < 0) j = procgrid[1]-1;
    procneigh[1][0] = grid2proc[ime][j][kme];
    j = jme+1;
    if (j == procgrid[1]) j = 0;
    procneigh[1][1] = grid2proc[ime][j][kme];

    k = kme-1;
    if (k < 0) k = procgrid[2]-1;
    procneigh[2][0] = grid2proc[ime][jme][k];
    k = kme+1;
    if (k == procgrid[2]) k = 0;
    procneigh[2][1] = grid2proc[ime][jme][k];
  }

  // set lamda box params after procs are assigned

  if (domain->triclinic) domain->set_lamda_box();

  if (me == 0) {
    if (screen) fprintf(screen,"  %d by %d by %d MPI processor grid\n",
			procgrid[0],procgrid[1],procgrid[2]);
    if (logfile) fprintf(logfile,"  %d by %d by %d MPI processor grid\n",
			 procgrid[0],procgrid[1],procgrid[2]);
  }

  // send my proc layout to another partition if requested

  if (send_to_partition >= 0) {
    if (me == 0) MPI_Send(procgrid,3,MPI_INT,
			  universe->root_proc[send_to_partition],0,
			  universe->uworld);
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
  // augment by velocity quantities if needed

  size_forward = atom->avec->size_forward;
  size_reverse = atom->avec->size_reverse;
  size_border = atom->avec->size_border;

  if (ghost_velocity) size_forward += atom->avec->size_velocity;
  if (ghost_velocity) size_border += atom->avec->size_velocity;

  // maxforward = # of datums in largest forward communication
  // maxreverse = # of datums in largest reverse communication
  // query pair,fix,compute,dump for their requirements

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

  // need = # of procs I need atoms from in each dim based on max cutoff
  // for 2d, don't communicate in z

  need[0] = static_cast<int> (cutghost[0] * procgrid[0] / prd[0]) + 1;
  need[1] = static_cast<int> (cutghost[1] * procgrid[1] / prd[1]) + 1;
  need[2] = static_cast<int> (cutghost[2] * procgrid[2] / prd[2]) + 1;
  if (domain->dimension == 2) need[2] = 0;

  // if non-periodic, do not communicate further than procgrid-1 away
  // this enables very large cutoffs in non-periodic systems

  int *periodicity = domain->periodicity;
  if (periodicity[0] == 0) need[0] = MIN(need[0],procgrid[0]-1);
  if (periodicity[1] == 0) need[1] = MIN(need[1],procgrid[1]-1);
  if (periodicity[2] == 0) need[2] = MIN(need[2],procgrid[2]-1);

  // allocate comm memory

  nswap = 2 * (need[0]+need[1]+need[2]);
  if (nswap > maxswap) grow_swap(nswap);

  // setup parameters for each exchange:
  // sendproc = proc to send to at each swap
  // recvproc = proc to recv from at each swap
  // for style SINGLE:
  //   slablo/slabhi = boundaries for slab of atoms to send at each swap
  //   use -BIG/midpt/BIG to insure all atoms included even if round-off occurs
  //   if round-off, atoms recvd across PBC can be < or > than subbox boundary
  //   note that borders() only loops over subset of atoms during each swap
  //   set slablo > slabhi for swaps across non-periodic boundaries
  //     this insures no atoms are swapped
  //     only for procs owning sub-box at non-periodic end of global box
  // for style MULTI:
  //   multilo/multihi is same as slablo/slabhi, only for each atom type
  // pbc_flag: 0 = nothing across a boundary, 1 = something across a boundary
  // pbc = -1/0/1 for PBC factor in each of 3/6 orthog/triclinic dirs
  // for triclinic, slablo/hi and pbc_border will be used in lamda (0-1) coords
  // 1st part of if statement is sending to the west/south/down
  // 2nd part of if statement is sending to the east/north/up

  int dim,ineed;

  int iswap = 0;
  for (dim = 0; dim < 3; dim++) {
    for (ineed = 0; ineed < 2*need[dim]; ineed++) {
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
	  if (periodicity[dim] == 0) {
	    if (style == SINGLE) slabhi[iswap] = slablo[iswap] - 1.0;
	    else 
	      for (i = 1; i <= ntypes; i++)
		multihi[iswap][i] = multilo[iswap][i] - 1.0;
	  } else {
	    pbc_flag[iswap] = 1;
	    pbc[iswap][dim] = 1;
	    if (triclinic) {
	      if (dim == 1) pbc[iswap][5] = 1;
	      else if (dim == 2) pbc[iswap][4] = pbc[iswap][3] = 1;
	    }
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
	  if (periodicity[dim] == 0) {
	    if (style == SINGLE) slabhi[iswap] = slablo[iswap] - 1.0;
	    else
	      for (i = 1; i <= ntypes; i++)
		multihi[iswap][i] = multilo[iswap][i] - 1.0;
	  } else {
	    pbc_flag[iswap] = 1;
	    pbc[iswap][dim] = -1;
	    if (triclinic) {
	      if (dim == 1) pbc[iswap][5] = -1;
	      else if (dim == 2) pbc[iswap][4] = pbc[iswap][3] = -1;
	    }
	  }
	}
      }
      
      iswap++;
    }
  }
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
	MPI_Irecv(buf,size_forward_recv[iswap],MPI_DOUBLE,
		  recvproc[iswap],0,world,&request);
	n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
			    buf_send,pbc_flag[iswap],pbc[iswap]);
	MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
	MPI_Wait(&request,&status);
      } else if (ghost_velocity) {
	MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,
		  recvproc[iswap],0,world,&request);
	n = avec->pack_comm_vel(sendnum[iswap],sendlist[iswap],
				buf_send,pbc_flag[iswap],pbc[iswap]);
	MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
	MPI_Wait(&request,&status);
	avec->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],buf_recv);
      } else {
	MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,
		  recvproc[iswap],0,world,&request);
	n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
			    buf_send,pbc_flag[iswap],pbc[iswap]);
	MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
	MPI_Wait(&request,&status);
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
	MPI_Irecv(buf_recv,size_reverse_recv[iswap],MPI_DOUBLE,
		  sendproc[iswap],0,world,&request);
	if (size_reverse_send[iswap]) buf = f[firstrecv[iswap]];
	else buf = NULL;
	MPI_Send(buf,size_reverse_send[iswap],MPI_DOUBLE,
		 recvproc[iswap],0,world);
	MPI_Wait(&request,&status);
      } else {
	MPI_Irecv(buf_recv,size_reverse_recv[iswap],MPI_DOUBLE,
		  sendproc[iswap],0,world,&request);
	n = avec->pack_reverse(recvnum[iswap],firstrecv[iswap],buf_send);
	MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,world);
	MPI_Wait(&request,&status);
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
  // new ghosts are created in borders()
  // map_set() is done at end of borders()

  if (map_style) atom->map_clear();

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
   this does equivalent of a communicate (so don't need to explicitly
     call communicate routine on reneighboring timestep)
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
------------------------------------------------------------------------- */

void Comm::borders()
{
  int i,n,itype,iswap,dim,ineed,maxneed,smax,rmax;
  int nsend,nrecv,nfirst,nlast,ngroup;
  double lo,hi;
  int *type;
  double **x;
  double *buf,*mlo,*mhi;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;

  // clear old ghosts and any ghost bonus data internal to AtomVec

  atom->nghost = 0;
  atom->avec->clear_bonus();

  // do swaps over all 3 dimensions

  iswap = 0;
  smax = rmax = 0;

  for (dim = 0; dim < 3; dim++) {
    nlast = 0;
    maxneed = 2*need[dim];
    for (ineed = 0; ineed < maxneed; ineed++) {

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

      // find send atoms according to SINGLE vs MULTI
      // all atoms eligible versus atoms in bordergroup
      // only need to limit loop to bordergroup for first sends (ineed < 2)
      // on these sends, break loop in two: owned (in group) and ghost

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

      // pack up list of border atoms

      if (nsend*size_border > maxsend)
	grow_send(nsend*size_border,0);
      if (ghost_velocity)
	n = avec->pack_border_vel(nsend,sendlist[iswap],buf_send,
				  pbc_flag[iswap],pbc[iswap]);
      else
	n = avec->pack_border(nsend,sendlist[iswap],buf_send,
			      pbc_flag[iswap],pbc[iswap]);
      
      // swap atoms with other proc
      // put incoming ghosts at end of my atom arrays
      // if swapping with self, simply copy, no messages

      if (sendproc[iswap] != me) {
	MPI_Sendrecv(&nsend,1,MPI_INT,sendproc[iswap],0,
		     &nrecv,1,MPI_INT,recvproc[iswap],0,world,&status);
	if (nrecv*size_border > maxrecv) 
	  grow_recv(nrecv*size_border);
	MPI_Irecv(buf_recv,nrecv*size_border,MPI_DOUBLE,
		  recvproc[iswap],0,world,&request);
	MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
	MPI_Wait(&request,&status);
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
      MPI_Irecv(buf_recv,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
		world,&request);
      MPI_Send(buf_send,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,world);
      MPI_Wait(&request,&status);
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
      MPI_Irecv(buf_recv,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
		world,&request);
      MPI_Send(buf_send,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,world);
      MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    pair->unpack_reverse_comm(sendnum[iswap],sendlist[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
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
      MPI_Irecv(buf_recv,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
		world,&request);
      MPI_Send(buf_send,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,world);
      MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    fix->unpack_comm(recvnum[iswap],firstrecv[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
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
      MPI_Irecv(buf_recv,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
		world,&request);
      MPI_Send(buf_send,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,world);
      MPI_Wait(&request,&status);
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
      MPI_Irecv(buf_recv,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
		world,&request);
      MPI_Send(buf_send,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,world);
      MPI_Wait(&request,&status);
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
      MPI_Irecv(buf_recv,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
		world,&request);
      MPI_Send(buf_send,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,world);
      MPI_Wait(&request,&status);
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
      MPI_Irecv(buf_recv,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
		world,&request);
      MPI_Send(buf_send,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,world);
      MPI_Wait(&request,&status);
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
      MPI_Irecv(buf_recv,n*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
		world,&request);
      MPI_Send(buf_send,n*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,world);
      MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    dump->unpack_reverse_comm(sendnum[iswap],sendlist[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   assign numprocs to 3d box so as to minimize surface area
   area = surface area of each of 3 faces of simulation box divided by sx,sy,sz
   for triclinic, area = cross product of 2 edge vectors stored in h matrix
   valid assignment will be factorization of numprocs = Px by Py by Pz
   user_factors = if non-zero, factor is specified by user
   sx,sy,sz = scale box xyz dimension vy dividing by sx,sy,sz
   return factors = # of procs assigned to each dimension
   return 1 if successully factor, 0 if not
------------------------------------------------------------------------- */

int Comm::procs2box(int numprocs, int user_factors[3], int factors[3], 
		    const int sx, const int sy, const int sz, int other)
{
  factors[0] = user_factors[0];
  factors[1] = user_factors[1];
  factors[2] = user_factors[2];

  // all 3 proc counts are specified

  if (factors[0] && factors[1] && factors[2]) return 1;

  // 2 out of 3 proc counts are specified

  if (factors[0] > 0 && factors[1] > 0) {
    factors[2] = nprocs/(factors[0]*factors[1]);
    return 1;
  } else if (factors[0] > 0 && factors[2] > 0) {
    factors[1] = nprocs/(factors[0]*factors[2]);
    return 1;
  } else if (factors[1] > 0 && factors[2] > 0) {
    factors[0] = nprocs/(factors[1]*factors[2]);
    return 1;
  } 

  // determine cross-sectional areas for orthogonal and triclinic boxes
  // area[0] = xy, area[1] = xz, area[2] = yz

  double area[3];
  if (domain->triclinic == 0) {
    area[0] = domain->xprd * domain->yprd / (sx * sy);
    area[1] = domain->xprd * domain->zprd / (sx * sz);
    area[2] = domain->yprd * domain->zprd / (sy * sz);
  } else {
    double *h = domain->h;
    double a[3],b[3],c[3];
    a[0] = h[0]; a[1] = 0.0; a[2] = 0.0;
    b[0] = h[5]; b[1] = h[1]; b[2] = 0.0;
    MathExtra::cross3(a,b,c);
    area[0] = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) / (sx * sy);
    a[0] = h[0]; a[1] = 0.0; a[2] = 0.0;
    b[0] = h[4]; b[1] = h[3]; b[2] = h[2];
    MathExtra::cross3(a,b,c);
    area[1] = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) / (sx * sz);
    a[0] = h[5]; a[1] = h[1]; a[2] = 0.0;
    b[0] = h[4]; b[1] = h[3]; b[2] = h[2];
    MathExtra::cross3(a,b,c);
    area[2] = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) / (sy * sz);
  }

  double bestsurf = 2.0 * (area[0]+area[1]+area[2]);

  // loop thru all possible factorizations of numprocs
  // only consider valid cases that match procgrid settings
  // surf = surface area of a proc sub-domain
  // only consider cases that match user_factors & other_procgrid settings
  // success = 1 if valid factoriztion is found
  // may not be if other constraint is enforced

  int ipx,ipy,ipz,valid;
  double surf;
  
  int success = 0;
  ipx = 1;
  while (ipx <= numprocs) {
    valid = 1;
    if (user_factors[0] && ipx != user_factors[0]) valid = 0;
    if (other == MULTIPLE && other_procgrid[0] % ipx) valid = 0;
    if (numprocs % ipx) valid = 0;

    if (!valid) {
      ipx++;
      continue;
    }

    ipy = 1;
    while (ipy <= numprocs/ipx) {
      valid = 1;
      if (user_factors[1] && ipy != user_factors[1]) valid = 0;
      if (other == MULTIPLE && other_procgrid[1] % ipy) valid = 0;
      if ((numprocs/ipx) % ipy) valid = 0;
      if (!valid) {
	ipy++;
	continue;
      }
      
      ipz = numprocs/ipx/ipy;
      valid = 1;
      if (user_factors[2] && ipz != user_factors[2]) valid = 0;
      if (other == MULTIPLE && other_procgrid[2] % ipz) valid = 0;
      if (domain->dimension == 2 && ipz != 1) valid = 0;
      if (!valid) {
	ipy++;
	continue;
      }
      
      surf = area[0]/ipx/ipy + area[1]/ipx/ipz + area[2]/ipy/ipz;
      if (surf < bestsurf) {
	success = 1;
	bestsurf = surf;
	factors[0] = ipx;
	factors[1] = ipy;
	factors[2] = ipz;
      }
      ipy++;
    }

    ipx++;
  }

  return success;
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA 
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void Comm::grow_send(int n, int flag)
{
  maxsend = static_cast<int> (BUFFACTOR * n);
  if (flag)
    memory->grow(buf_send,(maxsend+BUFEXTRA),"comm:buf_send");
  else {
    memory->destroy(buf_send);
    memory->create(buf_send,maxsend+BUFEXTRA,"comm:buf_send");
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
      cutghostuser = atof(arg[iarg+1]);
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
  else user_procgrid[0] = atoi(arg[0]);
  if (strcmp(arg[1],"*") == 0) user_procgrid[1] = 0;
  else user_procgrid[1] = atoi(arg[1]);
  if (strcmp(arg[2],"*") == 0) user_procgrid[2] = 0;
  else user_procgrid[2] = atoi(arg[2]);

  if (user_procgrid[0] < 0 || user_procgrid[1] < 0 || user_procgrid[2] < 0) 
    error->all(FLERR,"Illegal processors command");

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"part") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal processors command");
      if (universe->nworlds == 1)
	error->all(FLERR,
		   "Cannot use processors part command "
		   "without using partitions");
      int isend = atoi(arg[iarg+1]);
      int irecv = atoi(arg[iarg+2]);
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
      if (strcmp(arg[iarg+3],"multiple") == 0) {
	if (universe->iworld == irecv-1) other_partition_style = MULTIPLE;
      } else error->all(FLERR,"Illegal processors command");
      iarg += 4;

    } else if (strcmp(arg[iarg],"grid") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal processors command");
      if (strcmp(arg[iarg+1],"cart") == 0) layoutflag = CART;
      else if (strcmp(arg[iarg+1],"cart/reorder") == 0)
	layoutflag = CARTREORDER;
      else if (strcmp(arg[iarg+1],"xyz") == 0) layoutflag = XYZ;
      else if (strcmp(arg[iarg+1],"xzy") == 0) layoutflag = XZY;
      else if (strcmp(arg[iarg+1],"yxz") == 0) layoutflag = YXZ;
      else if (strcmp(arg[iarg+1],"yzx") == 0) layoutflag = YZX;
      else if (strcmp(arg[iarg+1],"zxy") == 0) layoutflag = ZXY;
      else if (strcmp(arg[iarg+1],"zyx") == 0) layoutflag = ZYX;
      else error->all(FLERR,"Illegal processors command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"numa") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal processors command");
      numa_nodes = atoi(arg[iarg+1]);
      if (numa_nodes < 0) error->all(FLERR,"Illegal processors command");
      iarg += 2;

    } else error->all(FLERR,"Illegal processors command");
  }
}

/* ----------------------------------------------------------------------
   setup 3d grid of procs based on box size, group neighbors by numa node
   return 1 if successful, return 0 if not
------------------------------------------------------------------------- */

int Comm::numa_set_proc_grid()
{
  // get names of all nodes

  int name_length;
  char node_name[MPI_MAX_PROCESSOR_NAME];
  char node_names[MPI_MAX_PROCESSOR_NAME*nprocs];
  MPI_Get_processor_name(node_name,&name_length);
  MPI_Allgather(&node_name,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,&node_names,
                MPI_MAX_PROCESSOR_NAME,MPI_CHAR,world);
  std::string node_string = std::string(node_name);
  
  // get number of procs per node

  std::map<std::string,int> name_map;
  std::map<std::string,int>::iterator np;
  for (int i = 0; i < nprocs; i++) {
    std::string i_string = std::string(&node_names[i*MPI_MAX_PROCESSOR_NAME]);
    np = name_map.find(i_string);
    if (np == name_map.end())
      name_map[i_string] = 1;
    else
      np->second++;
  }
  int procs_per_node = name_map.begin()->second;
  int procs_per_numa = procs_per_node / numa_nodes;
  
  // use regular mapping if:
  
  if (procs_per_numa < 4 ||               // 3 or less procs per numa node
      procs_per_node % numa_nodes != 0 || // Different # of procs per numa node
      nprocs % procs_per_numa != 0 ||     // Different # of procs per numa node
      nprocs <= procs_per_numa ||         // Only 1 numa node used
      user_procgrid[0] > 1 ||             // User specified grid dimension
      user_procgrid[1] > 1 ||             //    that is greater than 1
      user_procgrid[2] > 1) {             //    in any dimension
    if (me == 0) {
      if (screen) fprintf(screen,"  1 by 1 by 1 Node grid\n");
      if (logfile) fprintf(logfile,"  1 by 1 by 1 Node grid\n");
    }
    return 0;
  }
  
  // user settings for the factorization per numa node - currently always zero

  int user_numagrid[3];
  user_numagrid[0] = user_numagrid[1] = user_numagrid[2] = 0;
  
  // if user specifies 1 for a proc grid dimension,
  // also use 1 for the node grid dimension

  if (user_procgrid[0] == 1) user_numagrid[0] = 1;
  if (user_procgrid[1] == 1) user_numagrid[1] = 1;
  if (user_procgrid[2] == 1) user_numagrid[2] = 1;
  
  // get an initial factorization for each numa node,
  // if the user has not set the number of processors
  // can fail (on one partition) if constrained by other_partition_style

  int numagrid[3];
  int flag = procs2box(procs_per_numa,user_numagrid,numagrid,
		       1,1,1,other_partition_style);
  if (!flag) error->all(FLERR,"Could not layout grid of processors",1);
  if (numagrid[0] * numagrid[1] * numagrid[2] != procs_per_numa)
    error->all(FLERR,"Bad node grid of processors");
  
  // get a factorization for the grid of numa nodes
  // should never fail

  int node_count = nprocs / procs_per_numa;
  procs2box(node_count,user_procgrid,procgrid,
	    numagrid[0],numagrid[1],numagrid[2],0);
  if (procgrid[0] * procgrid[1] * procgrid[2] != node_count)
    error->all(FLERR,"Bad grid of processors");
  
  // repeat the numa node factorization using the subdomain sizes
  // this will refine the factorization if the user specified the node layout
  // should never fail

  procs2box(procs_per_numa,user_numagrid,numagrid,
	    procgrid[0],procgrid[1],procgrid[2],0);
  if (numagrid[0]*numagrid[1]*numagrid[2] != procs_per_numa)
    error->all(FLERR,"Bad Node grid of processors");
  if (domain->dimension == 2 && (procgrid[2] != 1 || numagrid[2] != 1))
    error->all(FLERR,"Processor count in z must be 1 for 2d simulation");
  
  // assign a unique id to each node

  int node_num = 0, node_id = 0;
  for (np = name_map.begin(); np != name_map.end(); ++np) {
    if (np->first == node_string) node_id = node_num;
    node_num++;
  }
  
  // setup a per node communicator and find rank within

  MPI_Comm node_comm;
  MPI_Comm_split(world,node_id,0,&node_comm);  
  int node_rank;
  MPI_Comm_rank(node_comm,&node_rank);
  
  // setup a per numa communicator and find rank within

  MPI_Comm numa_comm;
  int local_numa = node_rank / procs_per_numa;
  MPI_Comm_split(node_comm,local_numa,0,&numa_comm);     
  int numa_rank;
  MPI_Comm_rank(numa_comm,&numa_rank);
  
  // setup a communicator with the rank 0 procs from each numa node

  MPI_Comm numa_leaders;
  MPI_Comm_split(world,numa_rank,0,&numa_leaders);
  
  // use the MPI Cartesian routines to map the nodes to the grid
  // could implement layoutflag as in non-NUMA case

  int reorder = 0;
  int periods[3];
  periods[0] = periods[1] = periods[2] = 1;
  MPI_Comm cartesian;
  if (numa_rank == 0) {
    MPI_Cart_create(numa_leaders,3,procgrid,periods,reorder,&cartesian);
    MPI_Cart_get(cartesian,3,procgrid,periods,myloc);
  }
  
  // broadcast numa node location in grid to other procs in numa node

  MPI_Bcast(myloc,3,MPI_INT,0,numa_comm);
  
  // get storage for the process mapping

  if (grid2proc) memory->destroy(grid2proc);
  memory->create(grid2proc,procgrid[0]*numagrid[0],procgrid[1]*numagrid[1],
                 procgrid[2]*numagrid[2],"comm:grid2proc");
  
  // compute my location within the grid

  int z_offset = numa_rank / (numagrid[0] * numagrid[1]);
  int y_offset = (numa_rank % (numagrid[0] * numagrid[1]))/numagrid[0];
  int x_offset = numa_rank % numagrid[0];
  myloc[0] = myloc[0] * numagrid[0] + x_offset;
  myloc[1] = myloc[1] * numagrid[1] + y_offset;
  myloc[2] = myloc[2] * numagrid[2] + z_offset;
  procgrid[0] *= numagrid[0];
  procgrid[1] *= numagrid[1];
  procgrid[2] *= numagrid[2];
  
  // allgather of locations to fill grid2proc

  int **gridi;
  memory->create(gridi,nprocs,3,"comm:gridi");
  MPI_Allgather(&myloc,3,MPI_INT,gridi[0],3,MPI_INT,world);
  for (int i = 0; i < nprocs; i++)
    grid2proc[gridi[i][0]][gridi[i][1]][gridi[i][2]] = i;
  memory->destroy(gridi);
  
  // get my neighbors

  int minus, plus;
  for (int i = 0; i < 3; i++) {
    numa_shift(myloc[i],procgrid[i],minus,plus);
    procneigh[i][0] = minus;
    procneigh[i][1] = plus;
  }
  procneigh[0][0] = grid2proc[procneigh[0][0]][myloc[1]][myloc[2]];
  procneigh[0][1] = grid2proc[procneigh[0][1]][myloc[1]][myloc[2]];
  procneigh[1][0] = grid2proc[myloc[0]][procneigh[1][0]][myloc[2]];
  procneigh[1][1] = grid2proc[myloc[0]][procneigh[1][1]][myloc[2]];
  procneigh[2][0] = grid2proc[myloc[0]][myloc[1]][procneigh[2][0]];
  procneigh[2][1] = grid2proc[myloc[0]][myloc[1]][procneigh[2][1]];
  
  if (numa_rank == 0) MPI_Comm_free(&cartesian);
  MPI_Comm_free(&numa_leaders);
  MPI_Comm_free(&numa_comm);
  MPI_Comm_free(&node_comm);
  
  // set lamda box params after procs are assigned

  if (domain->triclinic) domain->set_lamda_box();
  
  if (me == 0) {
    if (screen) fprintf(screen,"  %d by %d by %d Node grid\n",
			numagrid[0],numagrid[1],numagrid[2]);
    if (logfile) fprintf(logfile,"  %d by %d by %d Node grid\n",
			 numagrid[0],numagrid[1],numagrid[2]);
    if (screen) fprintf(screen,"  %d by %d by %d processor grid\n",
			procgrid[0],procgrid[1],procgrid[2]);
    if (logfile) fprintf(logfile,"  %d by %d by %d processor grid\n",
			 procgrid[0],procgrid[1],procgrid[2]);
  }

  return 1;
}

/* ----------------------------------------------------------------------
   get the index to the neighboring processors in a dimension
------------------------------------------------------------------------- */

void Comm::numa_shift(int myloc, int num_procs, int &minus, int &plus) 
{
  minus = myloc - 1;
  if (minus < 0) minus = num_procs - 1;
  plus = myloc + 1;
  if (plus == num_procs) plus = 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory 
------------------------------------------------------------------------- */

bigint Comm::memory_usage()
{
  bigint bytes = 0;
  for (int i = 0; i < nswap; i++) 
    bytes += memory->usage(sendlist[i],maxsendlist[i]);
  bytes += memory->usage(buf_send,maxsend+BUFEXTRA);
  bytes += memory->usage(buf_recv,maxrecv);
  return bytes;
}
