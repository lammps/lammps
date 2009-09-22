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

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "comm.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "neighbor.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "compute.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define BIG 1.0e20

enum{SINGLE,MULTI};

/* ----------------------------------------------------------------------
   setup MPI and allocate buffer space 
------------------------------------------------------------------------- */

Comm::Comm(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  user_procgrid[0] = user_procgrid[1] = user_procgrid[2] = 0;
  grid2proc = NULL;

  bordergroup = 0;
  style = SINGLE;
  multilo = multihi = NULL;
  cutghostmulti = NULL;
  cutghostuser = 0.0;

  // initialize comm buffers & exchange memory

  maxsend = BUFMIN;
  buf_send = (double *) 
    memory->smalloc((maxsend+BUFEXTRA)*sizeof(double),"comm:buf_send");
  maxrecv = BUFMIN;
  buf_recv = (double *) 
    memory->smalloc(maxrecv*sizeof(double),"comm:buf_recv");

  maxswap = 6;
  allocate_swap(maxswap);

  sendlist = (int **) memory->smalloc(maxswap*sizeof(int *),"sendlist");
  maxsendlist = (int *) memory->smalloc(maxswap*sizeof(int),"maxsendlist");
  for (int i = 0; i < maxswap; i++) {
    maxsendlist[i] = BUFMIN;
    sendlist[i] = (int *) memory->smalloc(BUFMIN*sizeof(int),"sendlist[i]");
  }

  maxforward_fix = maxreverse_fix = 0;
  maxforward_pair = maxreverse_pair = 0;
}

/* ---------------------------------------------------------------------- */

Comm::~Comm()
{
  if (grid2proc) memory->destroy_3d_int_array(grid2proc);

  free_swap();
  if (style == MULTI) {
    free_multi();
    memory->destroy_2d_double_array(cutghostmulti);
  }

  memory->sfree(maxsendlist);
  if (sendlist) for (int i = 0; i < maxswap; i++) memory->sfree(sendlist[i]);
  memory->sfree(sendlist);

  memory->sfree(buf_send);
  memory->sfree(buf_recv);
}

/* ----------------------------------------------------------------------
   setup 3d grid of procs based on box size
------------------------------------------------------------------------- */

void Comm::set_procs()
{
  if (user_procgrid[0] == 0) procs2box();
  else {
    procgrid[0] = user_procgrid[0];
    procgrid[1] = user_procgrid[1];
    procgrid[2] = user_procgrid[2];
  }

  if (procgrid[0]*procgrid[1]*procgrid[2] != nprocs)
    error->all("Bad grid of processors");
  if (domain->dimension == 2 && procgrid[2] != 1)
    error->all("Proc grid in z != 1 for 2d simulation");

  if (grid2proc) memory->destroy_3d_int_array(grid2proc);
  grid2proc = memory->create_3d_int_array(procgrid[0],procgrid[1],procgrid[2],
					  "comm:grid2proc");

  // use MPI Cartesian routines to setup 3d grid of procs
  // grid2proc[i][j][k] = proc that owns i,j,k location in grid
  // let MPI compute it instead of LAMMPS in case it is machine optimized

  int reorder = 0;
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

  // set lamda box params after procs are assigned

  if (domain->triclinic) domain->set_lamda_box();

  if (me == 0) {
    if (screen) fprintf(screen,"  %d by %d by %d processor grid\n",
			procgrid[0],procgrid[1],procgrid[2]);
    if (logfile) fprintf(logfile,"  %d by %d by %d processor grid\n",
			 procgrid[0],procgrid[1],procgrid[2]);
  }
}

/* ---------------------------------------------------------------------- */

void Comm::init()
{
  triclinic = domain->triclinic;
  map_style = atom->map_style;

  // comm_only = 1 if only x,f are exchanged in forward/reverse comm

  comm_x_only = atom->avec->comm_x_only;
  comm_f_only = atom->avec->comm_f_only;

  // maxforward = # of datums in largest forward communication
  // maxreverse = # of datums in largest reverse communication
  // query pair,fix,compute for their requirements

  maxforward = MAX(atom->avec->size_comm,atom->avec->size_border);
  maxreverse = atom->avec->size_reverse;

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

  if (force->newton == 0) maxreverse = 0;

  // memory for multi-style communication

  if (style == MULTI && multilo == NULL) {
    allocate_multi(maxswap);
    cutghostmulti = 
      memory->create_2d_double_array(atom->ntypes+1,3,"comm:cutghostmulti");
  }
  if (style == SINGLE && multilo) {
    free_multi();
    memory->destroy_2d_double_array(cutghostmulti);
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
   communication of atom coords every timestep
   other stuff may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void Comm::communicate()
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
	if (size_comm_recv[iswap]) buf = x[firstrecv[iswap]];
	else buf = NULL;
	MPI_Irecv(buf,size_comm_recv[iswap],MPI_DOUBLE,
		  recvproc[iswap],0,world,&request);
	n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
			    buf_send,pbc_flag[iswap],pbc[iswap]);
	MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
	MPI_Wait(&request,&status);
      } else {
	MPI_Irecv(buf_recv,size_comm_recv[iswap],MPI_DOUBLE,
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
			      x[firstrecv[iswap]],pbc_flag[iswap],pbc[iswap]);
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
   other stuff can also be sent via pack/unpack routines
------------------------------------------------------------------------- */
      
void Comm::reverse_communicate()
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
	avec->copy(nlocal-1,i);
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
  int size_border = avec->size_border;

  // clear old ghosts

  atom->nghost = 0;

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

      avec->unpack_border(nrecv,atom->nlocal+atom->nghost,buf);

      // set all pointers & counters

      smax = MAX(smax,nsend);
      rmax = MAX(rmax,nrecv);
      sendnum[iswap] = nsend;
      recvnum[iswap] = nrecv;
      size_comm_recv[iswap] = nrecv * avec->size_comm;
      size_reverse_send[iswap] = nrecv * avec->size_reverse;
      size_reverse_recv[iswap] = nsend * avec->size_reverse;
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

void Comm::comm_pair(Pair *pair)
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

void Comm::comm_fix(Fix *fix)
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

void Comm::comm_compute(Compute *compute)
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
   communicate atoms to new owning procs via irregular communication
   can be used in place of exchange()
   unlike exchange(), allows atoms to have moved arbitrarily long distances
   first setup irregular comm pattern, then invoke it
   atoms must be remapped to be inside simulation box before this is called
   for triclinic:
     atoms must be in lamda coords (0-1) before this is called
------------------------------------------------------------------------- */

void Comm::irregular()
{
  // clear global->local map since atoms move to new procs
  // zero out ghosts so map_set() at end will operate only on local atoms
  // exchange() doesn't need to zero ghosts b/c borders()
  //   is called right after and it zeroes ghosts and calls map_set()

  if (map_style) atom->map_clear();
  atom->nghost = 0;

  // subbox bounds for orthogonal or triclinic

  double *sublo,*subhi;
  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // loop over atoms, flag any that are not in my sub-box
  // fill buffer with atoms leaving my box, using < and >=
  // assign which proc it belongs to via irregular_lookup()
  // if irregular_lookup() returns me, due to round-off
  //   in triclinic x2lamda(), then keep atom and don't send
  // when atom is deleted, fill it in with last atom

  AtomVec *avec = atom->avec;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  int nsend = 0;
  int nsendatom = 0;
  int *sizes = new int[nlocal];
  int *proclist = new int[nlocal];

  int i = 0;
  while (i < nlocal) {
    if (x[i][0] < sublo[0] || x[i][0] >= subhi[0] ||
	x[i][1] < sublo[1] || x[i][1] >= subhi[1] ||
	x[i][2] < sublo[2] || x[i][2] >= subhi[2]) {
      proclist[nsendatom] = irregular_lookup(x[i]);
      if (proclist[nsendatom] != me) {
	if (nsend > maxsend) grow_send(nsend,1);
	sizes[nsendatom] = avec->pack_exchange(i,&buf_send[nsend]);
	nsend += sizes[nsendatom];
	nsendatom++;
	avec->copy(nlocal-1,i);
	nlocal--;
      } else i++;
    } else i++;
  }
  atom->nlocal = nlocal;

  // create irregular communication plan, perform comm, destroy plan
  // returned nrecv = size of buffer needed for incoming atoms

  int nrecv;
  Plan *plan = irregular_create(nsendatom,sizes,proclist,&nrecv);
  if (nrecv > maxrecv) grow_recv(nrecv);
  irregular_perform(plan,buf_send,sizes,buf_recv);
  irregular_destroy(plan);

  delete [] sizes;
  delete [] proclist;

  // add received atoms to my list

  int m = 0;
  while (m < nrecv) m += avec->unpack_exchange(&buf_recv[m]);

  // reset global->local map

  if (map_style) atom->map_set();
}

/* ----------------------------------------------------------------------
   create an irregular communication plan
   n = # of atoms to send
   sizes = # of doubles for each atom
   proclist = proc to send each atom to (none to me)
   return nrecvsize = total # of doubles I will recv
------------------------------------------------------------------------- */

Comm::Plan *Comm::irregular_create(int n, int *sizes, int *proclist,
				   int *nrecvsize)
{
  int i;

  // allocate plan and work vectors

  Plan *plan = (struct Plan *) memory->smalloc(sizeof(Plan),"comm:plan");
  int *list = new int[nprocs];
  int *count = new int[nprocs];

  // nrecv = # of messages I receive

  for (i = 0; i < nprocs; i++) {
    list[i] = 0;
    count[i] = 1;
  }
  for (i = 0; i < n; i++) list[proclist[i]] = 1;

  int nrecv;
  MPI_Reduce_scatter(list,&nrecv,count,MPI_INT,MPI_SUM,world);

  // allocate receive arrays

  int *proc_recv = new int[nrecv];
  int *length_recv = new int[nrecv];
  MPI_Request *request = new MPI_Request[nrecv];
  MPI_Status *status = new MPI_Status[nrecv];

  // nsend = # of messages I send

  for (i = 0; i < nprocs; i++) list[i] = 0;
  for (i = 0; i < n; i++) list[proclist[i]] += sizes[i];

  int nsend = 0;
  for (i = 0; i < nprocs; i++)
    if (list[i]) nsend++;

  // allocate send arrays

  int *proc_send = new int[nsend];
  int *length_send = new int[nsend];
  int *num_send = new int[nsend];
  int *index_send = new int[n];
  int *offset_send = new int[n];

  // list still stores size of message for procs I send to
  // proc_send = procs I send to
  // length_send = total size of message I send to each proc
  // to balance pattern of send messages:
  //   each proc begins with iproc > me, continues until iproc = me
  // reset list to store which send message each proc corresponds to

  int iproc = me;
  int isend = 0;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (list[iproc] > 0) {
      proc_send[isend] = iproc;
      length_send[isend] = list[iproc];
      list[iproc] = isend;
      isend++;
    }
  }

  // num_send = # of datums I send to each proc

  for (i = 0; i < nsend; i++) num_send[i] = 0;
  for (i = 0; i < n; i++) {
    isend = list[proclist[i]];
    num_send[isend]++;
  }

  // count = offsets into n-length index_send for each proc I send to
  // index_send = list of which datums to send to each proc
  //   1st N1 values are datum indices for 1st proc,
  //   next N2 values are datum indices for 2nd proc, etc
  // offset_send = where each datum starts in send buffer

  count[0] = 0;
  for (i = 1; i < nsend; i++) count[i] = count[i-1] + num_send[i-1];

  for (i = 0; i < n; i++) {
    isend = list[proclist[i]];
    index_send[count[isend]++] = i;
    if (i) offset_send[i] = offset_send[i-1] + sizes[i-1];
    else offset_send[i] = 0;
  }

  // tell receivers how much data I send
  // sendmax = largest # of datums I send in a single message

  int sendmax = 0;
  for (i = 0; i < nsend; i++) {
    MPI_Send(&length_send[i],1,MPI_INT,proc_send[i],0,world);
    sendmax = MAX(sendmax,length_send[i]);
  }

  // receive incoming messages
  // proc_recv = procs I recv from
  // length_recv = total size of message each proc sends me
  // nrecvsize = total size of data I recv

  *nrecvsize = 0;
  for (i = 0; i < nrecv; i++) {
    MPI_Recv(&length_recv[i],1,MPI_INT,MPI_ANY_SOURCE,0,world,status);
    proc_recv[i] = status->MPI_SOURCE;
    *nrecvsize += length_recv[i];
  }

  // barrier to insure all MPI_ANY_SOURCE messages are received
  // else another proc could proceed to irregular_perform() and send to me

  MPI_Barrier(world);

  // free work vectors

  delete [] count;
  delete [] list;
    
  // initialize plan and return it

  plan->nsend = nsend;
  plan->nrecv = nrecv;
  plan->sendmax = sendmax;

  plan->proc_send = proc_send;
  plan->length_send = length_send;
  plan->num_send = num_send;
  plan->index_send = index_send;
  plan->offset_send = offset_send;
  plan->proc_recv = proc_recv;
  plan->length_recv = length_recv;
  plan->request = request;
  plan->status = status;

  return plan;
}

/* ----------------------------------------------------------------------
   perform irregular communication
   plan = previouly computed plan via irregular_create()
   sendbuf = list of atoms to send
   sizes = # of doubles for each atom
   recvbuf = received atoms
------------------------------------------------------------------------- */

void Comm::irregular_perform(Plan *plan, double *sendbuf, int *sizes,
			     double *recvbuf)
{
  int i,m,n,offset;

  // post all receives

  offset = 0;
  for (int irecv = 0; irecv < plan->nrecv; irecv++) {
    MPI_Irecv(&recvbuf[offset],plan->length_recv[irecv],MPI_DOUBLE,
	      plan->proc_recv[irecv],0,world,&plan->request[irecv]);
    offset += plan->length_recv[irecv];
  }

  // allocate buf for largest send

  double *buf = (double *) memory->smalloc(plan->sendmax*sizeof(double),
					   "comm::irregular");

  // send each message
  // pack buf with list of datums (datum = one atom)
  // m = index of datum in sendbuf

  n = 0;
  for (int isend = 0; isend < plan->nsend; isend++) {
    offset = 0;
    for (i = 0; i < plan->num_send[isend]; i++) {
      m = plan->index_send[n++];
      memcpy(&buf[offset],&sendbuf[plan->offset_send[m]],
	     sizes[m]*sizeof(double));
      offset += sizes[m];
    }
    MPI_Send(buf,plan->length_send[isend],MPI_DOUBLE,
	     plan->proc_send[isend],0,world);
  }       

  // free temporary send buffer

  memory->sfree(buf);

  // wait on all incoming messages

  if (plan->nrecv) MPI_Waitall(plan->nrecv,plan->request,plan->status);
}

/* ----------------------------------------------------------------------
   destroy an irregular communication plan
------------------------------------------------------------------------- */

void Comm::irregular_destroy(Plan *plan)
{
  delete [] plan->proc_send;
  delete [] plan->length_send;
  delete [] plan->num_send;
  delete [] plan->index_send;
  delete [] plan->offset_send;
  delete [] plan->proc_recv;
  delete [] plan->length_recv;
  delete [] plan->request;
  delete [] plan->status;
  memory->sfree(plan);
}

/* ----------------------------------------------------------------------
   determine which proc owns atom with x coord
   x will be in box (orthogonal) or lamda coords (triclinic)
------------------------------------------------------------------------- */

int Comm::irregular_lookup(double *x)
{
  int loc[3];
  if (triclinic == 0) {
    double *boxlo = domain->boxlo;
    double *boxhi = domain->boxhi;
    loc[0] = static_cast<int>
      (procgrid[0] * (x[0]-boxlo[0]) / (boxhi[0]-boxlo[0]));
    loc[1] = static_cast<int>
      (procgrid[1] * (x[1]-boxlo[1]) / (boxhi[1]-boxlo[1]));
    loc[2] = static_cast<int>
      (procgrid[2] * (x[2]-boxlo[2]) / (boxhi[2]-boxlo[2]));
  } else {
    loc[0] = static_cast<int> (procgrid[0] * x[0]);
    loc[1] = static_cast<int> (procgrid[1] * x[1]);
    loc[2] = static_cast<int> (procgrid[2] * x[2]);
  }

  if (loc[0] < 0) loc[0] = 0;
  if (loc[0] >= procgrid[0]) loc[0] = procgrid[0] - 1;
  if (loc[1] < 0) loc[1] = 0;
  if (loc[1] >= procgrid[1]) loc[1] = procgrid[1] - 1;
  if (loc[2] < 0) loc[2] = 0;
  if (loc[2] >= procgrid[2]) loc[2] = procgrid[2] - 1;

  return grid2proc[loc[0]][loc[1]][loc[2]];
}

/* ----------------------------------------------------------------------
   assign nprocs to 3d xprd,yprd,zprd box so as to minimize surface area 
   area = surface area of each of 3 faces of simulation box
   for triclinic, area = cross product of 2 edge vectors stored in h matrix
------------------------------------------------------------------------- */

void Comm::procs2box()
{
  double area[3];

  if (domain->triclinic == 0) {
    area[0] = domain->xprd * domain->yprd;
    area[1] = domain->xprd * domain->zprd;
    area[2] = domain->yprd * domain->zprd;
  } else {
    double *h = domain->h;
    double x,y,z;
    cross(h[0],0.0,0.0,h[5],h[1],0.0,x,y,z);
    area[0] = sqrt(x*x + y*y + z*z);
    cross(h[0],0.0,0.0,h[4],h[3],h[2],x,y,z);
    area[1] = sqrt(x*x + y*y + z*z);
    cross(h[5],h[1],0.0,h[4],h[3],h[2],x,y,z);
    area[2] = sqrt(x*x + y*y + z*z);
  }

  double bestsurf = 2.0 * (area[0]+area[1]+area[2]);

  // loop thru all possible factorizations of nprocs
  // surf = surface area of a proc sub-domain
  // for 2d, insure ipz = 1

  int ipx,ipy,ipz,nremain;
  double surf;

  ipx = 1;
  while (ipx <= nprocs) {
    if (nprocs % ipx == 0) {
      nremain = nprocs/ipx;
      ipy = 1;
      while (ipy <= nremain) {
	if (nremain % ipy == 0) {
	  ipz = nremain/ipy;
	  if (domain->dimension == 3 || ipz == 1) {
	    surf = area[0]/ipx/ipy + area[1]/ipx/ipz + area[2]/ipy/ipz;
	    if (surf < bestsurf) {
	      bestsurf = surf;
	      procgrid[0] = ipx;
	      procgrid[1] = ipy;
	      procgrid[2] = ipz;
	    }
	  }
	}
	ipy++;
      }
    }
    ipx++;
  }
}

/* ----------------------------------------------------------------------
   vector cross product: c = a x b
------------------------------------------------------------------------- */

void Comm::cross(double ax, double ay, double az, 
		 double bx, double by, double bz, 
		 double &cx, double &cy, double &cz)
{
  cx = ay*bz - az*by;
  cy = az*bx - ax*bz;
  cz = ax*by - ay*bx;
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
    buf_send = (double *) 
      memory->srealloc(buf_send,(maxsend+BUFEXTRA)*sizeof(double),
		       "comm:buf_send");
  else {
    memory->sfree(buf_send);
    buf_send = (double *) memory->smalloc((maxsend+BUFEXTRA)*sizeof(double),
					  "comm:buf_send");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR 
------------------------------------------------------------------------- */

void Comm::grow_recv(int n)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  memory->sfree(buf_recv);
  buf_recv = (double *) memory->smalloc(maxrecv*sizeof(double),
					"comm:buf_recv");
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR 
------------------------------------------------------------------------- */

void Comm::grow_list(int iswap, int n)
{
  maxsendlist[iswap] = static_cast<int> (BUFFACTOR * n);
  sendlist[iswap] = (int *) 
    memory->srealloc(sendlist[iswap],maxsendlist[iswap]*sizeof(int),
		     "comm:sendlist[iswap]");
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
  maxsendlist = (int *)
    memory->srealloc(maxsendlist,n*sizeof(int),"comm:maxsendlist");
  for (int i = maxswap; i < n; i++) {
    maxsendlist[i] = BUFMIN;
    sendlist[i] = (int *)
      memory->smalloc(BUFMIN*sizeof(int),"comm:sendlist[i]");
  }
  maxswap = n;
}

/* ----------------------------------------------------------------------
   allocation of swap info 
------------------------------------------------------------------------- */

void Comm::allocate_swap(int n)
{
  sendnum = (int *) memory->smalloc(n*sizeof(int),"comm:sendnum");
  recvnum = (int *) memory->smalloc(n*sizeof(int),"comm:recvnum");
  sendproc = (int *) memory->smalloc(n*sizeof(int),"comm:sendproc");
  recvproc = (int *) memory->smalloc(n*sizeof(int),"comm:recvproc");
  size_comm_recv = (int *) memory->smalloc(n*sizeof(int),"comm:size");
  size_reverse_send = (int *) memory->smalloc(n*sizeof(int),"comm:size");
  size_reverse_recv = (int *) memory->smalloc(n*sizeof(int),"comm:size");
  slablo = (double *) memory->smalloc(n*sizeof(double),"comm:slablo");
  slabhi = (double *) memory->smalloc(n*sizeof(double),"comm:slabhi");
  firstrecv = (int *) memory->smalloc(n*sizeof(int),"comm:firstrecv");
  pbc_flag = (int *) memory->smalloc(n*sizeof(int),"comm:pbc_flag");
  pbc = (int **) memory->create_2d_int_array(n,6,"comm:pbc");
}


/* ----------------------------------------------------------------------
   allocation of multi-type swap info
------------------------------------------------------------------------- */

void Comm::allocate_multi(int n)
{
  multilo = memory->create_2d_double_array(n,atom->ntypes+1,"comm:multilo");
  multihi = memory->create_2d_double_array(n,atom->ntypes+1,"comm:multihi");
}

/* ----------------------------------------------------------------------
   free memory for swaps 
------------------------------------------------------------------------- */

void Comm::free_swap()
{
  memory->sfree(sendnum);
  memory->sfree(recvnum);
  memory->sfree(sendproc);
  memory->sfree(recvproc);
  memory->sfree(size_comm_recv);
  memory->sfree(size_reverse_send);
  memory->sfree(size_reverse_recv);
  memory->sfree(slablo);
  memory->sfree(slabhi);
  memory->sfree(firstrecv);
  memory->sfree(pbc_flag);
  memory->destroy_2d_int_array(pbc);
}

/* ----------------------------------------------------------------------
   free memory for multi-type swaps
------------------------------------------------------------------------- */

void Comm::free_multi()
{
  memory->destroy_2d_double_array(multilo);
  memory->destroy_2d_double_array(multihi);
}

/* ----------------------------------------------------------------------
   set communication style
------------------------------------------------------------------------- */

void Comm::set(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal communicate command");

  if (strcmp(arg[0],"single") == 0) style = SINGLE;
  else if (strcmp(arg[0],"multi") == 0) style = MULTI;
  else error->all("Illegal communicate command");

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all("Illegal communicate command");
      bordergroup = group->find(arg[iarg+1]);
      if (bordergroup < 0)
	error->all("Invalid group in communicate command");
      if (bordergroup && (atom->firstgroupname == NULL || 
			  strcmp(arg[iarg+1],atom->firstgroupname) != 0))
	error->all("Communicate group != atom_modify first group");
      iarg += 2;
    } else if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+2 > narg) error->all("Illegal communicate command");
      cutghostuser = atof(arg[iarg+1]);
      if (cutghostuser < 0.0) 
	error->all("Invalid cutoff in communicate command");
      iarg += 2;
    } else error->all("Illegal communicate command");
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory 
------------------------------------------------------------------------- */

double Comm::memory_usage()
{
  double bytes = 0.0;

  for (int i = 0; i < nswap; i++) bytes += maxsendlist[i] * sizeof(int);
  bytes += maxsend * sizeof(double);
  bytes += maxrecv * sizeof(double);

  return bytes;
}
