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

#include "mpi.h"
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

/* ----------------------------------------------------------------------
   setup MPI and allocate buffer space 
------------------------------------------------------------------------- */

Comm::Comm(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  user_procgrid[0] = user_procgrid[1] = user_procgrid[2] = 0;

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

/* ----------------------------------------------------------------------
   destructor, free all memory 
------------------------------------------------------------------------- */

Comm::~Comm()
{
  free_swap();
  memory->sfree(maxsendlist);
  if (sendlist) for (int i = 0; i < maxswap; i++) memory->sfree(sendlist[i]);
  memory->sfree(sendlist);

  memory->sfree(buf_send);
  memory->sfree(buf_recv);
}

/* ---------------------------------------------------------------------- */

void Comm::init()
{
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
  if (force->dimension == 2 && procgrid[2] != 1)
    error->all("Proc grid in z != 1 for 2d simulation");

  int reorder = 0;
  int periods[3];
  periods[0] = periods[1] = periods[2] = 1;
  MPI_Comm cartesian;
      
  MPI_Cart_create(world,3,procgrid,periods,reorder,&cartesian);
  MPI_Cart_get(cartesian,3,procgrid,periods,myloc);
  MPI_Cart_shift(cartesian,0,1,&procneigh[0][0],&procneigh[0][1]);
  MPI_Cart_shift(cartesian,1,1,&procneigh[1][0],&procneigh[1][1]);
  MPI_Cart_shift(cartesian,2,1,&procneigh[2][0],&procneigh[2][1]);
  MPI_Comm_free(&cartesian);

  if (me == 0) {
    if (screen) fprintf(screen,"  %d by %d by %d processor grid\n",
			procgrid[0],procgrid[1],procgrid[2]);
    if (logfile) fprintf(logfile,"  %d by %d by %d processor grid\n",
			 procgrid[0],procgrid[1],procgrid[2]);
  }
}

/* ----------------------------------------------------------------------
   setup spatial-decomposition communication patterns 
   function of neighbor cutoff and current box size 
------------------------------------------------------------------------- */

void Comm::setup()
{
  // need = # of procs I need atoms from in each dim

  need[0] = static_cast<int> 
    (neighbor->cutneigh * procgrid[0] / domain->xprd) + 1;
  need[1] = static_cast<int> 
    (neighbor->cutneigh * procgrid[1] / domain->yprd) + 1;
  need[2] = static_cast<int> 
    (neighbor->cutneigh * procgrid[2] / domain->zprd) + 1;

  // for 2d, don't communicate in z

  if (force->dimension == 2) need[2] = 0;

  // if non-periodic, do not communicate further than procgrid-1 away
  // this enables very large cutoffs in non-periodic systems

  if (domain->xperiodic == 0) need[0] = MIN(need[0],procgrid[0]-1);
  if (domain->yperiodic == 0) need[1] = MIN(need[1],procgrid[1]-1);
  if (domain->zperiodic == 0) need[2] = MIN(need[2],procgrid[2]-1);

  // allocate comm memory

  nswap = 2 * (need[0]+need[1]+need[2]);
  if (nswap > maxswap) grow_swap();

  // setup parameters for each exchange:
  // sendproc = proc to send to at each swap
  // recvproc = proc to recv from at each swap
  // slablo/slabhi = boundaries for slab of atoms to send at each swap
  //   use -BIG/midpt/BIG to insure all atoms included even if round-off occurs
  //   if round-off, atoms recvd across PBC can be < or > than subbox boundary
  //   note that borders() only loops over subset of atoms during each swap
  // set slablo > slabhi for swaps across non-periodic boundaries
  //   this insures no atoms are swapped
  //   only for procs owning sub-box at non-periodic end of global box
  // pbc_flags = add-on factor for atoms sent across a periodic global boundary
  //    0 = not across a boundary
  //    1 = add box-length to coord when sending
  //   -1 = subtract box-length from coord when sending
  // 1st part of if statement is sending to the west/south/down
  // 2nd part of if statement is sending to the east/north/up

  int iswap = 0;
  for (int dim = 0; dim < 3; dim++) {
    for (int ineed = 0; ineed < 2*need[dim]; ineed++) {
      pbc_flags[iswap][0] = 0;
      pbc_flags[iswap][1] = 0;
      pbc_flags[iswap][2] = 0;
      pbc_flags[iswap][3] = 0;

      if (ineed % 2 == 0) {
	sendproc[iswap] = procneigh[dim][0];
	recvproc[iswap] = procneigh[dim][1];
	if (ineed < 2) slablo[iswap] = -BIG;
	else slablo[iswap] = 0.5 * (domain->sublo[dim] + domain->subhi[dim]);
	slabhi[iswap] = domain->sublo[dim] + neighbor->cutneigh;
	if (myloc[dim] == 0) {
	  if (domain->periodicity[dim] == 0)
	    slabhi[iswap] = slablo[iswap] - 1.0;
	  else {
	    pbc_flags[iswap][0] = 1;
	    pbc_flags[iswap][dim+1] = 1;
	  }
	}
      } else {
	sendproc[iswap] = procneigh[dim][1];
	recvproc[iswap] = procneigh[dim][0];
	slablo[iswap] = domain->subhi[dim] - neighbor->cutneigh;
	if (ineed < 2) slabhi[iswap] = BIG;
	else slabhi[iswap] = 0.5 * (domain->sublo[dim] + domain->subhi[dim]);
	if (myloc[dim] == procgrid[dim]-1) {
	  if (domain->periodicity[dim] == 0)
	    slabhi[iswap] = slablo[iswap] - 1.0;
	  else {
	    pbc_flags[iswap][0] = 1;
	    pbc_flags[iswap][dim+1] = -1;
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
			    buf_send,pbc_flags[iswap]);
	MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
	MPI_Wait(&request,&status);
      } else {
	MPI_Irecv(buf_recv,size_comm_recv[iswap],MPI_DOUBLE,
		  recvproc[iswap],0,world,&request);
	n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
			    buf_send,pbc_flags[iswap]);
	MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
	MPI_Wait(&request,&status);
	avec->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_recv);
      }

    } else {
      if (comm_x_only) {
	if (sendnum[iswap])
	  n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
			      x[firstrecv[iswap]],pbc_flags[iswap]);
      } else {
	n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
			    buf_send,pbc_flags[iswap]);
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
   exchange:
   move atoms to correct processors
   atoms exchanged with all 6 stencil neighbors
   send out atoms that have left my box, receive ones entering my box
   atoms will be lost if not inside some proc's box
     can happen if atom moves outside of non-periodic bounary
     or if atom moves more than one proc away
   this routine called before every reneighboring
------------------------------------------------------------------------- */

void Comm::exchange()
{
  int i,m,nsend,nrecv,nrecv1,nrecv2,nlocal;
  double lo,hi,value;
  double **x;
  double *buf;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;

  // clear global->local map since atoms move & new ghosts are created

  if (map_style) atom->map_clear();

  // loop over dimensions

  for (int dim = 0; dim < 3; dim++) {

    // fill buffer with atoms leaving my box, using < and >=
    // when atom is deleted, fill it in with last atom

    x = atom->x;
    lo = domain->sublo[dim];
    hi = domain->subhi[dim];
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
}

/* ----------------------------------------------------------------------
   borders:
   make lists of nearby atoms to send to neighboring procs at every timestep
   one list is created for every swap that will be made
   as list is made, actually do swaps
   this does equivalent of a communicate (so don't need to explicitly
     call communicate routine on reneighboring timestep)
   this routine is called before every reneighboring
------------------------------------------------------------------------- */

void Comm::borders()
{
  int i,n,iswap,dim,ineed,maxneed,nsend,nrecv,nfirst,nlast,smax,rmax;
  double lo,hi;
  double **x;
  double *buf;
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

      lo = slablo[iswap];
      hi = slabhi[iswap];
      x = atom->x;
      if (ineed % 2 == 0) {
	nfirst = nlast;
	nlast = atom->nlocal + atom->nghost;
      }

      nsend = 0;
      for (i = nfirst; i < nlast; i++)
	if (x[i][dim] >= lo && x[i][dim] <= hi) {
	  if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
	  sendlist[iswap][nsend++] = i;
	}
      
      // pack up list of border atoms

      if (nsend*size_border > maxsend)
	grow_send(nsend*size_border,0);
      n = avec->pack_border(nsend,sendlist[iswap],buf_send,pbc_flags[iswap]);

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
			buf_send,pbc_flags[iswap]);

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
		       buf_send,pbc_flags[iswap]);

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
			   buf_send,pbc_flags[iswap]);

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
   assign nprocs to 3d xprd,yprd,zprd box so as to minimize surface area 
------------------------------------------------------------------------- */

void Comm::procs2box()
{
  int ipx,ipy,ipz,nremain;
  double boxx,boxy,boxz,surf;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double bestsurf = 2.0 * (xprd*yprd + yprd*zprd + zprd*xprd);

  // loop thru all possible factorizations of nprocs
  // surf = surface area of a proc sub-domain
  // for 2d, insure ipz = 1

  ipx = 1;
  while (ipx <= nprocs) {
    if (nprocs % ipx == 0) {
      nremain = nprocs/ipx;
      ipy = 1;
      while (ipy <= nremain) {
	if (nremain % ipy == 0) {
	  ipz = nremain/ipy;
	  if (force->dimension == 3 || ipz == 1) {
	    boxx = xprd/ipx;
	    boxy = yprd/ipy;
	    boxz = zprd/ipz;
	    surf = boxx*boxy + boxy*boxz + boxz*boxx;
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

void Comm::grow_swap()
{
  free_swap();
  allocate_swap(nswap);

  sendlist = (int **) memory->srealloc(sendlist,nswap*sizeof(int *),
				       "comm:sendlist");
  maxsendlist = (int *) memory->srealloc(maxsendlist,nswap*sizeof(int),
					 "comm:maxsendlist");
  for (int i = maxswap; i < nswap; i++) {
    maxsendlist[i] = BUFMIN;
    sendlist[i] = (int *) memory->smalloc(BUFMIN*sizeof(int),
					  "comm:sendlist[i]");
  }
  maxswap = nswap;
}

/* ----------------------------------------------------------------------
   initial allocation of swap info 
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
  pbc_flags = (int **) memory->create_2d_int_array(n,4,"comm:pbc_flags");
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
  memory->destroy_2d_int_array(pbc_flags);
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory 
------------------------------------------------------------------------- */

int Comm::memory_usage()
{
  int bytes = 0;

  for (int i = 0; i < nswap; i++) bytes += maxsendlist[i] * sizeof(int);
  bytes += maxsend * sizeof(double);
  bytes += maxrecv * sizeof(double);

  return bytes;
}
