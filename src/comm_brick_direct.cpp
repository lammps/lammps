// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "comm_brick_direct.h"

#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "neighbor.h"

// NOTES:
// allocate requests to length of nrecv_direct
// what do lengths of send/recv bufs need to be

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CommBrickDirect::CommBrickDirect(LAMMPS *lmp) : CommBrick(lmp)
{
  style = Comm::BRICK_DIRECT;
}

/* ---------------------------------------------------------------------- */

CommBrickDirect::~CommBrickDirect()
{
}

/* ---------------------------------------------------------------------- */

CommBrickDirect::CommBrickDirect(LAMMPS *lmp, Comm *oldcomm) : CommBrick(lmp, oldcomm)
{
  if (oldcomm->layout == Comm::LAYOUT_TILED)
    error->all(FLERR,"Cannot change to comm_style brick/direct from tiled layout");

  style = Comm::BRICK_DIRECT;
}

/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

void CommBrickDirect::setup()
{
  CommBrick::setup();

}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommBrickDirect::forward_comm(int /*dummy*/)
{
  int n,iswap;
  AtomVec *avec = atom->avec;
  double **x = atom->x;
  double *buf;

  // exchange atoms directly with all neighbor procs, not via 6-way stencil
  // if other proc is self, just copy
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  // post all receives for ghost atoms
  // except for self copies

  int irecv = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (sendproc_direct[iswap] == me) continue;
    if (comm_x_only) {
      if (size_forward_recv_direct[iswap]) {
        buf = x[firstrecv_direct[iswap]];
        MPI_Irecv(buf,size_forward_recv_direct[iswap],MPI_DOUBLE,
                  recvproc_direct[iswap],0,world,&requests[irecv++]);
      }
    } else {
      if (size_forward_recv_direct[iswap])
        MPI_Irecv(buf_recv_direct[iswap],size_forward_recv_direct[iswap],MPI_DOUBLE,
                  recvproc_direct[iswap],0,world,&requests[irecv++]);
    }
  }
                    
  // send all owned atoms to receiving procs
  // except for self copies

  for (iswap = 0; iswap < ndirect; iswap++) {
    if (sendproc_direct[iswap] == me) continue;
    if (comm_x_only) {
      n = avec->pack_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                          pbc_flag_direct[iswap],pbc_direct[iswap]);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,sendproc_direct[iswap],0,world);
    } else if (ghost_velocity) {
      n = avec->pack_comm_vel(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                              pbc_flag_direct[iswap],pbc_direct[iswap]);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,sendproc_direct[iswap],0,world);
    } else{
      n = avec->pack_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                          pbc_flag_direct[iswap],pbc_direct[iswap]);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,sendproc_direct[iswap],0,world);
    }
  }

  // copy atoms to self via pack and unpack

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_direct[iself];
    if (comm_x_only) {
      if (sendnum_direct[iswap])
        avec->pack_comm(sendnum_direct[iswap],sendlist_direct[iswap],
                        x[firstrecv_direct[iswap]],pbc_flag_direct[iswap],pbc_direct[iswap]);
    } else if (ghost_velocity) {
      avec->pack_comm_vel(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                          pbc_flag_direct[iswap],pbc_direct[iswap]);
      avec->unpack_comm_vel(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    } else {
      avec->pack_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                      pbc_flag_direct[iswap],pbc_direct[iswap]);
      avec->unpack_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    }
  }

  // wait on incoming messages with ghost atoms
  // unpack each message as it arrives

  if (nrecv_direct) {
    if (comm_x_only) {
      MPI_Waitall(nrecv_direct,requests,MPI_STATUS_IGNORE);
    } else if (ghost_velocity) {
      for (int irecv; irecv < nrecv_direct; irecv++) {
        MPI_Waitany(nrecv_direct,requests,&iswap,MPI_STATUS_IGNORE);
        avec->unpack_comm_vel(recvnum_direct[iswap],firstrecv_direct[iswap],buf_recv_direct[iswap]);
      }
    } else {
      for (int irecv; irecv < nrecv_direct; irecv++) {
        MPI_Waitany(nrecv_direct,requests,&iswap,MPI_STATUS_IGNORE);
        avec->unpack_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_recv_direct[iswap]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication of forces on atoms every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommBrickDirect::reverse_comm()
{
  int n,iswap;
  MPI_Request request;
  AtomVec *avec = atom->avec;
  double **f = atom->f;
  double *buf;

  // exchange atoms directly with all neighbor procs, not via 6-way stencil
  // if other proc is self, just copy
  // if comm_f_only set, exchange or copy directly from f, don't pack

  // post all receives for owned atoms
  // except for self copies + sums

  int irecv = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (sendproc_direct[iswap] == me) continue;
    if (size_reverse_recv_direct[iswap]) {
      MPI_Irecv(buf_recv_direct[iswap],size_reverse_recv_direct[iswap],MPI_DOUBLE,
                sendproc_direct[iswap],0,world,&requests[irecv++]);
    }
  }

  // send all ghost atoms to receiving procs
  // except for self copies + sums

  for (iswap = 0; iswap < ndirect; iswap++) {
    if (sendproc_direct[iswap] == me) continue;
    if (comm_f_only) {
      if (size_reverse_send_direct[iswap]) {
        buf = f[firstrecv_direct[iswap]];
        MPI_Send(buf,size_reverse_send_direct[iswap],MPI_DOUBLE,recvproc_direct[iswap],0,world);
      }
    } else {
      n = avec->pack_reverse(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,recvproc_direct[iswap],0,world);
    }
  }

  // copy + sum atoms to self via pack and unpack

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_direct[iself];
    if (comm_f_only) {
      if (sendnum_direct[iswap])
        avec->unpack_reverse(sendnum_direct[iswap],sendlist_direct[iswap],
                             f[firstrecv_direct[iswap]]);
    } else {
      avec->pack_reverse(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
      avec->unpack_reverse(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct);
    }
  }

  // wait on incoming messages with owned atoms
  // unpack each message as it arrives

  if (nrecv_direct) {
    for (int irecv; irecv < nrecv_direct; irecv++) {
      MPI_Waitany(nrecv_direct,requests,&iswap,MPI_STATUS_IGNORE);
      avec->unpack_reverse(sendnum_direct[iswap],sendlist_direct[iswap],buf_recv_direct[iswap]);
    }
  }
}

/* ----------------------------------------------------------------------
   borders: list nearby atoms to send to neighboring procs at every timestep
   one list is created for every swap that will be made
   as list is made, actually do swaps
   this does equivalent of a forward_comm(), so don't need to explicitly
     call forward_comm() on reneighboring timestep
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
------------------------------------------------------------------------- */

void CommBrickDirect::borders()
{
  /*
  int i,n,itype,icollection,iswap,dim,ineed,twoneed;
  int nsend,nrecv,sendflag,nfirst,nlast,ngroup,nprior;
  double lo,hi;
  int *type;
  int *collection;
  double **x;
  double *buf,*mlo,*mhi;
  MPI_Request request;
  AtomVec *avec = atom->avec;

  // After exchanging/sorting, need to reconstruct collection array for border communication
  if (mode == Comm::MULTI) neighbor->build_collection(0);

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
      // store sent atom indices in sendlist for use in future timesteps

      x = atom->x;
      if (mode == Comm::SINGLE) {
        lo = slablo[iswap];
        hi = slabhi[iswap];
      } else if (mode == Comm::MULTI) {
        collection = neighbor->collection;
        mlo = multilo[iswap];
        mhi = multihi[iswap];
      } else {
        type = atom->type;
        mlo = multioldlo[iswap];
        mhi = multioldhi[iswap];
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
      // all atoms eligible versus only atoms in bordergroup
      // can only limit loop to bordergroup for first sends (ineed < 2)
      // on these sends, break loop in two: owned (in group) and ghost

      if (sendflag) {
        if (!bordergroup || ineed >= 2) {
          if (mode == Comm::SINGLE) {
            for (i = nfirst; i < nlast; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
          } else if (mode == Comm::MULTI) {
            for (i = nfirst; i < nlast; i++) {
              icollection = collection[i];
              if (x[i][dim] >= mlo[icollection] && x[i][dim] <= mhi[icollection]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
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
          if (mode == Comm::SINGLE) {
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
          } else if (mode == Comm::MULTI) {
            ngroup = atom->nfirst;
            for (i = 0; i < ngroup; i++) {
              icollection = collection[i];
              if (x[i][dim] >= mlo[icollection] && x[i][dim] <= mhi[icollection]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
            for (i = atom->nlocal; i < nlast; i++) {
              icollection = collection[i];
              if (x[i][dim] >= mlo[icollection] && x[i][dim] <= mhi[icollection]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
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
        n = avec->pack_border_vel(nsend,sendlist[iswap],buf_send,pbc_flag[iswap],pbc[iswap]);
      else
        n = avec->pack_border(nsend,sendlist[iswap],buf_send,pbc_flag[iswap],pbc[iswap]);

      // swap atoms with other proc
      // no MPI calls except SendRecv if nsend/nrecv = 0
      // put incoming ghosts at end of my atom arrays
      // if swapping with self, simply copy, no messages

      if (sendproc[iswap] != me) {
        MPI_Sendrecv(&nsend,1,MPI_INT,sendproc[iswap],0,
                     &nrecv,1,MPI_INT,recvproc[iswap],0,world,MPI_STATUS_IGNORE);
        if (nrecv*size_border > maxrecv) grow_recv(nrecv*size_border);
        if (nrecv) MPI_Irecv(buf_recv,nrecv*size_border,MPI_DOUBLE,
                             recvproc[iswap],0,world,&request);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        if (nrecv) MPI_Wait(&request,MPI_STATUS_IGNORE);
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
      nprior = atom->nlocal + atom->nghost;
      atom->nghost += nrecv;
      if (neighbor->style == Neighbor::MULTI) neighbor->build_collection(nprior);

      iswap++;
    }
  }

  // For molecular systems we lose some bits for local atom indices due
  // to encoding of special pairs in neighbor lists. Check for overflows.

  if ((atom->molecular != Atom::ATOMIC)
      && ((atom->nlocal + atom->nghost) > NEIGHMASK))
    error->one(FLERR,"Per-processor number of atoms is too large for "
               "molecular neighbor lists");

  // ensure send/recv buffers are long enough for all forward & reverse comm

  int max = MAX(maxforward*smax,maxreverse*rmax);
  if (max > maxsend) grow_send(max,0);
  max = MAX(maxforward*rmax,maxreverse*smax);
  if (max > maxrecv) grow_recv(max);

  // reset global->local map

  if (map_style != Atom::MAP_NONE) atom->map_set();
  */
}

