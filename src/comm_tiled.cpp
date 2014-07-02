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

#include "lmptype.h"
#include "comm_tiled.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "pair.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "dump.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BUFFACTOR 1.5

enum{SINGLE,MULTI};               // same as in Comm

/* ---------------------------------------------------------------------- */

CommTiled::CommTiled(LAMMPS *lmp) : Comm(lmp)
{
  style = 1;
  layout = 0;

  error->all(FLERR,"Comm_style tiled is not yet supported");
}

/* ---------------------------------------------------------------------- */

CommTiled::~CommTiled()
{
}

/* ---------------------------------------------------------------------- */

void CommTiled::init()
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
}

/* ----------------------------------------------------------------------
   setup spatial-decomposition communication patterns
   function of neighbor cutoff(s) & cutghostuser & current box size
   single mode sets slab boundaries (slablo,slabhi) based on max cutoff
   multi mode sets type-dependent slab boundaries (multilo,multihi)
------------------------------------------------------------------------- */

void CommTiled::setup()
{
  // error on triclinic or multi?
  // set nswap = 2*dim
  // setup neighbor proc info for exchange()
  // setup nsendproc and nrecvproc bounts
  // setup sendproc and recvproc lists
  // setup sendboxes
  // reallocate requests and statuses

  // check that cutoff is <= 1/2 of periodic box len?

  // loop over dims
  //  left:
  //    construct ghost boxes
  //      differnet in x,y,z
  //      account for ghost borders in y,z
  //      account for PBC by shifting
  //      split into multiple boxes if straddles PBC
  //    drop boxes down RCB tree
  //      count unique procs they cover
  //      what about self if crosses PBC
  //      for each proc they cover:
  //        compute box I send it to left
  //        is a message I will recv from right (don't care about box)
  //    for ghost-extended boxes
  //      do not count procs that do not overlap my owned box at all
  //      only touching edge of my owned box does not count
  //      in this case list I send to and recv from may be different?
  //  same thing to right

  // what need from decomp (RCB):
  // dropbox: return list of procs with overlap and overlapping boxes
  //   return n, proclist, boxlist
  // otherbox: bbox of another proc
  // dropatom: return what proc owns the atom coord
}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommTiled::forward_comm(int dummy)
{
  int i,irecv,n;
  MPI_Status status;
  AtomVec *avec = atom->avec;
  double **x = atom->x;

  // exchange data with another set of procs in each swap
  // if first proc in set is self, then is just self across PBC, just copy
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap][0] != me) {
      if (comm_x_only) {
        for (i = 0; i < nrecvproc[iswap]; i++)
          MPI_Irecv(x[firstrecv[iswap][i]],size_forward_recv[iswap][i],
                    MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        for (i = 0; i < nsendproc[iswap]; i++) {
          n = avec->pack_comm(sendnum[iswap][i],sendlist[iswap][i],
                              buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
        }
        MPI_Waitall(nrecvproc[iswap],requests,statuses);

      } else if (ghost_velocity) {
        for (i = 0; i < nrecvproc[iswap]; i++)
          MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                    size_forward_recv[iswap][i],
                    MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        for (i = 0; i < nsendproc[iswap]; i++) {
          n = avec->pack_comm_vel(sendnum[iswap][i],sendlist[iswap][i],
                                  buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
        }
        for (i = 0; i < nrecvproc[iswap]; i++) {
          MPI_Waitany(nrecvproc[iswap],requests,&irecv,&status);
          avec->unpack_comm_vel(recvnum[iswap][i],firstrecv[iswap][irecv],
                                &buf_recv[forward_recv_offset[iswap][irecv]]);
        }

      } else {
        for (i = 0; i < nrecvproc[iswap]; i++)
          MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                    size_forward_recv[iswap][i],
                    MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        for (i = 0; i < nsendproc[iswap]; i++) {
          n = avec->pack_comm(sendnum[iswap][i],sendlist[iswap][i],
                              buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
        }
        for (i = 0; i < nrecvproc[iswap]; i++) {
          MPI_Waitany(nrecvproc[iswap],requests,&irecv,&status);
          avec->unpack_comm(recvnum[iswap][i],firstrecv[iswap][irecv],
                            &buf_recv[forward_recv_offset[iswap][irecv]]);
        }
      }

    } else {
      if (comm_x_only) {
        if (sendnum[iswap][0])
          n = avec->pack_comm(sendnum[iswap][0],sendlist[iswap][0],
                              x[firstrecv[iswap][0]],pbc_flag[iswap][0],
                              pbc[iswap][0]);
      } else if (ghost_velocity) {
        n = avec->pack_comm_vel(sendnum[iswap][0],sendlist[iswap][0],
                                buf_send,pbc_flag[iswap][0],pbc[iswap][0]);
        avec->unpack_comm_vel(recvnum[iswap][0],firstrecv[iswap][0],buf_send);
      } else {
        n = avec->pack_comm(sendnum[iswap][0],sendlist[iswap][0],
                            buf_send,pbc_flag[iswap][0],pbc[iswap][0]);
        avec->unpack_comm(recvnum[iswap][0],firstrecv[iswap][0],buf_send);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication of forces on atoms every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommTiled::reverse_comm()
{
  int i,irecv,n;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;
  double **f = atom->f;

  // exchange data with another set of procs in each swap
  // if first proc in set is self, then is just self across PBC, just copy
  // if comm_f_only set, exchange or copy directly from f, don't pack

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (sendproc[iswap][0] != me) {
      if (comm_f_only) {
        for (i = 0; i < nsendproc[iswap]; i++)
          MPI_Irecv(&buf_recv[reverse_recv_offset[iswap][i]],
                    size_reverse_recv[iswap][i],MPI_DOUBLE,
                    sendproc[iswap][i],0,world,&requests[i]);
        for (i = 0; i < nrecvproc[iswap]; i++)
          MPI_Send(f[firstrecv[iswap][i]],size_reverse_send[iswap][i],
                   MPI_DOUBLE,recvproc[iswap][i],0,world);
        for (i = 0; i < nsendproc[iswap]; i++) {
          MPI_Waitany(nsendproc[iswap],requests,&irecv,&status);
          avec->unpack_reverse(sendnum[iswap][irecv],sendlist[iswap][irecv],
                               &buf_recv[reverse_recv_offset[iswap][irecv]]);
        }

      } else {
        for (i = 0; i < nsendproc[iswap]; i++)
          MPI_Irecv(&buf_recv[reverse_recv_offset[iswap][i]],
                    size_reverse_recv[iswap][i],MPI_DOUBLE,
                    sendproc[iswap][i],0,world,&requests[i]);
        for (i = 0; i < nrecvproc[iswap]; i++) {
          n = avec->pack_reverse(recvnum[iswap][i],firstrecv[iswap][i],
                                 buf_send);
          MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap][i],0,world);
        }
        for (i = 0; i < nsendproc[iswap]; i++) {
          MPI_Waitany(nsendproc[iswap],requests,&irecv,&status);
          avec->unpack_reverse(sendnum[iswap][irecv],sendlist[iswap][irecv],
                               &buf_recv[reverse_recv_offset[iswap][irecv]]);
        }
      }

    } else {
      if (comm_f_only) {
        if (sendnum[iswap][0])
          avec->unpack_reverse(sendnum[iswap][0],sendlist[iswap][0],
                               f[firstrecv[iswap][0]]);
      } else {
        n = avec->pack_reverse(recvnum[iswap][0],firstrecv[iswap][0],buf_send);
        avec->unpack_reverse(sendnum[iswap][0],sendlist[iswap][0],buf_send);
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

void CommTiled::exchange()
{
  // loop over atoms
  //   if not outside my box, continue
  //   find which proc it is in
  //   find which one of my touching procs it is, else lost
  //   make sure all atoms are "lost" that should be (e.g. outside non-PBC)
  //   add to list to send to that proc
  // loop over touching procs
  //   send buffer to them
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

void CommTiled::borders()
{
  int i,n,irecv,ngroup,nlast,nsend,rmaxswap;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double *bbox;
  double **x;
  MPI_Status status;
  AtomVec *avec = atom->avec;

  // smax = max size of single send in a swap/proc
  // rmax = max size of recvs from all procs for a swap

  int smax = 0;
  int rmax = 0;

  // loop over all swaps in all dimensions

  for (int iswap = 0; iswap < nswap; iswap++) {

    // find atoms within rectangles using <= and >=
    // for x-dim swaps, check owned atoms
    // for yz-dim swaps, check owned and ghost atoms
    // store sent atom indices in list for use in future timesteps
    // NOTE: assume SINGLE mode, add back in logic for MULTI mode later

    x = atom->x;

    for (i = 0; i < nsendproc[iswap]; i++) {
      bbox = sendbox[iswap][i];
      xlo = bbox[0]; xhi = bbox[1];
      ylo = bbox[2]; yhi = bbox[3];
      zlo = bbox[4]; zhi = bbox[5];

      ngroup = atom->nfirst;
      if (iswap < 2) nlast = atom->nlocal;
      else nlast = atom->nlocal + atom->nghost;

      nsend = 0;
      for (i = 0; i < ngroup; i++)
        if (x[i][0] >= xlo && x[i][0] <= xhi &&
            x[i][1] >= ylo && x[i][1] <= yhi &&
            x[i][2] >= zlo && x[i][2] <= zhi) {
          if (nsend == maxsendlist[iswap][i]) grow_list(iswap,i,nsend);
          sendlist[iswap][i][nsend++] = i;
        }
      for (i = atom->nlocal; i < nlast; i++)
        if (x[i][0] >= xlo && x[i][0] <= xhi &&
            x[i][1] >= ylo && x[i][1] <= yhi &&
            x[i][2] >= zlo && x[i][2] <= zhi) {
          if (nsend == maxsendlist[iswap][i]) grow_list(iswap,i,nsend);
          sendlist[iswap][i][nsend++] = i;
        }
      sendnum[iswap][i] = nsend;
      smax = MAX(smax,nsend);
    }

    // send sendnum counts to procs who recv from me

    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nrecvproc[iswap]; i++)
        MPI_Irecv(&recvnum[iswap][i],1,MPI_INT,
                  recvproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nsendproc[iswap]; i++)
        MPI_Send(&sendnum[iswap][i],1,MPI_INT,sendproc[iswap][i],0,world);
      MPI_Waitall(nrecvproc[iswap],requests,statuses);

    } else recvnum[iswap][0] = sendnum[iswap][0];

    // setup other per swap/proc values from sendnum and recvnum

    rmaxswap = 0;
    for (i = 0; i < nrecvproc[iswap]; i++) {
      rmaxswap += recvnum[iswap][i];
      size_forward_recv[iswap][i] = recvnum[iswap][i]*size_forward;
      size_reverse_send[iswap][i] = recvnum[iswap][i]*size_reverse;
      size_reverse_recv[iswap][i] = sendnum[iswap][i]*size_reverse;
      if (i == 0) {
        firstrecv[iswap][0] = atom->nlocal + atom->nghost;
        forward_recv_offset[iswap][0] = 0;
      } else {
        firstrecv[iswap][i] = firstrecv[iswap][i-1] + recvnum[iswap][i-1];
        forward_recv_offset[iswap][i] = 
          forward_recv_offset[iswap][i-1] + recvnum[iswap][i-1];
      }
    }
    rmax = MAX(rmax,rmaxswap);

    // insure send/recv buffers are large enough for border comm

    if (smax*size_border > maxsend) grow_send(smax*size_border,0);
    if (rmax*size_border > maxrecv) grow_recv(rmax*size_border);

    // swap atoms with other procs using pack_border(), unpack_border()

    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nsendproc[iswap]; i++) {
        if (ghost_velocity) {
          for (i = 0; i < nrecvproc[iswap]; i++)
            MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                      recvnum[iswap][i]*size_border,
                      MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
          for (i = 0; i < nsendproc[iswap]; i++) {
            n = avec->pack_border_vel(sendnum[iswap][i],sendlist[iswap][i],
                                      buf_send,pbc_flag[iswap][i],
                                      pbc[iswap][i]);
            MPI_Send(buf_send,n*size_border,MPI_DOUBLE,
                     sendproc[iswap][i],0,world);
          }
          for (i = 0; i < nrecvproc[iswap]; i++) {
            MPI_Waitany(nrecvproc[iswap],requests,&irecv,&status);
            avec->unpack_border(recvnum[iswap][i],firstrecv[iswap][irecv],
                                &buf_recv[forward_recv_offset[iswap][irecv]]);
          }
          
        } else {
          for (i = 0; i < nrecvproc[iswap]; i++)
            MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                      recvnum[iswap][i]*size_border,
                      MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
          for (i = 0; i < nsendproc[iswap]; i++) {
            n = avec->pack_border(sendnum[iswap][i],sendlist[iswap][i],
                                  buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
            MPI_Send(buf_send,n*size_border,MPI_DOUBLE,
                     sendproc[iswap][i],0,world);
          }
          for (i = 0; i < nrecvproc[iswap]; i++) {
            MPI_Waitany(nrecvproc[iswap],requests,&irecv,&status);
            avec->unpack_border(recvnum[iswap][i],firstrecv[iswap][irecv],
                                &buf_recv[forward_recv_offset[iswap][irecv]]);
          }
        }
      }

    } else {
      if (ghost_velocity) {
        n = avec->pack_border_vel(sendnum[iswap][0],sendlist[iswap][0],
                                  buf_send,pbc_flag[iswap][0],pbc[iswap][0]);
        avec->unpack_border_vel(recvnum[iswap][0],firstrecv[iswap][0],buf_send);
      } else {
        n = avec->pack_border(sendnum[iswap][0],sendlist[iswap][0],
                              buf_send,pbc_flag[iswap][0],pbc[iswap][0]);
        avec->unpack_border(recvnum[iswap][0],firstrecv[iswap][0],buf_send);
      }
    }

    // increment ghost atoms

    n = nrecvproc[iswap];
    atom->nghost += forward_recv_offset[iswap][n-1] + recvnum[iswap][n-1];
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
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::forward_comm_pair(Pair *pair)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nrecvproc[iswap]; i++)
        MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                  size_forward_recv[iswap][i],
                  MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nsendproc[iswap]; i++) {
        n = pair->pack_comm(sendnum[iswap][i],sendlist[iswap][i],
                            buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
        MPI_Send(buf_send,n*sendnum[iswap][i],MPI_DOUBLE,
                 sendproc[iswap][i],0,world);
      }
      for (i = 0; i < nrecvproc[iswap]; i++) {
        MPI_Waitany(nrecvproc[iswap],requests,&irecv,&status);
        pair->unpack_comm(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                          &buf_recv[forward_recv_offset[iswap][irecv]]);
      }

    } else {
      n = pair->pack_comm(sendnum[iswap][0],sendlist[iswap][0],
                          buf_send,pbc_flag[iswap][0],pbc[iswap][0]);
      pair->unpack_comm(recvnum[iswap][0],firstrecv[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Pair
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_pair(Pair *pair)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nsendproc[iswap]; i++)
        MPI_Irecv(&buf_recv[reverse_recv_offset[iswap][i]],
                  size_reverse_recv[iswap][i],MPI_DOUBLE,
                  sendproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nrecvproc[iswap]; i++) {
        n = pair->pack_reverse_comm(recvnum[iswap][i],firstrecv[iswap][i],
                                    buf_send);
        MPI_Send(buf_send,n*recvnum[iswap][i],MPI_DOUBLE,
                 recvproc[iswap][i],0,world);
      }
      for (i = 0; i < nsendproc[iswap]; i++) {
        MPI_Waitany(nsendproc[iswap],requests,&irecv,&status);
        pair->unpack_reverse_comm(sendnum[iswap][irecv],sendlist[iswap][irecv],
                                  &buf_recv[reverse_recv_offset[iswap][irecv]]);
      }

    } else {
      n = pair->pack_reverse_comm(recvnum[iswap][0],firstrecv[iswap][0],
                                  buf_send);
      pair->unpack_reverse_comm(sendnum[iswap][0],sendlist[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::forward_comm_fix(Fix *fix)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nrecvproc[iswap]; i++)
        MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                  size_forward_recv[iswap][i],
                  MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nsendproc[iswap]; i++) {
        n = fix->pack_comm(sendnum[iswap][i],sendlist[iswap][i],
                           buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
        MPI_Send(buf_send,n*sendnum[iswap][i],MPI_DOUBLE,
                 sendproc[iswap][i],0,world);
      }
      for (i = 0; i < nrecvproc[iswap]; i++) {
        MPI_Waitany(nrecvproc[iswap],requests,&irecv,&status);
        fix->unpack_comm(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                         &buf_recv[forward_recv_offset[iswap][irecv]]);
      }

    } else {
      n = fix->pack_comm(sendnum[iswap][0],sendlist[iswap][0],
                         buf_send,pbc_flag[iswap][0],pbc[iswap][0]);
      fix->unpack_comm(recvnum[iswap][0],firstrecv[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_fix(Fix *fix)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nsendproc[iswap]; i++)
        MPI_Irecv(&buf_recv[reverse_recv_offset[iswap][i]],
                  size_reverse_recv[iswap][i],MPI_DOUBLE,
                  sendproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nrecvproc[iswap]; i++) {
        n = fix->pack_reverse_comm(recvnum[iswap][i],firstrecv[iswap][i],
                                   buf_send);
        MPI_Send(buf_send,n*recvnum[iswap][i],MPI_DOUBLE,
                 recvproc[iswap][i],0,world);
      }
      for (i = 0; i < nsendproc[iswap]; i++) {
        MPI_Waitany(nsendproc[iswap],requests,&irecv,&status);
        fix->unpack_reverse_comm(sendnum[iswap][irecv],sendlist[iswap][irecv],
                                 &buf_recv[reverse_recv_offset[iswap][irecv]]);
      }

    } else {
      n = fix->pack_reverse_comm(recvnum[iswap][0],firstrecv[iswap][0],
                                 buf_send);
      fix->unpack_reverse_comm(sendnum[iswap][0],sendlist[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   n = total datums for all atoms, allows for variable number/atom
   NOTE: complicated b/c don't know # to recv a priori
------------------------------------------------------------------------- */

void CommTiled::forward_comm_variable_fix(Fix *fix)
{
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   n = total datums for all atoms, allows for variable number/atom
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_variable_fix(Fix *fix)
{
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Compute
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::forward_comm_compute(Compute *compute)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nrecvproc[iswap]; i++)
        MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                  size_forward_recv[iswap][i],
                  MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nsendproc[iswap]; i++) {
        n = compute->pack_comm(sendnum[iswap][i],sendlist[iswap][i],
                               buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
        MPI_Send(buf_send,n*sendnum[iswap][i],MPI_DOUBLE,
                 sendproc[iswap][i],0,world);
      }
      for (i = 0; i < nrecvproc[iswap]; i++) {
        MPI_Waitany(nrecvproc[iswap],requests,&irecv,&status);
        compute->unpack_comm(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                             &buf_recv[forward_recv_offset[iswap][irecv]]);
      }

    } else {
      n = compute->pack_comm(sendnum[iswap][0],sendlist[iswap][0],
                             buf_send,pbc_flag[iswap][0],pbc[iswap][0]);
      compute->unpack_comm(recvnum[iswap][0],firstrecv[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Compute
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_compute(Compute *compute)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nsendproc[iswap]; i++)
        MPI_Irecv(&buf_recv[reverse_recv_offset[iswap][i]],
                  size_reverse_recv[iswap][i],MPI_DOUBLE,
                  sendproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nrecvproc[iswap]; i++) {
        n = compute->pack_reverse_comm(recvnum[iswap][i],firstrecv[iswap][i],
                                   buf_send);
        MPI_Send(buf_send,n*recvnum[iswap][i],MPI_DOUBLE,
                 recvproc[iswap][i],0,world);
      }
      for (i = 0; i < nsendproc[iswap]; i++) {
        MPI_Waitany(nsendproc[iswap],requests,&irecv,&status);
        compute->
          unpack_reverse_comm(sendnum[iswap][irecv],sendlist[iswap][irecv],
                              &buf_recv[reverse_recv_offset[iswap][irecv]]);
      }

    } else {
      n = compute->pack_reverse_comm(recvnum[iswap][0],firstrecv[iswap][0],
                                 buf_send);
      compute->unpack_reverse_comm(sendnum[iswap][0],sendlist[iswap][0],
                                   buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Dump
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::forward_comm_dump(Dump *dump)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nrecvproc[iswap]; i++)
        MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                  size_forward_recv[iswap][i],
                  MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nsendproc[iswap]; i++) {
        n = dump->pack_comm(sendnum[iswap][i],sendlist[iswap][i],
                            buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
        MPI_Send(buf_send,n*sendnum[iswap][i],MPI_DOUBLE,
                 sendproc[iswap][i],0,world);
      }
      for (i = 0; i < nrecvproc[iswap]; i++) {
        MPI_Waitany(nrecvproc[iswap],requests,&irecv,&status);
        dump->unpack_comm(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                          &buf_recv[forward_recv_offset[iswap][irecv]]);
      }

    } else {
      n = dump->pack_comm(sendnum[iswap][0],sendlist[iswap][0],
                          buf_send,pbc_flag[iswap][0],pbc[iswap][0]);
      dump->unpack_comm(recvnum[iswap][0],firstrecv[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Dump
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_dump(Dump *dump)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nsendproc[iswap]; i++)
        MPI_Irecv(&buf_recv[reverse_recv_offset[iswap][i]],
                  size_reverse_recv[iswap][i],MPI_DOUBLE,
                  sendproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nrecvproc[iswap]; i++) {
        n = dump->pack_reverse_comm(recvnum[iswap][i],firstrecv[iswap][i],
                                    buf_send);
        MPI_Send(buf_send,n*recvnum[iswap][i],MPI_DOUBLE,
                 recvproc[iswap][i],0,world);
      }
      for (i = 0; i < nsendproc[iswap]; i++) {
        MPI_Waitany(nsendproc[iswap],requests,&irecv,&status);
        dump->unpack_reverse_comm(sendnum[iswap][irecv],sendlist[iswap][irecv],
                                  &buf_recv[reverse_recv_offset[iswap][irecv]]);
      }

    } else {
      n = dump->pack_reverse_comm(recvnum[iswap][0],firstrecv[iswap][0],
                                  buf_send);
      dump->unpack_reverse_comm(sendnum[iswap][0],sendlist[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   forward communication of N values in array
------------------------------------------------------------------------- */

void CommTiled::forward_comm_array(int n, double **array)
{
}

/* ----------------------------------------------------------------------
   exchange info provided with all 6 stencil neighbors
------------------------------------------------------------------------- */

int CommTiled::exchange_variable(int n, double *inbuf, double *&outbuf)
{
  int nrecv = n;
  return nrecv;
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR and bufextra
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void CommTiled::grow_send(int n, int flag)
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

void CommTiled::grow_recv(int n)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  memory->destroy(buf_recv);
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommTiled::grow_list(int iswap, int iwhich, int n)
{
  maxsendlist[iswap][iwhich] = static_cast<int> (BUFFACTOR * n);
  memory->grow(sendlist[iswap][iwhich],maxsendlist[iswap][iwhich],
               "comm:sendlist[iswap]");
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint CommTiled::memory_usage()
{
  bigint bytes = 0;
  return bytes;
}
