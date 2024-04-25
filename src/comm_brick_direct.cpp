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
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "neighbor.h"

// NOTES:
// allocate requests to length of nrecv_direct
// what do lengths of send/recv bufs need to be
// do not allow MULTI with brick/direct
// do not allow bordergroup
// how to order dswap

using namespace LAMMPS_NS;

static constexpr double BUFFACTOR = 1.5;
static constexpr int BUFMIN = 1024;

/* ---------------------------------------------------------------------- */

CommBrickDirect::CommBrickDirect(LAMMPS *lmp) : CommBrick(lmp)
{
  style = Comm::BRICK_DIRECT;

  dswap = nullptr;
  requests = nullptr;
  maxdirect = 0;
}

/* ---------------------------------------------------------------------- */

CommBrickDirect::~CommBrickDirect()
{
  delete [] dswap;
  delete [] requests;
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

  // use recvneed to create logical 3d grid of procs to perform direct comm with
  // stored in dswap = list of DirectSwaps
  
  int dim = domain->dimension;

  double *prd,*sublo,*subhi;

  if (triclinic == 0) {
    prd = domain->prd;
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    prd = domain->prd_lamda;
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }
  
  // ndirect = # of direct swaps with other procs, including self copies
  // subtract 1 for myself in center of 3d grid of surrounding procs
  // ijk lo/hi = bounds of stencil around my proc at center
  
  int ilo = -recvneed[0][0];
  int ihi = recvneed[0][1];
  int jlo = -recvneed[1][0];
  int jhi = recvneed[1][1];
  int klo = -recvneed[2][0];
  int khi = recvneed[2][1];

  ndirect = (ihi-ilo+1) * (jhi-jlo+1) * (khi-klo+1) - 1;

  if (ndirect > maxdirect) {
    delete [] dswap;
    dswap = new DirectSwap[ndirect];
    delete [] requests;
    requests = new MPI_Request[ndirect];
  }

  // loop over stencil and define each direct swap
  // NOTE: need to order the direct swaps as desired
  
  DirectSwap *ds;
  int ix,iy,iz;
  int igrid,jgrid,kgrid;
  int xpbc,ypbc,zpbc;
  
  int iswap = 0;

  for (iz = klo; iz <= khi; iz++) {
    for (iy = jlo; iy <= jhi; iy++) {
      for (ix = ilo; ix <= ihi; ix++) {

        // skip center of stencil = my subdomain
        
        if (ix == 0 && iy == 0 && iz == 0) continue;
  
        ds = &dswap[iswap];
        xpbc = ypbc = zpbc = 0;
        
        igrid = myloc[0] + ix;
        while (igrid < 0) {
          igrid += procgrid[0];
          xpbc++;
        }
        while (igrid >= procgrid[0]) {
          igrid -= procgrid[0];
          xpbc--;
        }

        jgrid = myloc[1] + iy;
        while (jgrid < 0) {
          jgrid += procgrid[1];
          ypbc++;
        }
        while (jgrid >= procgrid[1]) {
          jgrid -= procgrid[1];
          ypbc--;
        }

        kgrid = myloc[2] + iz;
        while (kgrid < 0) {
          kgrid += procgrid[2];
          zpbc++;
        }
        while (kgrid >= procgrid[2]) {
          kgrid -= procgrid[2];
          zpbc--;
        }

        ds->proc = grid2proc[igrid][jgrid][kgrid];

        // NOTE: for multiple procs in stencil, cutghost needs to
        //       have width of inbetween subdomains subtracted, via xyzsplit
        //       for orthogonal or triclinic

        if (ix > ilo && ix < ihi) ds->xcheck = 0;
        else {
          ds->xcheck = 1;
          if (ix == ilo) {
            ds->xlo = sublo[0];
            ds->xhi = sublo[0] + cutghost[0];
          } else if (ix == ihi) {
            ds->xlo = subhi[0] - cutghost[0];
            ds->xhi = subhi[0];
          }
        }

        if (iy > jlo && iy < jhi) ds->ycheck = 0;
        else {
          ds->ycheck = 1;
          if (iy == jlo) {
            ds->ylo = sublo[1];
            ds->yhi = sublo[1] + cutghost[1];
          } else if (iy == jhi) {
            ds->ylo = subhi[1] - cutghost[1];
            ds->yhi = subhi[1];
          }
        }

        if (dim == 2) ds->zcheck = 0;
        else if (iz > klo && iz < khi) ds->zcheck = 0;
        else {
          ds->zcheck = 1;
          if (iz == klo) {
            ds->zlo = sublo[2];
            ds->zhi = sublo[2] + cutghost[2];
          } else if (iz == khi) {
            ds->zlo = subhi[2] - cutghost[2];
            ds->zhi = subhi[2];
          }
        }
        
        if (!ds->xcheck and !ds->ycheck && !ds->zcheck) ds->allflag = 1;
        else ds->allflag = 0;

        ds->pbc_flag = 0;
        ds->pbc[0] = ds->pbc[1] = ds->pbc[2] = ds->pbc[3] = ds->pbc[4] = ds->pbc[5] = 0;
        if (xpbc || !ypbc || zpbc) {
          ds->pbc[0] = xpbc;
          ds->pbc[1] = ypbc;
          ds->pbc[2] = zpbc;
          if (triclinic) {
            pbc[5] = pbc[1];
            pbc[4] = pbc[3] = pbc[2];
          }
        }
        
        iswap++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
   exchange owned atoms directly with all neighbor procs,
     not via CommBrick 6-way stencil
------------------------------------------------------------------------- */

void CommBrickDirect::forward_comm(int /*dummy*/)
{
  int n,iswap,irecv,nrecv;
  AtomVec *avec = atom->avec;
  double **x = atom->x;
  double *buf;

  // post all receives for ghost atoms
  // except for self copies

  nrecv = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (comm_x_only) {
      if (size_forward_recv_direct[iswap]) {
        buf = x[firstrecv_direct[iswap]];
        MPI_Irecv(buf,size_forward_recv_direct[iswap],MPI_DOUBLE,
                  proc_direct[iswap],0,world,&requests[nrecv++]);
      }
    } else {
      if (size_forward_recv_direct[iswap]) {
        MPI_Irecv(buf_recv_direct[iswap],size_forward_recv_direct[iswap],MPI_DOUBLE,
                  proc_direct[iswap],0,world,&requests[nrecv++]);
      }
    }
  }
                    
  // send all owned atoms to receiving procs
  // except for self copies

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (ghost_velocity) {
      n = avec->pack_comm_vel(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                              pbc_flag_direct[iswap],pbc_direct[iswap]);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],0,world);
    } else {
      n = avec->pack_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                          pbc_flag_direct[iswap],pbc_direct[iswap]);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],0,world);
    }
  }

  // copy atoms to self via pack and unpack

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    if (comm_x_only) {
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

  if (nrecv == 0) return;
  
  if (comm_x_only) {
    MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);
  } else if (ghost_velocity) {
    for (int i = 0; i < nrecv; i++) {
      MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
      iswap = recv_indices_direct[irecv];
      avec->unpack_comm_vel(recvnum_direct[iswap],firstrecv_direct[iswap],buf_recv_direct[iswap]);
    }
  } else {
    for (int i = 0; i < nrecv; i++) {
      MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
      iswap = recv_indices_direct[irecv];
      avec->unpack_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_recv_direct[iswap]);
    }
  }
}
  
/* ----------------------------------------------------------------------
   reverse communication of forces on atoms every timestep
   other per-atom attributes may also be sent via pack/unpack routines
   exchange ghost atoms directly with all neighbor procs,
     not via CommBrick 6-way stencil
------------------------------------------------------------------------- */

void CommBrickDirect::reverse_comm()
{
  int n,iswap,irecv,nrecv;
  MPI_Request request;
  AtomVec *avec = atom->avec;
  double **f = atom->f;
  double *buf;

  // post all receives for owned atoms
  // except for self copy/sums

  nrecv = 0;
  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (recvproc[iswap] == me) continue;
    if (size_reverse_recv_direct[iswap])
      MPI_Irecv(buf_recv_direct[iswap],size_reverse_recv_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],0,world,&requests[nrecv++]);
  }

  // send all ghost atoms to receiving procs
  // except for self copy/sums

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (recvproc[iswap] == me) continue;
    if (comm_f_only) {
      if (size_reverse_send_direct[iswap]) {
        buf = f[firstrecv_direct[iswap]];
        MPI_Send(buf,size_reverse_send_direct[iswap],MPI_DOUBLE,proc_direct[iswap],0,world);
      }
    } else {
      n = avec->pack_reverse(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],0,world);
    }
  }

  // copy/sum atoms to self via pack and unpack

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    if (comm_f_only) {
      avec->unpack_reverse(sendnum_direct[iswap],sendlist_direct[iswap],
                           f[firstrecv_direct[iswap]]);
    } else {
      avec->pack_reverse(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
      avec->unpack_reverse(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct);
    }
  }

  // wait on incoming messages with owned atoms
  // unpack each message as it arrives

  if (nrecv == 0) return;
  
  for (int i; i < nrecv; i++) {
    MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = send_indices_direct[irecv];
    avec->unpack_reverse(sendnum_direct[iswap],sendlist_direct[iswap],buf_recv_direct[iswap]);
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
  // loop over conventional 6-way BRICK swaps in 3 dimensions
  // construct BRICK_DIRECT swaps from them
  // unlike borders() in CommBrick, cannot perform borders comm until end
  // this is b/c the swaps take place simultaneously in all dimensions
  //   and thus cannot contain ghost atoms in the forward comm
------------------------------------------------------------------------- */

void CommBrickDirect::borders()
{
  int i,n,iswap,irecv,nrecv;
  
  AtomVec *avec = atom->avec;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  // setup lists of atoms to send in each direct swap

  DirectSwap *ds;
  int nsend,allflag,xcheck,ycheck,zcheck;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  
  for (iswap = 0; iswap < ndirect; iswap++) {
    ds = &dswap[iswap];
    nsend = 0;
    maxsend = maxsendlist_direct[iswap];
      
    allflag = ds->allflag;

    // NOTE: have another option for this send of all atoms?
    
    if (allflag) {
      for (i = 0; i < nlocal; i++) {
        if (nsend == maxsend) grow_list_direct(iswap,nsend);
        sendlist_direct[iswap][nsend++] = i;
      }
      
    } else {
      xcheck = ds->xcheck;
      ycheck = ds->ycheck;
      zcheck = ds->zcheck;

      xlo = ds->xlo;
      xhi = ds->xlo;
      ylo = ds->ylo;
      yhi = ds->ylo;
      zlo = ds->zlo;
      zhi = ds->zlo;
      
      for (i = 0; i < nlocal; i++) {
        if (xcheck && (x[i][0] < xlo || x[i][0] > xhi)) continue;
        if (ycheck && (x[i][1] < ylo || x[i][1] > yhi)) continue;
        if (zcheck && (x[i][2] < zlo || x[i][2] > zhi)) continue;
        if (nsend == maxsend) grow_list(iswap,nsend);
        sendlist_direct[iswap][nsend++] = i;
      }
    }

    sendnum_direct[iswap] = nsend;
    proc_direct[iswap] = ds->proc;
  }

  // send value of nsend for each swap to each receiving proc
  // post receives, perform sends, copy to self, wait for all incoming messages

  // NOTE: how to distinguish multiple messages between same 2 procs - use MSG type ?
  //       how will both sender and receiver agree on MSG type ?
  
  nrecv = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    MPI_Irecv(&recvnum_direct[iswap],1,MPI_INT,
              proc_direct[iswap],0,world,&requests[nrecv++]);
  }

  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    MPI_Send(&sendnum_direct[iswap],1,MPI_INT,proc_direct[iswap],0,world);
  }

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    recvnum_direct[iswap] = sendnum_direct[iswap];
  }
  
  MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);


  // NOTE: set all per-swap header values to correct counts
  // NOTE: be sure to allocate all bufs to sufficient size, using nrecv*size_border


  
  // perform border comm via direct swaps
  // post receives, perform pack+sends, copy to self, wait for and unpack all incoming messages

  nrecv = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (size_forward_recv_direct[iswap]) {
      MPI_Irecv(buf_recv_direct[iswap],size_forward_recv_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],0,world,&requests[nrecv++]);
    }
  }
                    
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (ghost_velocity) {
      n = avec->pack_border_vel(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                              pbc_flag_direct[iswap],pbc_direct[iswap]);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],0,world);
    } else {
      n = avec->pack_border(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                          pbc_flag_direct[iswap],pbc_direct[iswap]);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],0,world);
    }
  }

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    if (ghost_velocity) {
      avec->pack_border_vel(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                          pbc_flag_direct[iswap],pbc_direct[iswap]);
      avec->unpack_border_vel(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    } else {
      avec->pack_border(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                      pbc_flag_direct[iswap],pbc_direct[iswap]);
      avec->unpack_border(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    }
  }

  if (nrecv) {
    if (ghost_velocity) {
      for (int i = 0; i < nrecv; i++) {
        MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
        iswap = recv_indices_direct[irecv];
        avec->unpack_border_vel(recvnum_direct[iswap],firstrecv_direct[iswap],buf_recv_direct[iswap]);
      }
    } else {
      for (int i = 0; i < nrecv; i++) {
        MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
        iswap = recv_indices_direct[irecv];
        avec->unpack_border(recvnum_direct[iswap],firstrecv_direct[iswap],buf_recv_direct[iswap]);
      }
    }
  }

  // set all pointers & counters

  /*
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
  */
  
  // for molecular systems some bits are lost for local atom indices
  //   due to encoding of special pairs in neighbor lists
  // check for overflow

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
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist_direct as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommBrickDirect::grow_list_direct(int iswap, int n)
{
  maxsendlist_direct[iswap] = static_cast<int> (BUFFACTOR * n);
  memory->grow(sendlist_direct[iswap],maxsendlist_direct[iswap],"comm:sendlist_direct[iswap]");
}
