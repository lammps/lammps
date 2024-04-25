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
// make sure all data structs are set by end of borders()
// allocate requests to length of nrecv_direct
// what do lengths of send/recv bufs need to be
// do not allow MULTI with brick/direct or bordergroup
// how to order dswap by shells within full stencil
// test msg tags with individual procs as muliple neighbors
// doc msg tags
// test for cutoffs >> box length
// change 512, 1025 to defines, error check 512 not exceeded
// reorg atom lists to just have 8 unique for 2d and 26 for 3d plus "all" list
// add an init() to check for some disallowed optons ?
// for outer shell of stencil procs, need to compute send proc cutoff
//    can do this for each of 6 directions, use xyzsplit for nonuniform bricks
//    for orthongonal or triclinic
// decide on PBC and tags being in struct or separate vecs

using namespace LAMMPS_NS;

static constexpr double BUFFACTOR = 1.5;
static constexpr int BUFMIN = 1024;
static constexpr int STENCIL_HALF = 512;
static constexpr int STENCIL_FULL = 1025;

/* ---------------------------------------------------------------------- */

CommBrickDirect::CommBrickDirect(LAMMPS *lmp) : CommBrick(lmp)
{
  style = Comm::BRICK_DIRECT;

  dswap = nullptr;
  requests = nullptr;
  maxdirect = 0;

  proc_direct = nullptr;
  pbc_flag_direct = nullptr;
  pbc_direct = nullptr;
  
  self_indices_direct = nullptr;
}

/* ---------------------------------------------------------------------- */

CommBrickDirect::~CommBrickDirect()
{
  delete [] dswap;
  delete [] requests;
  delete [] self_indices_direct;
}

/* ---------------------------------------------------------------------- */

CommBrickDirect::CommBrickDirect(LAMMPS *lmp, Comm *oldcomm) : CommBrick(lmp, oldcomm)
{
  if (oldcomm->layout == Comm::LAYOUT_TILED)
    error->all(FLERR,"Cannot change to comm_style brick/direct from tiled layout");

  style = Comm::BRICK_DIRECT;
}

/* ----------------------------------------------------------------------
   first perform CommBrick::setup() for 6-way stencil
   use its recvneed to create logical 3d grid of procs with me in center
   this is the stencil of direct swaps I make with each proc in stencil
   this proc will perform direct swaps (send and recv) with each proc in stencil
     same proc can appear multiple times in stencil, self proc can also appear
   create dswap = list of DirectSwap data structs
     same stencil and dswap info is used for both forward and reverse comm
------------------------------------------------------------------------- */

void CommBrickDirect::setup()
{
  CommBrick::setup();

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
  
  // ijk lo/hi = bounds of stencil around my proc at center
  // ndirect = # of direct swaps with other procs, including self copies
  // subtract 1 for myself in center of 3d grid of surrounding procs
  
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
    delete [] proc_direct;
    proc_direct = new int[ndirect];
    delete [] pbc_flag_direct;
    pbc_flag_direct = new int[ndirect];
    memory->destroy(pbc_direct);
    memory->create(pbc_direct,ndirect,3,"comm:pbc_direct");
  }

  // loop over stencil and define each direct swap
  
  int ix,iy,iz;
  int igrid,jgrid,kgrid;
  int xpbc,ypbc,zpbc;
  DirectSwap *ds;
  
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

        proc_direct[iswap] = grid2proc[igrid][jgrid][kgrid];

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

        pbc_flag[iswap] = 0;
        pbc[iswap][0] = pbc[iswap][1] = pbc[iswap][2] =
          pbc[iswap][3] = pbc[iswap][4] = pbc[iswap][5] = 0;
        if (xpbc || !ypbc || zpbc) {
          pbc[iswap][0] = xpbc;
          pbc[iswap][1] = ypbc;
          pbc[iswap][2] = zpbc;
          if (triclinic) {
            pbc[iswap][5] = pbc[iswap][1];
            pbc[iswap][4] = pbc[iswap][3] = pbc[iswap][2];
          }
        }

        // MPI tag is based on 3d offset between 2 procs from receiver's perspective
        
        sendtag[iswap] = STENCIL_FULL*STENCIL_FULL*(-klo+STENCIL_HALF) +
          STENCIL_FULL*(-jlo+STENCIL_HALF) + (-ilo+STENCIL_HALF);
        recvtag[iswap] = STENCIL_FULL*STENCIL_FULL*(klo+STENCIL_HALF) +
          STENCIL_FULL*(jlo+STENCIL_HALF) + (ilo+STENCIL_HALF);
        
        iswap++;
      }
    }
  }

  ndirect = iswap;

  // set nself, self_indices_direct

  nself_direct = 0;
  for (iswap = 0; iswap < ndirect; iswap++)
    if (dswap[iswap].proc == me) nself_direct++;

  delete [] self_indices_direct;
  self_indices_direct = new int[nself_direct];

  nself_direct = 0;
  for (iswap = 0; iswap < ndirect; iswap++)
    if (dswap[iswap].proc == me) self_indices_direct[nself_direct++] = iswap;
  
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
                  proc_direct[iswap],recvtag[iswap],world,&requests[nrecv++]);
      }
    } else {
      if (size_forward_recv_direct[iswap]) {
        MPI_Irecv(buf_recv_direct[iswap],size_forward_recv_direct[iswap],MPI_DOUBLE,
                  proc_direct[iswap],recvtag[iswap],world,&requests[nrecv++]);
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
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],sendtag[iswap],world);
    } else {
      n = avec->pack_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                          pbc_flag_direct[iswap],pbc_direct[iswap]);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],sendtag[iswap],world);
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
                proc_direct[iswap],sendtag[iswap],world,&requests[nrecv++]);
  }

  // send all ghost atoms to receiving procs
  // except for self copy/sums

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (recvproc[iswap] == me) continue;
    if (comm_f_only) {
      if (size_reverse_send_direct[iswap]) {
        buf = f[firstrecv_direct[iswap]];
        MPI_Send(buf,size_reverse_send_direct[iswap],MPI_DOUBLE,
                 proc_direct[iswap],recvtag[iswap],world);
      }
    } else {
      n = avec->pack_reverse(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],recvtag[iswap],world);
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

  nrecv = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    MPI_Irecv(&recvnum_direct[iswap],1,MPI_INT,
              proc_direct[iswap],ds->recvtag,world,&requests[nrecv++]);
  }

  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    MPI_Send(&sendnum_direct[iswap],1,MPI_INT,proc_direct[iswap],ds->sendtag,world);
  }

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    recvnum_direct[iswap] = sendnum_direct[iswap];
  }
  
  MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);

  // perform border comm via direct swaps
  // post receives, perform pack+sends, copy to self, wait for and unpack all incoming messages

  nrecv = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (size_forward_recv_direct[iswap]) {
      MPI_Irecv(buf_recv_direct[iswap],size_forward_recv_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],ds->recvtag,world,&requests[nrecv++]);
    }
  }
                    
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (ghost_velocity) {
      n = avec->pack_border_vel(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                              pbc_flag_direct[iswap],pbc_direct[iswap]);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],ds->sendtag,world);
    } else {
      n = avec->pack_border(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                          pbc_flag_direct[iswap],pbc_direct[iswap]);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],ds->sendtag,world);
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
  atom->nghost += nrecv;
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
