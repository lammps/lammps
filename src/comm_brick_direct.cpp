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
// do not allow MULTI with brick/direct or bordergroup
// how to order dswap by shells within full stencil
// test msg tags with individual procs as muliple neighbors
// doc msg tag logic in code
// test for cutoffs >> box length
// error check 512 not exceeded for tags and stencil
// reorg atom lists to just have 8 unique for 2d and 26 for 3d plus "all" list
// add an init() to check for some disallowed options ?
// for outer shell of stencil procs, need to compute send proc cutoff
//    can do this for each of 6 directions, use xyzsplit for nonuniform bricks
//    for orthongonal or triclinic

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
  proc_direct = nullptr;
  pbc_flag_direct = nullptr;
  pbc_direct = nullptr;
  sendtag = nullptr;
  recvtag = nullptr;
  send_indices_direct = nullptr;
  recv_indices_direct = nullptr;
  self_indices_direct = nullptr;
  sendnum_direct = nullptr;
  recvnum_direct = nullptr;
  size_forward_recv_direct = nullptr;
  size_reverse_send_direct = nullptr;
  size_reverse_recv_direct = nullptr;
  firstrecv_direct = nullptr;
  maxsendlist_direct = nullptr;
  sendlist_direct = nullptr;
  recv_offset_forward = nullptr;
  recv_offset_reverse = nullptr;
  recv_offset_border = nullptr;

  ndirect = maxdirect = 0;
  
  init_buffers_direct();
}

/* ---------------------------------------------------------------------- */

CommBrickDirect::~CommBrickDirect()
{
  deallocate_direct();

  memory->destroy(buf_send_direct);
  memory->destroy(buf_recv_direct);
}

/* ---------------------------------------------------------------------- */

CommBrickDirect::CommBrickDirect(LAMMPS *lmp, Comm *oldcomm) : CommBrick(lmp, oldcomm)
{
  if (oldcomm->layout == Comm::LAYOUT_TILED)
    error->all(FLERR,"Cannot change to comm_style brick/direct from tiled layout");

  style = Comm::BRICK_DIRECT;
  layout = oldcomm->layout;
  Comm::copy_arrays(oldcomm);
  init_buffers_direct();
}

/* ----------------------------------------------------------------------
   initialize comm buffers and other data structs local to CommBrickDirect
------------------------------------------------------------------------- */

void CommBrickDirect::init_buffers_direct()
{
  buf_send_direct = buf_recv_direct = nullptr;
  maxsend_direct = maxrecv_direct = BUFMIN;
  grow_send_direct(maxsend_direct,2);
  memory->create(buf_recv_direct,maxrecv_direct,"comm:buf_recv_direct");

  ndirect = 0;
  maxdirect = 26;
  allocate_direct();
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
    deallocate_direct();
    maxdirect = ndirect;
    allocate_direct();
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

        pbc_flag_direct[iswap] = 0;
        pbc_direct[iswap][0] = pbc_direct[iswap][1] = pbc_direct[iswap][2] =
          pbc_direct[iswap][3] = pbc_direct[iswap][4] = pbc_direct[iswap][5] = 0;
        
        if (xpbc || ypbc || zpbc) {
          pbc_flag_direct[iswap] = 1;
          pbc_direct[iswap][0] = xpbc;
          pbc_direct[iswap][1] = ypbc;
          pbc_direct[iswap][2] = zpbc;
          if (triclinic) {
            pbc_direct[iswap][5] = pbc_direct[iswap][1];
            pbc_direct[iswap][4] = pbc_direct[iswap][3] = pbc_direct[iswap][2];
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

  // set nself_direct and self_indices_direst

  nself_direct = 0;
  for (iswap = 0; iswap < ndirect; iswap++)
    if (proc_direct[iswap] == me) self_indices_direct[nself_direct++] = iswap;
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

  int offset;
  
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
        offset = recv_offset_forward[iswap];
        MPI_Irecv(&buf_recv_direct[offset],size_forward_recv_direct[iswap],MPI_DOUBLE,
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
      offset = recv_offset_forward[iswap];
      avec->unpack_comm_vel(recvnum_direct[iswap],firstrecv_direct[iswap],&buf_recv_direct[offset]);
    }
  } else {
    for (int i = 0; i < nrecv; i++) {
      MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
      iswap = recv_indices_direct[irecv];
      offset = recv_offset_forward[iswap];
      avec->unpack_comm(recvnum_direct[iswap],firstrecv_direct[iswap],&buf_recv_direct[offset]);
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

  int offset;
  
  nrecv = 0;
  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (recvproc[iswap] == me) continue;
    if (size_reverse_recv_direct[iswap]) {
      offset = recv_offset_forward[iswap];
      MPI_Irecv(&buf_recv_direct[offset],size_reverse_recv_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],sendtag[iswap],world,&requests[nrecv++]);
    }
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
    offset = recv_offset_reverse[iswap];
    avec->unpack_reverse(sendnum_direct[iswap],sendlist_direct[iswap],&buf_recv_direct[offset]);
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
  int i,n,iswap,isend,irecv,nsend,nrecv;
  
  AtomVec *avec = atom->avec;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  // setup lists of atoms to send in each direct swap

  DirectSwap *ds;
  int allflag,xcheck,ycheck,zcheck;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  
  for (iswap = 0; iswap < ndirect; iswap++) {
    ds = &dswap[iswap];
    nsend = 0;
    maxsend = maxsendlist_direct[iswap];
      
    allflag = ds->allflag;

    if (allflag) {
      for (i = 0; i < nlocal; i++) {
        if (nsend == maxsend) {
          grow_list_direct(iswap,nsend);
          maxsend = maxsendlist_direct[iswap];
        }
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
        if (nsend == maxsend) {
          grow_list_direct(iswap,nsend);
          maxsend = maxsendlist_direct[iswap];
        }
        sendlist_direct[iswap][nsend++] = i;
      }
    }

    sendnum_direct[iswap] = nsend;
  }

  // recvnum_direct = number of ghost atoms recieved in each swap
  // acquire by sending value of nsend for each swap to each receiving proc
  // post receives, perform sends, copy to self, wait for all incoming messages

  nrecv = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    MPI_Irecv(&recvnum_direct[iswap],1,MPI_INT,
              proc_direct[iswap],recvtag[iswap],world,&requests[nrecv++]);
  }

  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    MPI_Send(&sendnum_direct[iswap],1,MPI_INT,proc_direct[iswap],sendtag[iswap],world);
  }

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    recvnum_direct[iswap] = sendnum_direct[iswap];
  }
  
  MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);

  // set nghost = sum of recnum_direct over swaps
  // set firstrecv_direct = index to 1st ghost atom in each swap receive
  // set size_forward_recv_direct and size_reverse_send/recv_direct
  // set send_indices_direct and recv_indices_direct for non-empty swaps with other procs

  int nghost = 0;
  isend = irecv = 0;
  int offset_forward = 0;
  int offset_reverse = 0;
  int offset_border = 0;
  int smax = 0;
  int rmax = 0;
  int ssum = 0;
  int rsum = 0;
  
  for (iswap = 0; iswap < ndirect; iswap++) {
    nsend = sendnum_direct[iswap];
    nrecv = sendnum_direct[iswap];
    firstrecv_direct[iswap] = nlocal + nghost;
    nghost += nrecv;
    
    size_forward_recv_direct[iswap] = size_forward * nrecv;
    size_reverse_send_direct[iswap] = size_reverse * nrecv;
    size_reverse_recv_direct[iswap] = size_reverse * nsend;

    recv_offset_forward[iswap] = offset_forward;
    recv_offset_reverse[iswap] = offset_reverse;
    recv_offset_border[iswap] = offset_border;
    offset_forward += size_forward_recv_direct[iswap];
    offset_reverse += size_reverse_recv_direct[iswap];
    offset_border += size_border * nrecv;
    
    if (nsend) send_indices_direct[isend++] = iswap;
    if (nrecv) recv_indices_direct[irecv++] = iswap;
    smax = MAX(smax,nsend);
    rmax = MAX(rmax,nrecv);
    ssum += nsend;
    rsum += nrecv;
  }

  atom->nghost = nghost;

  // ensure send/recv buffers are large enough for all border & forward & reverse comm
  // size of send buf is for a single swap
  // size of recv buf is for all swaps
  
  int max = size_border * smax;
  max = MAX(max,maxforward*smax);
  max = MAX(max,maxreverse*rmax);
  if (max > maxsend_direct) grow_send_direct(max,0);

  max = size_border * rsum;
  max = MAX(max,maxforward*rsum);
  max = MAX(max,maxreverse*ssum);
  if (max > maxrecv_direct) grow_recv_direct(max);
  
  // perform border comm via direct swaps
  // use pack/unpack border and pack/unpack border_vel
  // post receives, perform sends, copy to self, wait for all incoming messages

  int offset;
  
  nrecv = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (size_forward_recv_direct[iswap]) {
      offset = recv_offset_border[iswap];
      MPI_Irecv(&buf_recv_direct[offset],size_forward_recv_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],recvtag[iswap],world,&requests[nrecv++]);
    }
  }
                    
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (ghost_velocity) {
      n = avec->pack_border_vel(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                              pbc_flag_direct[iswap],pbc_direct[iswap]);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],sendtag[iswap],world);
    } else {
      n = avec->pack_border(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                          pbc_flag_direct[iswap],pbc_direct[iswap]);
      if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],sendtag[iswap],world);
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
        offset = recv_offset_border[iswap];
        avec->unpack_border_vel(recvnum_direct[iswap],firstrecv_direct[iswap],
                                &buf_recv_direct[offset]);
      }
    } else {
      for (int i = 0; i < nrecv; i++) {
        MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
        iswap = recv_indices_direct[irecv];
        offset = recv_offset_border[iswap];
        avec->unpack_border(recvnum_direct[iswap],firstrecv_direct[iswap],
                            &buf_recv_direct[offset]);
      }
    }
  }
  
  // for molecular systems some bits are lost for local atom indices
  //   due to encoding of special pairs in neighbor lists
  // check for overflow

  if ((atom->molecular != Atom::ATOMIC)
      && ((atom->nlocal + atom->nghost) > NEIGHMASK))
    error->one(FLERR,"Per-processor number of atoms is too large for "
               "molecular neighbor lists");

  // reset global->local map

  if (map_style != Atom::MAP_NONE) atom->map_set();
}

/* ----------------------------------------------------------------------
   allocate all vectors that depend on maxdirect = size of direct swap stencil
------------------------------------------------------------------------- */

void CommBrickDirect::allocate_direct()
{
  dswap = (DirectSwap *) memory->smalloc(maxdirect*sizeof(DirectSwap),"comm:dswap");
  requests = (MPI_Request *) memory->smalloc(maxdirect*sizeof(MPI_Request),"comm:requests");
  memory->create(proc_direct,maxdirect,"comm:proc_direct");
  memory->create(pbc_flag_direct,maxdirect,"comm:pbc_flag_direct");
  memory->create(pbc_direct,maxdirect,6,"comm:pbc_direct");
  memory->create(sendtag,maxdirect,"comm:sendtag");
  memory->create(recvtag,maxdirect,"comm:recvtag");
  memory->create(send_indices_direct,maxdirect,"comm:self_indices_direct");
  memory->create(recv_indices_direct,maxdirect,"comm:self_indices_direct");
  memory->create(self_indices_direct,maxdirect,"comm:self_indices_direct");
  memory->create(sendnum_direct,maxdirect,"comm:sendnum_direct");
  memory->create(recvnum_direct,maxdirect,"comm:recvnum_direct");
  memory->create(size_forward_recv_direct,maxdirect,"comm:size_forward_recv_direct");
  memory->create(size_reverse_send_direct,maxdirect,"comm:size_reverse_send_direct");
  memory->create(size_reverse_recv_direct,maxdirect,"comm:size_reverse_recv_direct");
  memory->create(firstrecv_direct,maxdirect,"comm:recvnum_direct");
  memory->create(recv_offset_forward,maxdirect,"comm:recv_offset_forward");
  memory->create(recv_offset_reverse,maxdirect,"comm:recv_offset_reverse");
  memory->create(recv_offset_border,maxdirect,"comm:recv_offset_border");

  memory->create(maxsendlist_direct,maxdirect,"comm:maxsendlist_direct");
  sendlist_direct = (int **) memory->smalloc(maxdirect*sizeof(int *),"comm:sendlist_direct");
  for (int iswap = 0; iswap < maxdirect; iswap++) {
    maxsendlist_direct[iswap] = BUFMIN;
    memory->create(sendlist_direct[iswap],BUFMIN,"comm:sendlist_direct[iswap]");
  }
}

/* ----------------------------------------------------------------------
   deallocate all vectors that depend on maxdirect = size of direct swap stencil
------------------------------------------------------------------------- */

void CommBrickDirect::deallocate_direct()
{
  memory->sfree(dswap);
  memory->sfree(requests);
  memory->destroy(proc_direct);
  memory->destroy(pbc_flag_direct);
  memory->destroy(pbc_direct);
  memory->destroy(sendtag);
  memory->destroy(recvtag);
  memory->destroy(send_indices_direct);
  memory->destroy(recv_indices_direct);
  memory->destroy(self_indices_direct);
  memory->destroy(sendnum_direct);
  memory->destroy(recvnum_direct);
  memory->destroy(size_forward_recv_direct);
  memory->destroy(size_reverse_send_direct);
  memory->destroy(size_reverse_recv_direct);
  memory->destroy(firstrecv_direct);
  memory->destroy(recv_offset_forward);
  memory->destroy(recv_offset_reverse);
  memory->destroy(recv_offset_border);

  memory->destroy(maxsendlist_direct);
  for (int iswap = 0; iswap < maxdirect; iswap++)
    memory->destroy(sendlist_direct[iswap]);
  memory->sfree(sendlist_direct);
}

/* ----------------------------------------------------------------------
   realloc the size of the send_direct buffer as needed with BUFFACTOR
   do not use bufextra as in CommBrick, b/c not using buf_send_direct for exchange()
   flag = 0, don't need to realloc with copy, just free/malloc w/ BUFFACTOR
   flag = 1, realloc with BUFFACTOR
   flag = 2, free/malloc w/out BUFFACTOR
------------------------------------------------------------------------- */

void CommBrickDirect::grow_send_direct(int n, int flag)
{
  if (flag == 0) {
    maxsend = static_cast<int> (BUFFACTOR * n);
    memory->destroy(buf_send_direct);
    memory->create(buf_send_direct,maxsend,"comm:buf_send_direct");
  } else if (flag == 1) {
    maxsend = static_cast<int> (BUFFACTOR * n);
    memory->grow(buf_send_direct,maxsend,"comm:buf_send_direct");
  } else {
    memory->destroy(buf_send_direct);
    memory->grow(buf_send_direct,maxsend,"comm:buf_send_direct");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv_direct buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommBrickDirect::grow_recv_direct(int n)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  memory->destroy(buf_recv_direct);
  memory->create(buf_recv_direct,maxrecv,"comm:buf_recv_direct");
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist_direct as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommBrickDirect::grow_list_direct(int iswap, int n)
{
  maxsendlist_direct[iswap] = static_cast<int> (BUFFACTOR * n);
  memory->grow(sendlist_direct[iswap],maxsendlist_direct[iswap],"comm:sendlist_direct[iswap]");
}
