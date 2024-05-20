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
#include "compute.h"
#include "domain.h"
#include "dump.h"
#include "error.h"
#include "fix.h"
#include "memory.h"
#include "neighbor.h"
#include "pair.h"

// NOTES:
// still need cutoff calculation for nonuniform layout
// need forward_comm_array to test molecular systems
// test msg tags with individual procs as multiple neighbors via big stencil
// test when cutoffs >> box length
// test with triclinic
// doc msg tag logic in code
// doc stencil data structs and logic in code
// CommBrick could use local maxsend in its borders() check for sendlist realloc
//   instead of indexing the swap for each atom

using namespace LAMMPS_NS;

static constexpr double BUFFACTOR = 1.5;
static constexpr int BUFMIN = 1024;

/* ---------------------------------------------------------------------- */

CommBrickDirect::CommBrickDirect(LAMMPS *lmp) : CommBrick(lmp)
{
  style = Comm::BRICK_DIRECT;
  init_pointers();
  init_buffers_direct();
}

/* ---------------------------------------------------------------------- */

CommBrickDirect::~CommBrickDirect()
{
  deallocate_direct();
  deallocate_lists(maxlist);

  memory->destroy(buf_send_direct);
  memory->destroy(buf_recv_direct);
}

/* ----------------------------------------------------------------------
   initialize comm pointers to nullptr
------------------------------------------------------------------------- */

void CommBrickDirect::init_pointers()
{
  swaporder = nullptr;
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
  size_border_recv_direct = nullptr;
  swap2list = nullptr;
  sendlist_direct = nullptr;
  firstrecv_direct = nullptr;
  recv_offset_forward_direct = nullptr;
  recv_offset_reverse_direct = nullptr;
  recv_offset_border_direct = nullptr;
  recv_offset_forward_atoms = nullptr;
  recv_offset_reverse_atoms = nullptr;
  requests = nullptr;

  active_list = nullptr;
  check_list = nullptr;
  bounds_list = nullptr;
  sendnum_list = nullptr;
  sendatoms_list = nullptr;
  maxsendatoms_list = nullptr;
}

/* ---------------------------------------------------------------------- */
//IMPORTANT: we *MUST* pass "*oldcomm" to the Comm initializer here, as
//           the code below *requires* that the (implicit) copy constructor
//           for Comm is run and thus creating a shallow copy of "oldcomm".
//           The call to Comm::copy_arrays() then converts the shallow copy
//           into a deep copy of the class with the new layout.

CommBrickDirect::CommBrickDirect(LAMMPS *lmp, Comm *oldcomm) : CommBrick(lmp, oldcomm)
{
  if (oldcomm->layout == Comm::LAYOUT_TILED)
    error->all(FLERR,"Cannot change to comm_style brick/direct from tiled layout");

  style = Comm::BRICK_DIRECT;
  layout = oldcomm->layout;
  Comm::copy_arrays(oldcomm);
  init_pointers();
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
  if (domain->dimension == 2) {
    maxdirect = 8;
    maxlist = 9;
  } else {
    maxdirect = 26;
    maxlist = 27;
  }

  allocate_direct();
  allocate_lists();
}

/* ---------------------------------------------------------------------- */

void CommBrickDirect::init()
{
  CommBrick::init();

  // disallow options not supported by CommBrickDirect

  if (mode == Comm::MULTI || mode == Comm::MULTIOLD)
    error->all(FLERR,
               "Comm brick/direct does not yet support multi or multi/old");

  if (bordergroup)
    error->all(FLERR,
               "Comm brick/direct does not yet support comm_modify group");

  // allocate lists of atoms to send for first time if necessary
  // do now b/c domain->dimension may have changed since construction of this class
  //   if comm_style command specified before dimension command

  int oldlist = maxlist;
  if (domain->dimension == 2) maxlist = 9;
  else maxlist = 27;

  if (maxlist > oldlist) {
    deallocate_lists(oldlist);
    allocate_lists();
  }
}

/* ----------------------------------------------------------------------
   create stencil of direct swaps this procs make with each proc in stencil
   direct swap = send and recv
     same proc can appear multiple times in stencil, self proc can also appear
   stencil is used for border and forward and reverse comm
------------------------------------------------------------------------- */

void CommBrickDirect::setup()
{
  // first perform CommBrick::setup() for 6-way stencil
  // will use its recvneed to create logical 3d stencil of procs around me

  CommBrick::setup();

  // pointers for orthogonal or triclinic domains

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

  // ijk lo/hi = bounds of stencil around this proc at center

  int ilo = -recvneed[0][0];
  int ihi = recvneed[0][1];
  int jlo = -recvneed[1][0];
  int jhi = recvneed[1][1];
  int klo = -recvneed[2][0];
  int khi = recvneed[2][1];

  // calculate max stencil extent in each dim for any proc
  // max stencil extents are used for MPI tag calculations

  int stencil_half_local[3];
  stencil_half_local[0] = MAX(-ilo,ihi);
  stencil_half_local[1] = MAX(-jlo,jhi);
  stencil_half_local[2] = MAX(-klo,khi);

  int stencil_half[3];
  MPI_Allreduce(&stencil_half_local,&stencil_half,3,MPI_INT,MPI_MAX,world);

  int stencil_full[3];
  stencil_full[0] = 2*stencil_half[0] + 1;
  stencil_full[1] = 2*stencil_half[1] + 1;
  stencil_full[2] = 2*stencil_half[2] + 1;

  // ensure max possible tag does not exceed MPI limit

  bigint maxtag = (bigint) stencil_full[0] * stencil_full[1] * stencil_full[2];

  void *maxtag_mpi_ptr;
  int tmp;
  MPI_Comm_get_attr(world,MPI_TAG_UB,&maxtag_mpi_ptr,&tmp);
  int maxtag_mpi = *((int *) maxtag_mpi_ptr);

  if (maxtag > maxtag_mpi)
    error->all(FLERR,"Comm brick/direct stencil is too large");

  // ndirect = # of direct swaps this proc makes with other procs, including self copies
  // subtract 1 for self in center of 3d stencil of surrounding procs

  ndirect = (ihi-ilo+1) * (jhi-jlo+1) * (khi-klo+1) - 1;

  //printf("NDIRECT %d ijk lo/hi %d %d: %d %d: %d %d proc %d\n",
  //       ndirect,ilo,ihi,jlo,jhi,klo,khi,me);

  if (ndirect > maxdirect) {
    deallocate_direct();
    maxdirect = ndirect;
    allocate_direct();
  }

  // create swaporder = ordering of swaps within 3d stencil brick
  // each entry in swaporder is ijk indices of the swap within the brick

  order_swaps(ilo,ihi,jlo,jhi,klo,khi);

  // calculate 6 cutoffs within my subdomain
  // for sending owned atoms to procs on 6 faces of my stencil
  // cutxlo = upper x-cutoff within my subdomain to send owned atoms
  //          to proc on lower x-face of stencil

  if (layout == Comm::LAYOUT_UNIFORM) {
    int nbetween;

    double xwidth = prd[0] / procgrid[0];
    nbetween = -ilo - 1;
    cutxlo = cutghost[0] - nbetween*xwidth;
    nbetween = ihi - 1;
    cutxhi = cutghost[0] - nbetween*xwidth;

    double ywidth = prd[1] / procgrid[1];
    nbetween = -jlo - 1;
    cutylo = cutghost[1] - nbetween*ywidth;
    nbetween = jhi - 1;
    cutyhi = cutghost[1] - nbetween*ywidth;

    if (dim == 3) {
      double zwidth = prd[2] / procgrid[2];
      nbetween = -klo - 1;
      cutzlo = cutghost[2] - nbetween*zwidth;
      nbetween = khi - 1;
      cutzhi = cutghost[2] - nbetween*zwidth;
    }

  } else if (layout == Comm::LAYOUT_NONUNIFORM) {
    // NOTE: still needs to be coded
  }

  // calculate check_list and bounds_list for each list
  // used when building atom lists in borders()
  // unsetting of check_list is when a send list to a proc on a stencil face does not
  //   require a cutoff due to stencil being truncated by a non-PBC boundary

  int ix,iy,iz;

  for (int ilist = 0; ilist < maxlist; ilist++) {
    ix = ilist % 3;
    iy = (ilist/3) % 3;
    iz = ilist / 9;
    if (dim == 2) iz = 1;

    check_list[ilist][0] = 0;
    bounds_list[ilist][0][0] = bounds_list[ilist][0][1] = 0.0;
    if (ix == 0) {
      check_list[ilist][0] = 1;
      bounds_list[ilist][0][0] = sublo[0];
      bounds_list[ilist][0][1] = sublo[0] + cutxlo;
      if (bounds_list[ilist][0][1] >= subhi[0]) check_list[ilist][0] = 0;
    }
    if (ix == 2) {
      check_list[ilist][0] = 1;
      bounds_list[ilist][0][0] = subhi[0] - cutxhi;
      bounds_list[ilist][0][1] = subhi[0];
      if (bounds_list[ilist][0][0] <= sublo[0]) check_list[ilist][0] = 0;
    }

    check_list[ilist][1] = 0;
    bounds_list[ilist][1][0] = bounds_list[ilist][1][1] = 0.0;
    if (iy == 0) {
      check_list[ilist][1] = 1;
      bounds_list[ilist][1][0] = sublo[1];
      bounds_list[ilist][1][1] = sublo[1] + cutylo;
      if (bounds_list[ilist][1][1] >= subhi[1]) check_list[ilist][1] = 0;
    }
    if (iy == 2) {
      check_list[ilist][1] = 1;
      bounds_list[ilist][1][0] = subhi[1] - cutyhi;
      bounds_list[ilist][1][1] = subhi[1];
      if (bounds_list[ilist][1][0] <= sublo[1]) check_list[ilist][1] = 0;
    }

    check_list[ilist][2] = 0;
    bounds_list[ilist][2][0] = bounds_list[ilist][2][1] = 0.0;
    if (iz == 0) {
      check_list[ilist][2] = 1;
      bounds_list[ilist][2][0] = sublo[2];
      bounds_list[ilist][2][1] = sublo[2] + cutzlo;
      if (bounds_list[ilist][2][1] >= subhi[2]) check_list[ilist][2] = 0;
    }
    if (iz == 2) {
      check_list[ilist][2] = 1;
      bounds_list[ilist][2][0] = subhi[2] - cutzhi;
      bounds_list[ilist][2][1] = subhi[2];
      if (bounds_list[ilist][2][0] <= sublo[2]) check_list[ilist][2] = 0;
    }
  }

  // active_list = which lists of atoms are used by swaps
  // zero it before looping over swaps to set it

  for (int ilist = 0; ilist < maxlist; ilist++)
    active_list[ilist] = 0;

  // loop over stencil and define params for each direct swap

  int xpbc,ypbc,zpbc;
  int igrid,jgrid,kgrid;
  int ilistx,ilisty,ilistz,ilist;

  for (int iswap = 0; iswap < ndirect; iswap++) {
    ix = swaporder[iswap][0];
    iy = swaporder[iswap][1];
    iz = swaporder[iswap][2];

    // identify proc to swap with and atom coord PBC shifts required

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

    // identify which atom list this swap uses, based on ix,iy,iz
    // set swap2list and increment active_list for that list
    // put check for ixyz = 0 before lo or hi so that if lo or hi = 0
    //    then there will be no additional atom lists created for lo or hi

    if (ix == 0) ilistx = 1;
    else if (ix == ilo) ilistx = 0;
    else if (ix == ihi) ilistx = 2;
    else ilistx = 1;

    if (iy == 0) ilisty = 1;
    else if (iy == jlo) ilisty = 0;
    else if (iy == jhi) ilisty = 2;
    else ilisty = 1;

    if (iz == 0) ilistz = 1;
    else if (iz == klo) ilistz = 0;
    else if (iz == khi) ilistz = 2;
    else ilistz = 1;

    if (dim == 2) ilist = 3*ilisty + ilistx;
    if (dim == 3) ilist = 9*ilistz + 3*ilisty + ilistx;

    swap2list[iswap] = ilist;
    active_list[ilist]++;

    // set MPI tags based on 3d offset between 2 procs from receiver's perspective
    // this ensures for each swap, MPI Send/Recv on different procs will use same tag
    // necessary if multiple swaps are performed between same 2 procs
    //   so that receiver can identify which swap the received data is for

    sendtag[iswap] = stencil_full[1]*stencil_full[0]*(-iz+stencil_half[2]) +
      stencil_full[0] *(-iy+stencil_half[1]) + (-ix+stencil_half[0]) + 1;
    recvtag[iswap] = stencil_full[1]*stencil_full[0]*(iz+stencil_half[2]) +
      stencil_full[0]*(iy+stencil_half[1]) + (ix+stencil_half[0]) + 1;
  }

  // set nself_direct and self_indices_direct

  nself_direct = 0;
  for (int iswap = 0; iswap < ndirect; iswap++)
    if (proc_direct[iswap] == me) self_indices_direct[nself_direct++] = iswap;
}

/* ----------------------------------------------------------------------
   order the swaps within the 3d stencil of swaps
   swaporder[I][012] = 3 ijk indices within stencil of the Ith swap
   order the swaps by their stencil distance from the center (my proc)
   for ties, the swaps are stored in loop order (x first, y next, z last)
------------------------------------------------------------------------- */

void CommBrickDirect::order_swaps(int ilo, int ihi, int jlo, int jhi, int klo, int khi)
{
  // center of stencil: ix = iy = iz = 0
  // sdist = distance bewteen center of stencil (me) and another stencil proc
  //   sdist = abs(ix) + abs(iy) + abs(iz)
  // maxdistance = max distance of any corner point in stencil from center point
  // ixyz loop can include ceneter pt, b/c distance = 0, so not added to swaporder

  int imax = MAX(-ilo,ihi);
  int jmax = MAX(-jlo,jhi);
  int kmax = MAX(-klo,khi);
  int maxdistance = imax + jmax + kmax;

  int ix,iy,iz;
  int sdist;
  int idirect = 0;

  for (int distance = 1; distance <= maxdistance; distance++) {
    for (iz = klo; iz <= khi; iz++)
      for (iy = jlo; iy <= jhi; iy++)
        for (ix = ilo; ix <= ihi; ix++) {
          sdist = abs(ix) + abs(iy) + abs(iz);
          if (sdist == distance) {
            swaporder[idirect][0] = ix;
            swaporder[idirect][1] = iy;
            swaporder[idirect][2] = iz;
            idirect++;
          }
        }
  }

  if (idirect != ndirect) error->all(FLERR,"Mistake in stencil ordering");

  // move self to front out of place
}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
   exchange owned atoms directly with all neighbor procs,
     not via CommBrick 6-way stencil
------------------------------------------------------------------------- */

void CommBrickDirect::forward_comm(int /*dummy*/)
{
  int n,iswap,irecv;
  AtomVec *avec = atom->avec;
  double **x = atom->x;
  double *buf;

  // post all receives for ghost atoms
  // except for self copies

  int offset;

  int npost = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (comm_x_only) {
      if (recvnum_direct[iswap]) {
        buf = x[firstrecv_direct[iswap]];
        MPI_Irecv(buf,size_forward_recv_direct[iswap],MPI_DOUBLE,
                  proc_direct[iswap],recvtag[iswap],world,&requests[npost++]);
      }
    } else {
      if (recvnum_direct[iswap]) {
        offset = recv_offset_forward_direct[iswap];
        MPI_Irecv(&buf_recv_direct[offset],size_forward_recv_direct[iswap],MPI_DOUBLE,
                  proc_direct[iswap],recvtag[iswap],world,&requests[npost++]);
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

  if (npost == 0) return;

  if (comm_x_only) {
    MPI_Waitall(npost,requests,MPI_STATUS_IGNORE);
  } else if (ghost_velocity) {
    for (int ipost = 0; ipost < npost; ipost++) {
      MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
      iswap = recv_indices_direct[irecv];
      offset = recv_offset_forward_direct[iswap];
      avec->unpack_comm_vel(recvnum_direct[iswap],firstrecv_direct[iswap],&buf_recv_direct[offset]);
    }
  } else {
    for (int ipost = 0; ipost < npost; ipost++) {
      MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
      iswap = recv_indices_direct[irecv];
      offset = recv_offset_forward_direct[iswap];
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
  int n,iswap,irecv;
  AtomVec *avec = atom->avec;
  double **f = atom->f;
  double *buf;

  // post all receives for owned atoms
  // except for self copy/sums

  int offset;

  int npost = 0;
  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap]) {
      offset = recv_offset_reverse_direct[iswap];
      MPI_Irecv(&buf_recv_direct[offset],size_reverse_recv_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],sendtag[iswap],world,&requests[npost++]);
    }
  }

  // send all ghost atoms to receiving procs
  // except for self copy/sums

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
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

  if (npost == 0) return;

  for (int i; i < npost; i++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = send_indices_direct[irecv];
    offset = recv_offset_reverse_direct[iswap];
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
  int i,n,iswap,ilist,nsend,nrecv;

  // setup lists of atoms to send in each direct swap
  // only maxlist possible lists (27 in 3d, 9 in 2d) regardless of stencil size
  // skip non-active lists as flagged in setup()

  AtomVec *avec = atom->avec;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int dim = domain->dimension;

  int allflag;
  int xcheck,ycheck,zcheck;
  double xlo,xhi,ylo,yhi,zlo,zhi;

  for (ilist = 0; ilist < maxlist; ilist++) {
    if (!active_list[ilist]) continue;

    xcheck = check_list[ilist][0];
    ycheck = check_list[ilist][1];
    zcheck = check_list[ilist][2];
    xlo = bounds_list[ilist][0][0];
    xhi = bounds_list[ilist][0][1];
    ylo = bounds_list[ilist][1][0];
    yhi = bounds_list[ilist][1][1];
    zlo = bounds_list[ilist][2][0];
    zhi = bounds_list[ilist][2][1];

    nsend = 0;
    maxsend = maxsendatoms_list[ilist];

    if (!xcheck && !ycheck && !zcheck) {
      for (i = 0; i < nlocal; i++) {
        if (nsend == maxsend) {
          grow_list_direct(ilist,nsend);
          maxsend = maxsendatoms_list[ilist];
        }
        sendatoms_list[ilist][nsend++] = i;
      }
    } else {
      for (i = 0; i < nlocal; i++) {
        if (xcheck && (x[i][0] < xlo || x[i][0] > xhi)) continue;
        if (ycheck && (x[i][1] < ylo || x[i][1] > yhi)) continue;
        if (zcheck && (x[i][2] < zlo || x[i][2] > zhi)) continue;
        if (nsend == maxsend) {
          grow_list_direct(ilist,nsend);
          maxsend = maxsendatoms_list[ilist];
        }
        sendatoms_list[ilist][nsend++] = i;
      }
    }

    sendnum_list[ilist] = nsend;
  }

  // set sendnum_direct and sendlist_direct for all swaps from per-list data

  for (iswap = 0; iswap < ndirect; iswap++) {
    ilist = swap2list[iswap];
    sendnum_direct[iswap] = sendnum_list[ilist];
    sendlist_direct[iswap] = sendatoms_list[ilist];
  }

  // recvnum_direct = number of ghost atoms recieved in each swap
  // acquire by sending value of nsend for each swap to each receiving proc
  // post receives, perform sends, copy to self, wait for all incoming messages

  int npost = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    MPI_Irecv(&recvnum_direct[iswap],1,MPI_INT,
              proc_direct[iswap],recvtag[iswap],world,&requests[npost++]);
  }

  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    MPI_Send(&sendnum_direct[iswap],1,MPI_INT,proc_direct[iswap],sendtag[iswap],world);
  }

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    recvnum_direct[iswap] = sendnum_direct[iswap];
  }

  MPI_Waitall(npost,requests,MPI_STATUS_IGNORE);

  // set nghost = sum of recnum_direct over swaps
  // set firstrecv_direct = index to 1st ghost atom in each swap receive
  // set size_forward_recv_direct and size_reverse_send/recv_direct
  // set send_indices_direct and recv_indices_direct for non-empty swaps with other procs

  int nghost = 0;
  int isend = 0;
  int irecv = 0;
  int offset_forward = 0;
  int offset_reverse = 0;
  int offset_border = 0;
  int offset_forward_atoms = 0;
  int offset_reverse_atoms = 0;

  smax_direct = 0;
  rmax_direct = 0;
  ssum_direct = 0;
  rsum_direct = 0;

  for (iswap = 0; iswap < ndirect; iswap++) {
    nsend = sendnum_direct[iswap];
    nrecv = recvnum_direct[iswap];
    firstrecv_direct[iswap] = nlocal + nghost;
    nghost += nrecv;

    size_forward_recv_direct[iswap] = size_forward * nrecv;
    size_reverse_send_direct[iswap] = size_reverse * nrecv;
    size_reverse_recv_direct[iswap] = size_reverse * nsend;
    size_border_recv_direct[iswap] = size_border * nrecv;

    recv_offset_forward_direct[iswap] = offset_forward;
    recv_offset_reverse_direct[iswap] = offset_reverse;
    recv_offset_border_direct[iswap] = offset_border;
    recv_offset_forward_atoms[iswap] = offset_forward_atoms;
    recv_offset_reverse_atoms[iswap] = offset_reverse_atoms;

    offset_forward += size_forward_recv_direct[iswap];
    offset_reverse += size_reverse_recv_direct[iswap];
    offset_border += size_border * nrecv;
    offset_forward_atoms += nrecv;
    offset_reverse_atoms += nsend;

    if (proc_direct[iswap] != me) {
      if (nsend) send_indices_direct[isend++] = iswap;
      if (nrecv) recv_indices_direct[irecv++] = iswap;
    }

    smax_direct = MAX(smax_direct,nsend);
    rmax_direct = MAX(rmax_direct,nrecv);
    ssum_direct += nsend;
    rsum_direct += nrecv;
  }

  atom->nghost = nghost;

  // ensure send/recv buffers are large enough for all border & forward & reverse comm

  check_buffer_sizes();

  // perform border comm via direct swaps
  // use pack/unpack border and pack/unpack border_vel
  // post receives, perform sends, copy to self, wait for all incoming messages

  int offset;

  npost = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (size_border_recv_direct[iswap]) {
      offset = recv_offset_border_direct[iswap];
      MPI_Irecv(&buf_recv_direct[offset],size_border_recv_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],recvtag[iswap],world,&requests[npost++]);
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

  if (npost) {
    if (ghost_velocity) {
      for (int ipost = 0; ipost < npost; ipost++) {
        MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
        iswap = recv_indices_direct[irecv];
        offset = recv_offset_border_direct[iswap];
        avec->unpack_border_vel(recvnum_direct[iswap],firstrecv_direct[iswap],
                                &buf_recv_direct[offset]);
      }
    } else {
      for (int ipost = 0; ipost < npost; ipost++) {
        MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
        iswap = recv_indices_direct[irecv];
        offset = recv_offset_border_direct[iswap];
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
   forward communication invoked by a Pair
   nsize used to set recv buffer offsets and limits
   NOTE: this could be memory inefficient if Pair style has multiple forward comms ?
         leaving gaps in buf_recv_direct
   SOLUTION: add isize arg like forward_comm() for fixes
------------------------------------------------------------------------- */

void CommBrickDirect::forward_comm(Pair *pair)
{
  int n,iswap,irecv;
  double *buf;

  int nsize = pair->comm_forward;

  // post all receives for ghost atoms
  // except for self copies

  int offset;

  int npost = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (recvnum_direct[iswap]) {
      offset = nsize * recv_offset_forward_atoms[iswap];
      MPI_Irecv(&buf_recv_direct[offset],nsize*recvnum_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],recvtag[iswap],world,&requests[npost++]);
    }
  }

  // send all owned atoms to receiving procs
  // except for self copies

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    n = pair->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                                pbc_flag_direct[iswap],pbc_direct[iswap]);
    if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],sendtag[iswap],world);
  }

  // copy atoms to self via pack and unpack

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    pair->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                            pbc_flag_direct[iswap],pbc_direct[iswap]);
    pair->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
  }

  // wait on incoming messages with ghost atoms
  // unpack each message as it arrives

  if (npost == 0) return;

  for (int ipost = 0; ipost < npost; ipost++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = recv_indices_direct[irecv];
    offset = nsize * recv_offset_forward_atoms[iswap];
    pair->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],&buf_recv_direct[offset]);
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Pair
   nsize used to set recv buffer offsets and limits
------------------------------------------------------------------------- */

void CommBrickDirect::reverse_comm(Pair *pair)
{
  int n,iswap,irecv;
  double *buf;

  int nsize = MAX(pair->comm_reverse,pair->comm_reverse_off);

  // post all receives for owned atoms
  // except for self copy/sums

  int offset;

  int npost = 0;
  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap]) {
      offset = recv_offset_reverse_atoms[iswap];
      MPI_Irecv(&buf_recv_direct[offset],size_reverse_recv_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],sendtag[iswap],world,&requests[npost++]);
    }
  }

  // send all ghost atoms to receiving procs
  // except for self copy/sums

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    n = pair->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],recvtag[iswap],world);
  }

  // copy/sum atoms to self via pack and unpack

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    pair->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    pair->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct);
  }

  // wait on incoming messages with owned atoms
  // unpack each message as it arrives

  if (npost == 0) return;

  for (int i; i < npost; i++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = send_indices_direct[irecv];
    offset = recv_offset_reverse_atoms[iswap];
    pair->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],&buf_recv_direct[offset]);
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Bond
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommBrickDirect::forward_comm(Bond *bond)
{
  error->all(FLERR,"Comm_style brick/direct forward_comm for "
             "bond styles has not yet been implemented");
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Bond
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommBrickDirect::reverse_comm(Bond *bond)
{
  error->all(FLERR,"Comm_style brick/direct reverse_comm for "
             "bond styles has not yet been implemented");
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   size/nsize used only to set recv buffer offsets and limits
   size = 0 (default) -> use comm_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
     some are smaller than max stored in its comm_forward
------------------------------------------------------------------------- */

void CommBrickDirect::forward_comm(Fix *fix, int size)
{
  int n,iswap,irecv;
  double *buf;

  int nsize;
  if (size) nsize = size;
  else nsize = fix->comm_forward;

  // post all receives for ghost atoms
  // except for self copies

  int offset;

  int npost = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (recvnum_direct[iswap]) {
      offset = nsize * recv_offset_forward_atoms[iswap];
      MPI_Irecv(&buf_recv_direct[offset],nsize*recvnum_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],recvtag[iswap],world,&requests[npost++]);
    }
  }

  // send all owned atoms to receiving procs
  // except for self copies

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    n = fix->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                               pbc_flag_direct[iswap],pbc_direct[iswap]);
    if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],sendtag[iswap],world);
  }

  // copy atoms to self via pack and unpack

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    fix->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                           pbc_flag_direct[iswap],pbc_direct[iswap]);
    fix->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
  }

  // wait on incoming messages with ghost atoms
  // unpack each message as it arrives

  if (npost == 0) return;

  for (int ipost = 0; ipost < npost; ipost++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = recv_indices_direct[irecv];
    offset = nsize * recv_offset_forward_atoms[iswap];
    fix->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],&buf_recv_direct[offset]);
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   size/nsize used only to set recv buffer offsets and limits
   size = 0 (default) -> use comm_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
     some are smaller than max stored in its comm_forward
------------------------------------------------------------------------- */

void CommBrickDirect::reverse_comm(Fix *fix, int size)
{
  int n,iswap,irecv;
  double *buf;

  int nsize;
  if (size) nsize = size;
  else nsize = fix->comm_reverse;

  // post all receives for owned atoms
  // except for self copy/sums

  int offset;

  int npost = 0;
  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap]) {
      offset = recv_offset_reverse_atoms[iswap];
      MPI_Irecv(&buf_recv_direct[offset],size_reverse_recv_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],sendtag[iswap],world,&requests[npost++]);
    }
  }

  // send all ghost atoms to receiving procs
  // except for self copy/sums

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    n = fix->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],recvtag[iswap],world);
  }

  // copy/sum atoms to self via pack and unpack

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    fix->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    fix->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct);
  }

  // wait on incoming messages with owned atoms
  // unpack each message as it arrives

  if (npost == 0) return;

  for (int i; i < npost; i++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = send_indices_direct[irecv];
    offset = recv_offset_reverse_atoms[iswap];
    fix->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],&buf_recv_direct[offset]);
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix with variable size data
   query fix for pack size to ensure buf_send is big enough
   handshake sizes before each Irecv/Send to ensure buf_recv is big enough
------------------------------------------------------------------------- */

void CommBrickDirect::reverse_comm_variable(Fix *fix)
{
  error->all(FLERR,"Comm_style brick/direct reverse_comm for "
             "variables has not yet been implemented");
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Compute
   nsize used only to set recv buffer offsets and limits
------------------------------------------------------------------------- */

void CommBrickDirect::forward_comm(Compute *compute)
{
  int n,iswap,irecv;
  double *buf;

  int nsize = compute->comm_forward;

  // post all receives for ghost atoms
  // except for self copies

  int offset;

  int npost = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (recvnum_direct[iswap]) {
      offset = nsize * recv_offset_forward_atoms[iswap];
      MPI_Irecv(&buf_recv_direct[offset],nsize*recvnum_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],recvtag[iswap],world,&requests[npost++]);
    }
  }

  // send all owned atoms to receiving procs
  // except for self copies

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    n = compute->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                                pbc_flag_direct[iswap],pbc_direct[iswap]);
    if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],sendtag[iswap],world);
  }

  // copy atoms to self via pack and unpack

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    compute->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                            pbc_flag_direct[iswap],pbc_direct[iswap]);
    compute->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
  }

  // wait on incoming messages with ghost atoms
  // unpack each message as it arrives

  if (npost == 0) return;

  for (int ipost = 0; ipost < npost; ipost++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = recv_indices_direct[irecv];
    offset = nsize * recv_offset_forward_atoms[iswap];
    compute->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],&buf_recv_direct[offset]);
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Compute
   nsize used only to set recv buffer offsets and limits
------------------------------------------------------------------------- */

void CommBrickDirect::reverse_comm(Compute *compute)
{
  int n,iswap,irecv;
  double *buf;

  int nsize = compute->comm_reverse;

  // post all receives for owned atoms
  // except for self copy/sums

  int offset;

  int npost = 0;
  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap]) {
      offset = recv_offset_reverse_atoms[iswap];
      MPI_Irecv(&buf_recv_direct[offset],size_reverse_recv_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],sendtag[iswap],world,&requests[npost++]);
    }
  }

  // send all ghost atoms to receiving procs
  // except for self copy/sums

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    n = compute->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],recvtag[iswap],world);
  }

  // copy/sum atoms to self via pack and unpack

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    compute->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    compute->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct);
  }

  // wait on incoming messages with owned atoms
  // unpack each message as it arrives

  if (npost == 0) return;

  for (int i; i < npost; i++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = send_indices_direct[irecv];
    offset = recv_offset_reverse_atoms[iswap];
    compute->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],&buf_recv_direct[offset]);
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Dump
   nsize used only to set recv buffer offsets and limits
------------------------------------------------------------------------- */

void CommBrickDirect::forward_comm(Dump *dump)
{
  int n,iswap,irecv;
  double *buf;

  int nsize = dump->comm_forward;

  // post all receives for ghost atoms
  // except for self copies

  int offset;

  int npost = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (recvnum_direct[iswap]) {
      offset = nsize * recv_offset_forward_atoms[iswap];
      MPI_Irecv(&buf_recv_direct[offset],nsize*recvnum_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],recvtag[iswap],world,&requests[npost++]);
    }
  }

  // send all owned atoms to receiving procs
  // except for self copies

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    n = dump->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                                pbc_flag_direct[iswap],pbc_direct[iswap]);
    if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],sendtag[iswap],world);
  }

  // copy atoms to self via pack and unpack

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    dump->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                            pbc_flag_direct[iswap],pbc_direct[iswap]);
    dump->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
  }

  // wait on incoming messages with ghost atoms
  // unpack each message as it arrives

  if (npost == 0) return;

  for (int ipost = 0; ipost < npost; ipost++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = recv_indices_direct[irecv];
    offset = nsize * recv_offset_forward_atoms[iswap];
    dump->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],&buf_recv_direct[offset]);
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Dump
   nsize used only to set recv buffer offsets and limits
------------------------------------------------------------------------- */

void CommBrickDirect::reverse_comm(Dump *dump)
{
  int n,iswap,irecv;
  double *buf;

  int nsize = dump->comm_reverse;

  // post all receives for owned atoms
  // except for self copy/sums

  int offset;

  int npost = 0;
  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap]) {
      offset = recv_offset_reverse_atoms[iswap];
      MPI_Irecv(&buf_recv_direct[offset],size_reverse_recv_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],sendtag[iswap],world,&requests[npost++]);
    }
  }

  // send all ghost atoms to receiving procs
  // except for self copy/sums

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    n = dump->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    if (n) MPI_Send(buf_send_direct,n,MPI_DOUBLE,proc_direct[iswap],recvtag[iswap],world);
  }

  // copy/sum atoms to self via pack and unpack

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    dump->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    dump->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct);
  }

  // wait on incoming messages with owned atoms
  // unpack each message as it arrives

  if (npost == 0) return;

  for (int i; i < npost; i++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = send_indices_direct[irecv];
    offset = recv_offset_reverse_atoms[iswap];
    dump->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],&buf_recv_direct[offset]);
  }
}

/* ----------------------------------------------------------------------
   forward communication of N values in per-atom array
------------------------------------------------------------------------- */

void CommBrickDirect::forward_comm_array(int nsize, double **array)
{
  int i,j,k,m,iswap,irecv,offset,last;
  double *buf;

  // ensure send/recv bufs are big enough for nsize
  // based on smax/rmax/ssum/rsum from most recent borders() invocation

  if (nsize > maxforward) {
    maxforward = nsize;
    check_buffer_sizes();
  }

  // post all receives for ghost atoms
  // except for self copies

  int npost = 0;
  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (recvnum_direct[iswap]) {
      offset = nsize * recv_offset_forward_atoms[iswap];
      MPI_Irecv(&buf_recv_direct[offset],nsize*recvnum_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],recvtag[iswap],world,&requests[npost++]);
    }
  }

  // send all owned atoms to receiving procs
  // except for self copies

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    m = 0;
    for (i = 0; i < sendnum_direct[iswap]; i++) {
      j = sendlist_direct[iswap][i];
      for (k = 0; k < nsize; k++)
        buf_send_direct[m++] = array[j][k];
    }
    if (m) MPI_Send(buf_send_direct,m,MPI_DOUBLE,proc_direct[iswap],sendtag[iswap],world);
  }

  // copy atoms to self via pack and unpack

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    m = firstrecv_direct[iswap];
    for (i = 0; i < sendnum_direct[iswap]; i++) {
      j = sendlist_direct[iswap][i];
      for (k = 0; k < nsize; k++)
        array[m][k] = array[j][k];
      m++;
    }
  }

  // wait on incoming messages with ghost atoms
  // unpack each message as it arrives

  if (npost == 0) return;

  for (int ipost = 0; ipost < npost; ipost++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = recv_indices_direct[irecv];
    offset = nsize * recv_offset_forward_atoms[iswap];
    buf = &buf_recv_direct[offset];

    m = 0;
    last = firstrecv_direct[iswap] + recvnum_direct[iswap];
    for (i = firstrecv_direct[iswap]; i < last; i++) {
      for (k = 0; k < nsize; k++)
        array[i][k] = buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all vectors that depend on maxdirect = size of direct swap stencil
------------------------------------------------------------------------- */

void CommBrickDirect::allocate_direct()
{
  memory->create(swaporder,maxdirect,3,"comm:swaporder");
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
  memory->create(size_border_recv_direct,maxdirect,"comm:size_border_recv_direct");
  memory->create(swap2list,maxdirect,"comm:swap2list");
  sendlist_direct = (int **) memory->smalloc(maxdirect*sizeof(int *),"comm:sendlist_direct");
  memory->create(firstrecv_direct,maxdirect,"comm:recvnum_direct");
  memory->create(recv_offset_forward_direct,maxdirect,"comm:recv_offset_forward_direct");
  memory->create(recv_offset_reverse_direct,maxdirect,"comm:recv_offset_reverse_direct");
  memory->create(recv_offset_border_direct,maxdirect,"comm:recv_offset_border_direct");
  memory->create(recv_offset_forward_atoms,maxdirect,"comm:recv_offset_forward_atoms");
  memory->create(recv_offset_reverse_atoms,maxdirect,"comm:recv_offset_reverse_atoms");
  requests = (MPI_Request *) memory->smalloc(maxdirect*sizeof(MPI_Request),"comm:requests");
}

/* ----------------------------------------------------------------------
   allocate all send lists of atom indices
------------------------------------------------------------------------- */

void CommBrickDirect::allocate_lists()
{
  memory->create(active_list,maxlist,"comm:active_list");
  memory->create(check_list,maxlist,3,"comm:check_list");
  memory->create(bounds_list,maxlist,3,2,"comm:bounds_list");
  memory->create(sendnum_list,maxlist,"comm:sendnum_list");
  memory->create(maxsendatoms_list,maxlist,"comm:maxsendatoms_list");
  sendatoms_list = (int **) memory->smalloc(maxlist*sizeof(int *),"comm:sendatoms_list");
  for (int ilist = 0; ilist < maxlist; ilist++) {
    maxsendatoms_list[ilist] = BUFMIN;
    memory->create(sendatoms_list[ilist],BUFMIN,"comm:sendatoms_list[ilist]");
  }
}

/* ----------------------------------------------------------------------
   deallocate all vectors that depend on maxdirect = size of direct swap stencil
------------------------------------------------------------------------- */

void CommBrickDirect::deallocate_direct()
{
  memory->destroy(swaporder);
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
  memory->destroy(size_border_recv_direct);
  memory->destroy(swap2list);
  memory->sfree(sendlist_direct);
  memory->destroy(firstrecv_direct);
  memory->destroy(recv_offset_forward_direct);
  memory->destroy(recv_offset_reverse_direct);
  memory->destroy(recv_offset_border_direct);
  memory->destroy(recv_offset_forward_atoms);
  memory->destroy(recv_offset_reverse_atoms);
  memory->sfree(requests);
}

/* ----------------------------------------------------------------------
   deallocate all send lists of atom indices
------------------------------------------------------------------------- */

void CommBrickDirect::deallocate_lists(int nlist)
{
  memory->destroy(active_list);
  memory->destroy(check_list);
  memory->destroy(bounds_list);
  memory->destroy(sendnum_list);
  for (int ilist = 0; ilist < nlist; ilist++)
    memory->destroy(sendatoms_list[ilist]);
  memory->sfree(sendatoms_list);
  memory->destroy(maxsendatoms_list);
}

/* ----------------------------------------------------------------------
   ensure send/recv buffers are large enough for all border & forward & reverse comm
   size of send buf is for a single swap
   size of recv buf is for all swaps
------------------------------------------------------------------------- */

void CommBrickDirect::check_buffer_sizes()
{
  int max = size_border * smax_direct;
  max = MAX(max,maxforward*smax_direct);
  max = MAX(max,maxreverse*rmax_direct);
  if (max > maxsend_direct) grow_send_direct(max,0);

  max = size_border * rsum_direct;
  max = MAX(max,maxforward*rsum_direct);
  max = MAX(max,maxreverse*ssum_direct);
  if (max > maxrecv_direct) grow_recv_direct(max);
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
    maxsend_direct = static_cast<int> (BUFFACTOR * n);
    memory->destroy(buf_send_direct);
    memory->create(buf_send_direct,maxsend_direct,"comm:buf_send_direct");
  } else if (flag == 1) {
    maxsend_direct = static_cast<int> (BUFFACTOR * n);
    memory->grow(buf_send_direct,maxsend_direct,"comm:buf_send_direct");
  } else {
    memory->destroy(buf_send_direct);
    memory->grow(buf_send_direct,maxsend_direct,"comm:buf_send_direct");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv_direct buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommBrickDirect::grow_recv_direct(int n)
{
  maxrecv_direct = static_cast<int> (BUFFACTOR * n);
  memory->destroy(buf_recv_direct);
  memory->create(buf_recv_direct,maxrecv_direct,"comm:buf_recv_direct");
}

/* ----------------------------------------------------------------------
   realloc the size of the ilist entry in sendatoms_list as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommBrickDirect::grow_list_direct(int ilist, int n)
{
  maxsendatoms_list[ilist] = static_cast<int> (BUFFACTOR * n);
  memory->grow(sendatoms_list[ilist],maxsendatoms_list[ilist],"comm:sendatoms_list[ilist]");
}
