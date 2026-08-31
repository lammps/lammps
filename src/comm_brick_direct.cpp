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
//   until then init() rejects a non-uniform processor grid
// still need forward/reverse comm invoked by a Bond, and variable-size
//   reverse comm invoked by a Fix
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
  send_requests = nullptr;

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
}

/* ----------------------------------------------------------------------
   initialize comm buffers and other data structs local to CommBrickDirect
------------------------------------------------------------------------- */

void CommBrickDirect::init_buffers_direct()
{
  buf_send_direct = buf_recv_direct = nullptr;
  maxsend_direct = maxrecv_direct = BUFMIN;
  grow_send_direct(maxsend_direct,2);
  grow_recv_direct(maxrecv_direct);

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

  // allocate once: init() runs before every run, and init_buffers_direct()
  //   drops buf_send_direct/buf_recv_direct and calls allocate_direct() and
  //   allocate_lists() over pointers that are already live, so calling it
  //   again leaks all of them.  the oldcomm constructor does not call it, so
  //   the first init() is where those allocations happen for comm_style.

  if (!sendatoms_list) init_buffers_direct();

  // disallow options not supported by CommBrickDirect

  if (mode == Comm::MULTI)
    error->all(FLERR,"Comm brick/direct does not yet support comm_modify mode multi");

  if (bordergroup)
    error->all(FLERR,
               "Comm brick/direct does not yet support comm_modify group");

  // setup() only computes the 6 subdomain cutoffs for a uniform decomposition,
  //   so anything else (e.g. after a balance command) would build the send
  //   lists from cutoffs that were never set

  if (layout != Comm::LAYOUT_UNIFORM)
    error->all(FLERR,"Comm brick/direct requires a uniform processor grid");

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

  // copy atoms to self via pack and unpack
  // done before the sends below, b/c those use buf_send_direct from offset 0
  //   and their data must stay intact until the sends complete

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

  // send all owned atoms to receiving procs
  // except for self copies
  // each swap packs into its own region of buf_send_direct and is sent
  //   with a non-blocking send, so that packing a swap does not have to
  //   wait on the previous swap's message to be picked up by its receiver

  int nsendpost = 0;
  int send_offset = 0;

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap] == 0) continue;
    if (ghost_velocity) {
      n = avec->pack_comm_vel(sendnum_direct[iswap],sendlist_direct[iswap],
                              &buf_send_direct[send_offset],
                              pbc_flag_direct[iswap],pbc_direct[iswap]);
    } else {
      n = avec->pack_comm(sendnum_direct[iswap],sendlist_direct[iswap],
                          &buf_send_direct[send_offset],
                          pbc_flag_direct[iswap],pbc_direct[iswap]);
    }
    if (n) {
      MPI_Isend(&buf_send_direct[send_offset],n,MPI_DOUBLE,proc_direct[iswap],
                sendtag[iswap],world,&send_requests[nsendpost++]);
      send_offset += n;
    }
  }

  // wait on incoming messages with ghost atoms
  // unpack each message as it arrives

  if (npost == 0) {
    if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
    return;
  }

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

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
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
  // comm_f_only sends straight out of the ghost region of f, which the
  //   unpack below never touches (it sums into owned atoms), so the sends
  //   can be left in flight; otherwise each swap packs into its own region
  //   of buf_send_direct so that packing never waits on the previous send

  // copy/sum atoms to self first, since the packed path below uses
  //   buf_send_direct from offset 0 and its data must stay intact until
  //   the sends complete.  this only reads the ghost region of f and sums
  //   into owned atoms, so it does not disturb the sends either

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

  int nsendpost = 0;
  int send_offset = 0;

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (comm_f_only) {
      if (size_reverse_send_direct[iswap]) {
        buf = f[firstrecv_direct[iswap]];
        MPI_Isend(buf,size_reverse_send_direct[iswap],MPI_DOUBLE,
                  proc_direct[iswap],recvtag[iswap],world,&send_requests[nsendpost++]);
      }
    } else {
      n = avec->pack_reverse(recvnum_direct[iswap],firstrecv_direct[iswap],
                             &buf_send_direct[send_offset]);
      if (n) {
        MPI_Isend(&buf_send_direct[send_offset],n,MPI_DOUBLE,proc_direct[iswap],
                  recvtag[iswap],world,&send_requests[nsendpost++]);
        send_offset += n;
      }
    }
  }

  // wait on incoming messages with owned atoms
  // unpack each message as it arrives

  if (npost == 0) {
    if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
    return;
  }

  for (int i = 0; i < npost; i++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = send_indices_direct[irecv];
    offset = recv_offset_reverse_direct[iswap];
    avec->unpack_reverse(sendnum_direct[iswap],sendlist_direct[iswap],&buf_recv_direct[offset]);
  }

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
}

/* ----------------------------------------------------------------------
   build the list of owned atoms to send in each swap
------------------------------------------------------------------------- */

void CommBrickDirect::build_lists()
{
  int i,ilist,nsend;

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int dim = domain->dimension;

  // build every send list in a single pass over the owned atoms
  // a list is identified by its per-dim state: 0 = lo slab, 1 = unrestricted,
  //   2 = hi slab.  check_list/bounds_list for a dim depend only on that dim's
  //   state, so the bounds tests are hoisted out of the per-list loop and each
  //   atom's coords are read once instead of once per list

  int dcheck[3][3];
  double dlo[3][3],dhi[3][3];

  for (int idim = 0; idim < 3; idim++)
    for (int istate = 0; istate < 3; istate++) {
      dcheck[idim][istate] = 0;
      dlo[idim][istate] = dhi[idim][istate] = 0.0;
    }

  int nactive = 0;
  int alist[27],astate[27][3];
  int center_active = 0;

  for (ilist = 0; ilist < maxlist; ilist++) {
    if (!active_list[ilist]) continue;

    int state[3];
    state[0] = ilist % 3;
    state[1] = (ilist / 3) % 3;
    state[2] = (dim == 3) ? ilist / 9 : 1;

    alist[nactive] = ilist;
    astate[nactive][0] = state[0];
    astate[nactive][1] = state[1];
    astate[nactive][2] = state[2];
    nactive++;

    for (int idim = 0; idim < 3; idim++) {
      dcheck[idim][state[idim]] = check_list[ilist][idim];
      dlo[idim][state[idim]] = bounds_list[ilist][idim][0];
      dhi[idim][state[idim]] = bounds_list[ilist][idim][1];
    }

    // a stencil deeper than one proc can map an interior swap to the
    //   unrestricted list in every dim, which then holds all owned atoms

    if ((state[0] == 1) && (state[1] == 1) && (state[2] == 1)) center_active = 1;

    sendnum_list[ilist] = 0;
  }

  int q[3][3];
  q[0][1] = q[1][1] = q[2][1] = 1;

  for (i = 0; i < nlocal; i++) {
    const double xi = x[i][0];
    const double yi = x[i][1];
    const double zi = x[i][2];

    q[0][0] = !dcheck[0][0] || ((xi >= dlo[0][0]) && (xi <= dhi[0][0]));
    q[0][2] = !dcheck[0][2] || ((xi >= dlo[0][2]) && (xi <= dhi[0][2]));
    q[1][0] = !dcheck[1][0] || ((yi >= dlo[1][0]) && (yi <= dhi[1][0]));
    q[1][2] = !dcheck[1][2] || ((yi >= dlo[1][2]) && (yi <= dhi[1][2]));
    q[2][0] = !dcheck[2][0] || ((zi >= dlo[2][0]) && (zi <= dhi[2][0]));
    q[2][2] = !dcheck[2][2] || ((zi >= dlo[2][2]) && (zi <= dhi[2][2]));

    // an atom in no lo/hi slab in any dim can only belong to the unrestricted
    //   list, so skip it outright unless that list is in use

    if (!center_active && !q[0][0] && !q[0][2] && !q[1][0] &&
        !q[1][2] && !q[2][0] && !q[2][2]) continue;

    for (int a = 0; a < nactive; a++) {
      if (!q[0][astate[a][0]]) continue;
      if (!q[1][astate[a][1]]) continue;
      if (!q[2][astate[a][2]]) continue;
      ilist = alist[a];
      nsend = sendnum_list[ilist];
      if (nsend == maxsendatoms_list[ilist]) grow_list_direct(ilist,nsend);
      sendatoms_list[ilist][nsend] = i;
      sendnum_list[ilist] = nsend + 1;
    }
  }
}

/* ----------------------------------------------------------------------
   exchange border data for every swap
------------------------------------------------------------------------- */

void CommBrickDirect::borders_comm()
{
  int n,iswap,irecv,npost;
  AtomVec *avec = atom->avec;

  // perform border comm via direct swaps
  // use pack/unpack border and pack/unpack border_vel
  // post receives, copy to self, perform sends, wait for all incoming messages

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

  // copy to self first, since the sends below use buf_send_direct from
  //   offset 0 and their data must stay intact until they complete

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

  // each swap packs into its own region of buf_send_direct and is sent with
  //   a non-blocking send, so packing a swap does not wait on the previous
  //   swap's message being picked up by its receiver

  int nsendpost = 0;
  int send_offset = 0;

  for (iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap] == 0) continue;
    if (ghost_velocity) {
      n = avec->pack_border_vel(sendnum_direct[iswap],sendlist_direct[iswap],
                                &buf_send_direct[send_offset],
                                pbc_flag_direct[iswap],pbc_direct[iswap]);
    } else {
      n = avec->pack_border(sendnum_direct[iswap],sendlist_direct[iswap],
                            &buf_send_direct[send_offset],
                            pbc_flag_direct[iswap],pbc_direct[iswap]);
    }
    if (n) {
      MPI_Isend(&buf_send_direct[send_offset],n,MPI_DOUBLE,proc_direct[iswap],
                sendtag[iswap],world,&send_requests[nsendpost++]);
      send_offset += n;
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

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
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
  int iswap,ilist,nsend,nrecv;

  int nlocal = atom->nlocal;

  // setup lists of atoms to send in each direct swap
  // only maxlist possible lists (27 in 3d, 9 in 2d) regardless of stencil size
  // skip non-active lists as flagged in setup()

  AtomVec *avec = atom->avec;
  int dim = domain->dimension;

  build_lists();

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

  borders_comm();

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

void CommBrickDirect::forward_comm(Pair *pair, int size)
{
  int n,iswap,irecv;
  double *buf;

  int nsize;
  if (size) nsize = size;
  else nsize = pair->comm_forward;

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

  // copy atoms to self via pack and unpack
  // done before the sends below, b/c those use buf_send_direct from offset 0
  //   and their data must stay intact until the sends complete

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    pair->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                            pbc_flag_direct[iswap],pbc_direct[iswap]);
    pair->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
  }

  // send all owned atoms to receiving procs
  // except for self copies
  // each swap packs into its own region of buf_send_direct and is sent with
  //   a non-blocking send, so packing does not wait on the previous message

  int nsendpost = 0;
  int send_offset = 0;

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap] == 0) continue;
    n = pair->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],
                                &buf_send_direct[send_offset],
                                pbc_flag_direct[iswap],pbc_direct[iswap]);
    if (n) {
      MPI_Isend(&buf_send_direct[send_offset],n,MPI_DOUBLE,proc_direct[iswap],
                sendtag[iswap],world,&send_requests[nsendpost++]);
      send_offset += n;
    }
  }

  // wait on incoming messages with ghost atoms
  // unpack each message as it arrives

  if (npost == 0) {
    if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
    return;
  }

  for (int ipost = 0; ipost < npost; ipost++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = recv_indices_direct[irecv];
    offset = nsize * recv_offset_forward_atoms[iswap];
    pair->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],&buf_recv_direct[offset]);
  }

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Pair
   nsize used to set recv buffer offsets and limits
------------------------------------------------------------------------- */

void CommBrickDirect::reverse_comm(Pair *pair, int size)
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
      offset = nsize * recv_offset_reverse_atoms[iswap];
      MPI_Irecv(&buf_recv_direct[offset],nsize*sendnum_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],sendtag[iswap],world,&requests[npost++]);
    }
  }

  // copy/sum atoms to self via pack and unpack
  // done before the sends below, b/c those use buf_send_direct from offset 0
  //   and their data must stay intact until the sends complete

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    pair->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    pair->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct);
  }

  // send all ghost atoms to receiving procs
  // except for self copy/sums
  // each swap packs into its own region of buf_send_direct and is sent with
  //   a non-blocking send, so packing does not wait on the previous message

  int nsendpost = 0;
  int send_offset = 0;

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (recvnum_direct[iswap] == 0) continue;
    n = pair->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],
                                &buf_send_direct[send_offset]);
    if (n) {
      MPI_Isend(&buf_send_direct[send_offset],n,MPI_DOUBLE,proc_direct[iswap],
                recvtag[iswap],world,&send_requests[nsendpost++]);
      send_offset += n;
    }
  }

  // wait on incoming messages with owned atoms
  // unpack each message as it arrives

  if (npost == 0) {
    if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
    return;
  }

  for (int i = 0; i < npost; i++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = send_indices_direct[irecv];
    offset = nsize * recv_offset_reverse_atoms[iswap];
    pair->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],&buf_recv_direct[offset]);
  }

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Bond
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommBrickDirect::forward_comm(Bond * /*bond*/, int /*size*/)
{
  error->all(FLERR,"Comm_style brick/direct forward_comm for "
             "bond styles has not yet been implemented");
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Bond
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommBrickDirect::reverse_comm(Bond * /*bond*/, int /*size*/)
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

  // copy atoms to self via pack and unpack
  // done before the sends below, b/c those use buf_send_direct from offset 0
  //   and their data must stay intact until the sends complete

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    fix->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                           pbc_flag_direct[iswap],pbc_direct[iswap]);
    fix->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
  }

  // send all owned atoms to receiving procs
  // except for self copies
  // each swap packs into its own region of buf_send_direct and is sent with
  //   a non-blocking send, so packing does not wait on the previous message

  int nsendpost = 0;
  int send_offset = 0;

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap] == 0) continue;
    n = fix->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],
                               &buf_send_direct[send_offset],
                               pbc_flag_direct[iswap],pbc_direct[iswap]);
    if (n) {
      MPI_Isend(&buf_send_direct[send_offset],n,MPI_DOUBLE,proc_direct[iswap],
                sendtag[iswap],world,&send_requests[nsendpost++]);
      send_offset += n;
    }
  }

  // wait on incoming messages with ghost atoms
  // unpack each message as it arrives

  if (npost == 0) {
    if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
    return;
  }

  for (int ipost = 0; ipost < npost; ipost++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = recv_indices_direct[irecv];
    offset = nsize * recv_offset_forward_atoms[iswap];
    fix->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],&buf_recv_direct[offset]);
  }

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
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
      offset = nsize * recv_offset_reverse_atoms[iswap];
      MPI_Irecv(&buf_recv_direct[offset],nsize*sendnum_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],sendtag[iswap],world,&requests[npost++]);
    }
  }

  // copy/sum atoms to self via pack and unpack
  // done before the sends below, b/c those use buf_send_direct from offset 0
  //   and their data must stay intact until the sends complete

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    fix->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    fix->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct);
  }

  // send all ghost atoms to receiving procs
  // except for self copy/sums
  // each swap packs into its own region of buf_send_direct and is sent with
  //   a non-blocking send, so packing does not wait on the previous message

  int nsendpost = 0;
  int send_offset = 0;

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (recvnum_direct[iswap] == 0) continue;
    n = fix->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],
                               &buf_send_direct[send_offset]);
    if (n) {
      MPI_Isend(&buf_send_direct[send_offset],n,MPI_DOUBLE,proc_direct[iswap],
                recvtag[iswap],world,&send_requests[nsendpost++]);
      send_offset += n;
    }
  }

  // wait on incoming messages with owned atoms
  // unpack each message as it arrives

  if (npost == 0) {
    if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
    return;
  }

  for (int i = 0; i < npost; i++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = send_indices_direct[irecv];
    offset = nsize * recv_offset_reverse_atoms[iswap];
    fix->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],&buf_recv_direct[offset]);
  }

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix with variable size data
   query fix for pack size to ensure buf_send is big enough
   handshake sizes before each Irecv/Send to ensure buf_recv is big enough
------------------------------------------------------------------------- */

void CommBrickDirect::reverse_comm_variable(Fix * /*fix*/)
{
  error->all(FLERR,"Comm_style brick/direct reverse_comm for "
             "variables has not yet been implemented");
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Compute
   nsize used only to set recv buffer offsets and limits
------------------------------------------------------------------------- */

void CommBrickDirect::forward_comm(Compute *compute, int size)
{
  int n,iswap,irecv;
  double *buf;

  int nsize;
  if (size) nsize = size;
  else nsize = compute->comm_forward;

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

  // copy atoms to self via pack and unpack
  // done before the sends below, b/c those use buf_send_direct from offset 0
  //   and their data must stay intact until the sends complete

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    compute->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                               pbc_flag_direct[iswap],pbc_direct[iswap]);
    compute->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
  }

  // send all owned atoms to receiving procs
  // except for self copies
  // each swap packs into its own region of buf_send_direct and is sent with
  //   a non-blocking send, so packing does not wait on the previous message

  int nsendpost = 0;
  int send_offset = 0;

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap] == 0) continue;
    n = compute->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],
                                   &buf_send_direct[send_offset],
                                   pbc_flag_direct[iswap],pbc_direct[iswap]);
    if (n) {
      MPI_Isend(&buf_send_direct[send_offset],n,MPI_DOUBLE,proc_direct[iswap],
                sendtag[iswap],world,&send_requests[nsendpost++]);
      send_offset += n;
    }
  }

  // wait on incoming messages with ghost atoms
  // unpack each message as it arrives

  if (npost == 0) {
    if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
    return;
  }

  for (int ipost = 0; ipost < npost; ipost++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = recv_indices_direct[irecv];
    offset = nsize * recv_offset_forward_atoms[iswap];
    compute->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],&buf_recv_direct[offset]);
  }

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Compute
   nsize used only to set recv buffer offsets and limits
------------------------------------------------------------------------- */

void CommBrickDirect::reverse_comm(Compute *compute, int size)
{
  int n,iswap,irecv;
  double *buf;

  int nsize;
  if (size) nsize = size;
  else nsize = compute->comm_reverse;

  // post all receives for owned atoms
  // except for self copy/sums

  int offset;

  int npost = 0;
  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap]) {
      offset = nsize * recv_offset_reverse_atoms[iswap];
      MPI_Irecv(&buf_recv_direct[offset],nsize*sendnum_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],sendtag[iswap],world,&requests[npost++]);
    }
  }

  // copy/sum atoms to self via pack and unpack
  // done before the sends below, b/c those use buf_send_direct from offset 0
  //   and their data must stay intact until the sends complete

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    compute->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    compute->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct);
  }

  // send all ghost atoms to receiving procs
  // except for self copy/sums
  // each swap packs into its own region of buf_send_direct and is sent with
  //   a non-blocking send, so packing does not wait on the previous message

  int nsendpost = 0;
  int send_offset = 0;

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (recvnum_direct[iswap] == 0) continue;
    n = compute->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],
                                   &buf_send_direct[send_offset]);
    if (n) {
      MPI_Isend(&buf_send_direct[send_offset],n,MPI_DOUBLE,proc_direct[iswap],
                recvtag[iswap],world,&send_requests[nsendpost++]);
      send_offset += n;
    }
  }

  // wait on incoming messages with owned atoms
  // unpack each message as it arrives

  if (npost == 0) {
    if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
    return;
  }

  for (int i = 0; i < npost; i++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = send_indices_direct[irecv];
    offset = nsize * recv_offset_reverse_atoms[iswap];
    compute->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],&buf_recv_direct[offset]);
  }

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Dump
   nsize used only to set recv buffer offsets and limits
------------------------------------------------------------------------- */

void CommBrickDirect::forward_comm(Dump *dump, int size)
{
  int n,iswap,irecv;
  double *buf;

  int nsize;
  if (size) nsize = size;
  else nsize = dump->comm_forward;

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

  // copy atoms to self via pack and unpack
  // done before the sends below, b/c those use buf_send_direct from offset 0
  //   and their data must stay intact until the sends complete

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    dump->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct,
                            pbc_flag_direct[iswap],pbc_direct[iswap]);
    dump->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
  }

  // send all owned atoms to receiving procs
  // except for self copies
  // each swap packs into its own region of buf_send_direct and is sent with
  //   a non-blocking send, so packing does not wait on the previous message

  int nsendpost = 0;
  int send_offset = 0;

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap] == 0) continue;
    n = dump->pack_forward_comm(sendnum_direct[iswap],sendlist_direct[iswap],
                                &buf_send_direct[send_offset],
                                pbc_flag_direct[iswap],pbc_direct[iswap]);
    if (n) {
      MPI_Isend(&buf_send_direct[send_offset],n,MPI_DOUBLE,proc_direct[iswap],
                sendtag[iswap],world,&send_requests[nsendpost++]);
      send_offset += n;
    }
  }

  // wait on incoming messages with ghost atoms
  // unpack each message as it arrives

  if (npost == 0) {
    if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
    return;
  }

  for (int ipost = 0; ipost < npost; ipost++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = recv_indices_direct[irecv];
    offset = nsize * recv_offset_forward_atoms[iswap];
    dump->unpack_forward_comm(recvnum_direct[iswap],firstrecv_direct[iswap],&buf_recv_direct[offset]);
  }

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Dump
   nsize used only to set recv buffer offsets and limits
------------------------------------------------------------------------- */

void CommBrickDirect::reverse_comm(Dump *dump, int size)
{
  int n,iswap,irecv;
  double *buf;

  int nsize;
  if (size) nsize = size;
  else nsize = dump->comm_reverse;

  // post all receives for owned atoms
  // except for self copy/sums

  int offset;

  int npost = 0;
  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap]) {
      offset = nsize * recv_offset_reverse_atoms[iswap];
      MPI_Irecv(&buf_recv_direct[offset],nsize*sendnum_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],sendtag[iswap],world,&requests[npost++]);
    }
  }

  // copy/sum atoms to self via pack and unpack
  // done before the sends below, b/c those use buf_send_direct from offset 0
  //   and their data must stay intact until the sends complete

  for (int iself = 0; iself < nself_direct; iself++) {
    iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    dump->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],buf_send_direct);
    dump->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],buf_send_direct);
  }

  // send all ghost atoms to receiving procs
  // except for self copy/sums
  // each swap packs into its own region of buf_send_direct and is sent with
  //   a non-blocking send, so packing does not wait on the previous message

  int nsendpost = 0;
  int send_offset = 0;

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (recvnum_direct[iswap] == 0) continue;
    n = dump->pack_reverse_comm(recvnum_direct[iswap],firstrecv_direct[iswap],
                                &buf_send_direct[send_offset]);
    if (n) {
      MPI_Isend(&buf_send_direct[send_offset],n,MPI_DOUBLE,proc_direct[iswap],
                recvtag[iswap],world,&send_requests[nsendpost++]);
      send_offset += n;
    }
  }

  // wait on incoming messages with owned atoms
  // unpack each message as it arrives

  if (npost == 0) {
    if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
    return;
  }

  for (int i = 0; i < npost; i++) {
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    iswap = send_indices_direct[irecv];
    offset = nsize * recv_offset_reverse_atoms[iswap];
    dump->unpack_reverse_comm(sendnum_direct[iswap],sendlist_direct[iswap],&buf_recv_direct[offset]);
  }

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
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

  // copy atoms to self
  // done before the sends below, b/c those use buf_send_direct from offset 0
  //   and their data must stay intact until the sends complete

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

  // send all owned atoms to receiving procs
  // except for self copies
  // each swap packs into its own region of buf_send_direct and is sent with
  //   a non-blocking send, so packing does not wait on the previous message

  int nsendpost = 0;
  int send_offset = 0;

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap] == 0) continue;
    m = send_offset;
    for (i = 0; i < sendnum_direct[iswap]; i++) {
      j = sendlist_direct[iswap][i];
      for (k = 0; k < nsize; k++)
        buf_send_direct[m++] = array[j][k];
    }
    if (m > send_offset) {
      MPI_Isend(&buf_send_direct[send_offset],m-send_offset,MPI_DOUBLE,
                proc_direct[iswap],sendtag[iswap],world,&send_requests[nsendpost++]);
      send_offset = m;
    }
  }

  // wait on incoming messages with ghost atoms
  // unpack each message as it arrives

  if (npost == 0) {
    if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
    return;
  }

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

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);
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
  send_requests = (MPI_Request *)
    memory->smalloc(maxdirect*sizeof(MPI_Request),"comm:send_requests");
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
  memory->sfree(send_requests);
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
  if (sendatoms_list) {
    for (int ilist = 0; ilist < nlist; ilist++)
      memory->destroy(sendatoms_list[ilist]);
    memory->sfree(sendatoms_list);
    sendatoms_list = nullptr;
  }
  memory->destroy(maxsendatoms_list);
}

/* ----------------------------------------------------------------------
   ensure send/recv buffers are large enough for all border & forward & reverse comm
   size of send buf is for a single swap
   size of recv buf is for all swaps
------------------------------------------------------------------------- */

void CommBrickDirect::check_buffer_sizes()
{
  int max = size_border * ssum_direct;      // borders() packs all swaps at once
  max = MAX(max,maxforward*ssum_direct);   // forward_comm() packs all swaps at once
  max = MAX(max,maxreverse*rsum_direct);   // reverse_comm() packs all swaps at once
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
