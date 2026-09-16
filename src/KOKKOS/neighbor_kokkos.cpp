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

#include "neighbor_kokkos.h"

#include "angle.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "bond.h"
#include "domain.h"
#include "dihedral.h"
#include "error.h"
#include "force.h"
#include "improper.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "npair.h"
#include "style_nbin.h"
#include "style_npair.h"
#include "style_nstencil.h"
#include "style_ntopo.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NeighborKokkos::NeighborKokkos(LAMMPS *lmp) : Neighbor(lmp),
  neighbond_host(lmp),neighbond_device(lmp)
{
  device_flag = 0;
  bondlist = nullptr;
  anglelist = nullptr;
  dihedrallist = nullptr;
  improperlist = nullptr;
}

/* ---------------------------------------------------------------------- */

NeighborKokkos::~NeighborKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_cutneighsq,cutneighsq);
    cutneighsq = nullptr;

    memoryKK->destroy_kokkos(k_cutneighghostsq,cutneighghostsq);
    cutneighghostsq = nullptr;

    memoryKK->destroy_kokkos(k_ex_type,ex_type);
    memoryKK->destroy_kokkos(k_ex1_type,ex1_type);
    memoryKK->destroy_kokkos(k_ex2_type,ex2_type);
    memoryKK->destroy_kokkos(k_ex_mol_group,ex_mol_group);
    memoryKK->destroy_kokkos(k_ex1_bit,ex1_bit);
    memoryKK->destroy_kokkos(k_ex2_bit,ex2_bit);
    memoryKK->destroy_kokkos(k_ex_mol_bit,ex_mol_bit);
    memoryKK->destroy_kokkos(k_ex_mol_intra,ex_mol_intra);
  }
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init()
{
  atomKK = (AtomKokkos *) atom;
  Neighbor::init();

  // the pairwise neighbor list build of the KOKKOS package looks up special
  // bonds in the per-atom special list only.  with a molecule template that
  // list does not exist, so all special bonds would be silently ignored.
  // atom styles using a molecule template have no KOKKOS version (yet) and
  // are already rejected by AtomKokkos::new_avec(), but check here as well,
  // so that adding one cannot make the neighbor lists silently incorrect

  if (atom->molecular == Atom::TEMPLATE)
    error->all(FLERR,Error::NOLASTLINE,
               "KOKKOS package does not support atom styles with a molecule template");

  // Neighbor::init() allocates the host-side xhold array, but KOKKOS stores
  // the positions of the last build in its own view of the same name and
  // never fills the host array.  free it, so that Neighbor::get_xhold()
  // returns a null pointer instead of an array that was never written to

  memory->destroy(Neighbor::xhold);
  Neighbor::xhold = nullptr;

  // 1st time allocation of xhold

  if (dist_check)
      xhold = DAT::ttransform_kkfloat_1d_3_lr("neigh:xhold",maxhold);
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_cutneighsq_kokkos(int n)
{
  memoryKK->create_kokkos(k_cutneighsq,cutneighsq,n+1,n+1,"neigh:cutneighsq");
  k_cutneighsq.modify_host();
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_cutneighghostsq_kokkos(int n)
{
  memoryKK->create_kokkos(k_cutneighghostsq,cutneighghostsq,n+1,n+1,"neigh:cutneighghostsq");
  k_cutneighghostsq.modify_host();
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::create_kokkos_list(int i)
{
  if (style != Neighbor::BIN)
    error->all(FLERR,"KOKKOS package only supports 'bin' neighbor lists");

  if (requests[i]->kokkos_device) {
    lists[i] = new NeighListKokkos<LMPDeviceType>(lmp);
    device_flag = 1;
  } else if (requests[i]->kokkos_host)
    lists[i] = new NeighListKokkos<LMPHostType>(lmp);
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_ex_type_kokkos(int n)
{
  memoryKK->create_kokkos(k_ex_type,ex_type,n+1,n+1,"neigh:ex_type");
  k_ex_type.modify_host();
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_ex_bit_kokkos()
{
  memoryKK->create_kokkos(k_ex1_bit, ex1_bit, nex_group, "neigh:ex1_bit");
  k_ex1_bit.modify_host();
  memoryKK->create_kokkos(k_ex2_bit, ex2_bit, nex_group, "neigh:ex2_bit");
  k_ex2_bit.modify_host();
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_ex_mol_bit_kokkos()
{
  memoryKK->create_kokkos(k_ex_mol_bit, ex_mol_bit, nex_mol, "neigh:ex_mol_bit");
  k_ex_mol_bit.modify_host();
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::grow_ex_mol_intra_kokkos()
{
  memoryKK->grow_kokkos(k_ex_mol_intra, ex_mol_intra, maxex_mol, "neigh:ex_mol_intra");
  k_ex_mol_intra.modify_host();
}

/* ----------------------------------------------------------------------
   if any atom moved trigger distance (half of neighbor skin) return 1
   shrink trigger distance if box size has changed
   conservative shrink procedure:
     compute distance each of 8 corners of box has moved since last reneighbor
     reduce skin distance by sum of 2 largest of the 8 values
     if reduced skin distance is negative, set to zero
     new trigger = 1/2 of reduced skin distance
   for orthogonal box, only need 2 lo/hi corners
   for triclinic, need all 8 corners since deformations can displace all 8
------------------------------------------------------------------------- */

int NeighborKokkos::check_distance()
{
  if (device_flag)
    return check_distance_kokkos<LMPDeviceType>();
  else
    return check_distance_kokkos<LMPHostType>();
}

template<class DeviceType>
int NeighborKokkos::check_distance_kokkos()
{
  double delx,dely,delz;
  double delta,delta1,delta2;

  if (boxcheck) {
    if (triclinic == 0) {
      delx = bboxlo[0] - boxlo_hold[0];
      dely = bboxlo[1] - boxlo_hold[1];
      delz = bboxlo[2] - boxlo_hold[2];
      delta1 = sqrt(delx*delx + dely*dely + delz*delz);
      delx = bboxhi[0] - boxhi_hold[0];
      dely = bboxhi[1] - boxhi_hold[1];
      delz = bboxhi[2] - boxhi_hold[2];
      delta2 = sqrt(delx*delx + dely*dely + delz*delz);
      delta = 0.5 * (skin - (delta1+delta2));
      if (delta < 0.0) delta = 0.0;
      deltasq = delta*delta;
    } else {
      domain->box_corners();
      delta1 = delta2 = 0.0;
      for (int i = 0; i < 8; i++) {
        delx = corners[i][0] - corners_hold[i][0];
        dely = corners[i][1] - corners_hold[i][1];
        delz = corners[i][2] - corners_hold[i][2];
        delta = sqrt(delx*delx + dely*dely + delz*delz);
        if (delta > delta1) delta1 = delta;
        else if (delta > delta2) delta2 = delta;
      }
      delta = 0.5 * (skin - (delta1+delta2));
      if (delta < 0.0) delta = 0.0;
      deltasq = delta*delta;
    }
  } else deltasq = triggersq;

  atomKK->sync(ExecutionSpaceFromDevice<DeviceType>::space,X_MASK);
  x = atomKK->k_x;
  xhold.sync<DeviceType>();
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int flag = 0;
  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNeighborCheckDistance<DeviceType> >(0,nlocal),*this,flag);
  copymode = 0;

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall && ago == MAX(every,delay)) ndanger++;
  return flagall;
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void NeighborKokkos::operator()(TagNeighborCheckDistance<DeviceType>, const int &i, int &flag) const {
  const double delx = static_cast<double>(x.view<DeviceType>()(i,0) - xhold.view<DeviceType>()(i,0));
  const double dely = static_cast<double>(x.view<DeviceType>()(i,1) - xhold.view<DeviceType>()(i,1));
  const double delz = static_cast<double>(x.view<DeviceType>()(i,2) - xhold.view<DeviceType>()(i,2));
  const double rsq = delx*delx + dely*dely + delz*delz;
  if (rsq > deltasq) flag = 1;
}

/* ----------------------------------------------------------------------
   build perpetuals neighbor lists
   called at setup and every few timesteps during run or minimization
   topology lists also built if topoflag = 1, CUDA calls with topoflag = 0
------------------------------------------------------------------------- */


void NeighborKokkos::build(int topoflag)
{
  if (device_flag)
    build_kokkos<LMPDeviceType>(topoflag);
  else
    build_kokkos<LMPHostType>(topoflag);
}

template<class DeviceType>
void NeighborKokkos::build_kokkos(int topoflag)
{
  int i,m;

  ago = 0;
  ncalls++;
  lastcall = update->ntimestep;

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  // check that using special bond flags will not overflow neigh lists

  if (nall > NEIGHMASK)
    error->one(FLERR,Error::NOLASTLINE,"Too many local+ghost atoms for neighbor list");

  // store current atom positions and box size if needed

  if (dist_check) {
    atomKK->sync(ExecutionSpaceFromDevice<DeviceType>::space,X_MASK);
    x = atomKK->k_x;
    if (includegroup) nlocal = atom->nfirst;
    int maxhold_kokkos = xhold.view<DeviceType>().extent(0);
    if (atom->nmax > maxhold || maxhold_kokkos < maxhold) {
      maxhold = atom->nmax;
      xhold = DAT::ttransform_kkfloat_1d_3_lr("neigh:xhold",maxhold);
    }
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagNeighborXhold<DeviceType> >(0,nlocal),*this);
    copymode = 0;
    xhold.modify<DeviceType>();
    if (boxcheck) {
      if (triclinic == 0) {
        boxlo_hold[0] = bboxlo[0];
        boxlo_hold[1] = bboxlo[1];
        boxlo_hold[2] = bboxlo[2];
        boxhi_hold[0] = bboxhi[0];
        boxhi_hold[1] = bboxhi[1];
        boxhi_hold[2] = bboxhi[2];
      } else {
        domain->box_corners();
        corners = domain->corners;
        for (i = 0; i < 8; i++) {
          corners_hold[i][0] = corners[i][0];
          corners_hold[i][1] = corners[i][1];
          corners_hold[i][2] = corners[i][2];
        }
      }
    }
  }

  // bin atoms for all NBin instances
  // not just NBin associated with perpetual lists, also occasional lists
  // b/c cannot wait to bin occasional lists in build_one() call
  // if bin then, atoms may have moved outside of proc domain & bin extent,
  //   leading to errors or even a crash

  if (style != Neighbor::NSQ) {
    if (last_setup_bins < 0) setup_bins();
    for (int i = 0; i < nbin; i++) {
      if (!neigh_bin[i]->kokkos) atomKK->sync(Host,ALL_MASK);
      neigh_bin[i]->bin_atoms_setup(nall);
      neigh_bin[i]->bin_atoms();
    }
  }

  // build pairwise lists for all perpetual NPair/NeighList
  // grow() with nlocal/nall args so that only realloc if have to

  for (i = 0; i < npair_perpetual; i++) {
    m = plist[i];
    if (!lists[m]->kokkos) atomKK->sync(Host,ALL_MASK);
    if (!lists[m]->copy || lists[m]->trim || lists[m]->kk2cpu)
      lists[m]->grow(nlocal,nall);
    neigh_pair[m]->build_setup();
    neigh_pair[m]->build(lists[m]);
  }

  // build topology lists for bonds/angles/etc

  if ((atom->molecular != Atom::ATOMIC) && topoflag) build_topology();

  // reset last_build in all occasional lists
  // this will force them rebuild on next request
  // all occasional lists are now out-of-date b/c
  //   comm->exchange() occurred before neighbor->build()

  for (i = 0; i < npair_occasional; i++) {
    m = olist[i];
    neigh_pair[m]->last_build = -1;
  }
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void NeighborKokkos::operator()(TagNeighborXhold<DeviceType>, const int &i) const {
  xhold.view<DeviceType>()(i,0) = x.view<DeviceType>()(i,0);
  xhold.view<DeviceType>()(i,1) = x.view<DeviceType>()(i,1);
  xhold.view<DeviceType>()(i,2) = x.view<DeviceType>()(i,2);
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::modify_ex_type_grow_kokkos() {
  memoryKK->grow_kokkos(k_ex1_type,ex1_type,maxex_type,"neigh:ex1_type");
  k_ex1_type.modify_host();
  memoryKK->grow_kokkos(k_ex2_type,ex2_type,maxex_type,"neigh:ex2_type");
  k_ex2_type.modify_host();
}

/* ---------------------------------------------------------------------- */
void NeighborKokkos::modify_mol_group_grow_kokkos() {
  memoryKK->grow_kokkos(k_ex_mol_group,ex_mol_group,maxex_mol,"neigh:ex_mol_group");
  k_ex_mol_group.modify_host();
}

/* ---------------------------------------------------------------------- */
void NeighborKokkos::modify_mol_intra_grow_kokkos() {
  memoryKK->grow_kokkos(k_ex_mol_intra,ex_mol_intra,maxex_mol,"neigh:ex_mol_intra");
  k_ex_mol_intra.modify_host();
}

/* ---------------------------------------------------------------------- */
void NeighborKokkos::set_binsize_kokkos() {
  if (!binsizeflag && lmp->kokkos->ngpus > 0) {
    binsize_user = cutneighmax;
    binsizeflag = 1;
  }
}

/* ---------------------------------------------------------------------- */

void NeighborKokkos::init_topology() {
  if (device_flag) {
    neighbond_device.init_topology_kk();
  } else {
    neighbond_host.init_topology_kk();
  }
}

/* ----------------------------------------------------------------------
   build all topology neighbor lists every few timesteps
   normally built with pair lists, but CUDA separates them
------------------------------------------------------------------------- */

void NeighborKokkos::build_topology() {
  if (device_flag) {
    neighbond_device.build_topology_kk();

    k_bondlist = neighbond_device.k_bondlist;
    k_anglelist = neighbond_device.k_anglelist;
    k_dihedrallist = neighbond_device.k_dihedrallist;
    k_improperlist = neighbond_device.k_improperlist;

   } else {
    neighbond_host.build_topology_kk();

    k_bondlist = neighbond_host.k_bondlist;
    k_anglelist = neighbond_host.k_anglelist;
    k_dihedrallist = neighbond_host.k_dihedrallist;
    k_improperlist = neighbond_host.k_improperlist;
  }

  // transfer topology neighbor lists to the host for non-Kokkos styles,
  // which read them through the plain pointers in Neighbor

  if (force->bond && force->bond->execution_space == Host)
    k_bondlist.sync_host();
  if (force->angle && force->angle->execution_space == Host)
    k_anglelist.sync_host();
  if (force->dihedral && force->dihedral->execution_space == Host)
    k_dihedrallist.sync_host();
  if (force->improper && force->improper->execution_space == Host)
    k_improperlist.sync_host();
}
