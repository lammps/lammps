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

#include "fix_shake_kokkos.h"

#include "fix_rattle.h"
#include "atom_kokkos.h"
#include "atom_vec.h"
#include "molecule.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "domain.h"
#include "force.h"
#include "bond.h"
#include "angle.h"
#include "comm.h"
#include "group.h"
#include "fix_respa.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "atom_masks.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define RVOUS 1   // 0 for irregular, 1 for all2all

#define BIG 1.0e20
#define MASSDELTA 0.1

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixShakeKokkos<DeviceType>::FixShakeKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixShake(lmp, narg, arg)
{
  kokkosable = 1;
  forward_comm_device = exchange_comm_device = sort_device = 1;
  maxexchange = 9;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  memory->destroy(xshake);

  shake_flag_tmp = shake_flag;
  shake_atom_tmp = shake_atom;
  shake_type_tmp = shake_type;

  shake_flag = nullptr;
  shake_atom = nullptr;
  shake_type = nullptr;

  int nmax = atom->nmax;

  grow_arrays(nmax);

  for (int i = 0; i < atom->nlocal; i++) {
    k_shake_flag.h_view[i] = shake_flag_tmp[i];
    k_shake_atom.h_view(i,0) = shake_atom_tmp[i][0];
    k_shake_atom.h_view(i,1) = shake_atom_tmp[i][1];
    k_shake_atom.h_view(i,2) = shake_atom_tmp[i][2];
    k_shake_atom.h_view(i,3) = shake_atom_tmp[i][3];
    k_shake_type.h_view(i,0) = shake_type_tmp[i][0];
    k_shake_type.h_view(i,1) = shake_type_tmp[i][1];
    k_shake_type.h_view(i,2) = shake_type_tmp[i][2];
  }

  k_shake_flag.modify_host();
  k_shake_atom.modify_host();
  k_shake_type.modify_host();

  k_bond_distance = DAT::tdual_float_1d("fix_shake:bond_distance",atom->nbondtypes+1);
  k_angle_distance = DAT::tdual_float_1d("fix_shake:angle_distance",atom->nangletypes+1);

  d_bond_distance = k_bond_distance.view<DeviceType>();
  d_angle_distance = k_angle_distance.view<DeviceType>();

  // use 1D view for scalars to reduce GPU memory operations

  d_scalars = typename AT::t_int_1d("neighbor:scalars",2);
  h_scalars = HAT::t_int_1d("neighbor:scalars_mirror",2);

  d_error_flag = Kokkos::subview(d_scalars,0);
  d_nlist = Kokkos::subview(d_scalars,1);

  h_error_flag = Kokkos::subview(h_scalars,0);
  h_nlist = Kokkos::subview(h_scalars,1);

  d_count = typename AT::t_int_scalar("fix_shake:count");
  h_count = Kokkos::create_mirror_view(d_count);

  memory->destroy(shake_flag_tmp);
  memory->destroy(shake_atom_tmp);
  memory->destroy(shake_type_tmp);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixShakeKokkos<DeviceType>::~FixShakeKokkos()
{
  if (copymode) return;

  k_shake_flag.sync_host();
  k_shake_atom.sync_host();

  for (int i = 0; i < nlocal; i++) {
    if (shake_flag[i] == 0) continue;
    else if (shake_flag[i] == 1) {
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][1],1);
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][2],1);
      angletype_findset(i,shake_atom[i][1],shake_atom[i][2],1);
    } else if (shake_flag[i] == 2) {
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][1],1);
    } else if (shake_flag[i] == 3) {
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][1],1);
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][2],1);
    } else if (shake_flag[i] == 4) {
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][1],1);
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][2],1);
      bondtype_findset(i,shake_atom[i][0],shake_atom[i][3],1);
    }
  }

  memoryKK->destroy_kokkos(k_shake_flag,shake_flag);
  memoryKK->destroy_kokkos(k_shake_atom,shake_atom);
  memoryKK->destroy_kokkos(k_shake_type,shake_type);
  memoryKK->destroy_kokkos(k_xshake,xshake);
  memoryKK->destroy_kokkos(k_list,list);
  memoryKK->destroy_kokkos(k_closest_list,closest_list);

  memoryKK->destroy_kokkos(k_vatom,vatom);
}

/* ----------------------------------------------------------------------
   set bond and angle distances
   this init must happen after force->bond and force->angle inits
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::init()
{
  FixShake::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot yet use respa with Kokkos");

  if (rattle)
    error->all(FLERR,"Cannot yet use KOKKOS package with fix rattle");

  // set equilibrium bond distances

  for (int i = 1; i <= atom->nbondtypes; i++)
    k_bond_distance.h_view[i] = bond_distance[i];

  // set equilibrium angle distances

  for (int i = 1; i <= atom->nangletypes; i++)
    k_angle_distance.h_view[i] = angle_distance[i];

  k_bond_distance.modify_host();
  k_angle_distance.modify_host();

  k_bond_distance.sync<DeviceType>();
  k_angle_distance.sync<DeviceType>();
}


/* ----------------------------------------------------------------------
   run setup for minimization.
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::min_setup(int /*vflag*/)
{
  error->all(FLERR, "Cannot yet use fix {} during minimization with Kokkos", style);
}

/* ----------------------------------------------------------------------
   build list of SHAKE clusters to constrain
   if one or more atoms in cluster are on this proc,
     this proc lists the cluster exactly once
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::pre_neighbor()
{
  // local copies of atom quantities
  // used by SHAKE until next re-neighboring

  atomKK->sync(execution_space,X_MASK);

  ebond = 0.0;
  d_x = atomKK->k_x.view<DeviceType>();
  nlocal = atom->nlocal;

  map_style = atom->map_style;
  if (map_style == Atom::MAP_ARRAY) {
    k_map_array = atomKK->k_map_array;
    k_map_array.template sync<DeviceType>();
  } else if (map_style == Atom::MAP_HASH) {
    k_map_hash = atomKK->k_map_hash;
    k_map_hash.template sync<DeviceType>();
  }

  k_shake_flag.sync<DeviceType>();
  k_shake_atom.sync<DeviceType>();

  // extend size of SHAKE list if necessary

  if (nlocal > maxlist) {
    maxlist = atom->nmax;
    memoryKK->destroy_kokkos(k_list,list);
    memoryKK->create_kokkos(k_list,list,maxlist,"shake:list");
    d_list = k_list.view<DeviceType>();

    memoryKK->destroy_kokkos(k_closest_list,closest_list);
    memoryKK->create_kokkos(k_closest_list,closest_list,maxlist,4,"shake:closest_list");
    d_closest_list = k_closest_list.view<DeviceType>();
  }

  // Atom Map

  map_style = atom->map_style;

  if (map_style == Atom::MAP_ARRAY) {
    k_map_array = atomKK->k_map_array;
    k_map_array.template sync<DeviceType>();
  } else if (map_style == Atom::MAP_HASH) {
    k_map_hash = atomKK->k_map_hash;
    k_map_hash.template sync<DeviceType>();
  }

  atomKK->k_sametag.sync<DeviceType>();
  d_sametag = atomKK->k_sametag.view<DeviceType>();

  // build list of SHAKE clusters I compute

  Kokkos::deep_copy(d_scalars,0);

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixShakePreNeighbor>(0,nlocal),*this);

  copymode = 0;

  k_list.modify<DeviceType>();
  k_closest_list.modify<DeviceType>();

  Kokkos::deep_copy(h_scalars,d_scalars);
  nlist = h_nlist();

  if (h_error_flag() == 1) {
    error->one(FLERR,"Shake atoms missing on proc "
                                 "{} at step {}",comm->me,update->ntimestep);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixShakeKokkos<DeviceType>::operator()(TagFixShakePreNeighbor, const int &i) const {
  if (d_shake_flag[i]) {
    if (d_shake_flag[i] == 2) {
      int atom1 = AtomKokkos::map_kokkos<DeviceType>(d_shake_atom(i,0),map_style,k_map_array,k_map_hash);
      int atom2 = AtomKokkos::map_kokkos<DeviceType>(d_shake_atom(i,1),map_style,k_map_array,k_map_hash);
      if (atom1 == -1 || atom2 == -1) {
        d_error_flag() = 1;
      }
      atom1 = closest_image(i, atom1);
      atom2 = closest_image(i, atom2);
      if (i <= atom1 && i <= atom2) {
        const int nlist = Kokkos::atomic_fetch_add(&d_nlist(),1);
        d_list[nlist] = i;
        d_closest_list(nlist,0) = atom1;
        d_closest_list(nlist,1) = atom2;
      }
    } else if (d_shake_flag[i] % 2 == 1) {
      int atom1 = AtomKokkos::map_kokkos<DeviceType>(d_shake_atom(i,0),map_style,k_map_array,k_map_hash);
      int atom2 = AtomKokkos::map_kokkos<DeviceType>(d_shake_atom(i,1),map_style,k_map_array,k_map_hash);
      int atom3 = AtomKokkos::map_kokkos<DeviceType>(d_shake_atom(i,2),map_style,k_map_array,k_map_hash);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1)
        d_error_flag() = 1;
      atom1 = closest_image(i, atom1);
      atom2 = closest_image(i, atom2);
      atom3 = closest_image(i, atom3);
      if (i <= atom1 && i <= atom2 && i <= atom3) {
        const int nlist = Kokkos::atomic_fetch_add(&d_nlist(),1);
        d_list[nlist] = i;
        d_closest_list(nlist,0) = atom1;
        d_closest_list(nlist,1) = atom2;
        d_closest_list(nlist,2) = atom3;
      }
    } else {
      int atom1 = AtomKokkos::map_kokkos<DeviceType>(d_shake_atom(i,0),map_style,k_map_array,k_map_hash);
      int atom2 = AtomKokkos::map_kokkos<DeviceType>(d_shake_atom(i,1),map_style,k_map_array,k_map_hash);
      int atom3 = AtomKokkos::map_kokkos<DeviceType>(d_shake_atom(i,2),map_style,k_map_array,k_map_hash);
      int atom4 = AtomKokkos::map_kokkos<DeviceType>(d_shake_atom(i,3),map_style,k_map_array,k_map_hash);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
        d_error_flag() = 1;
      atom1 = closest_image(i, atom1);
      atom2 = closest_image(i, atom2);
      atom3 = closest_image(i, atom3);
      atom4 = closest_image(i, atom4);
      if (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4) {
        const int nlist = Kokkos::atomic_fetch_add(&d_nlist(),1);
        d_list[nlist] = i;
        d_closest_list(nlist,0) = atom1;
        d_closest_list(nlist,1) = atom2;
        d_closest_list(nlist,2) = atom3;
        d_closest_list(nlist,3) = atom4;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute the force adjustment for SHAKE constraint
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::post_force(int vflag)
{
  ebond = 0.0;
  copymode = 1;

  d_x = atomKK->k_x.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();
  d_mass = atomKK->k_mass.view<DeviceType>();
  nlocal = atomKK->nlocal;

  map_style = atom->map_style;
  if (map_style == Atom::MAP_ARRAY) {
    k_map_array = atomKK->k_map_array;
    k_map_array.template sync<DeviceType>();
  } else if (map_style == Atom::MAP_HASH) {
    k_map_hash = atomKK->k_map_hash;
    k_map_hash.template sync<DeviceType>();
  }

  if (d_rmass.data())
    atomKK->sync(execution_space,X_MASK|F_MASK|RMASS_MASK);
  else
    atomKK->sync(execution_space,X_MASK|F_MASK|TYPE_MASK);

  k_shake_flag.sync<DeviceType>();
  k_shake_type.sync<DeviceType>();
  k_list.sync<DeviceType>();
  k_closest_list.sync<DeviceType>();

  if (update->ntimestep == next_output)
    stats();

  // xshake = unconstrained move with current v,f
  // communicate results if necessary

  unconstrained_update();
  comm->forward_comm(this);
  k_xshake.sync<DeviceType>();

  // virial setup

  v_init(vflag);

  // reallocate per-atom arrays if necessary

  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"improper:vatom");
    d_vatom = k_vatom.template view<KKDeviceType>();
  }

  neighflag = lmp->kokkos->neighflag;

  // FULL neighlist still needs atomics in fix shake

  if (neighflag == FULL) {
    if (lmp->kokkos->nthreads > 1 || lmp->kokkos->ngpus > 0)
      neighflag = HALFTHREAD;
    else
      neighflag = HALF;
  }

  need_dup = 0;
  if (neighflag != HALF)
    need_dup = std::is_same<typename NeedDup<HALFTHREAD,DeviceType>::value,Kokkos::Experimental::ScatterDuplicated>::value;

  // allocate duplicated memory

  if (need_dup) {
    dup_f            = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_f);
    dup_vatom        = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f            = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_f);
    ndup_vatom        = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  Kokkos::deep_copy(d_error_flag,0);

  EV_FLOAT ev;

  // loop over clusters to add constraint forces

  if (neighflag == HALF) {
   if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixShakePostForce<HALF,1> >(0,nlist),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixShakePostForce<HALF,0> >(0,nlist),*this);
  } else {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixShakePostForce<HALFTHREAD,1> >(0,nlist),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixShakePostForce<HALFTHREAD,0> >(0,nlist),*this);
  }

  copymode = 0;

  Kokkos::deep_copy(h_error_flag,d_error_flag);

  if (h_error_flag() == 2)
    error->warning(FLERR,"Shake determinant < 0.0");
  else if (h_error_flag() == 3)
    error->one(FLERR,"Shake determinant = 0.0");

  // store vflag for coordinate_constraints_end_of_step()

  vflag_post_force = vflag;

  // reduction over duplicated memory

  if (need_dup)
    Kokkos::Experimental::contribute(d_f,dup_f);

  atomKK->modified(execution_space,F_MASK);

  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  // free duplicated memory

  if (need_dup) {
    dup_f = decltype(dup_f)();
    dup_vatom = decltype(dup_vatom)();
  }
}

/* ----------------------------------------------------------------------
   substitute shake constraints with very strong bonds
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::min_post_force(int vflag)
{
  // not yet ported to Kokkos

  atomKK->sync(Host,X_MASK | F_MASK);
  k_shake_flag.sync_host();
  k_shake_type.sync_host();
  k_list.sync_host();
  k_closest_list.sync_host();

  FixShake::min_post_force(vflag);

  atomKK->modified(Host,F_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void FixShakeKokkos<DeviceType>::operator()(TagFixShakePostForce<NEIGHFLAG,EVFLAG>, const int &i, EV_FLOAT& ev) const {
  const int m = d_list[i];
  if (d_shake_flag[m] == 2) shake<NEIGHFLAG,EVFLAG>(i,ev);
  else if (d_shake_flag[m] == 3) shake3<NEIGHFLAG,EVFLAG>(i,ev);
  else if (d_shake_flag[m] == 4) shake4<NEIGHFLAG,EVFLAG>(i,ev);
  else shake3angle<NEIGHFLAG,EVFLAG>(i,ev);
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void FixShakeKokkos<DeviceType>::operator()(TagFixShakePostForce<NEIGHFLAG,EVFLAG>, const int &i) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagFixShakePostForce<NEIGHFLAG,EVFLAG>(), i, ev);
}

/* ----------------------------------------------------------------------
   count # of degrees-of-freedom removed by SHAKE for atoms in igroup
------------------------------------------------------------------------- */

template<class DeviceType>
int FixShakeKokkos<DeviceType>::dof(int igroup)
{
  d_mask = atomKK->k_mask.view<DeviceType>();
  d_tag = atomKK->k_tag.view<DeviceType>();
  nlocal = atom->nlocal;

  atomKK->sync(execution_space,MASK_MASK|TAG_MASK);
  k_shake_flag.sync<DeviceType>();
  k_shake_atom.sync<DeviceType>();

  // count dof in a cluster if and only if
  // the central atom is in group and atom i is the central atom

  int n = 0;
  {
    // local variables for lambda capture

    auto d_shake_flag = this->d_shake_flag;
    auto d_shake_atom = this->d_shake_atom;
    auto tag = this->d_tag;
    auto mask = this->d_mask;
    auto groupbit = group->bitmask[igroup];

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,nlocal),
     LAMMPS_LAMBDA(const int& i, int& n) {
      if (!(mask[i] & groupbit)) return;
      if (d_shake_flag[i] == 0) return;
      if (d_shake_atom(i,0) != tag[i]) return;
      if (d_shake_flag[i] == 1) n += 3;
      else if (d_shake_flag[i] == 2) n += 1;
      else if (d_shake_flag[i] == 3) n += 2;
      else if (d_shake_flag[i] == 4) n += 3;
    },n);
  }

  int nall;
  MPI_Allreduce(&n,&nall,1,MPI_INT,MPI_SUM,world);
  return nall;
}


/* ----------------------------------------------------------------------
   assumes NVE update, seems to be accurate enough for NVT,NPT,NPH as well
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::unconstrained_update()
{
  d_x = atomKK->k_x.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();
  d_mass = atomKK->k_mass.view<DeviceType>();
  nlocal = atom->nlocal;

  if (d_rmass.data())
    atomKK->sync(execution_space,X_MASK|V_MASK|F_MASK|RMASS_MASK);
  else
    atomKK->sync(execution_space,X_MASK|V_MASK|F_MASK|TYPE_MASK);

  k_shake_flag.sync<DeviceType>();
  k_xshake.sync<DeviceType>();

  {
    // local variables for lambda capture

    auto d_shake_flag = this->d_shake_flag;
    auto d_xshake = this->d_xshake;
    auto x = this->d_x;
    auto v = this->d_v;
    auto f = this->d_f;
    auto dtfsq = this->dtfsq;
    auto dtv = this->dtv;

    if (d_rmass.data()) {

      auto rmass = this->d_rmass;

      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nlocal),
       LAMMPS_LAMBDA(const int& i) {
        if (d_shake_flag[i]) {
          const double dtfmsq = dtfsq / rmass[i];
          d_xshake(i,0) = x(i,0) + dtv*v(i,0) + dtfmsq*f(i,0);
          d_xshake(i,1) = x(i,1) + dtv*v(i,1) + dtfmsq*f(i,1);
          d_xshake(i,2) = x(i,2) + dtv*v(i,2) + dtfmsq*f(i,2);
        } else d_xshake(i,2) = d_xshake(i,1) = d_xshake(i,0) = 0.0;
      });
    } else {

      auto mass = this->d_mass;
      auto type = this->d_type;

      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nlocal),
       LAMMPS_LAMBDA(const int& i) {
        if (d_shake_flag[i]) {
          const double dtfmsq = dtfsq / mass[type[i]];
          d_xshake(i,0) = x(i,0) + dtv*v(i,0) + dtfmsq*f(i,0);
          d_xshake(i,1) = x(i,1) + dtv*v(i,1) + dtfmsq*f(i,1);
          d_xshake(i,2) = x(i,2) + dtv*v(i,2) + dtfmsq*f(i,2);
        } else d_xshake(i,2) = d_xshake(i,1) = d_xshake(i,0) = 0.0;
      });
    }
  }

  k_xshake.modify<DeviceType>();
}

/* ----------------------------------------------------------------------
   calculate SHAKE constraint forces for size 2 cluster = single bond
------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void FixShakeKokkos<DeviceType>::shake(int ilist, EV_FLOAT& ev) const
{

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  int atomlist[2];
  double v[6];
  double invmass0,invmass1;

  // local atom IDs and constraint distances

  int m = d_list[ilist];
  int i0 = d_closest_list(ilist,0);
  int i1 = d_closest_list(ilist,1);
  double bond1 = d_bond_distance[d_shake_type(m,0)];

  // r01 = distance vec between atoms

  double r01[3];
  r01[0] = d_x(i0,0) - d_x(i1,0);
  r01[1] = d_x(i0,1) - d_x(i1,1);
  r01[2] = d_x(i0,2) - d_x(i1,2);

  // s01 = distance vec after unconstrained update

  double s01[3];
  s01[0] = d_xshake(i0,0) - d_xshake(i1,0);
  s01[1] = d_xshake(i0,1) - d_xshake(i1,1);
  s01[2] = d_xshake(i0,2) - d_xshake(i1,2);

  // scalar distances between atoms

  double r01sq = r01[0]*r01[0] + r01[1]*r01[1] + r01[2]*r01[2];
  double s01sq = s01[0]*s01[0] + s01[1]*s01[1] + s01[2]*s01[2];

  // a,b,c = coeffs in quadratic equation for lamda

  if (d_rmass.data()) {
    invmass0 = 1.0/d_rmass[i0];
    invmass1 = 1.0/d_rmass[i1];
  } else {
    invmass0 = 1.0/d_mass[d_type[i0]];
    invmass1 = 1.0/d_mass[d_type[i1]];
  }

  double a = (invmass0+invmass1)*(invmass0+invmass1) * r01sq;
  double b = 2.0 * (invmass0+invmass1) *
    (s01[0]*r01[0] + s01[1]*r01[1] + s01[2]*r01[2]);
  double c = s01sq - bond1*bond1;

  // error check

  double determ = b*b - 4.0*a*c;
  if (determ < 0.0) {
    //error->warning(FLERR,"Shake determinant < 0.0",0);
    d_error_flag() = 2;
    determ = 0.0;
  }

  // exact quadratic solution for lamda

  double lamda,lamda1,lamda2;
  lamda1 = (-b+sqrt(determ)) / (2.0*a);
  lamda2 = (-b-sqrt(determ)) / (2.0*a);

  if (fabs(lamda1) <= fabs(lamda2)) lamda = lamda1;
  else lamda = lamda2;

  // update forces if atom is owned by this processor

  lamda /= dtfsq;

  if (i0 < nlocal) {
    a_f(i0,0) += lamda*r01[0];
    a_f(i0,1) += lamda*r01[1];
    a_f(i0,2) += lamda*r01[2];
  }

  if (i1 < nlocal) {
    a_f(i1,0) -= lamda*r01[0];
    a_f(i1,1) -= lamda*r01[1];
    a_f(i1,2) -= lamda*r01[2];
  }

  if (EVFLAG) {
    int count = 0;
    if (i0 < nlocal) atomlist[count++] = i0;
    if (i1 < nlocal) atomlist[count++] = i1;

    v[0] = lamda*r01[0]*r01[0];
    v[1] = lamda*r01[1]*r01[1];
    v[2] = lamda*r01[2]*r01[2];
    v[3] = lamda*r01[0]*r01[1];
    v[4] = lamda*r01[0]*r01[2];
    v[5] = lamda*r01[1]*r01[2];

    v_tally<NEIGHFLAG>(ev,count,atomlist,2.0,v);
  }
}

/* ----------------------------------------------------------------------
   calculate SHAKE constraint forces for size 3 cluster = two bonds
------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void FixShakeKokkos<DeviceType>::shake3(int ilist, EV_FLOAT& ev) const
{

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  int atomlist[3];
  double v[6];
  double invmass0,invmass1,invmass2;

  // local atom IDs and constraint distances

  int m = d_list[ilist];
  int i0 = d_closest_list(ilist,0);
  int i1 = d_closest_list(ilist,1);
  int i2 = d_closest_list(ilist,2);
  double bond1 = d_bond_distance[d_shake_type(m,0)];
  double bond2 = d_bond_distance[d_shake_type(m,1)];

  // r01,r02 = distance vec between atoms

  double r01[3];
  r01[0] = d_x(i0,0) - d_x(i1,0);
  r01[1] = d_x(i0,1) - d_x(i1,1);
  r01[2] = d_x(i0,2) - d_x(i1,2);

  double r02[3];
  r02[0] = d_x(i0,0) - d_x(i2,0);
  r02[1] = d_x(i0,1) - d_x(i2,1);
  r02[2] = d_x(i0,2) - d_x(i2,2);

  // s01,s02 = distance vec after unconstrained update

  double s01[3];
  s01[0] = d_xshake(i0,0) - d_xshake(i1,0);
  s01[1] = d_xshake(i0,1) - d_xshake(i1,1);
  s01[2] = d_xshake(i0,2) - d_xshake(i1,2);

  double s02[3];
  s02[0] = d_xshake(i0,0) - d_xshake(i2,0);
  s02[1] = d_xshake(i0,1) - d_xshake(i2,1);
  s02[2] = d_xshake(i0,2) - d_xshake(i2,2);

  // scalar distances between atoms

  double r01sq = r01[0]*r01[0] + r01[1]*r01[1] + r01[2]*r01[2];
  double r02sq = r02[0]*r02[0] + r02[1]*r02[1] + r02[2]*r02[2];
  double s01sq = s01[0]*s01[0] + s01[1]*s01[1] + s01[2]*s01[2];
  double s02sq = s02[0]*s02[0] + s02[1]*s02[1] + s02[2]*s02[2];

  // matrix coeffs and rhs for lamda equations

  if (d_rmass.data()) {
    invmass0 = 1.0/d_rmass[i0];
    invmass1 = 1.0/d_rmass[i1];
    invmass2 = 1.0/d_rmass[i2];
  } else {
    invmass0 = 1.0/d_mass[d_type[i0]];
    invmass1 = 1.0/d_mass[d_type[i1]];
    invmass2 = 1.0/d_mass[d_type[i2]];
  }

  double a11 = 2.0 * (invmass0+invmass1) *
    (s01[0]*r01[0] + s01[1]*r01[1] + s01[2]*r01[2]);
  double a12 = 2.0 * invmass0 *
    (s01[0]*r02[0] + s01[1]*r02[1] + s01[2]*r02[2]);
  double a21 = 2.0 * invmass0 *
    (s02[0]*r01[0] + s02[1]*r01[1] + s02[2]*r01[2]);
  double a22 = 2.0 * (invmass0+invmass2) *
    (s02[0]*r02[0] + s02[1]*r02[1] + s02[2]*r02[2]);

  // inverse of matrix

  double determ = a11*a22 - a12*a21;
  if (determ == 0.0) d_error_flag() = 3;
  //error->one(FLERR,"Shake determinant = 0.0");
  double determinv = 1.0/determ;

  double a11inv = a22*determinv;
  double a12inv = -a12*determinv;
  double a21inv = -a21*determinv;
  double a22inv = a11*determinv;

  // quadratic correction coeffs

  double r0102 = (r01[0]*r02[0] + r01[1]*r02[1] + r01[2]*r02[2]);

  double quad1_0101 = (invmass0+invmass1)*(invmass0+invmass1) * r01sq;
  double quad1_0202 = invmass0*invmass0 * r02sq;
  double quad1_0102 = 2.0 * (invmass0+invmass1)*invmass0 * r0102;

  double quad2_0202 = (invmass0+invmass2)*(invmass0+invmass2) * r02sq;
  double quad2_0101 = invmass0*invmass0 * r01sq;
  double quad2_0102 = 2.0 * (invmass0+invmass2)*invmass0 * r0102;

  // iterate until converged

  double lamda01 = 0.0;
  double lamda02 = 0.0;
  int niter = 0;
  int done = 0;

  double quad1,quad2,b1,b2,lamda01_new,lamda02_new;

  while (!done && niter < max_iter) {
    quad1 = quad1_0101 * lamda01*lamda01 + quad1_0202 * lamda02*lamda02 +
      quad1_0102 * lamda01*lamda02;
    quad2 = quad2_0101 * lamda01*lamda01 + quad2_0202 * lamda02*lamda02 +
      quad2_0102 * lamda01*lamda02;

    b1 = bond1*bond1 - s01sq - quad1;
    b2 = bond2*bond2 - s02sq - quad2;

    lamda01_new = a11inv*b1 + a12inv*b2;
    lamda02_new = a21inv*b1 + a22inv*b2;

    done = 1;
    if (fabs(lamda01_new-lamda01) > tolerance) done = 0;
    if (fabs(lamda02_new-lamda02) > tolerance) done = 0;

    lamda01 = lamda01_new;
    lamda02 = lamda02_new;

    // stop iterations before we have a floating point overflow
    // max double is < 1.0e308, so 1e150 is a reasonable cutoff

    if (fabs(lamda01) > 1e150 || fabs(lamda02) > 1e150) done = 1;

    niter++;
  }

  // update forces if atom is owned by this processor

  lamda01 = lamda01/dtfsq;
  lamda02 = lamda02/dtfsq;

  if (i0 < nlocal) {
    a_f(i0,0) += lamda01*r01[0] + lamda02*r02[0];
    a_f(i0,1) += lamda01*r01[1] + lamda02*r02[1];
    a_f(i0,2) += lamda01*r01[2] + lamda02*r02[2];
  }

  if (i1 < nlocal) {
    a_f(i1,0) -= lamda01*r01[0];
    a_f(i1,1) -= lamda01*r01[1];
    a_f(i1,2) -= lamda01*r01[2];
  }

  if (i2 < nlocal) {
    a_f(i2,0) -= lamda02*r02[0];
    a_f(i2,1) -= lamda02*r02[1];
    a_f(i2,2) -= lamda02*r02[2];
  }

  if (EVFLAG) {
    int count = 0;
    if (i0 < nlocal) atomlist[count++] = i0;
    if (i1 < nlocal) atomlist[count++] = i1;
    if (i2 < nlocal) atomlist[count++] = i2;

    v[0] = lamda01*r01[0]*r01[0] + lamda02*r02[0]*r02[0];
    v[1] = lamda01*r01[1]*r01[1] + lamda02*r02[1]*r02[1];
    v[2] = lamda01*r01[2]*r01[2] + lamda02*r02[2]*r02[2];
    v[3] = lamda01*r01[0]*r01[1] + lamda02*r02[0]*r02[1];
    v[4] = lamda01*r01[0]*r01[2] + lamda02*r02[0]*r02[2];
    v[5] = lamda01*r01[1]*r01[2] + lamda02*r02[1]*r02[2];

    v_tally<NEIGHFLAG>(ev,count,atomlist,3.0,v);
  }
}

/* ----------------------------------------------------------------------
   calculate SHAKE constraint forces for size 4 cluster = three bonds
------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void FixShakeKokkos<DeviceType>::shake4(int ilist, EV_FLOAT& ev) const
{

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  int atomlist[4];
  double v[6];
  double invmass0,invmass1,invmass2,invmass3;

  // local atom IDs and constraint distances

  int m = d_list[ilist];
  int i0 = d_closest_list(ilist,0);
  int i1 = d_closest_list(ilist,1);
  int i2 = d_closest_list(ilist,2);
  int i3 = d_closest_list(ilist,3);
  double bond1 = d_bond_distance[d_shake_type(m,0)];
  double bond2 = d_bond_distance[d_shake_type(m,1)];
  double bond3 = d_bond_distance[d_shake_type(m,2)];

  // r01,r02,r03 = distance vec between atoms

  double r01[3];
  r01[0] = d_x(i0,0) - d_x(i1,0);
  r01[1] = d_x(i0,1) - d_x(i1,1);
  r01[2] = d_x(i0,2) - d_x(i1,2);

  double r02[3];
  r02[0] = d_x(i0,0) - d_x(i2,0);
  r02[1] = d_x(i0,1) - d_x(i2,1);
  r02[2] = d_x(i0,2) - d_x(i2,2);

  double r03[3];
  r03[0] = d_x(i0,0) - d_x(i3,0);
  r03[1] = d_x(i0,1) - d_x(i3,1);
  r03[2] = d_x(i0,2) - d_x(i3,2);

  // s01,s02,s03 = distance vec after unconstrained update

  double s01[3];
  s01[0] = d_xshake(i0,0) - d_xshake(i1,0);
  s01[1] = d_xshake(i0,1) - d_xshake(i1,1);
  s01[2] = d_xshake(i0,2) - d_xshake(i1,2);

  double s02[3];
  s02[0] = d_xshake(i0,0) - d_xshake(i2,0);
  s02[1] = d_xshake(i0,1) - d_xshake(i2,1);
  s02[2] = d_xshake(i0,2) - d_xshake(i2,2);

  double s03[3];
  s03[0] = d_xshake(i0,0) - d_xshake(i3,0);
  s03[1] = d_xshake(i0,1) - d_xshake(i3,1);
  s03[2] = d_xshake(i0,2) - d_xshake(i3,2);

  // scalar distances between atoms

  double r01sq = r01[0]*r01[0] + r01[1]*r01[1] + r01[2]*r01[2];
  double r02sq = r02[0]*r02[0] + r02[1]*r02[1] + r02[2]*r02[2];
  double r03sq = r03[0]*r03[0] + r03[1]*r03[1] + r03[2]*r03[2];
  double s01sq = s01[0]*s01[0] + s01[1]*s01[1] + s01[2]*s01[2];
  double s02sq = s02[0]*s02[0] + s02[1]*s02[1] + s02[2]*s02[2];
  double s03sq = s03[0]*s03[0] + s03[1]*s03[1] + s03[2]*s03[2];

  // matrix coeffs and rhs for lamda equations

  if (d_rmass.data()) {
    invmass0 = 1.0/d_rmass[i0];
    invmass1 = 1.0/d_rmass[i1];
    invmass2 = 1.0/d_rmass[i2];
    invmass3 = 1.0/d_rmass[i3];
  } else {
    invmass0 = 1.0/d_mass[d_type[i0]];
    invmass1 = 1.0/d_mass[d_type[i1]];
    invmass2 = 1.0/d_mass[d_type[i2]];
    invmass3 = 1.0/d_mass[d_type[i3]];
  }

  double a11 = 2.0 * (invmass0+invmass1) *
    (s01[0]*r01[0] + s01[1]*r01[1] + s01[2]*r01[2]);
  double a12 = 2.0 * invmass0 *
    (s01[0]*r02[0] + s01[1]*r02[1] + s01[2]*r02[2]);
  double a13 = 2.0 * invmass0 *
    (s01[0]*r03[0] + s01[1]*r03[1] + s01[2]*r03[2]);
  double a21 = 2.0 * invmass0 *
    (s02[0]*r01[0] + s02[1]*r01[1] + s02[2]*r01[2]);
  double a22 = 2.0 * (invmass0+invmass2) *
    (s02[0]*r02[0] + s02[1]*r02[1] + s02[2]*r02[2]);
  double a23 = 2.0 * invmass0 *
    (s02[0]*r03[0] + s02[1]*r03[1] + s02[2]*r03[2]);
  double a31 = 2.0 * invmass0 *
    (s03[0]*r01[0] + s03[1]*r01[1] + s03[2]*r01[2]);
  double a32 = 2.0 * invmass0 *
    (s03[0]*r02[0] + s03[1]*r02[1] + s03[2]*r02[2]);
  double a33 = 2.0 * (invmass0+invmass3) *
    (s03[0]*r03[0] + s03[1]*r03[1] + s03[2]*r03[2]);

  // inverse of matrix;

  double determ = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 -
    a11*a23*a32 - a12*a21*a33 - a13*a22*a31;
  if (determ == 0.0) d_error_flag() = 3;
  //error->one(FLERR,"Shake determinant = 0.0");
  double determinv = 1.0/determ;

  double a11inv = determinv * (a22*a33 - a23*a32);
  double a12inv = -determinv * (a12*a33 - a13*a32);
  double a13inv = determinv * (a12*a23 - a13*a22);
  double a21inv = -determinv * (a21*a33 - a23*a31);
  double a22inv = determinv * (a11*a33 - a13*a31);
  double a23inv = -determinv * (a11*a23 - a13*a21);
  double a31inv = determinv * (a21*a32 - a22*a31);
  double a32inv = -determinv * (a11*a32 - a12*a31);
  double a33inv = determinv * (a11*a22 - a12*a21);

  // quadratic correction coeffs

  double r0102 = (r01[0]*r02[0] + r01[1]*r02[1] + r01[2]*r02[2]);
  double r0103 = (r01[0]*r03[0] + r01[1]*r03[1] + r01[2]*r03[2]);
  double r0203 = (r02[0]*r03[0] + r02[1]*r03[1] + r02[2]*r03[2]);

  double quad1_0101 = (invmass0+invmass1)*(invmass0+invmass1) * r01sq;
  double quad1_0202 = invmass0*invmass0 * r02sq;
  double quad1_0303 = invmass0*invmass0 * r03sq;
  double quad1_0102 = 2.0 * (invmass0+invmass1)*invmass0 * r0102;
  double quad1_0103 = 2.0 * (invmass0+invmass1)*invmass0 * r0103;
  double quad1_0203 = 2.0 * invmass0*invmass0 * r0203;

  double quad2_0101 = invmass0*invmass0 * r01sq;
  double quad2_0202 = (invmass0+invmass2)*(invmass0+invmass2) * r02sq;
  double quad2_0303 = invmass0*invmass0 * r03sq;
  double quad2_0102 = 2.0 * (invmass0+invmass2)*invmass0 * r0102;
  double quad2_0103 = 2.0 * invmass0*invmass0 * r0103;
  double quad2_0203 = 2.0 * (invmass0+invmass2)*invmass0 * r0203;

  double quad3_0101 = invmass0*invmass0 * r01sq;
  double quad3_0202 = invmass0*invmass0 * r02sq;
  double quad3_0303 = (invmass0+invmass3)*(invmass0+invmass3) * r03sq;
  double quad3_0102 = 2.0 * invmass0*invmass0 * r0102;
  double quad3_0103 = 2.0 * (invmass0+invmass3)*invmass0 * r0103;
  double quad3_0203 = 2.0 * (invmass0+invmass3)*invmass0 * r0203;

  // iterate until converged

  double lamda01 = 0.0;
  double lamda02 = 0.0;
  double lamda03 = 0.0;
  int niter = 0;
  int done = 0;

  double quad1,quad2,quad3,b1,b2,b3,lamda01_new,lamda02_new,lamda03_new;

  while (!done && niter < max_iter) {
    quad1 = quad1_0101 * lamda01*lamda01 +
      quad1_0202 * lamda02*lamda02 +
      quad1_0303 * lamda03*lamda03 +
      quad1_0102 * lamda01*lamda02 +
      quad1_0103 * lamda01*lamda03 +
      quad1_0203 * lamda02*lamda03;

    quad2 = quad2_0101 * lamda01*lamda01 +
      quad2_0202 * lamda02*lamda02 +
      quad2_0303 * lamda03*lamda03 +
      quad2_0102 * lamda01*lamda02 +
      quad2_0103 * lamda01*lamda03 +
      quad2_0203 * lamda02*lamda03;

    quad3 = quad3_0101 * lamda01*lamda01 +
      quad3_0202 * lamda02*lamda02 +
      quad3_0303 * lamda03*lamda03 +
      quad3_0102 * lamda01*lamda02 +
      quad3_0103 * lamda01*lamda03 +
      quad3_0203 * lamda02*lamda03;

    b1 = bond1*bond1 - s01sq - quad1;
    b2 = bond2*bond2 - s02sq - quad2;
    b3 = bond3*bond3 - s03sq - quad3;

    lamda01_new = a11inv*b1 + a12inv*b2 + a13inv*b3;
    lamda02_new = a21inv*b1 + a22inv*b2 + a23inv*b3;
    lamda03_new = a31inv*b1 + a32inv*b2 + a33inv*b3;

    done = 1;
    if (fabs(lamda01_new-lamda01) > tolerance) done = 0;
    if (fabs(lamda02_new-lamda02) > tolerance) done = 0;
    if (fabs(lamda03_new-lamda03) > tolerance) done = 0;

    lamda01 = lamda01_new;
    lamda02 = lamda02_new;
    lamda03 = lamda03_new;

    // stop iterations before we have a floating point overflow
    // max double is < 1.0e308, so 1e150 is a reasonable cutoff

    if (fabs(lamda01) > 1e150 || fabs(lamda02) > 1e150
        || fabs(lamda03) > 1e150) done = 1;

    niter++;
  }

  // update forces if atom is owned by this processor

  lamda01 = lamda01/dtfsq;
  lamda02 = lamda02/dtfsq;
  lamda03 = lamda03/dtfsq;

  if (i0 < nlocal) {
    a_f(i0,0) += lamda01*r01[0] + lamda02*r02[0] + lamda03*r03[0];
    a_f(i0,1) += lamda01*r01[1] + lamda02*r02[1] + lamda03*r03[1];
    a_f(i0,2) += lamda01*r01[2] + lamda02*r02[2] + lamda03*r03[2];
  }

  if (i1 < nlocal) {
    a_f(i1,0) -= lamda01*r01[0];
    a_f(i1,1) -= lamda01*r01[1];
    a_f(i1,2) -= lamda01*r01[2];
  }

  if (i2 < nlocal) {
    a_f(i2,0) -= lamda02*r02[0];
    a_f(i2,1) -= lamda02*r02[1];
    a_f(i2,2) -= lamda02*r02[2];
  }

  if (i3 < nlocal) {
    a_f(i3,0) -= lamda03*r03[0];
    a_f(i3,1) -= lamda03*r03[1];
    a_f(i3,2) -= lamda03*r03[2];
  }

  if (EVFLAG) {
    int count = 0;
    if (i0 < nlocal) atomlist[count++] = i0;
    if (i1 < nlocal) atomlist[count++] = i1;
    if (i2 < nlocal) atomlist[count++] = i2;
    if (i3 < nlocal) atomlist[count++] = i3;

    v[0] = lamda01*r01[0]*r01[0]+lamda02*r02[0]*r02[0]+lamda03*r03[0]*r03[0];
    v[1] = lamda01*r01[1]*r01[1]+lamda02*r02[1]*r02[1]+lamda03*r03[1]*r03[1];
    v[2] = lamda01*r01[2]*r01[2]+lamda02*r02[2]*r02[2]+lamda03*r03[2]*r03[2];
    v[3] = lamda01*r01[0]*r01[1]+lamda02*r02[0]*r02[1]+lamda03*r03[0]*r03[1];
    v[4] = lamda01*r01[0]*r01[2]+lamda02*r02[0]*r02[2]+lamda03*r03[0]*r03[2];
    v[5] = lamda01*r01[1]*r01[2]+lamda02*r02[1]*r02[2]+lamda03*r03[1]*r03[2];

    v_tally<NEIGHFLAG>(ev,count,atomlist,4.0,v);
  }
}

/* ----------------------------------------------------------------------
   calculate SHAKE constraint forces for size 3 cluster = two bonds + angle
------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void FixShakeKokkos<DeviceType>::shake3angle(int ilist, EV_FLOAT& ev) const
{

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  int atomlist[3];
  double v[6];
  double invmass0,invmass1,invmass2;

  // local atom IDs and constraint distances

  int m = d_list[ilist];
  int i0 = d_closest_list(ilist,0);
  int i1 = d_closest_list(ilist,1);
  int i2 = d_closest_list(ilist,2);
  double bond1 = d_bond_distance[d_shake_type(m,0)];
  double bond2 = d_bond_distance[d_shake_type(m,1)];
  double bond12 = d_angle_distance[d_shake_type(m,2)];

  // r01,r02,r12 = distance vec between atoms

  double r01[3];
  r01[0] = d_x(i0,0) - d_x(i1,0);
  r01[1] = d_x(i0,1) - d_x(i1,1);
  r01[2] = d_x(i0,2) - d_x(i1,2);

  double r02[3];
  r02[0] = d_x(i0,0) - d_x(i2,0);
  r02[1] = d_x(i0,1) - d_x(i2,1);
  r02[2] = d_x(i0,2) - d_x(i2,2);

  double r12[3];
  r12[0] = d_x(i1,0) - d_x(i2,0);
  r12[1] = d_x(i1,1) - d_x(i2,1);
  r12[2] = d_x(i1,2) - d_x(i2,2);

  // s01,s02,s12 = distance vec after unconstrained update

  double s01[3];
  s01[0] = d_xshake(i0,0) - d_xshake(i1,0);
  s01[1] = d_xshake(i0,1) - d_xshake(i1,1);
  s01[2] = d_xshake(i0,2) - d_xshake(i1,2);

  double s02[3];
  s02[0] = d_xshake(i0,0) - d_xshake(i2,0);
  s02[1] = d_xshake(i0,1) - d_xshake(i2,1);
  s02[2] = d_xshake(i0,2) - d_xshake(i2,2);

  double s12[3];
  s12[0] = d_xshake(i1,0) - d_xshake(i2,0);
  s12[1] = d_xshake(i1,1) - d_xshake(i2,1);
  s12[2] = d_xshake(i1,2) - d_xshake(i2,2);

  // scalar distances between atoms

  double r01sq = r01[0]*r01[0] + r01[1]*r01[1] + r01[2]*r01[2];
  double r02sq = r02[0]*r02[0] + r02[1]*r02[1] + r02[2]*r02[2];
  double r12sq = r12[0]*r12[0] + r12[1]*r12[1] + r12[2]*r12[2];
  double s01sq = s01[0]*s01[0] + s01[1]*s01[1] + s01[2]*s01[2];
  double s02sq = s02[0]*s02[0] + s02[1]*s02[1] + s02[2]*s02[2];
  double s12sq = s12[0]*s12[0] + s12[1]*s12[1] + s12[2]*s12[2];

  // matrix coeffs and rhs for lamda equations

  if (d_rmass.data()) {
    invmass0 = 1.0/d_rmass[i0];
    invmass1 = 1.0/d_rmass[i1];
    invmass2 = 1.0/d_rmass[i2];
  } else {
    invmass0 = 1.0/d_mass[d_type[i0]];
    invmass1 = 1.0/d_mass[d_type[i1]];
    invmass2 = 1.0/d_mass[d_type[i2]];
  }

  double a11 = 2.0 * (invmass0+invmass1) *
    (s01[0]*r01[0] + s01[1]*r01[1] + s01[2]*r01[2]);
  double a12 = 2.0 * invmass0 *
    (s01[0]*r02[0] + s01[1]*r02[1] + s01[2]*r02[2]);
  double a13 = - 2.0 * invmass1 *
    (s01[0]*r12[0] + s01[1]*r12[1] + s01[2]*r12[2]);
  double a21 = 2.0 * invmass0 *
    (s02[0]*r01[0] + s02[1]*r01[1] + s02[2]*r01[2]);
  double a22 = 2.0 * (invmass0+invmass2) *
    (s02[0]*r02[0] + s02[1]*r02[1] + s02[2]*r02[2]);
  double a23 = 2.0 * invmass2 *
    (s02[0]*r12[0] + s02[1]*r12[1] + s02[2]*r12[2]);
  double a31 = - 2.0 * invmass1 *
    (s12[0]*r01[0] + s12[1]*r01[1] + s12[2]*r01[2]);
  double a32 = 2.0 * invmass2 *
    (s12[0]*r02[0] + s12[1]*r02[1] + s12[2]*r02[2]);
  double a33 = 2.0 * (invmass1+invmass2) *
    (s12[0]*r12[0] + s12[1]*r12[1] + s12[2]*r12[2]);

  // inverse of matrix

  double determ = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 -
    a11*a23*a32 - a12*a21*a33 - a13*a22*a31;
  if (determ == 0.0) d_error_flag() = 3;
  //error->one(FLERR,"Shake determinant = 0.0");
  double determinv = 1.0/determ;

  double a11inv = determinv * (a22*a33 - a23*a32);
  double a12inv = -determinv * (a12*a33 - a13*a32);
  double a13inv = determinv * (a12*a23 - a13*a22);
  double a21inv = -determinv * (a21*a33 - a23*a31);
  double a22inv = determinv * (a11*a33 - a13*a31);
  double a23inv = -determinv * (a11*a23 - a13*a21);
  double a31inv = determinv * (a21*a32 - a22*a31);
  double a32inv = -determinv * (a11*a32 - a12*a31);
  double a33inv = determinv * (a11*a22 - a12*a21);

  // quadratic correction coeffs

  double r0102 = (r01[0]*r02[0] + r01[1]*r02[1] + r01[2]*r02[2]);
  double r0112 = (r01[0]*r12[0] + r01[1]*r12[1] + r01[2]*r12[2]);
  double r0212 = (r02[0]*r12[0] + r02[1]*r12[1] + r02[2]*r12[2]);

  double quad1_0101 = (invmass0+invmass1)*(invmass0+invmass1) * r01sq;
  double quad1_0202 = invmass0*invmass0 * r02sq;
  double quad1_1212 = invmass1*invmass1 * r12sq;
  double quad1_0102 = 2.0 * (invmass0+invmass1)*invmass0 * r0102;
  double quad1_0112 = - 2.0 * (invmass0+invmass1)*invmass1 * r0112;
  double quad1_0212 = - 2.0 * invmass0*invmass1 * r0212;

  double quad2_0101 = invmass0*invmass0 * r01sq;
  double quad2_0202 = (invmass0+invmass2)*(invmass0+invmass2) * r02sq;
  double quad2_1212 = invmass2*invmass2 * r12sq;
  double quad2_0102 = 2.0 * (invmass0+invmass2)*invmass0 * r0102;
  double quad2_0112 = 2.0 * invmass0*invmass2 * r0112;
  double quad2_0212 = 2.0 * (invmass0+invmass2)*invmass2 * r0212;

  double quad3_0101 = invmass1*invmass1 * r01sq;
  double quad3_0202 = invmass2*invmass2 * r02sq;
  double quad3_1212 = (invmass1+invmass2)*(invmass1+invmass2) * r12sq;
  double quad3_0102 = - 2.0 * invmass1*invmass2 * r0102;
  double quad3_0112 = - 2.0 * (invmass1+invmass2)*invmass1 * r0112;
  double quad3_0212 = 2.0 * (invmass1+invmass2)*invmass2 * r0212;

  // iterate until converged

  double lamda01 = 0.0;
  double lamda02 = 0.0;
  double lamda12 = 0.0;
  int niter = 0;
  int done = 0;

  double quad1,quad2,quad3,b1,b2,b3,lamda01_new,lamda02_new,lamda12_new;

  while (!done && niter < max_iter) {

    quad1 = quad1_0101 * lamda01*lamda01 +
      quad1_0202 * lamda02*lamda02 +
      quad1_1212 * lamda12*lamda12 +
      quad1_0102 * lamda01*lamda02 +
      quad1_0112 * lamda01*lamda12 +
      quad1_0212 * lamda02*lamda12;

    quad2 = quad2_0101 * lamda01*lamda01 +
      quad2_0202 * lamda02*lamda02 +
      quad2_1212 * lamda12*lamda12 +
      quad2_0102 * lamda01*lamda02 +
      quad2_0112 * lamda01*lamda12 +
      quad2_0212 * lamda02*lamda12;

    quad3 = quad3_0101 * lamda01*lamda01 +
      quad3_0202 * lamda02*lamda02 +
      quad3_1212 * lamda12*lamda12 +
      quad3_0102 * lamda01*lamda02 +
      quad3_0112 * lamda01*lamda12 +
      quad3_0212 * lamda02*lamda12;

    b1 = bond1*bond1 - s01sq - quad1;
    b2 = bond2*bond2 - s02sq - quad2;
    b3 = bond12*bond12 - s12sq - quad3;

    lamda01_new = a11inv*b1 + a12inv*b2 + a13inv*b3;
    lamda02_new = a21inv*b1 + a22inv*b2 + a23inv*b3;
    lamda12_new = a31inv*b1 + a32inv*b2 + a33inv*b3;

    done = 1;
    if (fabs(lamda01_new-lamda01) > tolerance) done = 0;
    if (fabs(lamda02_new-lamda02) > tolerance) done = 0;
    if (fabs(lamda12_new-lamda12) > tolerance) done = 0;

    lamda01 = lamda01_new;
    lamda02 = lamda02_new;
    lamda12 = lamda12_new;

    // stop iterations before we have a floating point overflow
    // max double is < 1.0e308, so 1e150 is a reasonable cutoff

    if (fabs(lamda01) > 1e150 || fabs(lamda02) > 1e150
        || fabs(lamda12) > 1e150) done = 1;

    niter++;
  }

  // update forces if atom is owned by this processor

  lamda01 = lamda01/dtfsq;
  lamda02 = lamda02/dtfsq;
  lamda12 = lamda12/dtfsq;

  if (i0 < nlocal) {
    a_f(i0,0) += lamda01*r01[0] + lamda02*r02[0];
    a_f(i0,1) += lamda01*r01[1] + lamda02*r02[1];
    a_f(i0,2) += lamda01*r01[2] + lamda02*r02[2];
  }

  if (i1 < nlocal) {
    a_f(i1,0) -= lamda01*r01[0] - lamda12*r12[0];
    a_f(i1,1) -= lamda01*r01[1] - lamda12*r12[1];
    a_f(i1,2) -= lamda01*r01[2] - lamda12*r12[2];
  }

  if (i2 < nlocal) {
    a_f(i2,0) -= lamda02*r02[0] + lamda12*r12[0];
    a_f(i2,1) -= lamda02*r02[1] + lamda12*r12[1];
    a_f(i2,2) -= lamda02*r02[2] + lamda12*r12[2];
  }

  if (EVFLAG) {
    int count = 0;
    if (i0 < nlocal) atomlist[count++] = i0;
    if (i1 < nlocal) atomlist[count++] = i1;
    if (i2 < nlocal) atomlist[count++] = i2;

    v[0] = lamda01*r01[0]*r01[0]+lamda02*r02[0]*r02[0]+lamda12*r12[0]*r12[0];
    v[1] = lamda01*r01[1]*r01[1]+lamda02*r02[1]*r02[1]+lamda12*r12[1]*r12[1];
    v[2] = lamda01*r01[2]*r01[2]+lamda02*r02[2]*r02[2]+lamda12*r12[2]*r12[2];
    v[3] = lamda01*r01[0]*r01[1]+lamda02*r02[0]*r02[1]+lamda12*r12[0]*r12[1];
    v[4] = lamda01*r01[0]*r01[2]+lamda02*r02[0]*r02[2]+lamda12*r12[0]*r12[2];
    v[5] = lamda01*r01[1]*r01[2]+lamda02*r02[1]*r02[2]+lamda12*r12[1]*r12[2];

    v_tally<NEIGHFLAG>(ev,count,atomlist,3.0,v);
  }
}

/* ----------------------------------------------------------------------
   print-out bond & angle statistics
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::stats()
{
  atomKK->sync(Host,X_MASK);
  k_shake_flag.sync_host();
  k_shake_atom.sync_host();
  k_shake_type.sync_host();
  k_list.sync_host();
  k_closest_list.sync_host();

  FixShake::stats();
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::grow_arrays(int nmax)
{
  memoryKK->grow_kokkos(k_shake_flag,shake_flag,nmax,"shake:shake_flag");
  memoryKK->grow_kokkos(k_shake_atom,shake_atom,nmax,4,"shake:shake_atom");
  memoryKK->grow_kokkos(k_shake_type,shake_type,nmax,3,"shake:shake_type");
  memoryKK->destroy_kokkos(k_xshake,xshake);
  memoryKK->create_kokkos(k_xshake,xshake,nmax,"shake:xshake");

  d_shake_flag = k_shake_flag.view<DeviceType>();
  d_shake_atom = k_shake_atom.view<DeviceType>();
  d_shake_type = k_shake_type.view<DeviceType>();
  d_xshake = k_xshake.view<DeviceType>();

  memory->destroy(ftmp);
  memory->create(ftmp,nmax,3,"shake:ftmp");
  memory->destroy(vtmp);
  memory->create(vtmp,nmax,3,"shake:vtmp");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::copy_arrays(int i, int j, int delflag)
{
  k_shake_flag.sync_host();
  k_shake_atom.sync_host();
  k_shake_type.sync_host();

  FixShake::copy_arrays(i,j,delflag);

  k_shake_flag.modify_host();
  k_shake_atom.modify_host();
  k_shake_type.modify_host();
}

/* ----------------------------------------------------------------------
   sort local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  // always sort on the device

  k_shake_flag.sync_device();
  k_shake_atom.sync_device();
  k_shake_type.sync_device();

  Sorter.sort(LMPDeviceType(), k_shake_flag.d_view);
  Sorter.sort(LMPDeviceType(), k_shake_atom.d_view);
  Sorter.sort(LMPDeviceType(), k_shake_type.d_view);

  k_shake_flag.modify_device();
  k_shake_atom.modify_device();
  k_shake_type.modify_device();
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::set_arrays(int i)
{
  k_shake_flag.sync_host();

  shake_flag[i] = 0;

  k_shake_flag.modify_host();
}

/* ----------------------------------------------------------------------
   update one atom's array values
   called when molecule is created from fix gcmc
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::update_arrays(int i, int atom_offset)
{
  k_shake_flag.sync_host();
  k_shake_atom.sync_host();

  FixShake::update_arrays(i,atom_offset);

  k_shake_flag.modify_host();
  k_shake_atom.modify_host();
}

/* ----------------------------------------------------------------------
   initialize a molecule inserted by another fix, e.g. deposit or pour
   called when molecule is created
   nlocalprev = # of atoms on this proc before molecule inserted
   tagprev = atom ID previous to new atoms in the molecule
   xgeom,vcm,quat ignored
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::set_molecule(int nlocalprev, tagint tagprev, int imol,
                            double * xgeom, double * vcm, double * quat)
{
  atomKK->sync(Host,TAG_MASK|MOLECULE_MASK);
  k_shake_flag.sync_host();
  k_shake_atom.sync_host();
  k_shake_type.sync_host();

  FixShake::set_molecule(nlocalprev,tagprev,imol,xgeom,vcm,quat);

  k_shake_flag.modify_host();
  k_shake_atom.modify_host();
  k_shake_type.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixShakeKokkos<DeviceType>::pack_exchange_item(const int &mysend, int &offset, const bool &final) const
{
  const int i = d_exchange_sendlist(mysend);
  int flag = d_shake_flag[i];

  if (!final) {
    if (flag == 1) offset += 7;
    else if (flag == 2) offset += 4;
    else if (flag == 3) offset += 6;
    else if (flag == 4) offset += 8;
    else offset++;
  } else {

    d_buf[mysend] = nsend + offset;
    int m = nsend + offset;
    d_buf[m++] = flag;
    if (flag == 1) {
      d_buf[m++] = d_shake_atom(i,0);
      d_buf[m++] = d_shake_atom(i,1);
      d_buf[m++] = d_shake_atom(i,2);
      d_buf[m++] = d_shake_type(i,0);
      d_buf[m++] = d_shake_type(i,1);
      d_buf[m++] = d_shake_type(i,2);
    } else if (flag == 2) {
      d_buf[m++] = d_shake_atom(i,0);
      d_buf[m++] = d_shake_atom(i,1);
      d_buf[m++] = d_shake_type(i,0);
    } else if (flag == 3) {
      d_buf[m++] = d_shake_atom(i,0);
      d_buf[m++] = d_shake_atom(i,1);
      d_buf[m++] = d_shake_atom(i,2);
      d_buf[m++] = d_shake_type(i,0);
      d_buf[m++] = d_shake_type(i,1);
    } else if (flag == 4) {
      d_buf[m++] = d_shake_atom(i,0);
      d_buf[m++] = d_shake_atom(i,1);
      d_buf[m++] = d_shake_atom(i,2);
      d_buf[m++] = d_shake_atom(i,3);
      d_buf[m++] = d_shake_type(i,0);
      d_buf[m++] = d_shake_type(i,1);
      d_buf[m++] = d_shake_type(i,2);
    }
    if (mysend == nsend-1) d_count() = m;
    offset = m - nsend;

    const int j = d_copylist(mysend);
    if (j > -1) {
      d_shake_flag[i] = d_shake_flag[j];
      int flag = d_shake_flag[i];
      if (flag == 1) {
        d_shake_atom(i,0) = d_shake_atom(j,0);
        d_shake_atom(i,1) = d_shake_atom(j,1);
        d_shake_atom(i,2) = d_shake_atom(j,2);
        d_shake_type(i,0) = d_shake_type(j,0);
        d_shake_type(i,1) = d_shake_type(j,1);
        d_shake_type(i,2) = d_shake_type(j,2);
      } else if (flag == 2) {
        d_shake_atom(i,0) = d_shake_atom(j,0);
        d_shake_atom(i,1) = d_shake_atom(j,1);
        d_shake_type(i,0) = d_shake_type(j,0);
      } else if (flag == 3) {
        d_shake_atom(i,0) = d_shake_atom(j,0);
        d_shake_atom(i,1) = d_shake_atom(j,1);
        d_shake_atom(i,2) = d_shake_atom(j,2);
        d_shake_type(i,0) = d_shake_type(j,0);
        d_shake_type(i,1) = d_shake_type(j,1);
      } else if (flag == 4) {
        d_shake_atom(i,0) = d_shake_atom(j,0);
        d_shake_atom(i,1) = d_shake_atom(j,1);
        d_shake_atom(i,2) = d_shake_atom(j,2);
        d_shake_atom(i,3) = d_shake_atom(j,3);
        d_shake_type(i,0) = d_shake_type(j,0);
        d_shake_type(i,1) = d_shake_type(j,1);
        d_shake_type(i,2) = d_shake_type(j,2);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixShakeKokkos<DeviceType>::pack_exchange_kokkos(
   const int &nsend, DAT::tdual_xfloat_2d &k_buf,
   DAT::tdual_int_1d k_exchange_sendlist, DAT::tdual_int_1d k_copylist,
   ExecutionSpace space)
{
  k_buf.sync<DeviceType>();
  k_copylist.sync<DeviceType>();
  k_exchange_sendlist.sync<DeviceType>();

  d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  d_copylist = k_copylist.view<DeviceType>();
  d_exchange_sendlist = k_exchange_sendlist.view<DeviceType>();
  this->nsend = nsend;

  k_shake_flag.template sync<DeviceType>();
  k_shake_atom.template sync<DeviceType>();
  k_shake_type.template sync<DeviceType>();

  Kokkos::deep_copy(d_count,0);

  copymode = 1;

  FixShakeKokkosPackExchangeFunctor<DeviceType> pack_exchange_functor(this);
  Kokkos::parallel_scan(nsend,pack_exchange_functor);

  copymode = 0;

  k_buf.modify<DeviceType>();

  if (space == Host) k_buf.sync<LMPHostType>();
  else k_buf.sync<LMPDeviceType>();

  k_shake_flag.template modify<DeviceType>();
  k_shake_atom.template modify<DeviceType>();
  k_shake_type.template modify<DeviceType>();

  Kokkos::deep_copy(h_count,d_count);

  return h_count();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixShakeKokkos<DeviceType>::operator()(TagFixShakeUnpackExchange, const int &i) const
{
  int index = d_indices(i);

  if (index > -1) {
    int m = d_buf[i];

    int flag = d_shake_flag[index] = static_cast<int> (d_buf[m++]);
    if (flag == 1) {
      d_shake_atom(index,0) = static_cast<tagint> (d_buf[m++]);
      d_shake_atom(index,1) = static_cast<tagint> (d_buf[m++]);
      d_shake_atom(index,2) = static_cast<tagint> (d_buf[m++]);
      d_shake_type(index,0) = static_cast<int> (d_buf[m++]);
      d_shake_type(index,1) = static_cast<int> (d_buf[m++]);
      d_shake_type(index,2) = static_cast<int> (d_buf[m++]);
    } else if (flag == 2) {
      d_shake_atom(index,0) = static_cast<tagint> (d_buf[m++]);
      d_shake_atom(index,1) = static_cast<tagint> (d_buf[m++]);
      d_shake_type(index,0) = static_cast<int> (d_buf[m++]);
    } else if (flag == 3) {
      d_shake_atom(index,0) = static_cast<tagint> (d_buf[m++]);
      d_shake_atom(index,1) = static_cast<tagint> (d_buf[m++]);
      d_shake_atom(index,2) = static_cast<tagint> (d_buf[m++]);
      d_shake_type(index,0) = static_cast<int> (d_buf[m++]);
      d_shake_type(index,1) = static_cast<int> (d_buf[m++]);
    } else if (flag == 4) {
      d_shake_atom(index,0) = static_cast<tagint> (d_buf[m++]);
      d_shake_atom(index,1) = static_cast<tagint> (d_buf[m++]);
      d_shake_atom(index,2) = static_cast<tagint> (d_buf[m++]);
      d_shake_atom(index,3) = static_cast<tagint> (d_buf[m++]);
      d_shake_type(index,0) = static_cast<int> (d_buf[m++]);
      d_shake_type(index,1) = static_cast<int> (d_buf[m++]);
      d_shake_type(index,2) = static_cast<int> (d_buf[m++]);
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixShakeKokkos<DeviceType>::unpack_exchange_kokkos(
  DAT::tdual_xfloat_2d &k_buf, DAT::tdual_int_1d &k_indices, int nrecv,
  ExecutionSpace /*space*/)
{
  k_buf.sync<DeviceType>();
  k_indices.sync<DeviceType>();

  d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  d_indices = k_indices.view<DeviceType>();

  k_shake_flag.template sync<DeviceType>();
  k_shake_atom.template sync<DeviceType>();
  k_shake_type.template sync<DeviceType>();

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixShakeUnpackExchange>(0,nrecv),*this);

  copymode = 0;

  k_shake_flag.template modify<DeviceType>();
  k_shake_atom.template modify<DeviceType>();
  k_shake_type.template modify<DeviceType>();
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixShakeKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_shake_flag.sync_host();
  k_shake_atom.sync_host();
  k_shake_type.sync_host();

  int m = FixShake::pack_exchange(i,buf);

  k_shake_flag.modify_host();
  k_shake_atom.modify_host();
  k_shake_type.modify_host();

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixShakeKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  k_shake_flag.sync_host();
  k_shake_atom.sync_host();
  k_shake_type.sync_host();

  int m = FixShake::unpack_exchange(nlocal,buf);

  k_shake_flag.modify_host();
  k_shake_atom.modify_host();
  k_shake_type.modify_host();

  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixShakeKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                                                        int iswap_in, DAT::tdual_xfloat_1d &k_buf,
                                                        int pbc_flag, int* pbc)
{
  d_sendlist = k_sendlist.view<DeviceType>();
  iswap = iswap_in;
  d_buf = k_buf.view<DeviceType>();

  if (domain->triclinic == 0) {
    dx = pbc[0]*domain->xprd;
    dy = pbc[1]*domain->yprd;
    dz = pbc[2]*domain->zprd;
  } else {
    dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
    dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
    dz = pbc[2]*domain->zprd;
  }

  if (pbc_flag)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixShakePackForwardComm<1> >(0,n),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixShakePackForwardComm<0> >(0,n),*this);
  return n*3;
}

template<class DeviceType>
template<int PBC_FLAG>
KOKKOS_INLINE_FUNCTION
void FixShakeKokkos<DeviceType>::operator()(TagFixShakePackForwardComm<PBC_FLAG>, const int &i) const {
  const int j = d_sendlist(iswap, i);

  if (PBC_FLAG == 0) {
    d_buf[3*i] = d_xshake(j,0);
    d_buf[3*i+1] = d_xshake(j,1);
    d_buf[3*i+2] = d_xshake(j,2);
  } else {
    d_buf[3*i] = d_xshake(j,0) + dx;
    d_buf[3*i+1] = d_xshake(j,1) + dy;
    d_buf[3*i+2] = d_xshake(j,2) + dz;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixShakeKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
                                int pbc_flag, int *pbc)
{
  k_xshake.sync_host();

  int m = FixShake::pack_forward_comm(n,list,buf,pbc_flag,pbc);

  k_xshake.modify_host();

  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first_in, DAT::tdual_xfloat_1d &buf)
{
  first = first_in;
  d_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixShakeUnpackForwardComm>(0,n),*this);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixShakeKokkos<DeviceType>::operator()(TagFixShakeUnpackForwardComm, const int &i) const {
  d_xshake(i + first,0) = d_buf[3*i];
  d_xshake(i + first,1) = d_buf[3*i+1];
  d_xshake(i + first,2) = d_buf[3*i+2];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  k_xshake.sync_host();

  FixShake::unpack_forward_comm(n,first,buf);

  k_xshake.modify_host();
}

/* ----------------------------------------------------------------------
   add coordinate constraining forces
   this method is called at the end of a timestep
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::shake_end_of_step(int vflag) {
  dtv     = update->dt;
  dtfsq   = 0.5 * update->dt * update->dt * force->ftm2v;
  FixShakeKokkos<DeviceType>::post_force(vflag);
  if (!rattle) dtfsq = update->dt * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   calculate constraining forces based on the current configuration
   change coordinates
------------------------------------------------------------------------- */

template<class DeviceType>
void FixShakeKokkos<DeviceType>::correct_coordinates(int vflag) {
  x = atom->x;
  v = atom->v;
  f = atom->f;
  mass = atom->mass;
  rmass = atom->rmass;
  type = atom->type;

  atomKK->sync(Host,X_MASK|V_MASK|F_MASK);

  // save current forces and velocities so that you
  // initialize them to zero such that FixShake::unconstrained_coordinate_update has no effect

  for (int j=0; j<nlocal; j++) {
    for (int k=0; k<3; k++) {

      // store current value of forces and velocities
      ftmp[j][k] = f[j][k];
      vtmp[j][k] = v[j][k];

      // set f and v to zero for SHAKE

      v[j][k] = 0;
      f[j][k] = 0;
    }
  }

  atomKK->modified(Host,V_MASK|F_MASK);

  // call SHAKE to correct the coordinates which were updated without constraints
  // IMPORTANT: use 1 as argument and thereby enforce velocity Verlet

  dtfsq   = 0.5 * update->dt * update->dt * force->ftm2v;
  FixShakeKokkos<DeviceType>::post_force(vflag);

  atomKK->sync(Host,X_MASK|F_MASK);

  // integrate coordinates: x' = xnp1 + dt^2/2m_i * f, where f is the constraining force
  // NOTE: After this command, the coordinates geometry of the molecules will be correct!

  double dtfmsq;
  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      dtfmsq = dtfsq/ rmass[i];
      x[i][0] = x[i][0] + dtfmsq*f[i][0];
      x[i][1] = x[i][1] + dtfmsq*f[i][1];
      x[i][2] = x[i][2] + dtfmsq*f[i][2];
    }
  }
  else {
    for (int i = 0; i < nlocal; i++) {
      dtfmsq = dtfsq / mass[type[i]];
      x[i][0] = x[i][0] + dtfmsq*f[i][0];
      x[i][1] = x[i][1] + dtfmsq*f[i][1];
      x[i][2] = x[i][2] + dtfmsq*f[i][2];
    }
  }

  // copy forces and velocities back

  for (int j=0; j<nlocal; j++) {
    for (int k=0; k<3; k++) {
      f[j][k] = ftmp[j][k];
      v[j][k] = vtmp[j][k];
    }
  }

  if (!rattle) dtfsq = update->dt * update->dt * force->ftm2v;

  // communicate changes
  // NOTE: for compatibility xshake is temporarily set to x, such that pack/unpack_forward
  //       can be used for communicating the coordinates.

  double **xtmp = xshake;
  xshake = x;
  forward_comm_device = 0;
  comm->forward_comm(this);
  forward_comm_device = 1;
  xshake = xtmp;

  atomKK->modified(Host,X_MASK|V_MASK|F_MASK);
}

/* ----------------------------------------------------------------------
   tally virial into global and per-atom accumulators
   n = # of local owned atoms involved, with local indices in list
   v = total virial for the interaction involving total atoms
   increment global virial by n/total fraction
   increment per-atom virial of each atom in list by 1/total fraction
   this method can be used when fix computes forces in post_force()
     e.g. fix shake, fix rigid: compute virial only on owned atoms
       whether newton_bond is on or off
     other procs will tally left-over fractions for atoms they own
------------------------------------------------------------------------- */
template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void FixShakeKokkos<DeviceType>::v_tally(EV_FLOAT &ev, int n, int *atomlist, double total,
     double *v) const
{
  int m;

  if (vflag_global) {
    double fraction = n/total;
    ev.v[0] += fraction*v[0];
    ev.v[1] += fraction*v[1];
    ev.v[2] += fraction*v[2];
    ev.v[3] += fraction*v[3];
    ev.v[4] += fraction*v[4];
    ev.v[5] += fraction*v[5];
  }

  if (vflag_atom) {
    double fraction = 1.0/total;
    for (int i = 0; i < n; i++) {
      auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
      auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
      m = atomlist[i];
      a_vatom(m,0) += fraction*v[0];
      a_vatom(m,1) += fraction*v[1];
      a_vatom(m,2) += fraction*v[2];
      a_vatom(m,3) += fraction*v[3];
      a_vatom(m,4) += fraction*v[4];
      a_vatom(m,5) += fraction*v[5];
    }
  }
}

/* ----------------------------------------------------------------------
   return local index of atom J or any of its images that is closest to atom I
   if J is not a valid index like -1, just return it
   copied from domain.cpp
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int FixShakeKokkos<DeviceType>::closest_image(const int i, int j) const
{
  if (j < 0) return j;

  const X_FLOAT xi0 = d_x(i,0);
  const X_FLOAT xi1 = d_x(i,1);
  const X_FLOAT xi2 = d_x(i,2);

  int closest = j;
  X_FLOAT delx = xi0 - d_x(j,0);
  X_FLOAT dely = xi1 - d_x(j,1);
  X_FLOAT delz = xi2 - d_x(j,2);
  X_FLOAT rsqmin = delx*delx + dely*dely + delz*delz;
  X_FLOAT rsq;

  while (d_sametag[j] >= 0) {
    j = d_sametag[j];
    delx = xi0 - d_x(j,0);
    dely = xi1 - d_x(j,1);
    delz = xi2 - d_x(j,2);
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < rsqmin) {
      rsqmin = rsq;
      closest = j;
    }
  }

  return closest;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixShakeKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixShakeKokkos<LMPHostType>;
#endif
}

