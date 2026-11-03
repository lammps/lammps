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

/* ----------------------------------------------------------------------
   Contributing authors: Trung Nguyen (U Chicago)
                         Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "fix_efield_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "error.h"
#include "input.h"
#include "kokkos_base.h"
#include "memory_kokkos.h"
#include "modify_kokkos.h"
#include "region.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixEfieldKokkos<DeviceType>::FixEfieldKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixEfield(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TORQUE_MASK | Q_MASK | MU_MASK | IMAGE_MASK | MASK_MASK;
  datamask_modify = F_MASK | TORQUE_MASK;

  memory->destroy(efield);
  memoryKK->create_kokkos(k_efield,efield,maxatom,4,"efield:efield");
  d_efield = k_efield.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixEfieldKokkos<DeviceType>::~FixEfieldKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_efield,efield);
  memoryKK->destroy_kokkos(k_vatom,vatom);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixEfieldKokkos<DeviceType>::init()
{
  FixEfield::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixEfieldKokkos<DeviceType>::post_force(int vflag)
{
  atomKK->sync(execution_space, datamask_read);

  d_x = atomKK->k_x.template view<DeviceType>();
  d_f = atomKK->k_f.template view<DeviceType>();
  d_q = atomKK->k_q.template view<DeviceType>();
  d_mu = atomKK->k_mu.template view<DeviceType>();
  d_torque = atomKK->k_torque.template view<DeviceType>();
  d_image = atomKK->k_image.template view<DeviceType>();
  d_mask = atomKK->k_mask.template view<DeviceType>();
  int nlocal = atomKK->nlocal;

  // virial setup

  // the per-atom virial is accumulated into a dual view, so the plain
  // base-class vatom array must not be allocated here (alloc = 0)

  v_init(vflag,0);

  // reallocate the per-atom virial dual view if necessary

  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"efield:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  // update region if necessary

  if (region) {
    KokkosBase* regionKKBase = dynamic_cast<KokkosBase*>(region);
    if (!regionKKBase)
      error->all(FLERR,"Cannot use fix efield/kk with region style {} that has no KOKKOS support",region->style);
    region->prematch();
    DAT::tdual_int_1d k_match = DAT::tdual_int_1d("efield:k_match",nlocal);
    regionKKBase->match_all_kokkos(groupbit,k_match);
    k_match.template sync<DeviceType>();
    d_match = k_match.template view<DeviceType>();
  }

  // reallocate efield array if necessary
  // an atom-style energy or potential variable writes into efield[i][3],
  // so the array is needed for those as well, not only for varflag == ATOM

  if (((varflag == ATOM) || (estyle == ATOM) || (pstyle == ATOM)) && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memoryKK->destroy_kokkos(k_efield,efield);
    memoryKK->create_kokkos(k_efield,efield,maxatom,4,"efield:efield");
    d_efield = k_efield.view<DeviceType>();
  }

  force_flag = 0;
  double result[10] = {0.0};

  if (varflag == CONSTANT) {

    prd = domain->prd;
    h = domain->h;
    triclinic = domain->triclinic;
    copymode = 1;

    if(qflag && muflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEfieldConstant<1,1> >(0,nlocal),*this,result);
    else if(qflag && !muflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEfieldConstant<1,0> >(0,nlocal),*this,result);
    else if(!qflag && muflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEfieldConstant<0,1> >(0,nlocal),*this,result);
    else
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEfieldConstant<0,0> >(0,nlocal),*this,result);

    copymode = 0;

  // variable force, wrap with clear/add

  } else {

    // the variables are evaluated on the host through atom->x, atom->v, ...
    // so the host copies have to be up to date before the base class runs

    atomKK->sync(Host, ALL_MASK);

    FixEfield::update_efield_variables();

    // atom-style variables are evaluated on the host, so the result has to be
    // copied to the device for the kernel below.  this is a real copy, not a
    // stale flag: it can only go away if variables are evaluated on the device.
    // an atom-style energy or potential variable writes into efield[i][3], so
    // it needs the copy as well

    if ((varflag == ATOM) || (estyle == ATOM) || (pstyle == ATOM)) {
      k_efield.modify_host();
      k_efield.sync<DeviceType>();
    }

    copymode = 1;

    if(qflag && muflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEfieldNonConstant<1,1> >(0,nlocal),*this,result);
    else if(qflag && !muflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEfieldNonConstant<1,0> >(0,nlocal),*this,result);
    else if(!qflag && muflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEfieldNonConstant<0,1> >(0,nlocal),*this,result);
    else
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEfieldNonConstant<0,0> >(0,nlocal),*this,result);

    copymode = 0;

  }

  atomKK->modified(execution_space, datamask_modify);

  fsum[0] = result[0];
  fsum[1] = result[1];
  fsum[2] = result[2];
  fsum[3] = result[3];

  if (vflag_global) {
    virial[0] += result[4];
    virial[1] += result[5];
    virial[2] += result[6];
    virial[3] += result[7];
    virial[4] += result[8];
    virial[5] += result[9];
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }
}

template<class DeviceType>
template<int QFLAG, int MUFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixEfieldKokkos<DeviceType>::operator()(TagFixEfieldConstant<QFLAG,MUFLAG>, const int &i, value_type result) const {
  if ( QFLAG && (d_mask(i) & groupbit)) {
    if (region && !d_match[i]) return;

    Few<double,3> x_i;
    x_i[0] = static_cast<double>(d_x(i,0));
    x_i[1] = static_cast<double>(d_x(i,1));
    x_i[2] = static_cast<double>(d_x(i,2));
    auto unwrapKK = DomainKokkos::unmap(prd,h,triclinic,x_i,d_image(i));
    const KK_FLOAT fx = d_q(i) * static_cast<KK_FLOAT>(ex);
    const KK_FLOAT fy = d_q(i) * static_cast<KK_FLOAT>(ey);
    const KK_FLOAT fz = d_q(i) * static_cast<KK_FLOAT>(ez);
    d_f(i,0) += static_cast<KK_ACC_FLOAT>(fx);
    d_f(i,1) += static_cast<KK_ACC_FLOAT>(fy);
    d_f(i,2) += static_cast<KK_ACC_FLOAT>(fz);
    result[0] -= static_cast<double>(fx) * unwrapKK[0] + static_cast<double>(fy) * unwrapKK[1] + static_cast<double>(fz) * unwrapKK[2];
    result[1] += static_cast<double>(fx);
    result[2] += static_cast<double>(fy);
    result[3] += static_cast<double>(fz);

    if (evflag) {
      KK_FLOAT v[6];
      v[0] = fx * static_cast<KK_FLOAT>(unwrapKK[0]);
      v[1] = fy * static_cast<KK_FLOAT>(unwrapKK[1]);
      v[2] = fz * static_cast<KK_FLOAT>(unwrapKK[2]);
      v[3] = fx * static_cast<KK_FLOAT>(unwrapKK[1]);
      v[4] = fx * static_cast<KK_FLOAT>(unwrapKK[2]);
      v[5] = fy * static_cast<KK_FLOAT>(unwrapKK[2]);
      v_tally(result, i, v);
    }

  }

  if (MUFLAG && (d_mask(i) & groupbit)) {
    if (region && !d_match[i]) return;
    d_torque(i,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(ez) * d_mu(i,1) - static_cast<KK_FLOAT>(ey) * d_mu(i,2));
    d_torque(i,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(ex) * d_mu(i,2) - static_cast<KK_FLOAT>(ez) * d_mu(i,0));
    d_torque(i,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(ey) * d_mu(i,0) - static_cast<KK_FLOAT>(ex) * d_mu(i,1));
    result[0] -= static_cast<double>(d_mu(i,0)) * ex + static_cast<double>(d_mu(i,1)) * ey + static_cast<double>(d_mu(i,2)) * ez;
  }
}

template<class DeviceType>
template<int QFLAG, int MUFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixEfieldKokkos<DeviceType>::operator()(TagFixEfieldNonConstant<QFLAG,MUFLAG>, const int &i, value_type result) const {
  if ( QFLAG && (d_mask(i) & groupbit)) {
    if (region && !d_match[i]) return;

    KK_FLOAT fx, fy, fz;

    if (xstyle == ATOM) fx = static_cast<KK_FLOAT>(qe2f) * d_q(i) * d_efield(i,0);
    else fx = d_q(i) * static_cast<KK_FLOAT>(ex);
    if (ystyle == ATOM) fy = static_cast<KK_FLOAT>(qe2f) * d_q(i) * d_efield(i,1);
    else fy = d_q(i) * static_cast<KK_FLOAT>(ey);
    if (zstyle == ATOM) fz = static_cast<KK_FLOAT>(qe2f) * d_q(i) * d_efield(i,2);
    else fz = d_q(i) * static_cast<KK_FLOAT>(ez);

    d_f(i,0) += static_cast<KK_ACC_FLOAT>(fx);
    d_f(i,1) += static_cast<KK_ACC_FLOAT>(fy);
    d_f(i,2) += static_cast<KK_ACC_FLOAT>(fz);
    result[1] += static_cast<double>(fx);
    result[2] += static_cast<double>(fy);
    result[3] += static_cast<double>(fz);

    if (pstyle == ATOM) result[0] += qe2f * static_cast<double>(d_q(i)) * static_cast<double>(d_efield(i,3));
    else if (estyle == ATOM) result[0] += static_cast<double>(d_efield(i,3));
  }

  if (MUFLAG && (d_mask(i) & groupbit)) {
    if (region && !d_match[i]) return;
    d_torque(i,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(ez) * d_mu(i,1) - static_cast<KK_FLOAT>(ey) * d_mu(i,2));
    d_torque(i,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(ex) * d_mu(i,2) - static_cast<KK_FLOAT>(ez) * d_mu(i,0));
    d_torque(i,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(ey) * d_mu(i,0) - static_cast<KK_FLOAT>(ex) * d_mu(i,1));

  }
}

/* ----------------------------------------------------------------------
   tally virial into global and per-atom accumulators
   i = local index of atom
   v = total virial for the interaction
   increment global virial by v
   increment per-atom virial by v
   this method can be used when fix computes forces in post_force()
   and the force depends on a distance to some external object
     e.g. fix wall/lj93: compute virial only on owned atoms
------------------------------------------------------------------------- */

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixEfieldKokkos<DeviceType>::v_tally(value_type result, int i, KK_FLOAT *v) const
{
  if (vflag_global) {
    result[4] += static_cast<double>(v[0]);
    result[5] += static_cast<double>(v[1]);
    result[6] += static_cast<double>(v[2]);
    result[7] += static_cast<double>(v[3]);
    result[8] += static_cast<double>(v[4]);
    result[9] += static_cast<double>(v[5]);
  }

  if (vflag_atom) {
    Kokkos::atomic_add(&(d_vatom(i,0)),static_cast<KK_ACC_FLOAT>(v[0]));
    Kokkos::atomic_add(&(d_vatom(i,1)),static_cast<KK_ACC_FLOAT>(v[1]));
    Kokkos::atomic_add(&(d_vatom(i,2)),static_cast<KK_ACC_FLOAT>(v[2]));
    Kokkos::atomic_add(&(d_vatom(i,3)),static_cast<KK_ACC_FLOAT>(v[3]));
    Kokkos::atomic_add(&(d_vatom(i,4)),static_cast<KK_ACC_FLOAT>(v[4]));
    Kokkos::atomic_add(&(d_vatom(i,5)),static_cast<KK_ACC_FLOAT>(v[5]));
  }
}

namespace LAMMPS_NS {
template class FixEfieldKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixEfieldKokkos<LMPHostType>;
#endif
}
