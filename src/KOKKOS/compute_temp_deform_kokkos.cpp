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
   Contributing authors: Emily Kahl (Uni. of QLD, e.kahl@uq.edu.au)
------------------------------------------------------------------------- */

#include "compute_temp_deform_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain_kokkos.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "math_extra.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeTempDeformKokkos<DeviceType>::ComputeTempDeformKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeTempDeform(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  domainKK = (DomainKokkos *) domain;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = EMPTY_MASK;

  maxbias = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::post_constructor() {
  if (tcomputeflag) {
    this->id_temp = utils::strdup(std::string(id) + "_temp");
    this->modify->add_compute(fmt::format("{} {} temp/kk", id_temp, this->group->names[igroup]));
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double ComputeTempDeformKokkos<DeviceType>::compute_scalar()
{
  atomKK->sync(execution_space,datamask_read);

  invoked_scalar = update->ntimestep;

  remove_deform_bias_all_kk();
  if (this->temperature->kokkosable)
    scalar = this->temperature->compute_scalar();
  else {
    atomKK->sync(this->temperature->execution_space,this->temperature->datamask_read);
    scalar = this->temperature->compute_scalar();
    atomKK->modified(this->temperature->execution_space,this->temperature->datamask_modify);
    atomKK->sync(this->execution_space,this->temperature->datamask_modify);
  }

  if (dynamic) dof = this->temperature->dof;
  restore_deform_bias_all_kk();

  return scalar;
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::compute_vector()
{
  atomKK->sync(execution_space,datamask_read);

  invoked_vector = update->ntimestep;

  remove_deform_bias_all_kk();
  if (this->temperature->kokkosable)
    this->temperature->compute_vector();
  else {
    atomKK->sync(this->temperature->execution_space,this->temperature->datamask_read);
    this->temperature->compute_vector();
    atomKK->modified(this->temperature->execution_space,this->temperature->datamask_modify);
    atomKK->sync(this->execution_space,this->temperature->datamask_modify);
  }

  vector = this->temperature->vector;
  if (dynamic) dof = this->temperature->dof;
  restore_deform_bias_all_kk();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::remove_bias_all()
{
  remove_bias_all_kk();
  atomKK->sync(Host,V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::remove_bias_all_kk()
{
  remove_deform_bias_all_kk();
  if (which == BIAS) {
    if (temperature->kokkosable) temperature->remove_bias_all_kk();
    else {
      atomKK->sync(this->temperature->execution_space,this->temperature->datamask_read);
      this->temperature->remove_bias_all();
      atomKK->modified(this->temperature->execution_space,this->temperature->datamask_modify);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::restore_bias_all()
{
  restore_bias_all_kk();
  atomKK->sync(Host,V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::restore_bias_all_kk()
{
  if (which == BIAS) {
    if (temperature->kokkosable) temperature->restore_bias_all();
    else {
      atomKK->sync(this->temperature->execution_space,this->temperature->datamask_read);
      this->temperature->restore_bias_all();
      atomKK->modified(this->temperature->execution_space,this->temperature->datamask_modify);
    }
  }
  restore_deform_bias_all_kk();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::remove_deform_bias_all()
{
  remove_deform_bias_all_kk();
  atomKK->sync(Host,V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::remove_deform_bias_all_kk()
{
  atomKK->sync(execution_space,X_MASK|V_MASK|MASK_MASK);
  v = atomKK->k_v.view<DeviceType>();
  x = atomKK->k_x.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  if (atom->nmax > maxbias) {
    maxbias = atom->nmax;
    vbiasall = typename AT::t_kkfloat_1d_3("temp/deform/kk:vbiasall", maxbias);
  }

  domainKK->x2lamda(nlocal);

  h_rate = domain->h_rate;
  h_ratelo = domain->h_ratelo;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempDeformRemoveBias >(0,nlocal),*this);
  copymode = 0;

  domainKK->lamda2x(nlocal);

  atomKK->modified(execution_space,V_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempDeformKokkos<DeviceType>::operator()(TagComputeTempDeformRemoveBias, const int &i) const {
  if (mask[i] & groupbit) {
    vbiasall(i,0) = h_rate[0]*x(i,0) + h_rate[5]*x(i,1) + h_rate[4]*x(i,2) + h_ratelo[0];
    vbiasall(i,1) = h_rate[1]*x(i,1) + h_rate[3]*x(i,2) + h_ratelo[1];
    vbiasall(i,2) = h_rate[2]*x(i,2) + h_ratelo[2];
    v(i,0) -= vbiasall(i,0);
    v(i,1) -= vbiasall(i,1);
    v(i,2) -= vbiasall(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::restore_deform_bias_all()
{
  restore_deform_bias_all_kk();
  atomKK->sync(Host,V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::restore_deform_bias_all_kk()
{
  atomKK->sync(execution_space,V_MASK|MASK_MASK);
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempDeformRestoreBias >(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempDeformKokkos<DeviceType>::operator()(TagComputeTempDeformRestoreBias, const int &i) const {
  if (mask[i] & groupbit) {
    v(i,0) += vbiasall(i,0);
    v(i,1) += vbiasall(i,1);
    v(i,2) += vbiasall(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::apply_deform_bias_all(double dtv)
{
  apply_deform_bias_all_kk(dtv);
  atomKK->sync(Host,V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::apply_deform_bias_all_kk(double dtv)
{
  atomKK->sync(execution_space,X_MASK | V_MASK | MASK_MASK);
  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  double grad_u[6];
  MathExtra::multiply_shape_shape(domain->h_rate,domain->h_inv,grad_u);
  d_grad_u = Few<double, 6>(grad_u);

  double xref[5];
  xref[0] = (domain->boxhi[0] + domain->boxlo[0])/2.;
  xref[1] = (domain->boxhi[1] + domain->boxlo[1])/2.;
  xref[2] = (domain->boxhi[2] + domain->boxlo[2])/2.;
  // if needed, integrate boxlo to account for box not being updated yet
  // xmid does not change
  xref[3] = xref[1] + (domain->boxlo[1] - xref[1])*exp(grad_u[1]*dtv);
  xref[4] = xref[2] + (domain->boxlo[2] - xref[2])*exp(grad_u[2]*dtv);
  d_xref = Few<double, 5>(xref);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempDeformApplyBias >(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempDeformKokkos<DeviceType>::operator()(TagComputeTempDeformApplyBias, const int &i) const {
  if (mask[i] & groupbit) {
    v(i,0) += (x(i,0) - d_xref[0]) * d_grad_u[0] + (x(i,1) - d_xref[3]) * d_grad_u[5] + (x(i,2) - d_xref[4]) * d_grad_u[4];
    v(i,1) += (x(i,1) - d_xref[1]) * d_grad_u[1] + (x(i,2) - d_xref[4]) * d_grad_u[3];
    v(i,2) += (x(i,2) - d_xref[2]) * d_grad_u[2];
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeTempDeformKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeTempDeformKokkos<LMPHostType>;
#endif
}
