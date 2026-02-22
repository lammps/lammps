// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/
   LAMMPS development team: developers@lammps.org, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Emily Kahl (Uni. of QLD, e.kahl@uq.edu.au)
------------------------------------------------------------------------- */

#include "fix_nvt_sllod_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "compute.h"
#include "compute_temp_deform.h"
#include "domain.h"
#include "error.h"
#include "fix_deform.h"
#include "group.h"
#include "kokkos_few.h"
#include "math_extra.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVTSllodKokkos<DeviceType>::FixNVTSllodKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNHKokkos<DeviceType>(lmp, narg, arg)
{
  atomKK = (AtomKokkos *) this->atom;
  this->kokkosable = 1;
  this->domainKK = (DomainKokkos *) this->domain;

  if (!this->tstat_flag)
    this->error->all(FLERR,"Temperature control must be used with fix nvt/sllod/kk");
  if (this->pstat_flag)
    this->error->all(FLERR,"Pressure control can not be used with fix nvt/sllod/kk");

  this->psllod_flag = 0;
  this->peculiar_flag = 0;
  this->integrator = REVERSIBLE;
  bool user_kick = false;
  if (this->mtchain_default_flag) this->mtchain = 1;

  // select SLLOD/p-SLLOD/g-SLLOD variant and velocity frame

  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"psllod") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix nvt/sllod/kk psllod", this->error);
      this->psllod_flag = utils::logical(FLERR,arg[iarg+1],false,this->lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"peculiar") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} peculiar", this->style), this->error);
      this->peculiar_flag = utils::logical(FLERR,arg[iarg+1],false,this->lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"kick") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} kick", this->style), this->error);
      this->kick_flag = utils::logical(FLERR,arg[iarg+1],false,this->lmp);
      user_kick = true;
      iarg += 2;
    } else if (strcmp(arg[iarg],"integrator") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix nvt/sllod integrator", this->error);
      if (strcmp(arg[iarg+1],"reversible") == 0)  this->integrator = REVERSIBLE;
      else if (strcmp(arg[iarg+1],"legacy") == 0) this->integrator = LEGACY;
      else this->error->all(FLERR, "Unknown fix {} integrator argument: {}", this->style, arg[iarg+1]);
      iarg += 2;
    } else iarg++;
  }

  // default to applying velocity kick in lab frame
  if (!user_kick) this->kick_flag = !this->peculiar_flag;

  if (this->integrator == LEGACY && this->peculiar_flag == 1)
    this->error->all(FLERR, "fix {} legacy integrator is not compatible with peculiar=yes", this->style);

  this->id_temp = utils::strdup(std::string(this->id)+"_temp");
  if (peculiar_flag) this->modify->add_compute(fmt::format("{} {} temp/kk",this->id_temp,this->group->names[this->igroup]));
  else this->modify->add_compute(fmt::format("{} {} temp/deform/kk",this->id_temp,this->group->names[this->igroup]));
  this->tcomputeflag = 1;
  this->nondeformbias = 0;

  this->execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  this->datamask_read =  EMPTY_MASK;
  this->datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVTSllodKokkos<DeviceType>::init()
{
  FixNHKokkos<DeviceType>::init();

  if (!this->peculiar_flag && !this->temperature->tempbias)
    this->error->all(FLERR,"Temperature for fix {} does not have a bias", this->style);

  this->nondeformbias = 0;
  if (!utils::strmatch(this->temperature->style,"^temp/deform")) {
    if (this->integrator == LEGACY) {
      this->nondeformbias = 1;
      if (this->kick_flag)
        this->error->all(FLERR, "fix {} with peculiar=no and kick=yes requires temperature bias "
                   "to be calculated by compute temp/deform", this->style);
    } else if (!this->peculiar_flag) {
      this->error->all(FLERR,"Fix {} used with lab-frame velocity and non-deform "
                     "temperature bias. For non-deform biases, either set peculiar = yes "
                     "or pass an explicit temp/deform with an extra bias", this->style);
    }
  }

  // check fix deform remap settings

  auto deform = this->modify->get_fix_by_style("^deform");
  if (deform.size() < 1)
    this->error->all(FLERR,"Using fix {} with no fix deform defined", this->style);

  for (auto &ifix : deform) {
    auto f = dynamic_cast<FixDeform *>(ifix);
    if ((peculiar_flag && f->remapflag != Domain::NO_REMAP) ||
        (!peculiar_flag && f->remapflag != Domain::V_REMAP))
      this->error->all(FLERR,"Using fix {} with inconsistent fix {} remap option", this->style, f->style);

    if (kick_flag) {
      // apply initial kick if velocity stored in lab frame
      // only kick once by default for correct dynamics with multiple run commands
      // make sure fix deform init happens first so h_rate is set
      if (!peculiar_flag) {
        f->init();
        if (this->comm->me == 0)
          utils::logmesg(this->lmp, "fix {} applying velocity profile kick.\n", this->style);
        dynamic_cast<ComputeTempDeform*>(this->temperature)->apply_deform_bias_all();
        kick_flag = 0;
      } else if (this->comm->me == 0) {
        this->error->warning(FLERR,"fix {} using peculiar frame velocity. "
                       "Ignoring kick flag.", this->style);
      }
    }

    break;
  }
}

/* ----------------------------------------------------------------------
   perform full-step update of positions with streaming velocity
   also perform sllod update reversibly
-----------------------------------------------------------------------*/

template<class DeviceType>
void FixNVTSllodKokkos<DeviceType>::nve_x()
{
  if (this->integrator == LEGACY) return FixNHKokkos<DeviceType>::nve_x();

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (this->igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  // x update by full step only for atoms in group
  // identical for SLLOD and p-SLLOD
  // velocity treated in peculiar frame relative to sllod streaming for
  //  reversibility, so need to manually account for change in streaming
  //  velocity

  double h_two[6];
  MathExtra::multiply_shape_shape(this->domain->h_rate,this->domain->h_inv,h_two);
  d_h_two = Few<double, 6>(h_two);

  double xfac[3];
  xfac[0] = exp(h_two[0]*this->dthalf);
  xfac[1] = exp(h_two[1]*this->dthalf);
  xfac[2] = exp(h_two[2]*this->dthalf);
  d_xfac = Few<double, 3>(xfac);
  double vfac[3];
  vfac[0] = exp(-h_two[0]*this->dthalf);
  vfac[1] = exp(-h_two[1]*this->dthalf);
  vfac[2] = exp(-h_two[2]*this->dthalf);
  d_vfac = Few<double, 3>(vfac);

  // fix deform uses the box center as origin for elongation,
  // and the lower corner for shear, so adjust for that
  // to avoid an apparent drift relative to the box and prevent
  // extra atom exchanges between MPI ranks
  double xmid[3];
  for (int i = 0; i < 3; ++i) {
    xmid[i] = (this->domain->boxhi[i] + this->domain->boxlo[i])/2.;
  }
  d_xmid = Few<double, 3>(xmid);

  double xlo[5];
  xlo[0] = this->domain->boxlo[0];
  xlo[1] = this->domain->boxlo[1];
  xlo[2] = this->domain->boxlo[2];
  // propagate boxlo to make second half step reversible
  // xmid does not change
  xlo[3] = xmid[1] + (xlo[1] - xmid[1])*xfac[1]*xfac[1];
  xlo[4] = xmid[2] + (xlo[2] - xmid[2])*xfac[2]*xfac[2];
  d_xlo = Few<double, 5>(xlo);

  atomKK->sync(this->execution_space, X_MASK | V_MASK | MASK_MASK);

  if (!peculiar_flag) {
    if (this->temperature->kokkosable)
      dynamic_cast<ComputeTempDeform*>(this->temperature)->remove_deform_bias_all_kk();
    else {
      atomKK->sync(this->temperature->execution_space,this->temperature->datamask_read);
      dynamic_cast<ComputeTempDeform*>(this->temperature)->remove_deform_bias_all();
      atomKK->modified(this->temperature->execution_space,this->temperature->datamask_modify);
      atomKK->sync(this->execution_space,this->temperature->datamask_modify);
    }
  }

  // parallel for
  this->copymode = 1;
  if (psllod_flag) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNVTSllod_nvex<true>>(0,nlocal),*this);
  else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNVTSllod_nvex<false>>(0,nlocal),*this);
  this->copymode = 0;

  atomKK->modified(this->execution_space, X_MASK | V_MASK);

  // x has changed, so can't just call restore_deform_bias_all
  // pass in dtv to account for update to box shape
  if (!peculiar_flag) {
    if (this->temperature->kokkosable)
      dynamic_cast<ComputeTempDeform*>(this->temperature)->apply_deform_bias_all_kk(this->dtv);
    else {
      atomKK->sync(this->temperature->execution_space,this->temperature->datamask_read);
      dynamic_cast<ComputeTempDeform*>(this->temperature)->apply_deform_bias_all(this->dtv);
      atomKK->modified(this->temperature->execution_space,this->temperature->datamask_modify);
      atomKK->sync(this->execution_space,this->temperature->datamask_modify);
    }
  }
}

template<class DeviceType>
template<bool PSLLOD>
KOKKOS_INLINE_FUNCTION
void FixNVTSllodKokkos<DeviceType>::operator()(TagFixNVTSllod_nvex<PSLLOD>, const int& i) const {
  if (mask[i] & this->groupbit) {
    // first half sllod update
    if (PSLLOD) {
      v(i,2) -= this->dthalf*this->d_h_two[2]*this->d_h_two[2]*x(i,2);
      v(i,1) -= this->dthalf*this->d_h_two[3]*v(i,2) + this->dthalf*this->d_h_two[1]*this->d_h_two[1]*x(i,1);
      v(i,0) -= this->dthalf*(this->d_h_two[5]*v(i,1) + this->d_h_two[4]*v(i,2))
                 + this->dthalf*this->d_h_two[0]*this->d_h_two[0]*x(i,0);
    } else {
      v(i,1) -= this->dthalf*this->d_h_two[3]*v(i,2);
      v(i,0) -= this->dthalf*(this->d_h_two[5]*v(i,1) + this->d_h_two[4]*v(i,2));
    }
    v(i,0) *= this->d_vfac[0];
    v(i,1) *= this->d_vfac[1];
    v(i,2) *= this->d_vfac[2];

    x(i,1) += this->dthalf * this->d_h_two[3]*(x(i,2) - this->d_xlo[2]);
    x(i,0) += this->dthalf * (this->d_h_two[5]*(x(i,1) - this->d_xlo[1]) + this->d_h_two[4]*(x(i,2) - this->d_xlo[2]));
    x(i,0) = this->d_xmid[0] + (x(i,0) - this->d_xmid[0])*this->d_xfac[0];
    x(i,1) = this->d_xmid[1] + (x(i,1) - this->d_xmid[1])*this->d_xfac[1];
    x(i,2) = this->d_xmid[2] + (x(i,2) - this->d_xmid[2])*this->d_xfac[2];

    // nve position update
    x(i,0) += this->dtv * v(i,0);
    x(i,1) += this->dtv * v(i,1);
    x(i,2) += this->dtv * v(i,2);

    // 2nd half sllod update
    x(i,0) = this->d_xmid[0] + (x(i,0) - this->d_xmid[0])*this->d_xfac[0];
    x(i,1) = this->d_xmid[1] + (x(i,1) - this->d_xmid[1])*this->d_xfac[1];
    x(i,2) = this->d_xmid[2] + (x(i,2) - this->d_xmid[2])*this->d_xfac[2];
    // d_xlo[3] is propagated xlo[1], d_xlo[4] is propagated xlo[2]
    x(i,0) += this->dthalf * (this->d_h_two[5]*(x(i,1) - this->d_xlo[3]) + this->d_h_two[4]*(x(i,2) - this->d_xlo[4]));
    x(i,1) += this->dthalf * this->d_h_two[3]*(x(i,2) - this->d_xlo[4]);

    // second half sllod velocity step
    // apply here so streaming component matches x when storing in lab frame
    v(i,0) *= this->d_vfac[0];
    v(i,1) *= this->d_vfac[1];
    v(i,2) *= this->d_vfac[2];
    if (PSLLOD) {
      v(i,0) -= this->dthalf*(this->d_h_two[5]*v(i,1) + this->d_h_two[4]*v(i,2))
                 + this->dthalf*this->d_h_two[0]*this->d_h_two[0]*x(i,0);
      v(i,1) -= this->dthalf*this->d_h_two[3]*v(i,2) + this->dthalf*this->d_h_two[1]*this->d_h_two[1]*x(i,1);
      v(i,2) -= this->dthalf*this->d_h_two[2]*this->d_h_two[2]*x(i,2);
    } else {
      v(i,0) -= this->dthalf*(this->d_h_two[5]*v(i,1) + this->d_h_two[4]*v(i,2));
      v(i,1) -= this->dthalf*this->d_h_two[3]*v(i,2);
    }
  }
}

/* ----------------------------------------------------------------------
   perform half-step scaling of velocities using legacy method
   NOTE: this style of integration is not time-reversible under mixed
         flows, and neglects the change in streaming velocity caused by
         the position update.
-----------------------------------------------------------------------*/

template<class DeviceType>
void FixNVTSllodKokkos<DeviceType>::nh_v_temp()
{
  if (this->integrator == REVERSIBLE) return FixNHKokkos<DeviceType>::nh_v_temp();

  // remove and restore bias = streaming velocity = Hrate*lamda + Hratelo
  // thermostat thermal velocity only
  // vdelu = SLLOD correction = Hrate*Hinv*vthermal
  // for non temp/deform BIAS:
  //   calculate temperature since some computes require temp
  //   computed on current nlocal atoms to remove bias

  if (nondeformbias) {
    atomKK->sync(this->temperature->execution_space,this->temperature->datamask_read);
    this->temperature->compute_scalar();
    atomKK->modified(this->temperature->execution_space,this->temperature->datamask_modify);
    atomKK->sync(this->execution_space,this->temperature->datamask_modify);
  }

  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (this->igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  double h_two[6];
  MathExtra::multiply_shape_shape(this->domain->h_rate,this->domain->h_inv,h_two);

  d_h_two = Few<double, 6>(h_two);

  if ((int)vdelu.extent(0) < atomKK->nmax)
    vdelu = typename AT::t_kkfloat_1d_3(Kokkos::NoInit("nvt/sllod/kk:vdelu"), atomKK->nmax);

  if (!this->psllod_flag) {
    if (this->temperature->kokkosable) this->temperature->remove_bias_all_kk();
    else {
      atomKK->sync(this->temperature->execution_space,this->temperature->datamask_read);
      this->temperature->remove_bias_all();
      atomKK->modified(this->temperature->execution_space,this->temperature->datamask_modify);
      atomKK->sync(this->execution_space,this->temperature->datamask_modify);
    }
  }

  atomKK->sync(this->execution_space,V_MASK | MASK_MASK);

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNVTSllod_temp1>(0,nlocal),*this);
  this->copymode = 0;

  if (this->psllod_flag) {
    if (this->temperature->kokkosable) this->temperature->remove_bias_all_kk();
    else {
      atomKK->sync(this->temperature->execution_space,this->temperature->datamask_read);
      this->temperature->remove_bias_all();
      atomKK->modified(this->temperature->execution_space,this->temperature->datamask_modify);
      atomKK->sync(this->execution_space,this->temperature->datamask_modify);
    }
  }

  atomKK->sync(this->execution_space,V_MASK | MASK_MASK);

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNVTSllod_temp2>(0,nlocal),*this);
  this->copymode = 0;

  atomKK->modified(this->execution_space,V_MASK);

  if (this->temperature->kokkosable) this->temperature->restore_bias_all_kk();
  else {
    atomKK->sync(this->temperature->execution_space,this->temperature->datamask_read);
    this->temperature->restore_bias_all();
    atomKK->modified(this->temperature->execution_space,this->temperature->datamask_modify);
    atomKK->sync(this->execution_space,this->temperature->datamask_modify);
  }
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixNVTSllodKokkos<DeviceType>::operator()(TagFixNVTSllod_temp1, const int &i) const {
  if (mask[i] & this->groupbit) {
    vdelu(i,0) = d_h_two[0]*v(i,0) + d_h_two[5]*v(i,1) + d_h_two[4]*v(i,2);
    vdelu(i,1) = d_h_two[1]*v(i,1) + d_h_two[3]*v(i,2);
    vdelu(i,2) = d_h_two[2]*v(i,2);
  }
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixNVTSllodKokkos<DeviceType>::operator()(TagFixNVTSllod_temp2, const int &i) const {
  if (mask[i] & this->groupbit) {
    v(i,0) = v(i,0)*this->factor_eta - this->dthalf*vdelu(i,0);
    v(i,1) = v(i,1)*this->factor_eta - this->dthalf*vdelu(i,1);
    v(i,2) = v(i,2)*this->factor_eta - this->dthalf*vdelu(i,2);
  }
}

/* ----------------------------------------------------------------------
    calculate the number of data to be packed
------------------------------------------------------------------------- */

template<class DeviceType>
int FixNVTSllodKokkos<DeviceType>::size_restart_global()
{
  return FixNH::size_restart_global() + 1;
}

/* ----------------------------------------------------------------------
   pack restart data
------------------------------------------------------------------------- */

template<class DeviceType>
int FixNVTSllodKokkos<DeviceType>::pack_restart_data(double *list)
{
  list[0] = this->kick_flag;
  int n = 1 + FixNHKokkos<DeviceType>::pack_restart_data(&list[1]);

  return n;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNVTSllodKokkos<DeviceType>::restart(char *buf)
{
  auto *list = (double *) buf;
  this->kick_flag = static_cast<int>(list[0]);
  FixNHKokkos<DeviceType>::restart(buf + sizeof(double));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixNVTSllodKokkos<DeviceType>::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"kick") == 0) {
    if (narg < 2) this->error->all(FLERR,"Illegal fix_modify command");
    kick_flag = utils::logical(FLERR,arg[1],false,this->lmp);
    return 2;
  }
  return FixNHKokkos<DeviceType>::modify_param(narg, arg);
}
namespace LAMMPS_NS {
template class FixNVTSllodKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNVTSllodKokkos<LMPHostType>;
#endif
}

