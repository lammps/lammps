// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

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

#include "atom.h"
#include "atom.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_deform_kokkos.h"
#include "group.h"
#include "kokkos_few.h"
#include "math_extra.h"
#include "memory_kokkos.h"
#include "modify.h"

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
    this->error->all(FLERR,"Temperature control must be used with fix nvt/kk");
  if (this->pstat_flag)
    this->error->all(FLERR,"Pressure control can not be used with fix nvt/kk");

  if (this->mtchain_default_flag) this->mtchain = 1;

  this->id_temp = utils::strdup(std::string(this->id)+"_temp");
  this->modify->add_compute(fmt::format("{} all temp/deform/kk",this->id_temp));
  this->tcomputeflag = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVTSllodKokkos<DeviceType>::init()
{
  FixNHKokkos<DeviceType>::init();

  if (!this->temperature->tempbias)
    this->error->all(FLERR,"Temperature for fix nvt/sllod does not have a bias");

  nondeformbias = 0;
  if (utils::strmatch(this->temperature->style,"^temp/deform")) nondeformbias = 1;

  // check fix deform remap settings

  int i;
  for (i = 0; i < this->modify->nfix; i++)
    if (utils::strmatch(this->modify->fix[i]->style,"^deform")) {
      if (((FixDeform *) this->modify->fix[i])->remapflag != Domain::V_REMAP)
        this->error->all(FLERR,"Using fix nvt/sllod with inconsistent fix deform remap option");
      break;
    }
  if (i == this->modify->nfix)
    this->error->all(FLERR,"Using fix nvt/sllod with no fix deform defined");
}

/* ----------------------------------------------------------------------
   perform half-step scaling of velocities
-----------------------------------------------------------------------*/

template<class DeviceType>
void FixNVTSllodKokkos<DeviceType>::nh_v_temp()
{
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
  }
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (this->igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  double h_two[6];
  MathExtra::multiply_shape_shape(this->domain->h_rate,this->domain->h_inv,h_two);

  d_h_two = Few<double, 6>(h_two);

  if (vdelu.extent(0) < atomKK->nmax)
    vdelu = typename AT::t_v_array(Kokkos::NoInit("nvt/sllod/kk:vdelu"), atomKK->nmax);

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNVTSllod_temp1>(0,nlocal),*this);
  this->copymode = 0;

  this->temperature->remove_bias_all();

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNVTSllod_temp2>(0,nlocal),*this);
  this->copymode = 0;

  this->temperature->restore_bias_all();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVTSllodKokkos<DeviceType>::operator()(TagFixNVTSllod_temp1, const int &i) const {
  if (mask[i] & this->groupbit) {
    vdelu(i,0) = d_h_two[0]*v(i,0) + d_h_two[5]*v(i,1) + d_h_two[4]*v(i,2);
    vdelu(i,1) = d_h_two[1]*v(i,1) + d_h_two[3]*v(i,2);
    vdelu(i,2) = d_h_two[2]*v(i,2);
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVTSllodKokkos<DeviceType>::operator()(TagFixNVTSllod_temp2, const int &i) const {
  if (mask[i] & this->groupbit) {
    v(i,0) = v(i,0)*this->factor_eta - this->dthalf*vdelu(i,0);
    v(i,1) = v(i,1)*this->factor_eta - this->dthalf*vdelu(i,1);
    v(i,2) = v(i,2)*this->factor_eta - this->dthalf*vdelu(i,2);
  }
}

namespace LAMMPS_NS {
template class FixNVTSllodKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNVTSllodKokkos<LMPHostType>;
#endif
}


