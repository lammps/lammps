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
   Contributing author: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "fix_rigid_nh_small_kokkos.h"

#include "atom_masks.h"
#include "error.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType, bool TSTAT, bool PSTAT>
FixRigidNHSmallKokkos<DeviceType,TSTAT,PSTAT>::
  FixRigidNHSmallKokkos(LAMMPS *lmp, int narg, char **arg) :
    FixRigidBaseKokkos<DeviceType, FixRigidNHSmall>(lmp, narg, arg)
{

  if constexpr (TSTAT || PSTAT) {
    this->restart_global = 1;
    this->extscalar = 1;
  }

  if constexpr (TSTAT) {
    if (tstat_flag == 0)
      error->all(FLERR,"Did not set temp for fix {}", style);
    if (t_start <= 0.0 || t_stop <= 0.0)
      error->all(FLERR,"Target temperature for fix {} cannot be 0.0", style);

  if (this->t_period <= 0.0)
    this->error->all(FLERR,"Fix {} period must be > 0.0", style);
  this->t_freq = 1.0 / this->t_period;

  if (this->t_chain < 1)
    this->error->all(FLERR,"Fix {} t_chain should not be less than 1", style);
  if (this->t_iter < 1)
    this->error->all(FLERR,"Fix {} t_iter should not be less than 1", style);
  if (this->t_order != 3 && this->t_order != 5)
    this->error->all(FLERR,"Fix {} t_order must be 3 or 5", style);

  }

  if constexpr (PSTAT) {
    if (pstat_flag == 0)
      error->all(FLERR,"Did not set press for fix {}", style);
    if (p_start[0] < 0.0 || p_start[1] < 0.0 || p_start[2] < 0.0 ||
        p_stop[0] < 0.0  || p_stop[1] < 0.0  || p_stop[2] < 0.0)
      error->all(FLERR,"Target pressure for fix {} cannot be 0.0", style);
    p_freq[0] = this->p_freq[1] = this->p_freq[2] = 0.0;
    if (this->p_flag[0]) this->p_freq[0] = 1.0 / this->p_period[0];
    if (this->p_flag[1]) this->p_freq[1] = 1.0 / this->p_period[1];
    if (this->p_flag[2]) this->p_freq[2] = 1.0 / this->p_period[2];
    this->id_temp = utils::strdup(std::string(this->id)+"_temp");
    this->modify->add_compute(fmt::format("{} all temp",this->id_temp));
  this->tcomputeflag = 1;

  this->id_press = utils::strdup(std::string(this->id)+"_press");
  this->modify->add_compute(fmt::format("{} all pressure {}",this->id_press,this->id_temp));
  this->pcomputeflag = 1;

  }

}

/* ----------------------------------------------------------------------
   FixRigidSmall PROTECTED METHODS
------------------------------------------------------------------------- */

template<class DeviceType, bool TSTAT, bool PSTAT>
void FixRigidNHSmallKokkos<DeviceType,TSTAT,PSTAT>::copy_arrays(int i, int j, int delflag)
{
  this->sync_host_base();
  FixRigidSmall::copy_arrays(i, j, delflag);
  this->modify_host_base();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, bool TSTAT, bool PSTAT>
void FixRigidNHSmallKokkos<DeviceType,TSTAT,PSTAT>::set_arrays(int i)
{
  this->sync_host_base();
  FixRigidSmall::set_arrays(i);
  this->modify_host_base();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, bool TSTAT, bool PSTAT>
void FixRigidNHSmallKokkos<DeviceType,TSTAT,PSTAT>::set_molecule(int nlocalprev, tagint tagprev,
                                                    int imol, double *xgeom,
                                                    double *vcm, double *quat)
{
  this->sync_host_base();
  FixRigidSmall::set_molecule(nlocalprev, tagprev, imol, xgeom, vcm, quat);
  this->modify_host_base();
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixRigidNHSmallKokkos<LMPDeviceType,false,false>; // NVE
template class FixRigidNHSmallKokkos<LMPDeviceType,true,false>;  // NVT
template class FixRigidNHSmallKokkos<LMPDeviceType,false,true>;  // NPH
template class FixRigidNHSmallKokkos<LMPDeviceType,true,true>;   // NPT
#ifdef LMP_KOKKOS_GPU
//template class FixRigidSmallKokkos<LMPHostType>;
#endif
}
