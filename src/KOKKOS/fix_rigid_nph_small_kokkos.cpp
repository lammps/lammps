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
   Contributing author: Mitch Murphy (alphataubio@gmail.com)
   references: Kamberaj et al., J. Chem. Phys. 122, 224114 (2005)
               Miller et al., J Chem Phys. 116, 8649-8659 (2002)
------------------------------------------------------------------------- */

#include "fix_rigid_nph_small_kokkos.h"

#include "error.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidNPHSmallKokkos<DeviceType>::FixRigidNPHSmallKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigidNHSmallKokkos<DeviceType>(lmp, narg, arg)
{
  this->scalar_flag = 1;
  this->restart_global = 1;
  this->extscalar = 1;

  if (this->pstat_flag == 0)
    this->error->all(FLERR,"Pressure control must be used with fix nph/small/kk");
  if (this->tstat_flag == 1)
    this->error->all(FLERR,"Temperature control must not be used with fix nph/small/kk");
  if (this->p_start[0] < 0.0 || this->p_start[1] < 0.0 || this->p_start[2] < 0.0 ||
      this->p_stop[0] < 0.0 || this->p_stop[1] < 0.0 || this->p_stop[2] < 0.0)
    this->error->all(FLERR,"Target pressure for fix rigid/nph/small/kk cannot be < 0.0");

  this->p_freq[0] = this->p_freq[1] = this->p_freq[2] = 0.0;
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

namespace LAMMPS_NS {
template class FixRigidNPHSmallKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRigidNPHSmallKokkos<LMPHostType>;
#endif
}
