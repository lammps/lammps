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
   Contributing author: Tony Sheh (U Michigan), Trung Dac Nguyen (U Michigan)
   references: Kamberaj et al., J. Chem. Phys. 122, 224114 (2005)
               Miller et al., J Chem Phys. 116, 8649-8659 (2002)
------------------------------------------------------------------------- */

#include "fix_rigid_npt.h"

#include "error.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixRigidNPT::FixRigidNPT(LAMMPS *lmp, int narg, char **arg) :
  FixRigidNH(lmp, narg, arg)
{
  // other setting are made by parent

  scalar_flag = 1;
  restart_global = 1;
  extscalar = 1;

  // error checks

  if (tstat_flag == 0 || pstat_flag == 0)
    error->all(FLERR,"Did not set temperature or pressure for fix rigid/npt");
  if (t_start <= 0.0 || t_stop <= 0.0)
    error->all(FLERR,"Target temperature for fix rigid/npt cannot be 0.0");
  if (t_period <= 0.0) error->all(FLERR,"Fix rigid/npt period must be > 0.0");

  // thermostat chain parameters

  if (t_chain < 1) error->all(FLERR,"Illegal fix rigid/npt command");
  if (t_iter < 1) error->all(FLERR,"Illegal fix rigid/npt command");
  if (t_order != 3 && t_order != 5)
    error->all(FLERR,"Fix rigid/npt temperature order must be 3 or 5");

  // convert input periods to frequency

  t_freq = 0.0;
  p_freq[0] = p_freq[1] = p_freq[2] = 0.0;

  t_freq = 1.0 / t_period;
  if (p_flag[0]) p_freq[0] = 1.0 / p_period[0];
  if (p_flag[1]) p_freq[1] = 1.0 / p_period[1];
  if (p_flag[2]) p_freq[2] = 1.0 / p_period[2];

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  //   and thus its KE/temperature contribution should use group all

  id_temp = utils::strdup(std::string(id)+"_temp");
  modify->add_compute(fmt::format("{} all temp",id_temp));
  tcomputeflag = 1;

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  id_press = utils::strdup(std::string(id)+"_press");
  modify->add_compute(fmt::format("{} all pressure {}",id_press,id_temp));
  pcomputeflag = 1;
}
