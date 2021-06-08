// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

#include "fix_rigid_nph.h"

#include "error.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixRigidNPH::FixRigidNPH(LAMMPS *lmp, int narg, char **arg) :
  FixRigidNH(lmp, narg, arg)
{
  // other setting are made by parent

  scalar_flag = 1;
  restart_global = 1;
  extscalar = 1;

  // error checks

  if (pstat_flag == 0)
    error->all(FLERR,"Did not set pressure for fix rigid/nph");
  if (tstat_flag == 1)
    error->all(FLERR,"Cannot set temperature for fix rigid/nph");

  // convert input periods to frequency

  p_freq[0] = p_freq[1] = p_freq[2] = 0.0;

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
