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

#include "fix_rigid_nvt.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixRigidNVT::FixRigidNVT(LAMMPS *lmp, int narg, char **arg) :
  FixRigidNH(lmp, narg, arg)
{
  // other settings are made by parent

  scalar_flag = 1;
  restart_global = 1;
  extscalar = 1;

  // error checking
  // convert input period to frequency

  if (tstat_flag == 0)
    error->all(FLERR,"Did not set temperature for fix rigid/nvt");
  if (t_start < 0.0 || t_stop <= 0.0)
    error->all(FLERR,"Target temperature for fix rigid/nvt cannot be 0.0");
  if (t_period <= 0.0) error->all(FLERR,"Fix rigid/nvt period must be > 0.0");
  t_freq = 1.0 / t_period;

  if (t_chain < 1) error->all(FLERR,"Illegal fix rigid/nvt command");
  if (t_iter < 1) error->all(FLERR,"Illegal fix rigid/nvt  command");
  if (t_order != 3 && t_order != 5)
    error->all(FLERR,"Fix rigid/nvt temperature order must be 3 or 5");
}
