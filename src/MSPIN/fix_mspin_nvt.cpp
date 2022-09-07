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
   Contributing author: Akhlak Mahmood

   Contact:
     Department of Materials Science and Engineering,
     North Carolina State University,
     Raleigh, NC, USA

     amahmoo3@ncsu.edu; mahmoodakhlak@gmail.com
------------------------------------------------------------------------- */

#include "fix_mspin_nvt.h"
#include "error.h"

using namespace LAMMPS_NS;

FixMspinNVT::FixMspinNVT(LAMMPS *lmp, int narg, char **arg) : FixMspinNH(lmp, narg, arg)
{
  // other settings are made by parent

  scalar_flag = 1;
  restart_global = 1;
  extscalar = 1;

  // error checking
  // convert input period to frequency

  if (tstat_flag == 0) error->all(FLERR, "Did not set temperature for fix mspin/nvt");
  if (t_start < 0.0 || t_stop <= 0.0)
    error->all(FLERR, "Target temperature for fix mspin/nvt cannot be 0.0");
  if (t_period <= 0.0) error->all(FLERR, "Fix mspin/nvt period must be > 0.0");
  t_freq = 1.0 / t_period;

  if (t_chain < 1) error->all(FLERR, "Illegal fix mspin/nvt command");
  if (t_iter < 1) error->all(FLERR, "Illegal fix mspin/nvt  command");
  if (t_order != 3 && t_order != 5)
    error->all(FLERR, "Fix mspin/nvt temperature order must be 3 or 5");
}
