/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
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

#include "fix_rigid_npt_omp.h"
#include <cstring>
#include "modify.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixRigidNPTOMP::FixRigidNPTOMP(LAMMPS *lmp, int narg, char **arg) :
  FixRigidNHOMP(lmp, narg, arg)
{
  // other setting are made by parent

  scalar_flag = 1;
  restart_global = 1;
  extscalar = 1;

  // error checks

  if (tstat_flag == 0 || pstat_flag == 0)
    error->all(FLERR,"Did not set temp or press for fix rigid/npt/omp");
  if (t_start <= 0.0 || t_stop <= 0.0)
    error->all(FLERR,"Target temperature for fix rigid/npt/omp cannot be 0.0");
  if (p_start[0] < 0.0 || p_start[1] < 0.0 || p_start[2] < 0.0 ||
      p_stop[0] < 0.0 || p_stop[1] < 0.0 || p_stop[2] < 0.0)
    error->all(FLERR,"Target pressure for fix rigid/npt/omp cannot be 0.0");

  if (t_period <= 0.0) error->all(FLERR,"Fix rigid/npt/omp period must be > 0.0");

  // thermostat chain parameters

  if (t_chain < 1) error->all(FLERR,"Illegal fix_modify command");
  if (t_iter < 1) error->all(FLERR,"Illegal fix_modify command");
  if (t_order != 3 && t_order != 5)
    error->all(FLERR,"Fix_modify order must be 3 or 5");

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

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "temp";
  modify->add_compute(3,newarg);
  delete [] newarg;
  tcomputeflag = 1;

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  n = strlen(id) + 7;
  id_press = new char[n];
  strcpy(id_press,id);
  strcat(id_press,"_press");

  newarg = new char*[4];
  newarg[0] = id_press;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure";
  newarg[3] = id_temp;
  modify->add_compute(4,newarg);
  delete [] newarg;
  pcomputeflag = 1;
}
