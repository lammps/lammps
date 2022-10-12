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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#include "fix_electrode_thermo.h"

#include "atom.h"
#include "error.h"
#include "fix_electrode_conp.h"
#include "force.h"
#include "input.h"
#include "random_mars.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

#define NUM_GROUPS 2
#define SMALL 0.00001

enum { CONST, EQUAL };

/* ----------------------------------------------------------------------- */

//     0        1      2                3    4
// fix fxupdate group1 electrode/thermo pot1 eta couple group2 pot2
FixElectrodeThermo::FixElectrodeThermo(LAMMPS *lmp, int narg, char **arg) :
    FixElectrodeConp(lmp, narg, arg)
{
  if (num_of_groups != NUM_GROUPS)
    error->all(FLERR, "Number of electrodes != two in electrode/thermo");
  if (group_psi_var_styles[0] != group_psi_var_styles[1])
    error->all(FLERR, "Potentials in electrode/thermo must have same style");
  if (symm) error->all(FLERR, "Keyword symm on not allowed in electrode/thermo");
  if (thermo_time < SMALL) error->all(FLERR, "Keyword temp not set or zero in electrode/thermo");

  thermo_random = new RanMars(lmp, thermo_init);
  if (group_psi_var_styles[0] == CONST) delta_psi_0 = group_psi[1] - group_psi[0];
}

/* ----------------------------------------------------------------------- */

FixElectrodeThermo::~FixElectrodeThermo()
{
  delete thermo_random;
}

/* ----------------------------------------------------------------------- */

void FixElectrodeThermo::compute_macro_matrices()
{
  FixElectrodeConp::compute_macro_matrices();
  vac_cap = (macro_capacitance[0][0] * macro_capacitance[1][1] -
             macro_capacitance[0][1] * macro_capacitance[0][1]) /
      (macro_capacitance[0][0] + macro_capacitance[1][1] + 2 * macro_capacitance[0][1]);
}

/* ----------------------------------------------------------------------- */

void FixElectrodeThermo::pre_update()
{
  // total electrode charges after last step, required for update psi
  int const nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *q = atom->q;
  for (int g = 0; g < NUM_GROUPS; g++) {
    group_q_old[g] = 0.;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & group_bits[g]) { group_q_old[g] += q[i]; }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &group_q_old, NUM_GROUPS, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------- */

void FixElectrodeThermo::update_psi()
{
  double const dt = update->dt;

  // group_q_eff is charge that corresponds to potential after previous step
  double group_q_eff[NUM_GROUPS] = {0., 0.};
  for (int g = 0; g < NUM_GROUPS; g++) { group_q_eff[g] = group_q_old[g] - sb_charges[g]; }
  double group_psi_old[NUM_GROUPS] = {0., 0.};
  for (int g = 0; g < NUM_GROUPS; g++) {
    double vtmp = 0;
    for (int h = 0; h < NUM_GROUPS; h++) { vtmp += macro_elastance[g][h] * group_q_eff[h]; }
    group_psi_old[g] = vtmp;
  }
  double const delta_psi = group_psi_old[1] - group_psi_old[0];

  // target potential difference from input parameters
  if (group_psi_var_styles[0] != CONST) {
    delta_psi_0 = input->variable->compute_equal(group_psi_var_ids[1]) -
        input->variable->compute_equal(group_psi_var_ids[0]);
  }

  double delta_charge = 0.5 * (group_q_old[1] - group_q_old[0]) -
      vac_cap * (delta_psi - delta_psi_0) * (1. - exp(-dt / thermo_time));
  delta_charge += sqrt((thermo_temp * vac_cap) * (1. - exp(-2. * dt / thermo_time))) *
      thermo_random->gaussian();

  double const group_remainder_q[NUM_GROUPS] = {-delta_charge - sb_charges[0],
                                                delta_charge - sb_charges[1]};

  for (int g = 0; g < NUM_GROUPS; g++) {
    double vtmp = 0;
    for (int h = 0; h < NUM_GROUPS; h++) { vtmp += macro_elastance[g][h] * group_remainder_q[h]; }
    group_psi[g] = vtmp;
  }
}
