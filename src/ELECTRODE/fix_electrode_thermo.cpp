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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (GU), Robert Meissner (Hereon, TUHH)
------------------------------------------------------------------------- */

#include "fix_electrode_thermo.h"

#include "atom.h"
#include "charge_solver.h"
#include "comm.h"
#include "error.h"
#include "fix_electrode_conp.h"
#include "input.h"
#include "random_mars.h"
#include "update.h"
#include "variable.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace LAMMPS_NS;
using namespace std;

static constexpr double SMALL = 1e-16;
static constexpr int NUM_GROUPS = 2;

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
  if (algo != Algo::MATRIX_INV) error->all(FLERR, "Algorithm not allowed in electrode/thermo");
  if (thermo_time < SMALL) error->all(FLERR, "Keyword temp not set or zero in electrode/thermo");

  thermo_random = new RanMars(lmp, thermo_init);
  if (group_psi_var_styles[0] == VarStyle::CONST)
    delta_v_0 = group_psi_const[1] - group_psi_const[0];
}

/* ----------------------------------------------------------------------- */

FixElectrodeThermo::~FixElectrodeThermo()
{
  delete thermo_random;
}

/* ----------------------------------------------------------------------
   configure charge solver with group charges
------------------------------------------------------------------------- */

void FixElectrodeThermo::update_psi_set_constraint()
{
  double const dt = update->dt;
  if (group_psi_var_styles[0] == VarStyle::EQUAL) {
    delta_v_0 = input->variable->compute_equal(group_psi_var_ids[1]) -
        input->variable->compute_equal(group_psi_var_ids[0]);
  }

  // calculate potential for current charges
  auto v_old = charge_solver->compute_potentials();
  assert(v_old.size() == NUM_GROUPS);

  // sums of group charges
  int const nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *q = atom->q;
  auto group_q_old = vector<double>(NUM_GROUPS, 0.);
  for (int g = 0; g < NUM_GROUPS; g++) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & group_bits[g]) { group_q_old[g] += q[i]; }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, group_q_old.data(), NUM_GROUPS, MPI_DOUBLE, MPI_SUM, world);

  // thermo-potentio-stat algorithm by Deissenbeck
  double const delta_v = v_old[1] - v_old[0];
  double const vac_cap = charge_solver->vacuum_capacitance();
  double delta_charge = 0.5 * (group_q_old[1] - group_q_old[0]) -
      vac_cap * (delta_v - delta_v_0) * (1. - exp(-dt / thermo_time));
  delta_charge += sqrt((thermo_temp * vac_cap) * (1. - exp(-2. * dt / thermo_time))) *
      thermo_random->gaussian();

  // configure solver with new group charges
  charge_solver->set_constraint({-delta_charge, delta_charge});
}
