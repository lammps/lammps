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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#include "fix_electrode_conq.h"

#include "comm.h"
#include "error.h"
#include "fix_electrode_conp.h"
#include "group.h"
#include "input.h"
#include "variable.h"

using namespace LAMMPS_NS;

FixElectrodeConq::FixElectrodeConq(LAMMPS *lmp, int narg, char **arg) :
    FixElectrodeConp(lmp, narg, arg)
{
  // copy const-style values across because update_psi will change group_psi
  group_q = group_psi;

  if (symm) {
    if (num_of_groups == 1)
      error->all(FLERR, "Keyword symm on not allowed in electrode/conq with only one electrode");
    if (comm->me == 0)
      error->warning(FLERR,
                     "Fix electrode/conq with keyword symm ignores the charge setting for the last "
                     "electrode listed");
    if (algo != Algo::MATRIX_INV) {
      double last_q = 0.;
      for (int g = 0; g < num_of_groups - 1; g++) last_q -= group_q[g];
      group_q.back() = last_q;    // needed for CG algos
    }
  }
}

void FixElectrodeConq::update_psi()
{
  int const numsymm = num_of_groups - ((symm) ? 1 : 0);
  bool symm_update_back = false;
  for (int g = 0; g < numsymm; g++) {
    if (group_psi_var_styles[g] == VarStyle::CONST) continue;
    group_q[g] = input->variable->compute_equal(group_psi_var_ids[g]);
    symm_update_back = true;
  }
  if (algo == Algo::MATRIX_INV) {
    std::vector<double> group_remainder_q(num_of_groups, 0.);
    for (int g = 0; g < numsymm; g++) { group_remainder_q[g] = group_q[g] - sb_charges[g]; }
    for (int g = 0; g < num_of_groups; g++) {
      double vtmp = 0;
      for (int h = 0; h < num_of_groups; h++) {
        vtmp += macro_elastance[g][h] * group_remainder_q[h];
      }
      group_psi[g] = vtmp;
    }
  } else {
    if (symm && symm_update_back) {    // needed for CG algos
      double last_q = 0.;
      for (int g = 0; g < num_of_groups - 1; g++) last_q -= group_q[g];
      group_q.back() = last_q;
    }
    for (double &g : group_psi) g = 0;
  }
}

/* ----------------------------------------------------------------------
   Correct charge of each electrode to target charge by adding a homogeneous charge
------------------------------------------------------------------------- */

std::vector<double> FixElectrodeConq::constraint_correction(std::vector<double> x)
{
  int const n = x.size();
  auto sums = std::vector<double>(num_of_groups, 0);
  for (int i = 0; i < n; i++) sums[iele_to_group_local[i]] += x[i];
  MPI_Allreduce(MPI_IN_PLACE, &sums.front(), num_of_groups, MPI_DOUBLE, MPI_SUM, world);
  for (int g = 0; g < num_of_groups; g++) {
    sums[g] -= group_q[g];
    sums[g] /= group->count(groups[g]);
  }
  for (int i = 0; i < n; i++) x[i] -= sums[iele_to_group_local[i]];
  return x;
}

/* ----------------------------------------------------------------------
   Project into direction that conserves charge of each electrode (cf. M. Shariff (1995))
------------------------------------------------------------------------- */

std::vector<double> FixElectrodeConq::constraint_projection(std::vector<double> x)
{
  int const n = x.size();
  auto sums = std::vector<double>(num_of_groups, 0);
  for (int i = 0; i < n; i++) sums[iele_to_group_local[i]] += x[i];
  MPI_Allreduce(MPI_IN_PLACE, &sums.front(), num_of_groups, MPI_DOUBLE, MPI_SUM, world);
  for (int g = 0; g < num_of_groups; g++) sums[g] /= group->count(groups[g]);
  for (int i = 0; i < n; i++) x[i] -= sums[iele_to_group_local[i]];
  return x;
}

/* ----------------------------------------------------------------------
   Recompute group potential as average for output if using cg algo
------------------------------------------------------------------------- */

void FixElectrodeConq::recompute_potential(std::vector<double> b, std::vector<double> q_local)
{
  int const n = b.size();
  auto a = ele_ele_interaction(q_local);
  auto psi_sums = std::vector<double>(num_of_groups, 0);
  for (int i = 0; i < n; i++) { psi_sums[iele_to_group_local[i]] += (a[i] + b[i]) / evscale; }
  MPI_Allreduce(MPI_IN_PLACE, &psi_sums.front(), num_of_groups, MPI_DOUBLE, MPI_SUM, world);
  for (int g = 0; g < num_of_groups; g++) group_psi[g] = psi_sums[g] / group->count(groups[g]);
}
