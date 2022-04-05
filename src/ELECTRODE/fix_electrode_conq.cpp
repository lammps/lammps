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

#include "fix_electrode_conq.h"

#include "error.h"
#include "fix_electrode_conp.h"
#include "input.h"
#include "variable.h"

using namespace LAMMPS_NS;

#define SMALL 0.00001

enum { CONST, EQUAL };

//     0        1      2              3    4
// fix fxupdate group1 electrode/conp pot1 eta couple group2 pot2
FixElectrodeConq::FixElectrodeConq(LAMMPS *lmp, int narg, char **arg) :
    FixElectrodeConp(lmp, narg, arg)
{
  // copy const-style values across because update_psi will change group_psi
  group_q = group_psi;

  if (symm) { error->all(FLERR, "Keyword symm on not allowed in electrode/conq"); }
}

void FixElectrodeConq::update_psi()
{
  // don't need MPI_Barrier because always preceded by MPI_Allreduce
  for (int g = 0; g < num_of_groups; g++) {
    if (group_psi_var_styles[g] == CONST) continue;
    group_q[g] = input->variable->compute_equal(group_psi_var_ids[g]);
  }

  std::vector<double> group_remainder_q(num_of_groups);
  for (int g = 0; g < num_of_groups; g++) { group_remainder_q[g] = group_q[g] - sb_charges[g]; }

  for (int g = 0; g < num_of_groups; g++) {
    double vtmp = 0;
    for (int h = 0; h < num_of_groups; h++) {
      vtmp += macro_elastance[g][h] * group_remainder_q[h];
    }
    group_psi[g] = vtmp;
  }
}
