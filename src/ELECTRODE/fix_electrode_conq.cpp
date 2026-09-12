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

#include "fix_electrode_conq.h"

#include "charge_solver.h"
#include "comm.h"
#include "error.h"
#include "fix_electrode_conp.h"
#include "input.h"
#include "variable.h"

using namespace LAMMPS_NS;

FixElectrodeConq::FixElectrodeConq(LAMMPS *lmp, int narg, char **arg) :
    FixElectrodeConp(lmp, narg, arg)
{
  if (qtotal_var_style != VarStyle::UNSET)
    error->all(FLERR, "qtotal and symm keyword are not available for {}", this->style);
  // copy const-style values across because update_psi will change group_psi
  group_q = group_psi_const;
}

/* ----------------------------------------------------------------------
   configure charge solver with group charges
------------------------------------------------------------------------- */

void FixElectrodeConq::update_psi_set_constraint()
{
  for (int g = 0; g < num_of_groups; g++) {
    if (group_psi_var_styles[g] == VarStyle::CONST) continue;
    group_q[g] = input->variable->compute_equal(group_psi_var_ids[g]);
  }
  charge_solver->set_constraint(group_q);
}

