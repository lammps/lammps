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
   Contributing author: Vaibhav Palkar (Kuksenok Lab, Clemson University)
   based on the pair srp code by Timothy Sirk (ARL)

This pair style interfaces the pair style srp with the stochastic reaction
fixes bond/break and bond/create by updating pseudo beads corresponding to
bonds as bond breaking and formation takes place. This is useful in
simulation of reactions in polymers with soft potentials such as DPD.

See the doc page for pair_style srp/react command for usage instructions.

There is an example script for this package in examples/PACKAGES/srp_react/.
------------------------------------------------------------------------- */

#include "fix_srp_react.h"

#include "modify.h"
#include "neighbor.h"

#include "fix_bond_break.h"
#include "fix_bond_create.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

// clang-format on
/* ---------------------------------------------------------------------- */

FixSRPREACT::FixSRPREACT(LAMMPS *lmp, int narg, char **arg) :
    FixSRP(lmp, narg, arg), f_bb(nullptr), idbreak(nullptr), f_bc(nullptr), idcreate(nullptr)
{
  pair_name = "srp/react";
}

/* ---------------------------------------------------------------------- */

FixSRPREACT::~FixSRPREACT()
{
  // free memory
  delete[] idbreak;
  delete[] idcreate;
}

/* ---------------------------------------------------------------------- */

int FixSRPREACT::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= PRE_EXCHANGE;
  mask |= POST_RUN;
  mask |= POST_NEIGHBOR;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSRPREACT::init()
{
  FixSRP::init();

  // find fix bond break
  if (idbreak) f_bb = (FixBondBreak *) modify->get_fix_by_id(idbreak);

  // find fix bond create
  if (idcreate) f_bc = (FixBondCreate *) modify->get_fix_by_id(idcreate);
}

/* ----------------------------------------------------------------------
   rebuild bond particle array
------------------------------------------------------------------------- */
void FixSRPREACT::post_neighbor()
{
  // store ncalls as it is reset in fix srp setup pre force
  int ncalls = neighbor->ncalls;

  if (f_bb) {
    if (f_bb->breakcount) {
      setup_pre_force(0);

      //reset break count before exiting
      // not reseting breakcount would lead to redundant rebuilds
      f_bb->breakcount = 0;

      // count additional call during setup_pre_force
      neighbor->ncalls = ncalls + 1;
    }
  }

  if (f_bc) {
    if (f_bc->createcount) {
      setup_pre_force(0);

      //reset create count before exiting
      // not reseting createcount would lead to redundant rebuilds
      f_bc->createcount = 0;

      // count additional call during setup_pre_force
      neighbor->ncalls = ncalls + 1;
    }
  }
}

/* ----------------------------------------------------------------------
   interface with pair class
   in addition to bond type and bond particle type,
     pair srp react sets id of either fix bond break or bond create
------------------------------------------------------------------------- */

int FixSRPREACT::modify_param(int /*narg*/, char **arg)
{
  if (strcmp(arg[0], "btype") == 0) {
    btype = utils::inumeric(FLERR, arg[1], false, lmp);
    return 2;
  }
  if (strcmp(arg[0], "bptype") == 0) {
    bptype = utils::inumeric(FLERR, arg[1], false, lmp);
    return 2;
  }
  if (strcmp(arg[0], "bond/break") == 0) {
    idbreak = utils::strdup(arg[1]);
    return 2;
  }
  if (strcmp(arg[0], "bond/create") == 0) {
    idcreate = utils::strdup(arg[1]);
    return 2;
  }
  return 0;
}
