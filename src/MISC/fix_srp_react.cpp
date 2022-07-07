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

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"

#include "fix_bond_break.h"
#include "fix_bond_create.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSRPREACT::FixSRPREACT(LAMMPS *lmp, int narg, char **arg) : FixSRP(lmp, narg, arg)
{
  // default idbreak and idcreate = NULL
  idbreak = nullptr;
  idcreate = nullptr;
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
  if (force->pair_match("hybrid",1) == nullptr && force->pair_match("hybrid/overlay",1) == nullptr)
    error->all(FLERR,"Cannot use pair srp without pair_style hybrid");

  int has_rigid = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (utils::strmatch(modify->fix[i]->style,"^rigid")) ++has_rigid;

  if (has_rigid > 0)
    error->all(FLERR,"Pair srp is not compatible with rigid fixes.");

  if ((bptype < 1) || (bptype > atom->ntypes))
    error->all(FLERR,"Illegal bond particle type");

  // this fix must come before any fix which migrates atoms in its pre_exchange()
  // because this fix's pre_exchange() creates per-atom data structure
  // that data must be current for atom migration to carry it along

  for (int i = 0; i < modify->nfix; i++) {
    if (modify->fix[i] == this) break;
    if (modify->fix[i]->pre_exchange_migrate)
      error->all(FLERR,"Fix SRP comes after a fix which "
                 "migrates atoms in pre_exchange");
  }

  // setup neigh exclusions for diff atom types
  // bond particles do not interact with other types
  // type bptype only interacts with itself

  for (int z = 1; z < atom->ntypes; z++) {
    if (z == bptype)
      continue;
    neighbor->modify_params(fmt::format("exclude type {} {}",z,bptype));
  }

  // find fix bond break
  if( idbreak != nullptr )
    f_bb = (FixBondBreak *) modify->get_fix_by_id(idbreak);

  // find fix bond create
  if( idcreate != nullptr )
    f_bc = (FixBondCreate *) modify->get_fix_by_id(idcreate);

  // free memory
  delete [] idbreak;
  delete [] idcreate;
}

/* ----------------------------------------------------------------------
   rebuild bond particle array
------------------------------------------------------------------------- */
void FixSRPREACT::post_neighbor()
{
  // store ncalls as it is reset in fix srp setup pre force
  int ncalls = neighbor->ncalls;

  if( idbreak != nullptr)
    if(f_bb->breakcount)
    {
      setup_pre_force(0);

      //reset break count before exiting
      // not reseting breakcount would lead to redundant rebuilds
      f_bb->breakcount=0;

      // count additional call during setup_pre_force
      neighbor->ncalls = ncalls+1;
    }
  if( idcreate != nullptr)
    if(f_bc->createcount)
    {
      setup_pre_force(0);

      //reset create count before exiting
      // not reseting createcount would lead to redundant rebuilds
      f_bc->createcount=0;

      // count additional call during setup_pre_force
      neighbor->ncalls = ncalls+1;
    }
}

/* ----------------------------------------------------------------------
   interface with pair class
   in addition to bond type and bond particle type,
     pair srp react sets id of either fix bond break or bond create
------------------------------------------------------------------------- */

int FixSRPREACT::modify_param(int /*narg*/, char **arg)
{
  if (strcmp(arg[0],"btype") == 0) {
    btype = utils::inumeric(FLERR,arg[1],false,lmp);
    return 2;
  }
  if (strcmp(arg[0],"bptype") == 0) {
    bptype = utils::inumeric(FLERR,arg[1],false,lmp);
    return 2;
  }
  if (strcmp(arg[0],"bond/break") == 0) {
    int n = strlen(arg[1]) + 1;
    idbreak = new char[n];
    strcpy(idbreak,arg[1]);
    return 2;
  }
  if (strcmp(arg[0],"bond/create") == 0) {
    int n = strlen(arg[1]) + 1;
    idcreate = new char[n];
    strcpy(idcreate,arg[1]);
    return 2;
  }
  return 0;
}
