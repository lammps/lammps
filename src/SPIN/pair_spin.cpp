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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "pair_spin.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix_nve_spin.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairSpin::PairSpin(LAMMPS *lmp) : Pair(lmp), emag(nullptr)
{
  hbar = force->hplanck/MY_2PI;
  single_enable = 0;
  respa_enable = 0;
  no_virial_fdotr_compute = 1;
  lattice_flag = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpin::settings(int narg, char **/*arg*/)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect number of args in pair_style pair/spin command");

  // pair spin/* need the metal unit style

  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Spin pair styles require metal units");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSpin::init_style()
{
  if (!atom->sp_flag)
    error->all(FLERR,"Pair spin requires atom/spin style");

  // checking if nve/spin or neb/spin is a listed fix

  if ((comm->me == 0) && ((modify->get_fix_by_style("^nve/spin").size()
                           + modify->get_fix_by_style("^neb/spin").size()) == 0))
    error->warning(FLERR,"Using spin pair style without nve/spin or neb/spin");

  // check if newton pair is on

  if ((force->newton_pair == 0) && (comm->me == 0))
    error->all(FLERR,"Pair style spin requires newton pair on");

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);

  // get the lattice_flag from nve/spin

  auto fixes = modify->get_fix_by_style("^nve/spin");
  if (fixes.size() == 1)
    lattice_flag = (dynamic_cast<FixNVESpin *>(fixes.front()))->lattice_flag;
  else if (fixes.size() > 1)
    error->warning(FLERR,"Using multiple instances of fix nve/spin or neb/spin");

  // init. size of energy stacking lists

  nlocal_max = atom->nlocal;
  memory->grow(emag,nlocal_max,"pair/spin:emag");
}
