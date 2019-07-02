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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include <cstring>
#include "pair_spin.h"
#include "atom.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "math_const.h"
#include "modify.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairSpin::PairSpin(LAMMPS *lmp) : Pair(lmp)
{
  hbar = force->hplanck/MY_2PI;
  single_enable = 0;
  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

PairSpin::~PairSpin() {}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpin::settings(int narg, char **/*arg*/)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect number of args in pair_style pair/spin command");

  // pair/spin need the metal unit style

  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"pair/spin style requires metal units");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSpin::init_style()
{
  if (!atom->sp_flag)
    error->all(FLERR,"Pair spin requires atom/spin style");

  // checking if nve/spin is a listed fix

  int ifix = 0;
  while (ifix < modify->nfix) {
    if (strcmp(modify->fix[ifix]->style,"nve/spin") == 0) break;
    ifix++;
  }
  if (ifix == modify->nfix)
    error->all(FLERR,"pair/spin style requires nve/spin");
}
