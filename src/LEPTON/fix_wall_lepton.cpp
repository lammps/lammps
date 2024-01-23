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

#include "fix_wall_lepton.h"
#include "atom.h"
#include "error.h"

#include "Lepton.h"
#include "lepton_utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixWallLepton::FixWallLepton(LAMMPS *lmp, int narg, char **arg) : FixWall(lmp, narg, arg)
{
  dynamic_group_allow = 1;
}

/* ---------------------------------------------------------------------- */

void FixWallLepton::post_constructor()
{
  // check validity of lepton expression

  for (int m = 0; m < nwall; ++m) {
    // remove whitespace and quotes from expression string and then
    // check if the expression can be parsed and evaluated without error
    std::string exp_one = LeptonUtils::condense(lstr[m]);
    try {
      auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(exp_one, lmp));
      auto wallpot = parsed.createCompiledExpression();
      auto wallforce = parsed.differentiate("r").createCompiledExpression();
      wallpot.getVariableReference("r") = 0.0;
      wallforce.getVariableReference("r") = 0.0;
      wallpot.evaluate();
      wallforce.evaluate();
    } catch (std::exception &e) {
      error->all(FLERR, e.what());
    }
  }
}

/* ----------------------------------------------------------------------
   compute the potential energy offset so it can be shifted to zero at the cutoff
------------------------------------------------------------------------- */

void FixWallLepton::precompute(int m)
{
  std::string exp_one = LeptonUtils::condense(lstr[m]);
  auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(exp_one, lmp));
  auto wallpot = parsed.createCompiledExpression();

  try {
    wallpot.getVariableReference("rc") = cutoff[m];
  } catch (std::exception &) {
    ;    // do nothing
  }

  wallpot.getVariableReference("r") = cutoff[m];
  offset[m] = wallpot.evaluate();
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

void FixWallLepton::wall_particle(int m, int which, double coord)
{
  std::string exp_one = LeptonUtils::condense(lstr[m]);
  auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(exp_one, lmp));
  auto wallpot = parsed.createCompiledExpression();
  auto wallforce = parsed.differentiate("r").createCompiledExpression();

  // set cutoff value, if used
  try {
    wallpot.getVariableReference("rc") = cutoff[m];
    wallforce.getVariableReference("rc") = cutoff[m];
  } catch (std::exception &) {
    ;    // do nothing
  }

  double delta, fwall, vn;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = which / 2;
  int side = which % 2;
  if (side == 0) side = -1;

  int onflag = 0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (side < 0)
        delta = x[i][dim] - coord;
      else
        delta = coord - x[i][dim];
      if (delta >= cutoff[m]) continue;
      if (delta <= 0.0) {
        onflag = 1;
        continue;
      }
      wallpot.getVariableReference("r") = delta;
      wallforce.getVariableReference("r") = delta;

      fwall = side * wallforce.evaluate();
      f[i][dim] += fwall;
      ewall[0] += wallpot.evaluate() - offset[m];
      ewall[m + 1] += fwall;

      if (evflag) {
        if (side < 0)
          vn = -fwall * delta;
        else
          vn = fwall * delta;
        v_tally(dim, i, vn);
      }
    }
  }
  if (onflag) error->one(FLERR, "Particle on or inside fix {} surface", style);
}
