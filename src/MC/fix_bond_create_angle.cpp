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
   Contributing authors:
   Joao Gregorio, Oliver Henrich (University of Strathclyde, Glasgow, UK)
------------------------------------------------------------------------- */

#include "fix_bond_create_angle.h"

#include "atom.h"
#include "domain.h"
#include "error.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

int FixBondCreateAngle::constrain(int i, int j, double amin, double amax)
{
  int **nspecial = atom->nspecial;

  // an atom without 1-2 neighbors imposes no constraint
  // the angle at each atom with 1-2 neighbors must be within amin and amax
  // the comparison also rejects a NaN angle

  if (nspecial[i][0] != 0) {
    double theta = bond_angle(i, j);
    if (!((amin <= theta) && (theta <= amax))) return 0;
  }

  if (nspecial[j][0] != 0) {
    double theta = bond_angle(j, i);
    if (!((amin <= theta) && (theta <= amax))) return 0;
  }

  return 1;
}

/* ----------------------------------------------------------------------
   angle at atom i between atom j and the first 1-2 neighbor of atom i
   j is a neighbor list entry of i or vice versa, so it is the closest image
   atom->map() may return any image of the 1-2 neighbor, use the closest one
------------------------------------------------------------------------- */

double FixBondCreateAngle::bond_angle(int i, int j)
{
  double **x = atom->x;

  int k = atom->map(atom->special[i][0]);
  if (k < 0)
    error->one(FLERR, Error::NOLASTLINE,
               "Fix {} needs ghost atoms from further away; use "
               "comm_modify cutoff to extend the communication range",
               style);
  k = domain->closest_image(i, k);

  double delx1 = x[i][0] - x[j][0];
  double dely1 = x[i][1] - x[j][1];
  double delz1 = x[i][2] - x[j][2];
  double delx2 = x[i][0] - x[k][0];
  double dely2 = x[i][1] - x[k][1];
  double delz2 = x[i][2] - x[k][2];

  double c = (delx1 * delx2 + dely1 * dely2 + delz1 * delz2) /
      sqrt((delx1 * delx1 + dely1 * dely1 + delz1 * delz1) *
           (delx2 * delx2 + dely2 * dely2 + delz2 * delz2));
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  return acos(c);
}
