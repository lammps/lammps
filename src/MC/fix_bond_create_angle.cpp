// clang-format off
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

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

int FixBondCreateAngle::constrain(int i, int j, double amin, double amax)
{
  double **x = atom->x;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  double v1x = 0.0;
  double v1y = 0.0;
  double v1z = 0.0;
  double v2x = 0.0;
  double v2y = 0.0;
  double v2z = 0.0;

  double angle1,angle2;

  int flag = 0;

  // pass if both atoms have no neighbors: bond is always formed

  if (nspecial[i][0] == 0 && nspecial[j][0] == 0) flag = 1;

  // pass if i has at least one neighbor and angle constraint is met

  if (nspecial[i][0] != 0 && nspecial[j][0] == 0) {

    // calculate first vector between i and j
    // calculate second vector between i and its first neighbor

    v1x = x[i][0] - x[j][0];
    v1y = x[i][1] - x[j][1];
    v1z = x[i][2] - x[j][2];
    v2x = x[i][0] - x[atom->map(special[i][0])][0];
    v2y = x[i][1] - x[atom->map(special[i][0])][1];
    v2z = x[i][2] - x[atom->map(special[i][0])][2];

    // calculate angle between both vectors
    // set flag if the angle constraint is met

    angle1 = acos((v1x*v2x + v1y*v2y + v1z*v2z)/
        (sqrt(v1x*v1x + v1y*v1y + v1z*v1z)*sqrt(v2x*v2x + v2y*v2y + v2z*v2z)));
    if (amin <= angle1 && angle1 <= amax) flag = 1;
  }

  // pass if j has at least one neighbor and angle constraint is met

  if (nspecial[j][0] != 0 && nspecial[i][0] == 0) {

    // calculate first vector between j and i
    // calculate second vector between j and its first neighbor

    v1x = x[j][0] - x[i][0];
    v1y = x[j][1] - x[i][1];
    v1z = x[j][2] - x[i][2];
    v2x = x[j][0] - x[atom->map(special[j][0])][0];
    v2y = x[j][1] - x[atom->map(special[j][0])][1];
    v2z = x[j][2] - x[atom->map(special[j][0])][2];

    // calculate angle between both vectors
    // set flag if angle constraint is met

    angle1 = acos((v1x*v2x + v1y*v2y + v1z*v2z) /
        (sqrt(v1x*v1x + v1y*v1y + v1z*v1z)*sqrt(v2x*v2x + v2y*v2y + v2z*v2z)));
    if (amin <= angle1 && angle1 <= amax) flag = 1;
  }

  // pass if both i and j have at least one neighbor and angle constraint is met

  if (nspecial[i][0] != 0 && nspecial[j][0] != 0) {

    // Calculate 1st angle
    // calculate first vector between i and j
    // calculate second vector between i and its first neighbor

    v1x = x[i][0] - x[j][0];
    v1y = x[i][1] - x[j][1];
    v1z = x[i][2] - x[j][2];
    v2x = x[i][0] - x[atom->map(special[i][0])][0];
    v2y = x[i][1] - x[atom->map(special[i][0])][1];
    v2z = x[i][2] - x[atom->map(special[i][0])][2];

    // calculate angle between both vectors

    angle1 = acos((v1x*v2x + v1y*v2y + v1z*v2z) /
        (sqrt(v1x*v1x + v1y*v1y + v1z*v1z)*sqrt(v2x*v2x + v2y*v2y + v2z*v2z)));

    // Calculate 2nd angle
    // calculate first vector between i and j
    // calculate second vector between i and its first neighbor

    v1x = x[j][0] - x[i][0];
    v1y = x[j][1] - x[i][1];
    v1z = x[j][2] - x[i][2];
    v2x = x[j][0] - x[atom->map(special[j][0])][0];
    v2y = x[j][1] - x[atom->map(special[j][0])][1];
    v2z = x[j][2] - x[atom->map(special[j][0])][2];

    // calculate angle between both vectors

    angle2 = acos((v1x*v2x + v1y*v2y + v1z*v2z) /
        (sqrt(v1x*v1x + v1y*v1y + v1z*v1z)*sqrt(v2x*v2x + v2y*v2y + v2z*v2z)));

    // set flag if the angle constraint is met

    if ((amin <= angle1 && angle1 <= amax) && (amin <= angle2 && angle2 <= amax))
      flag = 1;
  }

  return flag;
}
