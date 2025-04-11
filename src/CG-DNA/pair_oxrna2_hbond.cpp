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
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "pair_oxrna2_hbond.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairOxrna2Hbond::PairOxrna2Hbond(LAMMPS *lmp) : PairOxdnaHbond(lmp)
{
  single_enable = 0;
  writedata = 1;
  trim_flag = 0;

  // sequence-specific base-pairing strength
  // A:0 C:1 G:2 U:3, 5'- [i][j] -3'

  alpha_hb[0][0] = 1.00000;
  alpha_hb[0][1] = 1.00000;
  alpha_hb[0][2] = 1.00000;
  alpha_hb[0][3] = 0.94253;

  alpha_hb[1][0] = 1.00000;
  alpha_hb[1][1] = 1.00000;
  alpha_hb[1][2] = 1.22288;
  alpha_hb[1][3] = 1.00000;

  alpha_hb[2][0] = 1.00000;
  alpha_hb[2][1] = 1.22288;
  alpha_hb[2][2] = 1.00000;
  alpha_hb[2][3] = 0.58655;

  alpha_hb[3][0] = 0.94253;
  alpha_hb[3][1] = 1.00000;
  alpha_hb[3][2] = 0.58655;
  alpha_hb[3][3] = 1.00000;
}
