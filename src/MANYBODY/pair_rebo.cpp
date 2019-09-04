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

#include "pair_rebo.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairREBO::PairREBO(LAMMPS *lmp) : PairAIREBO(lmp) {
  variant = REBO_2;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairREBO::settings(int narg, char ** /* arg */)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

  cutlj = 0.0;
  ljflag = torflag = 0;
}

/* ----------------------------------------------------------------------
   initialize spline knot values
------------------------------------------------------------------------- */

void PairREBO::spline_init() {
  PairAIREBO::spline_init();

  PCCf[0][2] = 0.007860700254745;
  PCCf[0][3] = 0.016125364564267;
  PCCf[1][1] = 0.003026697473481;
  PCCf[1][2] = 0.006326248241119;
  PCCf[2][0] = 0.;
  PCCf[2][1] = 0.003179530830731;

  for (int nH = 0; nH < 4; nH++) {
    for (int nC = 0; nC < 4; nC++) {
      double y[4] = {0}, y1[4] = {0}, y2[4] = {0};
      y[0] = PCCf[nC][nH];
      y[1] = PCCf[nC][nH+1];
      y[2] = PCCf[nC+1][nH];
      y[3] = PCCf[nC+1][nH+1];
      Spbicubic_patch_coeffs(nC, nC+1, nH, nH+1, y, y1, y2, &pCC[nC][nH][0]);
      y[0] = PCHf[nC][nH];
      y[1] = PCHf[nC][nH+1];
      y[2] = PCHf[nC+1][nH];
      y[3] = PCHf[nC+1][nH+1];
      Spbicubic_patch_coeffs(nC, nC+1, nH, nH+1, y, y1, y2, &pCH[nC][nH][0]);
    }
  }
}
