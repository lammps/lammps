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

void PairREBO::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

  cutlj = 0.0;
  ljflag = torflag = 0;
}

/* ----------------------------------------------------------------------
   initialize spline knot values
------------------------------------------------------------------------- */

void PairREBO::spline_init()
{
  int i,j,k;

  for (i = 0; i < 5; i++) {
    for (j = 0; j < 5; j++) {
      PCCf[i][j] = 0.0;
      PCCdfdx[i][j] = 0.0;
      PCCdfdy[i][j] = 0.0;
      PCHf[i][j] = 0.0;
      PCHdfdx[i][j] = 0.0;
      PCHdfdy[i][j] = 0.0;
    }
  }

  PCCf[0][2] = 0.007860700254745;
  PCCf[0][3] = 0.016125364564267;
  PCCf[1][1] = 0.003026697473481;
  PCCf[1][2] = 0.006326248241119;
  PCCf[2][0] = 0.;
  PCCf[2][1] = 0.003179530830731;

  PCHf[0][1] = 0.2093367328250380;
  PCHf[0][2] = -0.064449615432525;
  PCHf[0][3] = -0.303927546346162;
  PCHf[1][0] = 0.010;
  PCHf[1][1] = -0.1251234006287090;
  PCHf[1][2] = -0.298905245783;
  PCHf[2][0] = -0.1220421462782555;
  PCHf[2][1] = -0.3005291724067579;
  PCHf[3][0] = -0.307584705066;

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

  for (i = 0; i < 5; i++) {
    for (j = 0; j < 5; j++) {
      for (k = 0; k < 10; k++) {
        piCCf[i][j][k] = 0.0;
        piCCdfdx[i][j][k] = 0.0;
        piCCdfdy[i][j][k] = 0.0;
        piCCdfdz[i][j][k] = 0.0;
        piCHf[i][j][k] = 0.0;
        piCHdfdx[i][j][k] = 0.0;
        piCHdfdy[i][j][k] = 0.0;
        piCHdfdz[i][j][k] = 0.0;
        piHHf[i][j][k] = 0.0;
        piHHdfdx[i][j][k] = 0.0;
        piHHdfdy[i][j][k] = 0.0;
        piHHdfdz[i][j][k] = 0.0;
        Tf[i][j][k] = 0.0;
        Tdfdx[i][j][k] = 0.0;
        Tdfdy[i][j][k] = 0.0;
        Tdfdz[i][j][k] = 0.0;
      }
    }
  }

  for (i = 3; i < 10; i++) piCCf[0][0][i] = 0.0049586079;
  piCCf[1][0][1] = 0.021693495;
  piCCf[0][1][1] = 0.021693495;
  for (i = 2; i < 10; i++) piCCf[1][0][i] = 0.0049586079;
  for (i = 2; i < 10; i++) piCCf[0][1][i] = 0.0049586079;
  piCCf[1][1][1] = 0.05250;
  piCCf[1][1][2] = -0.002088750;
  for (i = 3; i < 10; i++) piCCf[1][1][i] = -0.00804280;
  piCCf[2][0][1] = 0.024698831850;
  piCCf[0][2][1] = 0.024698831850;
  piCCf[2][0][2] = -0.00597133450;
  piCCf[0][2][2] = -0.00597133450;
  for (i = 3; i < 10; i++) piCCf[2][0][i] = 0.0049586079;
  for (i = 3; i < 10; i++) piCCf[0][2][i] = 0.0049586079;
  piCCf[2][1][1] = 0.00482478490;
  piCCf[1][2][1] = 0.00482478490;
  piCCf[2][1][2] = 0.0150;
  piCCf[1][2][2] = 0.0150;
  piCCf[2][1][3] = -0.010;
  piCCf[1][2][3] = -0.010;
  piCCf[2][1][4] = -0.01168893870;
  piCCf[1][2][4] = -0.01168893870;
  piCCf[2][1][5] = -0.013377877400;
  piCCf[1][2][5] = -0.013377877400;
  piCCf[2][1][6] = -0.015066816000;
  piCCf[1][2][6] = -0.015066816000;
  for (i = 7; i < 10; i++) piCCf[2][1][i] = -0.015066816000;
  for (i = 7; i < 10; i++) piCCf[1][2][i] = -0.015066816000;
  piCCf[2][2][1] = 0.0472247850;
  piCCf[2][2][2] = 0.0110;
  piCCf[2][2][3] = 0.0198529350;
  piCCf[2][2][4] = 0.01654411250;
  piCCf[2][2][5] = 0.013235290;
  piCCf[2][2][6] = 0.00992646749999 ;
  piCCf[2][2][7] = 0.006617644999;
  piCCf[2][2][8] = 0.00330882250;
  piCCf[3][0][1] = -0.05989946750;
  piCCf[0][3][1] = -0.05989946750;
  piCCf[3][0][2] = -0.05989946750;
  piCCf[0][3][2] = -0.05989946750;
  for (i = 3; i < 10; i++) piCCf[3][0][i] = 0.0049586079;
  for (i = 3; i < 10; i++) piCCf[0][3][i] = 0.0049586079;
  piCCf[3][1][2] = -0.0624183760;
  piCCf[1][3][2] = -0.0624183760;
  for (i = 3; i < 10; i++) piCCf[3][1][i] = -0.0624183760;
  for (i = 3; i < 10; i++) piCCf[1][3][i] = -0.0624183760;
  piCCf[3][2][1] = -0.02235469150;
  piCCf[2][3][1] = -0.02235469150;
  for (i = 2; i < 10; i++) piCCf[3][2][i] = -0.02235469150;
  for (i = 2; i < 10; i++) piCCf[2][3][i] = -0.02235469150;

  piCCdfdx[2][1][1] = -0.026250;
  piCCdfdx[2][1][5] = -0.0271880;
  piCCdfdx[2][1][6] = -0.0271880;
  for (i = 7; i < 10; i++) piCCdfdx[2][1][i] = -0.0271880;
  piCCdfdx[1][3][2] = 0.0187723882;
  for (i = 2; i < 10; i++) piCCdfdx[2][3][i] = 0.031209;

  piCCdfdy[1][2][1] = -0.026250;
  piCCdfdy[1][2][5] = -0.0271880;
  piCCdfdy[1][2][6] = -0.0271880;
  for (i = 7; i < 10; i++) piCCdfdy[1][2][i] = -0.0271880;
  piCCdfdy[3][1][2] = 0.0187723882;
  for (i = 2; i < 10; i++) piCCdfdy[3][2][i] = 0.031209;

  piCCdfdz[1][1][2] = -0.0302715;
  piCCdfdz[2][1][4] = -0.0100220;
  piCCdfdz[1][2][4] = -0.0100220;
  piCCdfdz[2][1][5] = -0.0100220;
  piCCdfdz[1][2][5] = -0.0100220;
  for (i = 4; i < 9; i++) piCCdfdz[2][2][i] = -0.0033090;

  //  make top end of piCC flat instead of zero
  i = 4;
  for (j = 0; j < 4; j++){
      for (k = 1; k < 11; k++){
          piCCf[i][j][k] = piCCf[i-1][j][k];
      }
  }
  for (i = 0; i < 4; i++){ // also enforces some symmetry
      for (j = i+1; j < 5; j++){
          for (k = 1; k < 11; k++){
              piCCf[i][j][k] = piCCf[j][i][k];
          }
      }
  }
  for (k = 1; k < 11; k++) piCCf[4][4][k] = piCCf[3][4][k];
  k = 10;
  for (i = 0; i < 5; i++){
      for (j = 0; j < 5; j++){
      piCCf[i][j][k] = piCCf[i][j][k-1];
      }
  }

  piCHf[1][1][1] = -0.050;
  piCHf[1][1][2] = -0.050;
  piCHf[1][1][3] = -0.30;
  for (i = 4; i < 10; i++) piCHf[1][1][i] = -0.050;
  for (i = 5; i < 10; i++) piCHf[2][0][i] = -0.004523893758064;
  for (i = 5; i < 10; i++) piCHf[0][2][i] = -0.004523893758064;
  piCHf[2][1][2] = -0.250;
  piCHf[1][2][2] = -0.250;
  piCHf[2][1][3] = -0.250;
  piCHf[1][2][3] = -0.250;
  piCHf[3][1][1] = -0.10;
  piCHf[1][3][1] = -0.10;
  piCHf[3][1][2] = -0.125;
  piCHf[1][3][2] = -0.125;
  piCHf[3][1][3] = -0.125;
  piCHf[1][3][3] = -0.125;
  for (i = 4; i < 10; i++) piCHf[3][1][i] = -0.10;
  for (i = 4; i < 10; i++) piCHf[1][3][i] = -0.10;

  // make top end of piCH flat instead of zero
 // also enforces some symmetry

  i = 4;
  for (j = 0; j < 4; j++){
      for (k = 1; k < 11; k++){
          piCHf[i][j][k] = piCHf[i-1][j][k];
      }
  }
  for (i = 0; i < 4; i++){
      for (j = i+1; j < 5; j++){
          for (k = 1; k < 11; k++){
              piCHf[i][j][k] = piCHf[j][i][k];
          }
      }
  }
  for (k = 1; k < 11; k++) piCHf[4][4][k] = piCHf[3][4][k];
  k = 10;
  for (i = 0; i < 5; i++){
      for (j = 0; j < 5; j++){
      piCHf[i][j][k] = piCHf[i][j][k-1];
      }
  }

  piHHf[1][1][1] = 0.124915958;

  Tf[2][2][1] = -0.035140;
  for (i = 2; i < 10; i++) Tf[2][2][i] = -0.0040480;

  for (int nH = 0; nH < 4; nH++) {
    for (int nC = 0; nC < 4; nC++) {
      // Note: Spline knot values exist up to "10", but are never used because
      // they are clamped down to 9.
      for (int nConj = 0; nConj < 9; nConj++) {
        double y[8] = {0}, y1[8] = {0}, y2[8] = {0}, y3[8] = {0};
        #define FILL_KNOTS_TRI(dest, src)      \
          dest[0] = src[nC+0][nH+0][nConj+0];  \
          dest[1] = src[nC+0][nH+0][nConj+1];  \
          dest[2] = src[nC+0][nH+1][nConj+0];  \
          dest[3] = src[nC+0][nH+1][nConj+1];  \
          dest[4] = src[nC+1][nH+0][nConj+0];  \
          dest[5] = src[nC+1][nH+0][nConj+1];  \
          dest[6] = src[nC+1][nH+1][nConj+0];  \
          dest[7] = src[nC+1][nH+1][nConj+1];
        FILL_KNOTS_TRI(y, piCCf)
        FILL_KNOTS_TRI(y1, piCCdfdx)
        FILL_KNOTS_TRI(y2, piCCdfdy)
        FILL_KNOTS_TRI(y3, piCCdfdz)
        Sptricubic_patch_coeffs(nC, nC+1, nH, nH+1, nConj, nConj+1, y, y1, y2, y3, &piCC[nC][nH][nConj][0]);
        FILL_KNOTS_TRI(y, piCHf)
        FILL_KNOTS_TRI(y1, piCHdfdx)
        FILL_KNOTS_TRI(y2, piCHdfdy)
        FILL_KNOTS_TRI(y3, piCHdfdz)
        Sptricubic_patch_coeffs(nC, nC+1, nH, nH+1, nConj, nConj+1, y, y1, y2, y3, &piCH[nC][nH][nConj][0]);
        FILL_KNOTS_TRI(y, piHHf)
        FILL_KNOTS_TRI(y1, piHHdfdx)
        FILL_KNOTS_TRI(y2, piHHdfdy)
        FILL_KNOTS_TRI(y3, piHHdfdz)
        Sptricubic_patch_coeffs(nC, nC+1, nH, nH+1, nConj, nConj+1, y, y1, y2, y3, &piHH[nC][nH][nConj][0]);
        FILL_KNOTS_TRI(y, Tf)
        FILL_KNOTS_TRI(y1, Tdfdx)
        FILL_KNOTS_TRI(y2, Tdfdy)
        FILL_KNOTS_TRI(y3, Tdfdz)
        Sptricubic_patch_coeffs(nC, nC+1, nH, nH+1, nConj, nConj+1, y, y1, y2, y3, &Tijc[nC][nH][nConj][0]);
        #undef FILL_KNOTS_TRI
      }
    }
  }
}
