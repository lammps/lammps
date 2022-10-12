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
#include "meam.h"

#include <cmath>

using namespace LAMMPS_NS;

template <typename TYPE, int maxi, int maxj>
static inline void setall2d(TYPE (&arr)[maxi][maxj], const TYPE v) {
  for (int i = 0; i < maxi; i++)
    for (int j = 0; j < maxj; j++)
      arr[i][j] = v;
}

template <typename TYPE, int maxi, int maxj, int maxk>
static inline void setall3d(TYPE (&arr)[maxi][maxj][maxk], const TYPE v) {
  for (int i = 0; i < maxi; i++)
    for (int j = 0; j < maxj; j++)
      for (int k = 0; k < maxk; k++)
        arr[i][j][k] = v;
}

void
MEAM::meam_setup_global(int nelt, lattice_t* lat, int* ielement, double* /*atwt*/, double* alpha,
                        double* b0, double* b1, double* b2, double* b3, double* alat, double* esub,
                        double* asub, double* t0, double* t1, double* t2, double* t3, double* rozero,
                        int* ibar)
{

  int i;
  double tmplat[maxelt];

  this->neltypes = nelt;

  for (i = 0; i < nelt; i++) {
    this->lattce_meam[i][i] = lat[i];

    this->ielt_meam[i] = ielement[i];
    this->alpha_meam[i][i] = alpha[i];
    this->beta0_meam[i] = b0[i];
    this->beta1_meam[i] = b1[i];
    this->beta2_meam[i] = b2[i];
    this->beta3_meam[i] = b3[i];
    tmplat[i] = alat[i];
    this->Ec_meam[i][i] = esub[i];
    this->A_meam[i] = asub[i];
    this->t0_meam[i] = t0[i];
    this->t1_meam[i] = t1[i];
    this->t2_meam[i] = t2[i];
    this->t3_meam[i] = t3[i];
    this->rho0_meam[i] = rozero[i];
    this->ibar_meam[i] = ibar[i];

    switch(this->lattce_meam[i][i]) {
      case FCC:
        this->re_meam[i][i] = tmplat[i] / sqrt(2.0);
        break;
      case BCC:
        this->re_meam[i][i] = tmplat[i] * sqrt(3.0) / 2.0;
        break;
      case HCP:
      case DIM:
      case CH4:
      case LIN:
      case ZIG:
      case TRI:
      case SC:
        this->re_meam[i][i] = tmplat[i];
        break;
      case DIA:
      case DIA3:
        this->re_meam[i][i] = tmplat[i] * sqrt(3.0) / 4.0;
        break;
      case B1:
      case B2:
      case C11:
      case L12:
        // do nothing
        break;
      default:
        ;
      //  error
    }
  }

  // Set some defaults
  this->rc_meam = 4.0;
  this->delr_meam = 0.1;
  setall2d(this->attrac_meam, 0.0);
  setall2d(this->repuls_meam, 0.0);
  setall3d(this->Cmax_meam, 2.8);
  setall3d(this->Cmin_meam, 2.0);
  setall2d(this->ebound_meam, (2.8*2.8) / (4.0 * (2.8 - 1.0)));
  setall2d(this->delta_meam, 0.0);
  setall2d(this->nn2_meam, 0);
  setall2d(this->zbl_meam, 1);
  this->gsmooth_factor = 99.0;
  this->augt1 = 1;
  this->ialloy = 0;
  this->mix_ref_t = 0;
  this->emb_lin_neg = 0;
  this->bkgd_dyn = 0;
  this->erose_form = 0;
  // for trimer, zigzag, line refernece structure, sungkwang
  setall2d(this->stheta_meam, 1.0); // stheta = sin(theta/2*pi/180) where theta is 180, so 1.0
  setall2d(this->ctheta_meam, 0.0); // stheta = cos(theta/2*pi/180) where theta is 180, so 0
}
