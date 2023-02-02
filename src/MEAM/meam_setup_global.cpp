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
#include "meam.h"

#include <cmath>

using namespace LAMMPS_NS;

template <typename TYPE, int maxi, int maxj>
static inline void setall2d(TYPE (&arr)[maxi][maxj], const TYPE v)
{
  for (int i = 0; i < maxi; i++)
    for (int j = 0; j < maxj; j++) arr[i][j] = v;
}

template <typename TYPE, int maxi, int maxj, int maxk>
static inline void setall3d(TYPE (&arr)[maxi][maxj][maxk], const TYPE v)
{
  for (int i = 0; i < maxi; i++)
    for (int j = 0; j < maxj; j++)
      for (int k = 0; k < maxk; k++) arr[i][j][k] = v;
}

void MEAM::meam_setup_global(int nelt, lattice_t *lat, int *ielement, double * /*atwt*/,
                             double *alpha, double *b0, double *b1, double *b2, double *b3,
                             double *alat, double *esub, double *asub, double *t0, double *t1,
                             double *t2, double *t3, double *rozero, int *ibar, double *b1m,
                             double *b2m, double *b3m, double *t1m, double *t2m, double *t3m)
{
  int i;
  double tmplat[maxelt];

  neltypes = nelt;

  for (i = 0; i < nelt; i++) {
    lattce_meam[i][i] = lat[i];

    ielt_meam[i] = ielement[i];
    alpha_meam[i][i] = alpha[i];
    beta0_meam[i] = b0[i];
    beta1_meam[i] = b1[i];
    beta2_meam[i] = b2[i];
    beta3_meam[i] = b3[i];
    if (msmeamflag) {
      beta1m_meam[i] = b1m[i];
      beta2m_meam[i] = b2m[i];
      beta3m_meam[i] = b3m[i];
    }
    tmplat[i] = alat[i];
    Ec_meam[i][i] = esub[i];
    A_meam[i] = asub[i];
    t0_meam[i] = t0[i];
    t1_meam[i] = t1[i];
    t2_meam[i] = t2[i];
    t3_meam[i] = t3[i];
    if (msmeamflag) {
      t1m_meam[i] = t1m[i];
      t2m_meam[i] = t2m[i];
      t3m_meam[i] = t3m[i];
    }
    rho0_meam[i] = rozero[i];
    ibar_meam[i] = ibar[i];

    switch (lattce_meam[i][i]) {
      case FCC:
        re_meam[i][i] = tmplat[i] / sqrt(2.0);
        break;
      case BCC:
        re_meam[i][i] = tmplat[i] * sqrt(3.0) / 2.0;
        break;
      case HCP:
      case DIM:
      case CH4:
      case LIN:
      case ZIG:
      case TRI:
      case SC:
        re_meam[i][i] = tmplat[i];
        break;
      case DIA:
      case DIA3:
        re_meam[i][i] = tmplat[i] * sqrt(3.0) / 4.0;
        break;
      case B1:
      case B2:
      case C11:
      case L12:
        // do nothing
        break;
      default:;
        //  error
    }
  }

  // Set some defaults
  rc_meam = 4.0;
  delr_meam = 0.1;
  setall2d(attrac_meam, 0.0);
  setall2d(repuls_meam, 0.0);
  setall3d(Cmax_meam, 2.8);
  setall3d(Cmin_meam, 2.0);
  setall2d(ebound_meam, (2.8 * 2.8) / (4.0 * (2.8 - 1.0)));
  setall2d(delta_meam, 0.0);
  setall2d(nn2_meam, 0);
  setall2d(zbl_meam, 1);
  gsmooth_factor = 99.0;
  augt1 = 1;
  ialloy = 0;
  mix_ref_t = 0;
  emb_lin_neg = 0;
  bkgd_dyn = 0;
  erose_form = 0;
  // for trimer, zigzag, line refernece structure, sungkwang
  setall2d(stheta_meam, 1.0);    // stheta = sin(theta/2*pi/180) where theta is 180, so 1.0
  setall2d(ctheta_meam, 0.0);    // stheta = cos(theta/2*pi/180) where theta is 180, so 0
}
