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

#include "math_const.h"

#include <cassert>
#include <cmath>

using namespace LAMMPS_NS;
using namespace MEAM_NS;
using MathConst::MY_PI;


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

/* ----------------------------------------------------------------------
   Validate index of "arg(i, j)"-type parameters
------------------------------------------------------------------------- */

static void check_index(int num, int lim, int nidx, int *idx /*idx(3)*/, int *ierr)
{
  //: idx[0..2]
  *ierr = 0;
  if (nidx < num) {
    *ierr = 2;
    return;
  }

  for (int i = 0; i < num; i++) {
    if ((idx[i] < 0) || (idx[i] >= lim)) {
      *ierr = 3;
      return;
    }
  }

  if (STRICT_IJ_ORDER && (num >= 2) && (idx[0] > idx[1])) {
    *ierr = 4;
    return;
  }
}

/* ----------------------------------------------------------------------
   Ingest element-specific parameters from the library file
------------------------------------------------------------------------- */

void MEAM::setup_library(int nelt, lattice_t *lat, int *ielement, double * /*atwt*/,
                         double *alpha, double *b0, double *b1, double *b2, double *b3,
                         double *alat, double *esub, double *asub, double *t0, double *t1,
                         double *t2, double *t3, double *rozero, int *ibar)
{
  int i;
  double tmplat[MAXELT];

  neltypes = nelt;

  for (i = 0; i < nelt; i++) {
    lattce_meam[i][i] = lat[i];

    ielt_meam[i] = ielement[i];
    alpha_meam[i][i] = alpha[i];
    beta0_meam[i] = b0[i];
    beta1_meam[i] = b1[i];
    beta2_meam[i] = b2[i];
    beta3_meam[i] = b3[i];
    tmplat[i] = alat[i];
    Ec_meam[i][i] = esub[i];
    A_meam[i] = asub[i];
    t0_meam[i] = t0[i];
    t1_meam[i] = t1[i];
    t2_meam[i] = t2[i];
    t3_meam[i] = t3[i];
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

/* ----------------------------------------------------------------------
   Ingest element-specific parameters from the library file for MS-MEAM
------------------------------------------------------------------------- */

void MEAM::setup_library_ms(int nelt, double *b1m, double *b2m, double *b3m, double *t1m, double *t2m, double *t3m)
{
  int i;

  assert(neltypes == nelt);
  assert(msmeamflag != 0);

  for (i = 0; i < nelt; i++) {
    beta1m_meam[i] = b1m[i];
    beta2m_meam[i] = b2m[i];
    beta3m_meam[i] = b3m[i];
    t1m_meam[i] = t1m[i];
    t2m_meam[i] = t2m[i];
    t3m_meam[i] = t3m[i];
  }

  ialloy = 1;
}

/* ----------------------------------------------------------------------
   Ingest pair parameters from the param file

   The "which" argument corresponds to the index of the "keyword" array
   in pair_meam.cpp:

   0 = Ec_meam
   1 = alpha_meam
   2 = rho0_meam
   3 = delta_meam
   4 = lattce_meam
   5 = attrac_meam
   6 = repuls_meam
   7 = nn2_meam
   8 = Cmin_meam
   9 = Cmax_meam
   10 = rc_meam
   11 = delr_meam
   12 = augt1
   13 = gsmooth_factor
   14 = re_meam
   15 = ialloy
   16 = mixture_ref_t
   17 = erose_form
   18 = zbl_meam
   19 = emb_lin_neg
   20 = bkgd_dyn
   21 = theta

   The returned errorflag has the following meanings:

   0 = no error
   1 = "which" out of range / invalid keyword
   2 = not enough indices given
   3 = an element index is out of range
------------------------------------------------------------------------- */

void MEAM::setup_param(int which, double value, int nindex, int *index /*index(3)*/, int *errorflag)
{
  //: index[0..2]
  int i1, i2;
  lattice_t vlat;
  *errorflag = 0;

  switch (which) {
    //     0 = Ec_meam
    case 0:
      check_index(2, neltypes, nindex, index, errorflag);
      if (*errorflag != 0) return;
      Ec_meam[index[0]][index[1]] = value;
      break;

    //     1 = alpha_meam
    case 1:
      check_index(2, neltypes, nindex, index, errorflag);
      if (*errorflag != 0) return;
      alpha_meam[index[0]][index[1]] = value;
      break;

    //     2 = rho0_meam
    case 2:
      check_index(1, neltypes, nindex, index, errorflag);
      if (*errorflag != 0) return;
      rho0_meam[index[0]] = value;
      break;

    //     3 = delta_meam
    case 3:
      check_index(2, neltypes, nindex, index, errorflag);
      if (*errorflag != 0) return;
      delta_meam[index[0]][index[1]] = value;
      break;

    //     4 = lattce_meam
    case 4:
      check_index(2, neltypes, nindex, index, errorflag);
      if (*errorflag != 0) return;
      vlat = (lattice_t) value;

      lattce_meam[index[0]][index[1]] = vlat;
      break;

    //     5 = attrac_meam
    case 5:
      check_index(2, neltypes, nindex, index, errorflag);
      if (*errorflag != 0) return;
      attrac_meam[index[0]][index[1]] = value;
      break;

    //     6 = repuls_meam
    case 6:
      check_index(2, neltypes, nindex, index, errorflag);
      if (*errorflag != 0) return;
      repuls_meam[index[0]][index[1]] = value;
      break;

    //     7 = nn2_meam
    case 7:
      check_index(2, neltypes, nindex, index, errorflag);
      if (*errorflag != 0) return;
      i1 = std::min(index[0], index[1]);
      i2 = std::max(index[0], index[1]);
      nn2_meam[i1][i2] = (int) value;
      break;

    //     8 = Cmin_meam
    case 8:
      check_index(3, neltypes, nindex, index, errorflag);
      if (*errorflag != 0) return;
      Cmin_meam[index[0]][index[1]][index[2]] = value;
      break;

    //     9 = Cmax_meam
    case 9:
      check_index(3, neltypes, nindex, index, errorflag);
      if (*errorflag != 0) return;
      Cmax_meam[index[0]][index[1]][index[2]] = value;
      break;

    //     10 = rc_meam
    case 10:
      rc_meam = value;
      break;

    //     11 = delr_meam
    case 11:
      delr_meam = value;
      break;

    //     12 = augt1
    case 12:
      augt1 = (int) value;
      break;

    //     13 = gsmooth
    case 13:
      gsmooth_factor = value;
      break;

    //     14 = re_meam
    case 14:
      check_index(2, neltypes, nindex, index, errorflag);
      if (*errorflag != 0) return;
      re_meam[index[0]][index[1]] = value;
      break;

    //     15 = ialloy
    case 15:
      ialloy = (int) value;
      break;

    //     16 = mixture_ref_t
    case 16:
      mix_ref_t = (int) value;
      break;

    //     17 = erose_form
    case 17:
      erose_form = (int) value;
      break;

    //     18 = zbl_meam
    case 18:
      check_index(2, neltypes, nindex, index, errorflag);
      if (*errorflag != 0) return;
      i1 = std::min(index[0], index[1]);
      i2 = std::max(index[0], index[1]);
      zbl_meam[i1][i2] = (int) value;
      break;

    //     19 = emb_lin_neg
    case 19:
      emb_lin_neg = (int) value;
      break;

    //     20 = bkgd_dyn
    case 20:
      bkgd_dyn = (int) value;
      break;

    //     21 = theta
    // see alloyparams(void) in meam_setup_done.cpp
    case 21:
      check_index(2, neltypes, nindex, index, errorflag);
      if (*errorflag != 0) return;
      i1 = std::min(index[0], index[1]);
      i2 = std::max(index[0], index[1]);
      // we don't use theta, instead stheta and ctheta
      stheta_meam[i1][i2] = sin(value / 2 * MY_PI / 180.0);
      ctheta_meam[i1][i2] = cos(value / 2 * MY_PI / 180.0);
      break;

    default:
      *errorflag = 1;
  }
}

/* ----------------------------------------------------------------------
   Finish setup, fill constant and index arrays
------------------------------------------------------------------------- */

void MEAM::setup_finish(double* cutmax)
{
  int nv2, nv3, m, n, p;

  //     Force cutoff
  cutforce = rc_meam;
  cutforcesq = cutforce * cutforce;

  //     Pass cutoff back to calling program
  *cutmax = cutforce;

  //     Augment t1 term
  for (int i = 0; i < MAXELT; i++)
    t1_meam[i] = t1_meam[i] + augt1 * 3.0 / 5.0 * t3_meam[i];

  //     Compute off-diagonal alloy parameters
  alloyparams();

  // indices and factors for Voight notation
  nv2 = 0;
  nv3 = 0;
  for (m = 0; m < 3; m++) {
    for (n = m; n < 3; n++) {
      vind2D[m][n] = nv2;
      vind2D[n][m] = nv2;
      nv2 = nv2 + 1;
      for (p = n; p < 3; p++) {
        vind3D[m][n][p] = nv3;
        vind3D[m][p][n] = nv3;
        vind3D[n][m][p] = nv3;
        vind3D[n][p][m] = nv3;
        vind3D[p][m][n] = nv3;
        vind3D[p][n][m] = nv3;
        nv3 = nv3 + 1;
      }
    }
  }

  v2D[0] = 1;
  v2D[1] = 2;
  v2D[2] = 2;
  v2D[3] = 1;
  v2D[4] = 2;
  v2D[5] = 1;

  v3D[0] = 1;
  v3D[1] = 3;
  v3D[2] = 3;
  v3D[3] = 3;
  v3D[4] = 6;
  v3D[5] = 3;
  v3D[6] = 1;
  v3D[7] = 3;
  v3D[8] = 3;
  v3D[9] = 1;

  nv2 = 0;
  for (m = 0; m < neltypes; m++) {
    for (n = m; n < neltypes; n++) {
      eltind[m][n] = nv2;
      eltind[n][m] = nv2;
      nv2 = nv2 + 1;
    }
  }
}

/* ----------------------------------------------------------------------
   Fill off-diagonal alloy parameters
------------------------------------------------------------------------- */

void MEAM::alloyparams()
{

  int i, j, k;
  double eb;

  // Loop over pairs
  for (i = 0; i < neltypes; i++) {
    for (j = 0; j < neltypes; j++) {
      // Treat off-diagonal pairs
      // If i>j, set all equal to i<j case (which has already been set,
      // here or in the input file)
      if (i > j) {
        re_meam[i][j] = re_meam[j][i];
        Ec_meam[i][j] = Ec_meam[j][i];
        alpha_meam[i][j] = alpha_meam[j][i];
        lattce_meam[i][j] = lattce_meam[j][i];
        nn2_meam[i][j] = nn2_meam[j][i];
        // theta for lin,tri,zig references
        stheta_meam[i][j] = stheta_meam[j][i];
        ctheta_meam[i][j] = ctheta_meam[j][i];
        // If i<j and term is unset, use default values (e.g. mean of i-i and
        // j-j)
      } else if (j > i) {
        if (iszero(Ec_meam[i][j])) {
          switch (lattce_meam[i][j]) {
            case L12:
              Ec_meam[i][j] = (3 * Ec_meam[i][i] + Ec_meam[j][j]) / 4.0 - delta_meam[i][j];
              break;
            case C11:
              if (lattce_meam[i][i] == DIA)
                Ec_meam[i][j] = (2 * Ec_meam[i][i] + Ec_meam[j][j]) / 3.0 - delta_meam[i][j];
              else
                Ec_meam[i][j] = (Ec_meam[i][i] + 2 * Ec_meam[j][j]) / 3.0 - delta_meam[i][j];
              break;
            default:
              Ec_meam[i][j] = (Ec_meam[i][i] + Ec_meam[j][j]) / 2.0 - delta_meam[i][j];
          }
        }
        if (iszero(alpha_meam[i][j]))
          alpha_meam[i][j] = (alpha_meam[i][i] + alpha_meam[j][j]) / 2.0;
        if (iszero(re_meam[i][j]))
          re_meam[i][j] = (re_meam[i][i] + re_meam[j][j]) / 2.0;
      }
    }
  }

  // Cmin[i][k][j] is symmetric in i-j, but not k.  For all triplets
  // where i>j, set equal to the i<j element.  Likewise for Cmax.
  for (i = 1; i < neltypes; i++) {
    for (j = 0; j < i; j++) {
      for (k = 0; k < neltypes; k++) {
        Cmin_meam[i][j][k] = Cmin_meam[j][i][k];
        Cmax_meam[i][j][k] = Cmax_meam[j][i][k];
      }
    }
  }

  // ebound gives the squared distance such that, for rik2 or rjk2>ebound,
  // atom k definitely lies outside the screening function ellipse (so
  // there is no need to calculate its effects).  Here, compute it for all
  // triplets [i][j][k] so that ebound[i][j] is the maximized over k
  for (i = 0; i < neltypes; i++) {
    for (j = 0; j < neltypes; j++) {
      for (k = 0; k < neltypes; k++) {
        eb = (Cmax_meam[i][j][k] * Cmax_meam[i][j][k]) / (4.0 * (Cmax_meam[i][j][k] - 1.0));
        ebound_meam[i][j] = std::max(ebound_meam[i][j], eb);
      }
    }
  }
}
