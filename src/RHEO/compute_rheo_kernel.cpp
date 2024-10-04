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
   Joel Clemmer (SNL), Thomas O'Connor (CMU)
----------------------------------------------------------------------- */

#include "compute_rheo_kernel.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_interface.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"
#include "utils.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace RHEO_NS;
using namespace MathConst;
using namespace MathExtra;

// max value of Mdim 1 + dim + dim * (dim + 1) / 2 with dim = 3
static constexpr int MAX_MDIM = 12;

// declare LAPACK functions

extern "C" {
  void dpotrf_(const char *uplo, const int *n, double *a, const int *lda, int *info);
  void dpotri_(const char *uplo, const int *n, double *a, const int *lda, int *info);
}

/* ---------------------------------------------------------------------- */

ComputeRHEOKernel::ComputeRHEOKernel(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), coordination(nullptr), fix_rheo(nullptr), C(nullptr), C0(nullptr),
    list(nullptr), compute_interface(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute rheo/kernel command");

  kernel_style = utils::inumeric(FLERR, arg[3], false, lmp);

  if (kernel_style == QUINTIC || kernel_style == WENDLANDC4) {
    correction_order = -1;
  } else if (kernel_style == RK0) {
    correction_order = 0;
  } else if (kernel_style == RK1) {
    correction_order = 1;
  } else if (kernel_style == RK2) {
    correction_order = 2;
  }

  dim = domain->dimension;

  comm_forward = 1;
  ncor = 0;
  Mdim = 0;
  if (kernel_style == RK1) {
    Mdim = 1 + dim;
    ncor = 1 + dim;
    comm_forward = ncor * Mdim;
  } else if (kernel_style == RK2) {
    //Polynomial basis size (up to quadratic order)
    Mdim = 1 + dim + dim * (dim + 1) / 2;
    //Number of sets of correction coefficients  (1 x y xx yy)  + z zz (3D)
    ncor = 1 + 2 * dim;
    comm_forward = ncor * Mdim;
  }

  comm_forward_save = comm_forward;
  corrections_calculated = 0;
  lapack_error_flag = 0;
}

/* ---------------------------------------------------------------------- */

ComputeRHEOKernel::~ComputeRHEOKernel()
{
  memory->destroy(coordination);
  memory->destroy(C);
  memory->destroy(C0);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOKernel::init()
{
  neighbor->add_request(this, NeighConst::REQ_FULL);

  interface_flag = fix_rheo->interface_flag;
  compute_interface = fix_rheo->compute_interface;

  zmin = fix_rheo->zmin_kernel;
  cut = fix_rheo->cut;
  cutsq = cut * cut;
  cutinv = 1.0 / cut;
  cutsqinv = cutinv * cutinv;

  if (kernel_style != WENDLANDC4) {
    if (dim == 3) {
      pre_w = 1.0 / (120.0 * MY_PI) * 27.0 * cutsqinv * cutinv;
      pre_wp = pre_w * 3.0 * cutinv;
    } else {
      pre_w = 7.0 / (478.0 * MY_PI) * 9 * cutsqinv;
      pre_wp = pre_w * 3.0 * cutinv;
    }
  } else {
    if (dim == 3) {
      pre_w = 495.0 / (32.0 * MY_PI * cutsq * cut);
      pre_wp = pre_w * cutinv;
    } else {
      pre_w = 9.0 / (MY_PI * cutsq);
      pre_wp = pre_w * cutinv;
    }
  }

  nmax_store = atom->nmax;
  memory->create(coordination, nmax_store, "rheo:coordination");
  if (kernel_style == RK0) {
    memory->create(C0, nmax_store, "rheo/kernel:C0");
  } else if (kernel_style == RK1) {
    memory->create(C, nmax_store, ncor, Mdim, "rheo/kernel:C");
  } else if (kernel_style == RK2) {
    memory->create(C, nmax_store, ncor, Mdim, "rheo/kernel:C");
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOKernel::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOKernel::check_corrections(int i)
{
  // Skip if there were lapack errors for this atom
  if (lapack_error_flag)
    if (lapack_error_tags.find(atom->tag[i]) != lapack_error_tags.end()) return 0;

  // Skip if undercoordinated
  if (coordination[i] < zmin) return 0;

  // Skip if corrections not yet calculated
  if (!corrections_calculated) return 0;

  return 1;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_w_self()
{
  if (kernel_style == WENDLANDC4)
    return calc_w_wendlandc4(0.0);
  else
    return calc_w_quintic(0.0);
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_w(int i, int j, double delx, double dely, double delz, double r)
{
  double w = 0.0;
  int corrections_i, corrections_j, corrections;

  if (kernel_style == WENDLANDC4) return calc_w_wendlandc4(r);

  if (kernel_style != QUINTIC) {
    corrections_i = check_corrections(i);
    corrections_j = check_corrections(j);
    corrections = corrections_i & corrections_j;
  } else {
    corrections = 0;
  }

  if (!corrections)
    w = calc_w_quintic(r);
  else if (kernel_style == RK0)
    w = calc_w_rk0(i, j, r);
  else if (kernel_style == RK1)
    w = calc_w_rk1(i, j, delx, dely, delz, r);
  else if (kernel_style == RK2)
    w = calc_w_rk2(i, j, delx, dely, delz, r);

  return w;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_dw(int i, int j, double delx, double dely, double delz, double r)
{
  double wp;
  int corrections_i, corrections_j;

  if (kernel_style == WENDLANDC4) return calc_dw_wendlandc4(delx, dely, delz, r, dWij, dWji);

  if (kernel_style != QUINTIC) {
    corrections_i = check_corrections(i);
    corrections_j = check_corrections(j);
  }

  // Calc wp and default dW's, a bit inefficient but can redo later
  wp = calc_dw_quintic(delx, dely, delz, r, dWij, dWji);

  // Overwrite if there are corrections
  if (kernel_style == RK1) {
    if (corrections_i) calc_dw_rk1(i, delx, dely, delz, r, dWij);
    if (corrections_j) calc_dw_rk1(j, -delx, -dely, -delz, r, dWji);
  } else if (kernel_style == RK2) {
    if (corrections_i) calc_dw_rk2(i, delx, dely, delz, r, dWij);
    if (corrections_j) calc_dw_rk2(j, -delx, -dely, -delz, r, dWji);
  }

  return wp;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_w_quintic(double r)
{
  double w, tmp1, tmp2, tmp3, tmp1sq, tmp2sq, tmp3sq, s;
  s = r * 3.0 * cutinv;

  if (s > 3.0) { w = 0.0; }

  if (s <= 3.0) {
    tmp3 = 3.0 - s;
    tmp3sq = tmp3 * tmp3;
    w = tmp3sq * tmp3sq * tmp3;
  }
  if (s <= 2.0) {
    tmp2 = 2.0 - s;
    tmp2sq = tmp2 * tmp2;
    w -= 6.0 * tmp2sq * tmp2sq * tmp2;
  }
  if (s <= 1.0) {
    tmp1 = 1.0 - s;
    tmp1sq = tmp1 * tmp1;
    w += 15.0 * tmp1sq * tmp1sq * tmp1;
  }

  w *= pre_w;

  Wij = w;
  Wji = w;

  return w;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_dw_quintic(double delx, double dely, double delz, double r,
                                          double *dW1, double *dW2)
{
  double wp, tmp1, tmp2, tmp3, tmp1sq, tmp2sq, tmp3sq, s, wprinv;

  s = r * 3.0 * cutinv;

  if (s > 3.0) { wp = 0.0; }
  if (s <= 3.0) {
    tmp3 = 3.0 - s;
    tmp3sq = tmp3 * tmp3;
    wp = -5.0 * tmp3sq * tmp3sq;
  }
  if (s <= 2.0) {
    tmp2 = 2.0 - s;
    tmp2sq = tmp2 * tmp2;
    wp += 30.0 * tmp2sq * tmp2sq;
  }
  if (s <= 1.0) {
    tmp1 = 1.0 - s;
    tmp1sq = tmp1 * tmp1;
    wp -= 75.0 * tmp1sq * tmp1sq;
  }

  wp *= pre_wp;
  wprinv = wp / r;
  dW1[0] = delx * wprinv;
  dW1[1] = dely * wprinv;
  dW1[2] = delz * wprinv;

  dW2[0] = -delx * wprinv;
  dW2[1] = -dely * wprinv;
  dW2[2] = -delz * wprinv;

  return wp;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_w_wendlandc4(double r)
{
  double w, tmp6, s;
  s = r * cutinv;

  if (s > 1.0) {
    w = 0.0;
  } else {
    tmp6 = (1.0 - s) * (1.0 - s);
    tmp6 *= tmp6 * tmp6;
    w = tmp6 * (1.0 + 6.0 * s + 35.0 * THIRD * s * s);
  }

  w *= pre_w;

  Wij = w;
  Wji = w;

  return w;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_dw_wendlandc4(double delx, double dely, double delz, double r,
                                             double *dW1, double *dW2)
{
  double wp, tmp1, tmp5, tmp6, s, wprinv;

  s = r * cutinv;

  if (s > 1.0) {
    wp = 0.0;
  } else {
    tmp1 = 1.0 - s;
    tmp5 = tmp1 * tmp1;
    tmp5 = tmp5 * tmp5 * tmp1;
    tmp6 = tmp5 * tmp1;
    wp = tmp6 * (6.0 + 70.0 * THIRD * s);
    wp -= 6 * tmp5 * (1.0 + 6.0 * s + 35.0 * THIRD * s * s);
  }

  wp *= pre_wp;
  wprinv = wp / r;
  dW1[0] = delx * wprinv;
  dW1[1] = dely * wprinv;
  dW1[2] = delz * wprinv;

  dW2[0] = -delx * wprinv;
  dW2[1] = -dely * wprinv;
  dW2[2] = -delz * wprinv;

  return wp;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_w_rk0(int i, int j, double r)
{
  double w;

  w = calc_w_quintic(r);

  Wij = C0[i] * w;
  Wji = C0[j] * w;

  return w;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_w_rk1(int i, int j, double delx, double dely, double delz, double r)
{
  int b;
  double w, dx[3], H[MAX_MDIM];

  dx[0] = delx;
  dx[1] = dely;
  dx[2] = delz;
  w = calc_w_quintic(r);

  if (dim == 2) {
    H[0] = 1.0;
    H[1] = dx[0] * cutinv;
    H[2] = dx[1] * cutinv;
  } else {
    H[0] = 1.0;
    H[1] = dx[0] * cutinv;
    H[2] = dx[1] * cutinv;
    H[3] = dx[2] * cutinv;
  }
  Wij = 0;
  for (b = 0; b < Mdim; b++) {
    Wij += C[i][0][b] * H[b];    // C columns: 1 x y (z) xx yy (zz)
  }
  Wij *= w;

  //Now compute Wji
  H[1] *= -1;
  H[2] *= -1;
  if (dim == 3) H[3] *= -1;

  Wji = 0;
  for (b = 0; b < Mdim; b++) {
    Wji += C[j][0][b] * H[b];    // C columns: 1 x y (z) xx yy (zz)
  }
  Wji *= w;

  return w;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_w_rk2(int i, int j, double delx, double dely, double delz, double r)
{
  int b;
  double w, dx[3], H[MAX_MDIM];
  dx[0] = delx;
  dx[1] = dely;
  dx[2] = delz;
  w = calc_w_quintic(r);

  if (dim == 2) {
    H[0] = 1.0;
    H[1] = dx[0] * cutinv;
    H[2] = dx[1] * cutinv;
    H[3] = 0.5 * dx[0] * dx[0] * cutsqinv;
    H[4] = 0.5 * dx[1] * dx[1] * cutsqinv;
    H[5] = dx[0] * dx[1] * cutsqinv;
  } else {
    H[0] = 1.0;
    H[1] = dx[0] * cutinv;
    H[2] = dx[1] * cutinv;
    H[3] = dx[2] * cutinv;
    H[4] = 0.5 * dx[0] * dx[0] * cutsqinv;
    H[5] = 0.5 * dx[1] * dx[1] * cutsqinv;
    H[6] = 0.5 * dx[2] * dx[2] * cutsqinv;
    H[7] = dx[0] * dx[1] * cutsqinv;
    H[8] = dx[0] * dx[2] * cutsqinv;
    H[9] = dx[1] * dx[2] * cutsqinv;
  }
  Wij = 0;
  for (b = 0; b < Mdim; b++) {
    Wij += C[i][0][b] * H[b];    // C columns: 1 x y (z) xx yy (zz)
  }
  Wij *= w;

  //Now compute Wji
  H[1] *= -1;
  H[2] *= -1;
  if (dim == 3) H[3] *= -1;

  Wji = 0;
  for (b = 0; b < Mdim; b++) {
    Wji += C[j][0][b] * H[b];    // C columns: 1 x y (z) xx yy (zz)
  }
  Wji *= w;

  return w;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOKernel::calc_dw_rk1(int i, double delx, double dely, double delz, double r,
                                    double *dW)
{
  int a, b;
  double w, dx[3], H[MAX_MDIM];
  dx[0] = delx;
  dx[1] = dely;
  dx[2] = delz;

  w = calc_w_quintic(r);

  //Populate correction basis
  if (dim == 2) {
    H[0] = 1.0;
    H[1] = dx[0] * cutinv;
    H[2] = dx[1] * cutinv;
  } else {
    H[0] = 1.0;
    H[1] = dx[0] * cutinv;
    H[2] = dx[1] * cutinv;
    H[3] = dx[2] * cutinv;
  }

  // dWij[] = dWx dWy (dWz)
  //compute derivative operators
  for (a = 0; a < dim; a++) {
    dW[a] = 0.0;
    for (b = 0; b < Mdim; b++) {
      //First derivative kernels
      dW[a] += C[i][1 + a][b] * H[b];    // C columns: 1 x y (z)
    }
    dW[a] *= w;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOKernel::calc_dw_rk2(int i, double delx, double dely, double delz, double r,
                                    double *dW)
{
  int a, b;
  double w, dx[3], H[MAX_MDIM];
  dx[0] = delx;
  dx[1] = dely;
  dx[2] = delz;

  w = calc_w_quintic(r);

  //Populate correction basis
  if (dim == 2) {
    H[0] = 1.0;
    H[1] = dx[0] * cutinv;
    H[2] = dx[1] * cutinv;
    H[3] = 0.5 * dx[0] * dx[0] * cutsqinv;
    H[4] = 0.5 * dx[1] * dx[1] * cutsqinv;
    H[5] = dx[0] * dx[1] * cutsqinv;
  } else {
    H[0] = 1.0;
    H[1] = dx[0] * cutinv;
    H[2] = dx[1] * cutinv;
    H[3] = dx[2] * cutinv;
    H[4] = 0.5 * dx[0] * dx[0] * cutsqinv;
    H[5] = 0.5 * dx[1] * dx[1] * cutsqinv;
    H[6] = 0.5 * dx[2] * dx[2] * cutsqinv;
    H[7] = dx[0] * dx[1] * cutsqinv;
    H[8] = dx[0] * dx[2] * cutsqinv;
    H[9] = dx[1] * dx[2] * cutsqinv;
  }

  // dWij[] = dWx dWy (dWz)
  //compute derivative operators
  for (a = 0; a < dim; a++) {
    dW[a] = 0.0;
    for (b = 0; b < Mdim; b++) {
      //First derivative kernels
      dW[a] += C[i][1 + a][b] * H[b];    // C columns: 1 x y (z) xx yy (zz)
    }
    dW[a] *= w;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOKernel::compute_peratom()
{
  lapack_error_flag = 0;
  lapack_error_tags.clear();

  if (kernel_style == QUINTIC) return;
  corrections_calculated = 1;

  int i, j, ii, jj, inum, jnum, a, b, lapack_error;
  double xtmp, ytmp, ztmp, r, rsq, w, vj, rhoj;
  double dx[3];

  double **x = atom->x;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *rho = atom->rho;
  int *status = atom->rheo_status;
  tagint *tag = atom->tag;

  int *ilist, *jlist, *numneigh, **firstneigh;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Grow arrays if necessary
  if (nmax_store < atom->nmax) grow_arrays(atom->nmax);

  if (kernel_style == RK0) {

    double M;
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      jlist = firstneigh[i];
      jnum = numneigh[i];

      //Initialize M to zero:
      M = 0;

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        dx[0] = xtmp - x[j][0];
        dx[1] = ytmp - x[j][1];
        dx[2] = ztmp - x[j][2];
        rsq = lensq3(dx);

        if (rsq < cutsq) {
          r = sqrt(rsq);
          w = calc_w_quintic(r);
          rhoj = rho[j];
          if (interface_flag)
            if (status[j] & PHASECHECK) rhoj = compute_interface->correct_rho(j);

          if (rmass)
            vj = rmass[j] / rhoj;
          else
            vj = mass[type[j]] / rhoj;
          M += w * vj;
        }
      }

      // Inverse of 1x1 matrix
      if (coordination[i] >= zmin) C0[i] = 1.0 / M;
    }
  } else if (correction_order > 0) {

    // Moment matrix M and polynomial basis vector cut (1d for LAPACK compatibility)
    double H[MAX_MDIM], M[MAX_MDIM * MAX_MDIM];

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      jlist = firstneigh[i];
      jnum = numneigh[i];

      // Zero upper-triangle M and cut (will be symmetric):
      for (a = 0; a < Mdim; a++) {
        for (b = a; b < Mdim; b++) {
          M[a * Mdim + b] = 0;
        }
      }

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        dx[0] = xtmp - x[j][0];
        dx[1] = ytmp - x[j][1];
        dx[2] = ztmp - x[j][2];

        rsq = lensq3(dx);

        if (rsq < cutsq) {
          r = sqrt(rsq);
          w = calc_w_quintic(r);

          rhoj = rho[j];
          if (interface_flag)
            if (status[j] & PHASECHECK) rhoj = compute_interface->correct_rho(j);

          if (rmass)
            vj = rmass[j] / rhoj;
          else
            vj = mass[type[j]] / rhoj;

          //Populate the H-vector of polynomials (2D)
          if (dim == 2) {
            H[0] = 1.0;
            H[1] = dx[0] * cutinv;
            H[2] = dx[1] * cutinv;
            if (kernel_style == RK2) {
              H[3] = 0.5 * dx[0] * dx[0] * cutsqinv;
              H[4] = 0.5 * dx[1] * dx[1] * cutsqinv;
              H[5] = dx[0] * dx[1] * cutsqinv;
            }
          } else {
            H[0] = 1.0;
            H[1] = dx[0] * cutinv;
            H[2] = dx[1] * cutinv;
            H[3] = dx[2] * cutinv;
            if (kernel_style == RK2) {
              H[4] = 0.5 * dx[0] * dx[0] * cutsqinv;
              H[5] = 0.5 * dx[1] * dx[1] * cutsqinv;
              H[6] = 0.5 * dx[2] * dx[2] * cutsqinv;
              H[7] = dx[0] * dx[1] * cutsqinv;
              H[8] = dx[0] * dx[2] * cutsqinv;
              H[9] = dx[1] * dx[2] * cutsqinv;
            }
          }

          // Populate the upper triangle
          for (a = 0; a < Mdim; a++) {
            for (b = a; b < Mdim; b++) {
              M[a * Mdim + b] += H[a] * H[b] * w * vj;
            }
          }
        }
      }

      // Populate the lower triangle from the symmetric entries of M:
      for (a = 0; a < Mdim; a++) {
        for (b = a; b < Mdim; b++) {
          M[b * Mdim + a] = M[a * Mdim + b];
        }
      }

      // Skip if undercoordinated
      if (coordination[i] < zmin) continue;

      // Use LAPACK to get Minv, use Cholesky decomposition since the
      // polynomials are independent, M is symmetrix & positive-definite
      const char uplo = 'U';
      dpotrf_(&uplo, &Mdim, M, &Mdim, &lapack_error);

      if (lapack_error) {
        // Revert to uncorrected SPH for this particle
        lapack_error_flag = 1;
        lapack_error_tags.insert(tag[i]);

        // check if not positive-definite
        if (lapack_error > 0)
          error->warning(FLERR, "Failed DPOTRF2 decomposition in rheo/kernel, info = {}",
                         lapack_error);

        continue;
      }

      // M is now M^-1
      dpotri_(&uplo, &Mdim, M, &Mdim, &lapack_error);

      // make result matrix symmetric
      for (a = 0; a < Mdim; a++) {
        for (b = a + 1; b < Mdim; b++) {
          M[a * Mdim + b] = M[b * Mdim + a];
        }
      }

      // Correction coefficients are columns of M^-1 multiplied by an appropriate coefficient
      // Solve the linear system several times to get coefficientns
      // M:    1   x   y  (z)  x^2  y^2 (z^2) xy   (xz)   (yz)
      // ----------------------------------------------------------
      //       0   1   2       3     4        5                 || 2D indexing
      //       0   1   2   3   4     5   6    7     8      9    || 3D indexing
      //  W    1   .   .   .    .    .   .    .     .      .
      // dWx   .  -1   .   .    .    .   .    .     .      .
      // dWy   .   .  -1   .    .    .   .    .     .      .
      // dWz   .   .   .  (-1)  .    .   .    .     .      .
      // d2Wx  .   .   .   .    2    .   .    .     .      .
      // d2Wy  .   .   .   .    .    2   .    .     .      .
      // d2Wz  .   .   .   .    .    .  (2)   .     .      .

      //0 1 2 3 4
      //0 1 2 3 4 5 6

      // Pack coefficients into C
      for (a = 0; a < Mdim; a++) {
        C[i][0][a] = M[a * Mdim + 0];    // all rows of column 0
        for (b = 0; b < dim; b++) {
          //First derivatives
          C[i][1 + b][a] = -M[a * Mdim + b + 1] * cutinv;
          // columns 1-2 (2D)  or 1-3 (3D)

          //Second derivatives
          if (kernel_style == RK2) C[i][1 + dim + b][a] = M[a * Mdim + b + 1 + dim] * cutsqinv;
          // columns 3-4 (2D) or 4-6 (3D)
        }
      }
    }
  }

  // communicate calculated quantities
  comm_stage = 1;
  comm_forward = comm_forward_save;
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOKernel::compute_coordination()
{
  int i, j, ii, jj, inum, jnum;
  double xtmp, ytmp, ztmp, rsq;
  double dx[3];

  double **x = atom->x;

  int *ilist, *jlist, *numneigh, **firstneigh;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Grow arrays if necessary
  if (nmax_store < atom->nmax) grow_arrays(atom->nmax);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    coordination[i] = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = lensq3(dx);

      if (rsq < cutsq) coordination[i] += 1;
    }
  }

  // communicate calculated quantities
  comm_stage = 0;
  comm_forward = 1;
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOKernel::grow_arrays(int nmax)
{
  memory->grow(coordination, nmax, "rheo:coordination");

  if (kernel_style == RK0) {
    memory->grow(C0, nmax, "rheo/kernel:C0");
  } else if (correction_order > 0) {
    memory->grow(C, nmax, ncor, Mdim, "rheo/kernel:C");
  }

  nmax_store = nmax;
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOKernel::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                         int * /*pbc*/)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    if (comm_stage == 0) {
      buf[m++] = coordination[j];
    } else {
      if (kernel_style == RK0) {
        buf[m++] = C0[j];
      } else {
        for (int a = 0; a < ncor; a++)
          for (int b = 0; b < Mdim; b++) buf[m++] = C[j][a][b];
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOKernel::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    if (comm_stage == 0) {
      coordination[i] = buf[m++];
    } else {
      if (kernel_style == RK0) {
        C0[i] = buf[m++];
      } else {
        for (int a = 0; a < ncor; a++)
          for (int b = 0; b < Mdim; b++) C[i][a][b] = buf[m++];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::memory_usage()
{
  double bytes = 0.0;
  bytes = (size_t) nmax_store * sizeof(int);

  if (kernel_style == RK0) {
    bytes += (size_t) nmax_store * sizeof(double);
  } else if (correction_order > 0) {
    bytes += (size_t) nmax_store * ncor * Mdim * sizeof(double);
  }
  return bytes;
}
