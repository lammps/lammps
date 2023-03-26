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

#include "compute_rheo_kernel.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_interface.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cblas.h>

using namespace LAMMPS_NS;
using namespace MathExtra;

enum {QUINTIC, CRK0, CRK1, CRK2};
#define DELTA 2000

Todo: convert delx dely delz to an array
Should vshift be using kernel quintic?
Move away from h notation, use cut?

/* ---------------------------------------------------------------------- */

ComputeRHEOKernel::ComputeRHEOKernel(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  C(nullptr), C0(nullptr), compute_interface(nullptr);
{
  if (narg != 3) error->all(FLERR,"Illegal compute rheo/kernel command");

  comm_forward = 1; // Always minimum for coordination
  solid_flag = 0;
  dim = domain->dimension;
}

/* ---------------------------------------------------------------------- */

ComputeRHEOKernel::~ComputeRHEOKernel()
{
  // Remove custom property if it exists
  int tmp1, tmp2, index;
  index = atom->find_custom("rheo_coordination", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index_coord, 1, 0);

  memory->destroy(C);
  memory->destroy(C0);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOKernel::init()
{
  neighbor->add_request(this, NeighConst::REQ_FULL);

  auto fixes = modify->get_fix_by_style("rheo");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use compute rheo/kernel");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);

  int icompute = modify->find_compute("rheo_interface");
  if (icompute != -1) {
    compute_interface = ((ComputeRHEOInterface *) modify->compute[icompute]);
    solid_flag = 1;
  }

  kernel_style = fix_rheo->kernel_style;
  zmin = fix_rheo->zmin_kernel;
  h = fix_rheo->h;
  hsq = h * h;
  hinv = 1.0 / h;
  hsqinv = hinv * hinv;

  if (kernel_style == FixRHEO::QUINTIC) {
    correction_order = -1;
  } else if (kernel_style == FixRHEO::CRK0) {
    correction_order = 0;
  } else if (kernel_style == FixRHEO::CRK1) {
    correction_order = 1;
  } else if (kernel_style == FixRHEO::CRK2) {
    correction_order = 2;
  }

  if (dim == 3) {
    pre_w = 0.002652582384864922 * 27.0 * ihsq * ih;
    pre_wp = pre_w * 3.0 * ih;
  } else {
    pre_w = 0.004661441847879780 * 9 * ihsq;
    pre_wp = pre_w * 3.0 * ih;
  }

  // Create coordination array if it doesn't already exist
  // Create a custom atom property so it works with compute property/atom
  // Do not create grow callback as there's no reason to copy/exchange data
  // Manually grow if nmax_old exceeded

  int tmp1, tmp2;
  int nmax = atom->nmax;
  index_coord = atom->find_custom("rheo_coordination", tmp1, tmp2);
  if (index_coord == -1) {
    index_coord = atom->add_custom("rheo_coordination", 0, 0);
    nmax_old = nmax;
  }

  comm_forward = 1;
  ncor = 0;
  Mdim = 0;
  if (kernel_type == CRK0) {
    memory->create(C0, nmax, "rheo/kernel:C0");
  } else if (kernel_type == CRK1) {
    Mdim = 1 + dim;
    ncor = 1 + dim;
    memory->create(C, nmax, ncor, Mdim, "rheo/kernel:C");
    comm_forward = ncor * Mdim;
  } else if (kernel_type == CRK2) {
    //Polynomial basis size (up to quadratic order)
    Mdim = 1 + dim + dim * (dim + 1) / 2;
    //Number of sets of correction coefficients  (1 x y xx yy)  + z zz (3D)
    ncor = 1 + 2 * dim;
    memory->create(C, nmax, ncor, Mdim, "rheo/kernel:C");
    comm_forward = ncor * Mdim;
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
  int corrections = 1;

  if (gsl_error_flag) {
    // If there were errors, check to see if it occured for this atom
    if (gsl_error_tags.find(atom->tag[i]) != gsl_error_tags.end())
      corrections = 0;
  }

  int *coordination = atom->ivector[index_coord];
  if (coordination[i] < zmin) corrections = 0;

  return corrections;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_w(int i, int j, double delx, double dely, double delz, double r)
{
  double w;

  int corrections_i = check_corrections(i);
  int corrections_j = check_corrections(j);
  int corrections = corrections_i & corrections_j;

  if (kernel_type == QUINTIC || !corrections) w = calc_w_quintic(i,j,delx,dely,delz,r);
  else if (kernel_type == CRK0) w = calc_w_crk0(i,j,delx,dely,delz,r);
  else if (kernel_type == CRK1) w = calc_w_crk1(i,j,delx,dely,delz,r);
  else if (kernel_type == CRK2) w = calc_w_crk2(i,j,delx,dely,delz,r);

  return w;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_dw(int i, int j, double delx, double dely, double delz, double r)
{
  double wp;

  int corrections_i = check_corrections(i);
  int corrections_j = check_corrections(j);

  // Calc wp and default dW's, a bit inefficient but can redo later
  wp = calc_dw_quintic(i,j,delx,dely,delz,r,dWij,dWji);
  if(kernel_type == CRK1) {
    //check if kernel correction calculated successfully. If not, revert to quintic
    if (corrections_i) calc_dw_crk1(i,j,delx,dely,delz,r,dWij);
    if (corrections_j) calc_dw_crk1(j,i,-delx,-dely,-delz,r,dWji);
  } else if(kernel_type == CRK2) {
    if (corrections_i) calc_dw_crk2(i,j,delx,dely,delz,r,dWij);
    if (corrections_j) calc_dw_crk2(j,i,-delx,-dely,-delz,r,dWji);
  }

  return wp;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_w_quintic(int i, int j, double delx, double dely, double delz, double r)
{
  double w, tmp1, tmp2, tmp3, tmp1sq, tmp2sq, tmp3sq, s;
  s = r * 3.0 * ih;

	if (s > 3.0) {
	  w = 0.0;
	}

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

double ComputeRHEOKernel::calc_dw_quintic(int i, int j, double delx, double dely, double delz, double r, double *dW1, double *dW2)
{
  double wp, tmp1, tmp2, tmp3, tmp1sq, tmp2sq, tmp3sq, s, wprinv;
  double *mass = atom->mass;
  int *type = atom->type;

  s = r * 3.0 * ih;

  if (s > 3.0) {
    wp = 0.0;
  }
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

double ComputeRHEOKernel::calc_w_crk0(int i, int j, double delx, double dely, double delz, double r)
{
  double w;

  w = calc_w_quintic(i,j,delx,dely,delz,r);

  Wij = C0[i] * w;
  Wji = C0[j] * w;

  return w;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_w_crk1(int i, int j, double delx, double dely, double delz, double r)
{
  int b;
  double w, wR, dx[3], H[Mdim];

  dx[0] = delx;
  dx[1] = dely;
  dx[2] = delz;
  w = calc_w_quintic(i,j,delx,dely,delz,r);

  if (dim == 2) {
    H[0] = 1.0;
    H[1] = dx[0] * hinv;
    H[2] = dx[1] * hinv;
  } else {
    H[0] = 1.0;
    H[1] = dx[0] * hinv;
    H[2] = dx[1] * hinv;
    H[3] = dx[2] * hinv;
  }
  Wij = 0;
  for (b = 0; b < Mdim; b++) {
    Wij += C[i][0][b] * H[b];  // C columns: 1 x y (z) xx yy (zz)
  }
  Wij *= w;

  //Now compute Wji
  H[1] *= -1;
  H[2] *= -1;
  if (dim == 3) H[3] *= -1;

  Wji = 0;
  for (b = 0; b < Mdim; b++) {
    Wji += C[j][0][b] * H[b];  // C columns: 1 x y (z) xx yy (zz)
  }
  Wji *= w;

  return w;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::calc_w_crk2(int i, int j, double delx, double dely, double delz, double r)
{
  int b;
  double w, wR, dx[3], H[Mdim];
  dx[0] = delx;
  dx[1] = dely;
  dx[2] = delz;
  w = calc_w_quintic(i,j,delx,dely,delz,r);

  if (dim == 2) {
    H[0] = 1.0;
    H[1] = dx[0] * hinv;
    H[2] = dx[1] * hinv;
    H[3] = 0.5 * dx[0] * dx[0] * hsqinv;
    H[4] = 0.5 * dx[1] * dx[1] * hsqinv;
    H[5] = dx[0] * dx[1] * hsqinv;
  } else {
    H[0] = 1.0;
    H[1] = dx[0] * hinv;
    H[2] = dx[1] * hinv;
    H[3] = dx[2] * hinv;
    H[4] = 0.5 * dx[0] * dx[0] * hsqinv;
    H[5] = 0.5 * dx[1] * dx[1] * hsqinv;
    H[6] = 0.5 * dx[2] * dx[2] * hsqinv;
    H[7] = dx[0] * dx[1] * hsqinv;
    H[8] = dx[0] * dx[2] * hsqinv;
    H[9] = dx[1] * dx[2] * hsqinv;
  }
  Wij = 0;
  for (b = 0; b < Mdim; b++) {
    Wij += C[i][0][b] * H[b];  // C columns: 1 x y (z) xx yy (zz)
  }
  Wij *= w;

  //Now compute Wji
  H[1] *= -1;
  H[2] *= -1;
  if (dim == 3) H[3] *= -1;

  Wji = 0;
  for (b = 0; b < Mdim; b++) {
    Wji += C[j][0][b] * H[b];  // C columns: 1 x y (z) xx yy (zz)
  }
  Wji *= w;

  return w;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOKernel::calc_dw_crk1(int i, int j, double delx, double dely, double delz, double r, double *dW)
{
  int a, b;
  double w, dx[3], H[Mdim];
  dx[0] = delx;
  dx[1] = dely;
  dx[2] = delz;

  w = calc_w_quintic(i,j,delx,dely,delz,r);

  //Populate correction basis
  if (dim == 2) {
    H[0] = 1.0;
    H[1] = dx[0] * hinv;
    H[2] = dx[1] * hinv;
  } else {
    H[0] = 1.0;
    H[1] = dx[0] * hinv;
    H[2] = dx[1] * hinv;
    H[3] = dx[2] * hinv;
  }

  // dWij[] = dWx dWy (dWz)
  //compute derivative operators
  for (a = 0; a < dim; a++) {
    dW[a] = 0.0;
    for (b = 0; b < Mdim; b++) {
      //First derivative kernels
      dW[a] += C[i][1 + a][b] * H[b]; // C columns: 1 x y (z)
    }
    dW[a] *= w;
  }
}


/* ---------------------------------------------------------------------- */

void ComputeRHEOKernel::calc_dw_crk2(int i, int j, double delx, double dely, double delz, double r, double *dW)
{
  int a, b;
  double w, dx[3], H[Mdim];
  dx[0] = delx;
  dx[1] = dely;
  dx[2] = delz;

  w = calc_w_quintic(i,j,delx,dely,delz,r);

  //Populate correction basis
  if (dim == 2) {
    H[0] = 1.0;
    H[1] = dx[0] * hinv;
    H[2] = dx[1] * hinv;
    H[3] = 0.5 * dx[0] * dx[0] * hsqinv;
    H[4] = 0.5 * dx[1] * dx[1] * hsqinv;
    H[5] = dx[0] * dx[1] * hsqinv;
  } else {
    H[0] = 1.0;
    H[1] = dx[0] * hinv;
    H[2] = dx[1] * hinv;
    H[3] = dx[2] * hinv;
    H[4] = 0.5 * dx[0] * dx[0] * hsqinv;
    H[5] = 0.5 * dx[1] * dx[1] * hsqinv;
    H[6] = 0.5 * dx[2] * dx[2] * hsqinv;
    H[7] = dx[0] * dx[1] * hsqinv;
    H[8] = dx[0] * dx[2] * hsqinv;
    H[9] = dx[1] * dx[2] * hsqinv;
  }

  // dWij[] = dWx dWy (dWz)
  //compute derivative operators
  for (a = 0; a < dim; a++) {
    dW[a] = 0.0;
    for (b = 0; b < Mdim; b++) {
      //First derivative kernels
      dW[a] += C[i][1 + a][b] * H[b]; // C columns: 1 x y (z) xx yy (zz)
    }
    dW[a] *= w;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOKernel::compute_peratom()
{
  gsl_error_flag = 0;
  gsl_error_tags.clear();

  int i, j, ii, jj, jnum, g, a, b, gsl_error;
  double xtmp, ytmp, ztmp, r, rsq, w, vj;
  double dx[3];
  gsl_matrix_view gM;

  // Turn off GSL error handler, revert RK to Quintic when insufficient neighbors
  gsl_set_error_handler_off();

  double **x = atom->x;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rho = atom->rho;
  int *status = atom->status;
  int *coordination = atom->ivector[index_coord];
  tagint *tag = atom->tag;

  int inum, *ilist, *jlist, *numneigh, **firstneigh;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Grow arrays if necessary
  int nmax = atom->nmax;
  if (nmax_old < nmax)
    memory->grow(coordination, nmax, "atom:rheo_coordination");

    if (kernel_type == FixRHEO::CRK0) {
      memory->grow(C0, nmax, "rheo/kernel:C0");
    } else if (correction_order > 0) {
      memory->grow(C, nmax, ncor, Mdim, "rheo/kernel:C");
    }

    nmax_old = atom->nmax;
  }

  if (kernel_type == FixRHEO::QUINTIC) {
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
        rsq = lensq(dx);

        if (rsq < hsq) {
          coordination[i] += 1;
        }
      }
    }
  } else if (kernel_type == FixRHEO::CRK0) {

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
      coordination[i] = 0;

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        dx[0] = xtmp - x[j][0];
        dx[1] = ytmp - x[j][1];
        dx[2] = ztmp - x[j][2];
        rsq = lensq(dx);

        if (rsq < hsq) {
          r = sqrt(rsq);
          w = calc_w_quintic(i,j,dx[0],dx[1],dx[2],r);
          if (!(status[j] & FixRHEO::STATUS_FLUID) && solid_flag) {
            vj = mass[type[j]] / compute_interface->correct_rho(j,i);
          } else vj = mass[type[j]] / rho[j];

          coordination[i] += 1;
          M += w * vj;
        }
      }

      // Inverse of 1x1 matrix
      if (coordination[i] >= zmin) C0[i] = 1.0 / M;
    }
  } else if (correction_order > 0) {

    // Moment matrix M and polynomial basis vector H (1d for gsl compatibility)
    double H[Mdim], M[Mdim * Mdim];

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      jlist = firstneigh[i];
      jnum = numneigh[i];
      itype = type[i];

      // Zero upper-triangle M and H (will be symmetric):
      for (a = 0; a < Mdim; a++) {
        for (b = a; b < Mdim; b++) {
          M[a * Mdim + b] = 0;
        }
      }
      coordination[i] = 0;

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        dx[0] = xtmp - x[j][0];
        dx[1] = ytmp - x[j][1];
        dx[2] = ztmp - x[j][2];

        rsq = lensq(dx);

        if (rsq < cutsq) {
          r = sqrt(rsq);
          w = calc_w_quintic(i,j,dx[0],dx[1],dx[2],r);

          if (status[j] > FixRHEO::FLUID_MAX && solid_flag)
            vj = mass[type[j]]/compute_interface->correct_rho(j,i);
          else vj = mass[type[j]]/rho[j];

          //Populate the H-vector of polynomials (2D)
          if (dim == 2) {
            H[0] = 1.0;
            H[1] = dx[0] * hinv;
            H[2] = dx[1] * hinv;
            if (kernel_type == FixRHEO::CRK2) {
              H[3] = 0.5 * dx[0] * dx[0] * hsqinv;
              H[4] = 0.5 * dx[1] * dx[1] * hsqinv;
              H[5] = dx[0] * dx[1] * hsqinv;
            }
          } else {
            H[0] = 1.0;
            H[1] = dx[0] * hinv;
            H[2] = dx[1] * hinv;
            H[3] = dx[2] * hinv;
            if (kernel_type == FixRHEO::CRK2) {
              H[4] = 0.5 * dx[0] * dx[0] * hsqinv;
              H[5] = 0.5 * dx[1] * dx[1] * hsqinv;
              H[6] = 0.5 * dx[2] * dx[2] * hsqinv;
              H[7] = dx[0] * dx[1] * hsqinv;
              H[8] = dx[0] * dx[2] * hsqinv;
              H[9] = dx[1] * dx[2] * hsqinv;
            }
          }

          coordination[i] += 1;

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
      if (coordination[i] < zmin) continue

      // Use gsl to get Minv, use Cholesky decomposition since the
      // polynomials are independent, M is symmetrix & positive-definite
      gM = gsl_matrix_view_array(M,Mdim,Mdim);
      gsl_error = gsl_linalg_cholesky_decomp(&gM.matrix);

      if (gsl_error) {
        //Revert to uncorrected SPH for this particle
        gsl_error_flag = 1;
        gsl_error_tags.insert(tag[i]);

        //check if not positive-definite
        if (gsl_error != GSL_EDOM)
          error->warn(FLERR, "Failed decomposition in rheo/kernel, gsl_error = {}", gsl_error);

        continue;
      }

      gsl_linalg_cholesky_invert(&gM.matrix);   //M is now M^-1

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
        C[i][0][a] = M[a * Mdim + 0]; // all rows of column 0
        for (b = 0; b < dim; b++) {
          //First derivatives
          C[i][1 + b][a] = -M[a * Mdim + b + 1] * hinv;
          // columns 1-2 (2D)  or 1-3 (3D)

          //Second derivatives
          if (kernel_type == FixRHEO::CRK2)
            C[i][1 + dim + b][a] = M[a * Mdim + b + 1 + dim] * hsqinv;
            // columns 3-4 (2D) or 4-6 (3D)
        }
      }
    }
  }

  // communicate calculated quantities
  comm->forward_comm_compute(this);
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOKernel::pack_forward_comm(int n, int *list, double *buf,
                                        int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,k,m,a,b;
  coordination = atom->ivector[index_coord];
  m = 0;
  if (correction_order > 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      for (a = 0; a < ncor; a++) {
        for (b = 0; b < Mdim; b++) {
          buf[m++] = C[j][a][b];
        }
      }
      buf[m++] = coordination[j];
    }
  } else if (kernel_type == FixRHEO::CRK0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = C0[j];
      buf[m++] = coordination[j];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = coordination[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOKernel::unpack_forward_comm(int n, int first, double *buf)
{
  int i, k, m, last,a,b;
  coordination = atom->ivector[index_coord];
  m = 0;
  last = first + n;
  if (correction_order > 0) {
    for (i = first; i < last; i++) {
      for (a = 0; a < ncor; a++) {
        for (b = 0; b < Mdim; b++) {
          C[i][a][b] = buf[m++];
        }
      }
      coordination[i] = buf[m++];
    }
  } else if (kernel_type == FixRHEO::CRK0) {
    for (i = first; i < last; i++) {
      C0[i] = buf[m++];
      coordination[i] = buf[m++];
    }
  } else {
    for (i = first; i < last; i++) {
      coordination[i] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOKernel::memory_usage()
{
  double bytes = 0.0;
  bytes = (size_t) atom->nmax * sizeof(int);

  if (kernel_type == FixRHEO::CRK0) {
    bytes += (size_t) atom->nmax * sizeof(double);
  } else if (correction_order > 0) {
    bytes += (size_t) atom->nmax * ncor * Mdim * sizeof(double);
  }
  return bytes;
}
