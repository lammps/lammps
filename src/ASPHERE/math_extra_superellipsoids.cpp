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
   Contributing authors: Jacopo Bilotto (EPFL), Jibril B. Coulibaly
------------------------------------------------------------------------- */

#include "math_extra_superellipsoids.h"

#include "math_const.h"
#include "math_extra.h"

#include <cmath>

namespace {
constexpr double TIKHONOV_SCALE = 1e-14;

constexpr int ITERMAX_NR = 100;
constexpr double TOL_NR_RES = 1e-10 * 1e-10;
constexpr double TOL_NR_POS = 1e-6 * 1e-6;

constexpr int ITERMAX_LS = 10;
constexpr double PARAMETER_LS = 1e-4;
constexpr double CUTBACK_LS = 0.5;

constexpr double TOL_OVERLAP = 1e-8;
constexpr unsigned int ITERMAX_OVERLAP = 20;
constexpr double MINSLOPE_OVERLAP = 1e-12;

constexpr double REGULARIZATION_EPSILON = 1e-12;
constexpr double MAX_B_FAST = 1e30;
}    // namespace

// clang-format off

namespace MathExtraSuperellipsoids {

/* ----------------------------------------------------------------------
   curvature of superellipsoid
   source https://en.wikipedia.org/wiki/Mean_curvature
------------------------------------------------------------------------- */

double mean_curvature_superellipsoid(const double *shape, const double *block, const int flag, const double R[3][3], const double *surf_global_point, const double *xc)
{
  // this code computes the mean curvature on the superellipsoid surface
  // for the given global point
  double hess[3][3], grad[3], normal[3];
  double xlocal[3], tmp_v[3];
  MathExtra::sub3(surf_global_point, xc, tmp_v); // here tmp_v is the vector from center to surface point
  MathExtra::transpose_matvec(R, tmp_v, xlocal);
  (void) shape_and_derivatives_local(xlocal, shape, block, flag, grad, hess); // computation of curvature is independent of local or global frame
  MathExtra::normalize3(grad, normal);
  MathExtra::matvec(hess, normal, tmp_v); // here tmp_v is intermediate product
  double F_mag = sqrt(MathExtra::dot3(grad, grad));
  double curvature = fabs(MathExtra::dot3(normal, tmp_v) - (hess[0][0] + hess[1][1] + hess[2][2])) / (2.0 * F_mag);
  return curvature;
}

double gaussian_curvature_superellipsoid(const double *shape, const double *block, const int flag, const double R[3][3], const double *surf_global_point, const double *xc)
{
  // this code computes the gaussian curvature coefficient
  // for the given global point
  double hess[3][3], grad[3], normal[3];
  double xlocal[3], tmp_v[3];
  MathExtra::sub3(surf_global_point, xc, tmp_v); // here tmp_v is the vector from center to surface point
  MathExtra::transpose_matvec(R, tmp_v, xlocal);
  (void) shape_and_derivatives_local(xlocal, shape, block, flag, grad, hess); // computation of curvature is independent of local or global frame
  MathExtra::normalize3(grad, normal);

  double temp[3];
  MathExtra::matvec(hess, normal, temp);
  double F_mag = sqrt(MathExtra::dot3(grad, grad));

  double fx = grad[0];
  double fy = grad[1];
  double fz = grad[2];

  double fxx = hess[0][0];
  double fxy = hess[0][1];
  double fxz = hess[0][2];

  double fyy = hess[1][1];
  double fyz = hess[1][2];

  double fzz = hess[2][2];

  double mat[4][4] = {
    {fxx, fxy, fxz, fx},
    {fxy, fyy, fyz, fy},
    {fxz, fyz, fzz, fz},
    {fx,  fy,  fz, 0.0}
  };

  double K = -det4_M44_zero(mat) / (F_mag*F_mag*F_mag*F_mag);
  double curvature =  sqrt(fabs(K));
  return curvature;
}


/* ----------------------------------------------------------------------
   express global (system level) to local (particle level) coordinates
------------------------------------------------------------------------- */

void global2local_vector(const double *v, const double *quat, double *local_v)
{
    double qc[4];
    MathExtra::qconjugate(const_cast<double*>(quat), qc);
    MathExtra::quatrotvec(qc, const_cast<double*>(v), local_v);
}

/* ----------------------------------------------------------------------
   Possible regularization for the shape functions
   Instead of F(x,y,z) = 0 we use (F(x,y,z)+1)^(1/n1) -1 = G(x,y,z) = 0
   We also scale G by the average radius to have better "midway" points
------------------------------------------------------------------------- */

void apply_regularization_shape_function(double n1, const double avg_radius, double *value, double *grad, double hess[3][3])
{
  // value is F - 1
  double base = fmax(*value + 1.0, REGULARIZATION_EPSILON);
  const double inv_F = 1.0 / base;
  const double inv_n1 = 1.0 / n1;

  // P = base^(1/n)
  const double F_pow_inv_n1 = pow(base, inv_n1);

  // Scale for Gradient: S1 = R * (1/n) * base^(1/n - 1)
  const double scale_grad = avg_radius * inv_n1 * F_pow_inv_n1 * inv_F;

  // Scale for Hessian addition: S2 = S1 * (1/n - 1) * base^-1
  const double scale_hess_add = scale_grad * (inv_n1 - 1.0) * inv_F;

  // H_new = scale_grad * H_old + scale_hess_add * (grad_old x grad_old^T)
  for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
          double grad_outer_prod = grad[i] * grad[j];
          hess[i][j] = (hess[i][j] * scale_grad) + (scale_hess_add * grad_outer_prod);
      }
  }

  // grad_new = scale_grad * grad_old
  for (int i = 0; i < 3; i++) {
      grad[i] *= scale_grad;
  }

  // G = R * (base^(1/n) - 1)
  *value = avg_radius * (F_pow_inv_n1 - 1.0);
};

/* ----------------------------------------------------------------------
   shape function computations for superellipsoids
------------------------------------------------------------------------- */

double shape_and_derivatives_local(const double* xlocal, const double* shape, const double* block, const int flag, double* grad, double hess[3][3])
{
  double shapefunc = 0.0;
  // TODO: Not sure how to make flag values more clear
  // Cannot forward declare the enum AtomVecEllipsoid::BlockType
  // Could use scoped (enum class) but no implicit conversion:
  //    must pass type `LAMMPS_NS::AtomVecEllipsoid::BlockType` instead of int,
  //    and/or static_cast the enum class to int, which is similar to current
  // Could define the enum in a dedicated header
  //    seems overkill just for one enum
  // I think the comment below making reference to the BlockType should be enough
  // Feel free to change to a better design
  switch (flag) { // LAMMPS_NS::AtomVecEllipsoid::BlockType
    case 0: {
      shapefunc = shape_and_derivatives_local_ellipsoid(xlocal, shape, grad, hess);
      break;
    }
    case 1: {
      shapefunc = shape_and_derivatives_local_n1equaln2(xlocal, shape, block[0], grad, hess);
      break;
    }
    case 2: {
      shapefunc = shape_and_derivatives_local_superquad(xlocal, shape, block, grad, hess);
      break;
    }
  }

  return shapefunc;
}

/* ----------------------------------------------------------------------
   General case for n1 != n2 > 2
------------------------------------------------------------------------- */

double shape_and_derivatives_local_superquad(const double* xlocal, const double* shape, const double* block, double* grad, double hess[3][3])
{
  double a_inv = 1.0 / shape[0];
  double b_inv = 1.0 / shape[1];
  double c_inv = 1.0 / shape[2];
  double x_a = fabs(xlocal[0] * a_inv);
  double y_b = fabs(xlocal[1] * b_inv);
  double z_c = fabs(xlocal[2] * c_inv);
  double n1 = block[0];
  double n2 = block[1];
  double x_a_pow_n2_m2 = pow(x_a, n2 - 2.0);
  double x_a_pow_n2_m1 = x_a_pow_n2_m2 * x_a;
  double y_b_pow_n2_m2 = pow(y_b, n2 - 2.0);
  double y_b_pow_n2_m1 = y_b_pow_n2_m2 * y_b;

  double nu = (x_a_pow_n2_m1 * x_a) + (y_b_pow_n2_m1 * y_b);
  double nu_pow_n1_n2_m2 = pow(nu, n1/n2 - 2.0);
  double nu_pow_n1_n2_m1 = nu_pow_n1_n2_m2 * nu;

  double z_c_pow_n1_m2 = pow(z_c, n1 -2.0);
  double z_c_pow_n1_m1 = z_c_pow_n1_m2 * z_c;

  // Equation (14)
  double signx = xlocal[0] > 0.0 ? 1.0 : -1.0;
  double signy = xlocal[1] > 0.0 ? 1.0 : -1.0;
  double signz = xlocal[2] > 0.0 ? 1.0 : -1.0;
  grad[0] = n1 * a_inv * x_a_pow_n2_m1 * nu_pow_n1_n2_m1 * signx;
  grad[1] = n1 * b_inv * y_b_pow_n2_m1 * nu_pow_n1_n2_m1 * signy;
  grad[2] = n1 * c_inv * z_c_pow_n1_m1 * signz;

  // Equation (15)
  double signxy = signx * signy;
  hess[0][0] = a_inv * a_inv * (n1 * (n2 - 1.0) * x_a_pow_n2_m2 * nu_pow_n1_n2_m1 +
                                (n1 - n2) * n1 * (x_a_pow_n2_m1 * x_a_pow_n2_m1) * nu_pow_n1_n2_m2);
  hess[1][1] = b_inv * b_inv * (n1 * (n2 - 1.0) * y_b_pow_n2_m2 * nu_pow_n1_n2_m1 +
                                (n1 - n2) * n1 * (y_b_pow_n2_m1 * y_b_pow_n2_m1) * nu_pow_n1_n2_m2);
  hess[0][1] = hess[1][0] = a_inv * b_inv * (n1 - n2) * n1 * x_a_pow_n2_m1 * y_b_pow_n2_m1 * nu_pow_n1_n2_m2 * signxy;
  hess[2][2] = c_inv * c_inv * n1 * (n1 - 1.0) * z_c_pow_n1_m2;
  hess[0][2] = hess[2][0] = hess[1][2] = hess[2][1] = 0.0;

  return (nu_pow_n1_n2_m1 * nu) + (z_c_pow_n1_m1 * z_c) - 1.0;
}

/* ----------------------------------------------------------------------
   Special case for n2 = n2 = n > 2
------------------------------------------------------------------------- */

double shape_and_derivatives_local_n1equaln2(const double* xlocal, const double* shape, const double n, double* grad, double hess[3][3])
{
  double a_inv = 1.0 / shape[0];
  double b_inv = 1.0 / shape[1];
  double c_inv = 1.0 / shape[2];
  double x_a = fabs(xlocal[0] * a_inv);
  double y_b = fabs(xlocal[1] * b_inv);
  double z_c = fabs(xlocal[2] * c_inv);
  double x_a_pow_n_m2 = pow(x_a, n - 2.0);
  double x_a_pow_n_m1 = x_a_pow_n_m2 * x_a;
  double y_b_pow_n_m2 = pow(y_b, n - 2.0);
  double y_b_pow_n_m1 = y_b_pow_n_m2 * y_b;
  double z_c_pow_n_m2 = pow(z_c, n - 2.0);
  double z_c_pow_n_m1 = z_c_pow_n_m2 * z_c;

  // Equation (14)
  double signx = xlocal[0] > 0.0 ? 1.0 : -1.0;
  double signy = xlocal[1] > 0.0 ? 1.0 : -1.0;
  double signz = xlocal[2] > 0.0 ? 1.0 : -1.0;
  grad[0] = n * a_inv * x_a_pow_n_m1 * signx;
  grad[1] = n * b_inv * y_b_pow_n_m1 * signy;
  grad[2] = n * c_inv * z_c_pow_n_m1 * signz;

  // Equation (15)
  hess[0][0] = a_inv * a_inv * n * (n - 1.0) * x_a_pow_n_m2;
  hess[1][1] = b_inv * b_inv * n * (n - 1.0) * y_b_pow_n_m2;
  hess[2][2] = c_inv * c_inv * n * (n - 1.0) * z_c_pow_n_m2;
  hess[0][1] = hess[1][0] = hess[0][2] = hess[2][0] = hess[1][2] = hess[2][1] = 0.0;

  return (x_a_pow_n_m1 * x_a) + (y_b_pow_n_m1 * y_b) + (z_c_pow_n_m1 * z_c) - 1.0;
}

/* ----------------------------------------------------------------------
   Special case for n1 = n2 = 2
------------------------------------------------------------------------- */

double shape_and_derivatives_local_ellipsoid(const double* xlocal, const double* shape, double* grad, double hess[3][3])
{
  double a = 2.0 / (shape[0] * shape[0]);
  double b = 2.0 / (shape[1] * shape[1]);
  double c = 2.0 / (shape[2] * shape[2]);

  // Equation (14) simplified for n1 = n2 = 2
  grad[0] = a * xlocal[0];
  grad[1] = b * xlocal[1];
  grad[2] = c * xlocal[2];

  // Equation (15)
  hess[0][0] = a;
  hess[1][1] = b;
  hess[2][2] = c;
  hess[0][1] = hess[1][0] = hess[0][2] = hess[2][0] = hess[1][2] = hess[2][1] = 0.0;

  return 0.5 * (grad[0]*xlocal[0] + grad[1]*xlocal[1] + grad[2]*xlocal[2]) - 1.0;
}

/* ---------------------------------------------------------------------- */

double shape_and_derivatives_global(const double* xc, const double R[3][3],
    const double* shape, const double* block, const int flag,
    const double* X0, double* grad, double hess[3][3],
    const int formulation, const double avg_radius)
{
  double xlocal[3], tmp_v[3], tmp_m[3][3];
  MathExtra::sub3(X0, xc, tmp_v);
  MathExtra::transpose_matvec(R, tmp_v, xlocal);
  double shapefunc = shape_and_derivatives_local(xlocal, shape, block, flag, tmp_v, hess);
  if (formulation == FORMULATION_GEOMETRIC) {
     apply_regularization_shape_function(block[0], avg_radius, &shapefunc, tmp_v, hess);
  }
  MathExtra::matvec(R, tmp_v, grad);
  MathExtra::times3_transpose(hess, R, tmp_m);
  MathExtra::times3(R, tmp_m, hess);

  return shapefunc;
}

// double compute_residual(const double shapefunci, const double* gradi_global, const double shapefuncj, const double* gradj_global, const double mu2, double* residual) {
//   // Equation (23)
//   MathExtra::scaleadd3(mu2, gradj_global, gradi_global, residual);
//   residual[3] = shapefunci - shapefuncj;
//   // Normalize residual Equation (23)
//   // shape functions and gradients dimensions are not homogeneous
//   // Gradient equality F1' + mu2 * F2' evaluated relative to magnitude of gradient ||F1'|| = ||mu2 * F2'||
//   // Shape function equality F1 - F2 evaluated relative to magnitude of shape function + 1
//   //    the shift f = polynomial - 1 is not necessary and cancels out in F1 - F2
//   // Last component homogeneous to shape function
//   return MathExtra::lensq3(residual) / MathExtra::lensq3(gradi_global) +
//          residual[3] * residual[3] / ((shapefunci + 1) * (shapefunci + 1));
// }

/* ---------------------------------------------------------------------- */

double compute_residual(const double shapefunci, const double* gradi_global,
                        const double shapefuncj, const double* gradj_global,
                        const double mu2, double* residual,
                        const int formulation, const double radius_scale)
{
  // Equation (23): Spatial residual (Gradient match)
  MathExtra::scaleadd3(mu2, gradj_global, gradi_global, residual);
  residual[3] = shapefunci - shapefuncj;

  // --- Spatial Normalization ---
  // Algebraic: Gradients are ~1/R. Dividing by lensq3 normalizes this.
  // We take average of gradient for polydisperse case
  // Geometric: Gradients are unit vectors. lensq3 is 1.0. This works for both.

  double gradi_global_mag = 0.5 * (MathExtra::lensq3(gradi_global) + MathExtra::lensq3(gradj_global));
  double spatial_norm = MathExtra::lensq3(residual) / gradi_global_mag;

  // --- Scalar Normalization ---
  double scalar_denom;

  if (formulation == FORMULATION_GEOMETRIC) {
      // GEOMETRIC: G is a distance (Length).
      scalar_denom = radius_scale;
  } else {
      // ALGEBRAIC: F is dimensionless (approx 0 at surface).
      scalar_denom = shapefunci + 1.0;
  }

  // Prevent division by zero in weird edge cases (e.g. very negative shape function)
  if (fabs(scalar_denom) < 1e-12) scalar_denom = 1.0;

  return spatial_norm + (residual[3] * residual[3]) / (scalar_denom * scalar_denom);
}

/* ---------------------------------------------------------------------- */

void compute_jacobian(const double* gradi_global, const double hessi_global[3][3],
                     const double* gradj_global, const double hessj_global[3][3], const double mu2, double* jacobian)
{
  // Jacobian (derivative of residual)
  for (int row = 0 ; row < 3 ; row++) {
    for (int col = 0 ; col < 3 ; col++) {
      jacobian[row*4 + col] = hessi_global[row][col] + mu2 * hessj_global[row][col];
    }
    jacobian[row*4 + 3] = gradj_global[row];
  }
  for (int col = 0 ; col < 3 ; col++) {
    jacobian[3*4 + col] = gradi_global[col] - gradj_global[col];
  }
  jacobian[15] = 0.0;
}

/* ---------------------------------------------------------------------- */

double compute_residual_and_jacobian(const double* xci, const double Ri[3][3], const double* shapei, const double* blocki, const int flagi,
                                     const double* xcj, const double Rj[3][3], const double* shapej, const double* blockj, const int flagj,
                                     const double* X, double* shapefunc, double* residual, double* jacobian,
                                     const int formulation, const double avg_radius_i, const double avg_radius_j)
{
  double gradi[3], hessi[3][3], gradj[3], hessj[3][3];
  shapefunc[0] = shape_and_derivatives_global(xci, Ri, shapei, blocki, flagi, X, gradi, hessi, formulation, avg_radius_i);
  shapefunc[1] = shape_and_derivatives_global(xcj, Rj, shapej, blockj, flagj, X, gradj, hessj, formulation, avg_radius_j);
  compute_jacobian(gradi, hessi, gradj, hessj, X[3], jacobian);
  return compute_residual(shapefunc[0], gradi, shapefunc[1], gradj, X[3], residual, formulation, (avg_radius_i + avg_radius_j) * 0.5);
}

/* ---------------------------------------------------------------------- */

int determine_contact_point(const double* xci, const double Ri[3][3], const double* shapei, const double* blocki, const int flagi,
                            const double* xcj, const double Rj[3][3], const double* shapej, const double* blockj, const int flagj,
                            double* X0, double* nij, int formulation)
{
  double norm, norm_old, shapefunc[2], residual[4], jacobian[16];
  double lsq = MathExtra::distsq3(xci, xcj);
  bool converged(false);

  double rhs_old[3];

  // avg radii for regularization if GEOMETRIC formulation
  double avg_radius_i = 1;
  double avg_radius_j = 1;
  double max_step = sqrt(lsq) * 0.2;
  if (formulation == FORMULATION_GEOMETRIC) {
    avg_radius_i = (shapei[0] + shapei[1] + shapei[2]) * LAMMPS_NS::MathConst::THIRD;
    avg_radius_j = (shapej[0] + shapej[1] + shapej[2]) * LAMMPS_NS::MathConst::THIRD;
  }

  norm = compute_residual_and_jacobian(xci, Ri, shapei, blocki, flagi, xcj, Rj, shapej, blockj, flagj, X0, shapefunc, residual, jacobian, formulation, avg_radius_i, avg_radius_j);
  // testing for convergence before attempting Newton's method.
  // the initial guess is the old X0, so with temporal coherence, it might still pass tolerance if deformation is slow!
  if (norm < TOL_NR_RES) {

    //  must compute the normal vector nij before returning since the Newton loop normally handles this upon convergence.
    double xilocal[3], tmp_v[3], gradi[3];

    // Transform global X0 to local frame of particle I
    MathExtra::sub3(X0, xci, tmp_v);
    MathExtra::transpose_matvec(Ri, tmp_v, xilocal);

    // Compute local gradient
    // Algebraic gradient is fine for direction even if we used Geometric for solving
    // TODO: might use a simpler function to simply compute the gradient, to
    // avoid computing quantities already computed in compute_residual_and_jacobian
    if (flagi <= 1)
      (void) shape_and_gradient_local_n1equaln2_surfacesearch(xilocal, shapei, blocki[0], tmp_v);
    else
      (void) shape_and_gradient_local_superquad_surfacesearch(xilocal, shapei, blocki, tmp_v);

    // Rotate gradient back to global frame to get normal
    MathExtra::matvec(Ri, tmp_v, gradi);
    MathExtra::normalize3(gradi, nij);

    // Return status
    if (shapefunc[0] > 0.0 || shapefunc[1] > 0.0)
      return 1; // Converged, but no contact (separated)

    return 0; // Converged and Contacting
  }


  for (int iter = 0 ; iter < ITERMAX_NR ; iter++) {
    norm_old = norm;

    double rhs[4];
    bool gauss_elim_solved = false;
    double A_fast[16];
    double b_fast[4];

    for(int i = 0; i < 16; ++i) {
        A_fast[i] = jacobian[i];
    }

    b_fast[0] = -residual[0]; b_fast[1] = -residual[1];
    b_fast[2] = -residual[2]; b_fast[3] = -residual[3];

    // Try Fast Solver
    gauss_elim_solved = MathExtraSuperellipsoids::solve_4x4_robust_unrolled(A_fast, b_fast);

    // check for divergence or numerical issues in the fast solver
    // and fall back to regularized solver if necessary
    bool fail0 = !std::isfinite(b_fast[0]) || (fabs(b_fast[0]) > MAX_B_FAST);
    bool fail1 = !std::isfinite(b_fast[1]) || (fabs(b_fast[1]) > MAX_B_FAST);
    bool fail2 = !std::isfinite(b_fast[2]) || (fabs(b_fast[2]) > MAX_B_FAST);
    bool fail3 = !std::isfinite(b_fast[3]) || (fabs(b_fast[3]) > MAX_B_FAST);
    if (fail0 || fail1 || fail2 || fail3) {
        gauss_elim_solved = false;
    }

    rhs[0] = b_fast[0]; rhs[1] = b_fast[1];
    rhs[2] = b_fast[2]; rhs[3] = b_fast[3];

    if (!gauss_elim_solved) {
      // restore matrix
      for(int i = 0; i < 16; ++i) {
        A_fast[i] = jacobian[i];
      }

      b_fast[0] = -residual[0]; b_fast[1] = -residual[1];
      b_fast[2] = -residual[2]; b_fast[3] = -residual[3];
       // enforce a minimum regularization to avoid zero pivots in edge cases (flat on flat)
      double trace = jacobian[0] + jacobian[5] + jacobian[10];
      double diag_weight = fmax(TIKHONOV_SCALE * trace, TIKHONOV_SCALE);
      A_fast[0]  += diag_weight;
      A_fast[5]  += diag_weight;
      A_fast[10] += diag_weight;

      if (MathExtraSuperellipsoids::solve_4x4_robust_unrolled(A_fast, b_fast)) {
          rhs[0] = b_fast[0]; rhs[1] = b_fast[1];
          rhs[2] = b_fast[2]; rhs[3] = b_fast[3];
          gauss_elim_solved = true;
      }
    }

    MathExtra::copy3(rhs, rhs_old);

    // Backtracking line search
    double X_line[4];
    int iter_ls;
    double a = 1.0;

    // Limit the max step size to avoid jumping too far
    // normalize residual vector if step was limited
    double spatial_residual_norm = sqrt(rhs[0]*rhs[0] + rhs[1]*rhs[1] + rhs[2]*rhs[2]);

    if (spatial_residual_norm > max_step) {
        double scale = max_step / spatial_residual_norm;
        rhs[0] *= scale;
        rhs[1] *= scale;
        rhs[2] *= scale;
    }

    for (iter_ls = 0 ; iter_ls < ITERMAX_LS ; iter_ls++) {
      X_line[0] = X0[0] + a * rhs[0];
      X_line[1] = X0[1] + a * rhs[1];
      X_line[2] = X0[2] + a * rhs[2];
      X_line[3] = X0[3] + a * rhs[3];

      // Line search iterates not selected for the next Newton iteration
      // do not need to compute the expensive Jacobian, only the residual.
      // We want to avoid calling `compute_residual_and_jacobian()` for each
      // line search iterate.
      // However, many intermediate variables that are costly to compute
      // are shared by the local gradient and local hessian calculations.
      // We want to avoid calling `compute_residual()` followed by `compute_jacobian()`
      // for the iterates that satisfy the descent condition.
      // To do so, we duplicate `compute_residual_and_jacobian()`, but only
      // build the global hessians if the descent condition is satisfied and
      // the iterate will be used in the next Newton step.
      // This leads to some code duplication, and still computes
      // the local hessians even when they are not necessary.
      // This seems to be an acceptable in-between of performance and clean code.
      // As most of the cost in the Hessian is in the 2 matrix products to
      // Compute the global matrix from the local one

      // One alternative would be to store the intermediate variables from
      // the local gradient calculation when calling `shape_and_gradient_local()`,
      // and re-use them during the local hessian calculation (function that
      // calculates only the Hessian from these intermediate values would need
      // to be implemented).
      // This seems a bit clunky just to save the few multiplications of the
      // local hessian calculation, that is why I did not do it. I am open to
      // other ideas and solutions.
      // Even then, we would have some code duplication with `compute_residual_and_jacobian()`
      // So maybe I am overthinking this...

      double xilocal[3], gradi[3], hessi[3][3], xjlocal[3], gradj[3], hessj[3][3], tmp_v[3];

      MathExtra::sub3(X_line, xci, tmp_v);
      MathExtra::transpose_matvec(Ri, tmp_v, xilocal);
      shapefunc[0] = shape_and_derivatives_local(xilocal, shapei, blocki, flagi, tmp_v, hessi);
      if (formulation == FORMULATION_GEOMETRIC) {
          apply_regularization_shape_function(blocki[0], avg_radius_i, &shapefunc[0], tmp_v, hessi);
      }
      MathExtra::matvec(Ri, tmp_v, gradi);

      MathExtra::sub3(X_line, xcj, tmp_v);
      MathExtra::transpose_matvec(Rj, tmp_v, xjlocal);
      shapefunc[1] = shape_and_derivatives_local(xjlocal, shapej, blockj, flagj, tmp_v, hessj);
      if (formulation == FORMULATION_GEOMETRIC) {
          apply_regularization_shape_function(blockj[0], avg_radius_j, &shapefunc[1], tmp_v, hessj);
      }
      MathExtra::matvec(Rj, tmp_v, gradj);

      norm = compute_residual(shapefunc[0], gradi, shapefunc[1], gradj, X_line[3], residual, formulation, (avg_radius_i + avg_radius_j)/2.0);

      if ((norm <= TOL_NR_RES) &&
          (MathExtra::lensq3(rhs) * a * a <= TOL_NR_POS * lsq)) {
        converged = true;

        MathExtra::normalize3(gradi, nij);
        break;
      } else if (norm > norm_old - PARAMETER_LS * a * norm_old) { // Armijo - Goldstein condition not met
        // Tested after convergence check because tiny values of norm and norm_old < TOL_NR
        // Can still fail the Armijo - Goldstein condition`
        a *= CUTBACK_LS;
      } else {
        // Only compute the jacobian if there is another Newton iteration to come
        double tmp_m[3][3];
        MathExtra::times3_transpose(hessi, Ri, tmp_m);
        MathExtra::times3(Ri, tmp_m, hessi);
        MathExtra::times3_transpose(hessj, Rj, tmp_m);
        MathExtra::times3(Rj, tmp_m, hessj);
        compute_jacobian(gradi, hessi, gradj, hessj, X_line[3], jacobian);
        break;
      }
    }
    // Take full step if no descent at the end of line search
    // Try to escape bad region
    if (iter_ls == ITERMAX_LS) {
      X0[0] += rhs[0];
      X0[1] += rhs[1];
      X0[2] += rhs[2];
      X0[3] += rhs[3];
      norm = compute_residual_and_jacobian(xci, Ri, shapei, blocki, flagi, xcj, Rj, shapej, blockj, flagj, X0, shapefunc, residual, jacobian, formulation, avg_radius_i, avg_radius_j);
      if (norm < TOL_NR_RES) {
        converged = true;
        // must re-compute the normal 'nij' for this final point
        double xilocal[3], tmp_v[3], gradi[3], hess_dummy[3][3];
        MathExtra::sub3(X0, xci, tmp_v);
        MathExtra::transpose_matvec(Ri, tmp_v, xilocal);

        // We only need the gradient for the normal
        shape_and_derivatives_local(xilocal, shapei, blocki, flagi, tmp_v, hess_dummy);
        if (formulation == FORMULATION_GEOMETRIC) {
            // If you use regularization, apply it here too for consistency
            apply_regularization_shape_function(blocki[0], avg_radius_i, &shapefunc[0], tmp_v, hess_dummy);
        }
        MathExtra::matvec(Ri, tmp_v, gradi);
        MathExtra::normalize3(gradi, nij);
      }

    } else {
      X0[0] = X_line[0];
      X0[1] = X_line[1];
      X0[2] = X_line[2];
      X0[3] = X_line[3];
    }

    if (converged) break;
  }

  if (!converged) {
    if (shapefunc[0] > 0.0 || shapefunc[1] > 0.0) return 1;
    return 2;
  } // not failing if not converged but shapefuncs positive (i.e., no contact)
              // might be risky to assume no contact if not converged, NR might have gone to a far away point
              // but no guarantee there is no contact
  if (shapefunc[0] > 0.0 || shapefunc[1] > 0.0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   Functions to compute shape function and gradient only when called for
     newton method to avoid computing hessian when not needed and having
     smoother landscape for the line search
   General case for n1 != n2 > 2
------------------------------------------------------------------------- */

double shape_and_gradient_local_superquad_surfacesearch(const double* xlocal, const double* shape, const double* block, double* grad)
{
  double a_inv = 1.0 / shape[0];
  double b_inv = 1.0 / shape[1];
  double c_inv = 1.0 / shape[2];
  double x_a = fabs(xlocal[0] * a_inv);
  double y_b = fabs(xlocal[1] * b_inv);
  double z_c = fabs(xlocal[2] * c_inv);
  double n1 = block[0];
  double n2 = block[1];
  double x_a_pow_n2_m2 = pow(x_a, n2 - 2.0);
  double x_a_pow_n2_m1 = x_a_pow_n2_m2 * x_a;
  double y_b_pow_n2_m2 = pow(y_b, n2 - 2.0);
  double y_b_pow_n2_m1 = y_b_pow_n2_m2 * y_b;

  double nu = (x_a_pow_n2_m1 * x_a) + (y_b_pow_n2_m1 * y_b);
  double nu_pow_n1_n2_m2 = pow(nu, n1/n2 - 2.0);
  double nu_pow_n1_n2_m1 = nu_pow_n1_n2_m2 * nu;

  double z_c_pow_n1_m2 = pow(z_c, n1 -2.0);
  double z_c_pow_n1_m1 = z_c_pow_n1_m2 * z_c;

  // Equation (14)
  double signx = xlocal[0] > 0.0 ? 1.0 : -1.0;
  double signy = xlocal[1] > 0.0 ? 1.0 : -1.0;
  double signz = xlocal[2] > 0.0 ? 1.0 : -1.0;
  grad[0] = n1 * a_inv * x_a_pow_n2_m1 * nu_pow_n1_n2_m1 * signx;
  grad[1] = n1 * b_inv * y_b_pow_n2_m1 * nu_pow_n1_n2_m1 * signy;
  grad[2] = n1 * c_inv * z_c_pow_n1_m1 * signz;

  double F = (nu_pow_n1_n2_m1 * nu) + (z_c_pow_n1_m1 * z_c);

  double scale_factor = pow(F, 1.0/n1 -1.0) / n1;

  grad[0] *= scale_factor;
  grad[1] *= scale_factor;
  grad[2] *= scale_factor;

  return pow(F, 1.0/n1) - 1.0;
}

/* ----------------------------------------------------------------------
   Special case for n2 = n2 = n > 2
------------------------------------------------------------------------- */

double shape_and_gradient_local_n1equaln2_surfacesearch(const double* xlocal, const double* shape, const double n, double* grad)
{
  double a_inv = 1.0 / shape[0];
  double b_inv = 1.0 / shape[1];
  double c_inv = 1.0 / shape[2];
  double x_a = fabs(xlocal[0] * a_inv);
  double y_b = fabs(xlocal[1] * b_inv);
  double z_c = fabs(xlocal[2] * c_inv);
  double x_a_pow_n_m2 = pow(x_a, n - 2.0);
  double x_a_pow_n_m1 = x_a_pow_n_m2 * x_a;
  double y_b_pow_n_m2 = pow(y_b, n - 2.0);
  double y_b_pow_n_m1 = y_b_pow_n_m2 * y_b;
  double z_c_pow_n_m2 = pow(z_c, n - 2.0);
  double z_c_pow_n_m1 = z_c_pow_n_m2 * z_c;

  // Equation (14)
  double signx = xlocal[0] > 0.0 ? 1.0 : -1.0;
  double signy = xlocal[1] > 0.0 ? 1.0 : -1.0;
  double signz = xlocal[2] > 0.0 ? 1.0 : -1.0;
  grad[0] = n * a_inv * x_a_pow_n_m1 * signx;
  grad[1] = n * b_inv * y_b_pow_n_m1 * signy;
  grad[2] = n * c_inv * z_c_pow_n_m1 * signz;

  double F = (x_a_pow_n_m1 * x_a) + (y_b_pow_n_m1 * y_b) + (z_c_pow_n_m1 * z_c);
  double scale_factor = pow(F, 1.0/n -1.0) / n;

  grad[0] *= scale_factor;
  grad[1] *= scale_factor;
  grad[2] *= scale_factor;

  return pow(F, 1.0/n) - 1.0;
}

/* ----------------------------------------------------------------------
   Newton Rapson method to find the overlap distance from the contact point given the normal
------------------------------------------------------------------------- */

double compute_overlap_distance(
  const double* shape, const double* block, const double Rot[3][3], const int flag,
  const double* global_point, const double* global_normal,
  const double* center)
{
  double local_point[3], local_normal[3];
  double del[3];
  double overlap;
  MathExtra::sub3(global_point, center, del);  // bring origin to 0.0
  MathExtra::transpose_matvec(Rot, del, local_point);
  MathExtra::transpose_matvec(Rot, global_normal, local_normal);

  double local_grad[3];

  // elliposid analytical solution, might need to double check the math
  // there is an easy way to find this by parametrizing the straight line as
  // X0 + t * n anf then substituting in the ellipsoid equation  for x, y, z
  // this results in a quadratic equation and we take the positive solution since
  // we are taking the outward facing normal for each grain

  if (flag == 0) {

    double a_inv2 = 1.0 / (shape[0] * shape[0]);
    double b_inv2 = 1.0 / (shape[1] * shape[1]);
    double c_inv2 = 1.0 / (shape[2] * shape[2]);

    // Coefficients for At^2 + Bt + C = 0
    double A = (local_normal[0] * local_normal[0] * a_inv2) +
               (local_normal[1] * local_normal[1] * b_inv2) +
               (local_normal[2] * local_normal[2] * c_inv2);

    double B = 2.0 * ( (local_point[0] * local_normal[0] * a_inv2) +
                     (local_point[1] * local_normal[1] * b_inv2) +
                     (local_point[2] * local_normal[2] * c_inv2) );

    double C = (local_point[0] * local_point[0] * a_inv2) +
               (local_point[1] * local_point[1] * b_inv2) +
               (local_point[2] * local_point[2] * c_inv2) - 1.0;

    // Discriminant
    double delta = B*B - 4.0*A*C;

    // Clamp delta to zero just in case numerical noise makes it negative
    if (delta < 0.0) delta = 0.0;
    overlap = (-B + sqrt(delta)) / (2.0 * A);
  } else {
    // --- Superquadric Case (Newton-Raphson on Distance Estimator) ---

    overlap = 0.0; // Distance along the normal
    double current_p[3];
    double val;
    for (unsigned int iter = 0; iter < ITERMAX_OVERLAP; iter++) {
      // Update current search position: P = Start + t * Normal
      current_p[0] = local_point[0] + overlap * local_normal[0];
      current_p[1] = local_point[1] + overlap * local_normal[1];
      current_p[2] = local_point[2] + overlap * local_normal[2];

      // Calculate Distance Estimator value and Gradient
      if (flag == 1) {
        val = shape_and_gradient_local_n1equaln2_surfacesearch(current_p, shape, block[0], local_grad);
      } else {
        val = shape_and_gradient_local_superquad_surfacesearch(current_p, shape, block, local_grad);
      }

      // Convergence Check
      if (fabs(val) < TOL_OVERLAP) break;

      // Newton Step
      double slope = local_grad[0] * local_normal[0] +
                     local_grad[1] * local_normal[1] +
                     local_grad[2] * local_normal[2];

      // Safety check to prevent divide-by-zero if ray grazes surface
      if (fabs(slope) < MINSLOPE_OVERLAP) break;

      overlap -= val / slope;
    }
  }
  return overlap;
}

} // namespace MathExtraSuperellipsoids
