/* -*- c++ -*- ----------------------------------------------------------
   Jacobi eigenvalue algorithm for small symmetric matrices.

   Computes all eigenvalues and eigenvectors of a real symmetric matrix
   using cyclic Jacobi rotations. Suitable for MSEVB Hamiltonians
   (typically < 20x20). No external dependencies (no LAPACK).

   All arrays are flat (row-major) double* for Kokkos compatibility.

   Usage:
     int n = 4;
     double H[16];       // symmetric matrix, row-major (destroyed on output)
     double evals[4];    // eigenvalues in ascending order
     double evecs[16];   // eigenvectors as columns, row-major

     int info = JacobiEigen::solve(n, H, evals, evecs);
     // info == 0: success
     // info >  0: did not converge in max_sweeps

   Reference: Golub & Van Loan, "Matrix Computations", Section 8.4
------------------------------------------------------------------------- */

#ifndef LMP_JACOBI_EIGEN_H
#define LMP_JACOBI_EIGEN_H

#include <cmath>

namespace LAMMPS_NS {
namespace JacobiEigen {

  // Access element (i,j) of row-major n×n matrix
  inline double &elem(double *A, int n, int i, int j)
  {
    return A[i * n + j];
  }
  inline double elem(const double *A, int n, int i, int j)
  {
    return A[i * n + j];
  }

  /* ---------------------------------------------------------------------- */

  // Solve eigenvalue problem for n×n symmetric matrix H (in-place destroyed).
  // Returns eigenvalues in ascending order in evals[n].
  // Returns eigenvectors as columns of evecs[n×n] (row-major).
  // Returns 0 on success, 1 if not converged.

  inline int solve(int n, double *H, double *evals, double *evecs, int max_sweeps = 100)
  {
    const double tol = 1.0e-15;

    // Initialize evecs to identity
    for (int i = 0; i < n * n; i++) evecs[i] = 0.0;
    for (int i = 0; i < n; i++) elem(evecs, n, i, i) = 1.0;

    // Cyclic Jacobi sweeps
    for (int sweep = 0; sweep < max_sweeps; sweep++) {

      // Check convergence: sum of squares of off-diagonal elements
      double off_diag = 0.0;
      for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) off_diag += elem(H, n, i, j) * elem(H, n, i, j);

      if (off_diag < tol) break;

      if (sweep == max_sweeps - 1) {
        // Copy diagonal to evals before returning failure
        for (int i = 0; i < n; i++) evals[i] = elem(H, n, i, i);
        return 1;
      }

      // Sweep over all upper-triangular pairs
      for (int p = 0; p < n - 1; p++) {
        for (int q = p + 1; q < n; q++) {
          double Hpq = elem(H, n, p, q);
          if (std::fabs(Hpq) < tol) continue;

          double Hpp = elem(H, n, p, p);
          double Hqq = elem(H, n, q, q);
          double theta = 0.5 * (Hqq - Hpp) / Hpq;

          // Compute t = sign(theta) / (|theta| + sqrt(theta^2 + 1))
          double t;
          if (std::fabs(theta) > 1.0e15) {
            t = 0.5 / theta;
          } else {
            t = 1.0 / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
            if (theta < 0.0) t = -t;
          }

          double c = 1.0 / std::sqrt(t * t + 1.0);
          double s = t * c;
          double tau = s / (1.0 + c);

          // Update H matrix
          elem(H, n, p, p) -= t * Hpq;
          elem(H, n, q, q) += t * Hpq;
          elem(H, n, p, q) = 0.0;
          elem(H, n, q, p) = 0.0;

          for (int r = 0; r < n; r++) {
            if (r == p || r == q) continue;
            double Hrp = elem(H, n, r, p);
            double Hrq = elem(H, n, r, q);
            elem(H, n, r, p) = Hrp - s * (Hrq + tau * Hrp);
            elem(H, n, p, r) = elem(H, n, r, p);
            elem(H, n, r, q) = Hrq + s * (Hrp - tau * Hrq);
            elem(H, n, q, r) = elem(H, n, r, q);
          }

          // Accumulate eigenvector rotations
          for (int r = 0; r < n; r++) {
            double Vrp = elem(evecs, n, r, p);
            double Vrq = elem(evecs, n, r, q);
            elem(evecs, n, r, p) = Vrp - s * (Vrq + tau * Vrp);
            elem(evecs, n, r, q) = Vrq + s * (Vrp - tau * Vrq);
          }
        }
      }
    }

    // Extract eigenvalues from diagonal
    for (int i = 0; i < n; i++) evals[i] = elem(H, n, i, i);

    // Sort eigenvalues in ascending order, rearrange eigenvectors to match
    for (int i = 0; i < n - 1; i++) {
      int k = i;
      for (int j = i + 1; j < n; j++)
        if (evals[j] < evals[k]) k = j;
      if (k != i) {
        double tmp = evals[i];
        evals[i] = evals[k];
        evals[k] = tmp;
        for (int r = 0; r < n; r++) {
          tmp = elem(evecs, n, r, i);
          elem(evecs, n, r, i) = elem(evecs, n, r, k);
          elem(evecs, n, r, k) = tmp;
        }
      }
    }

    return 0;
  }

}    // namespace JacobiEigen
}    // namespace LAMMPS_NS

#endif
