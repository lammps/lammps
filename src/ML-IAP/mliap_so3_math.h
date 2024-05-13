/*
 * MLIAP_SO3_math.h
 *
 *  Created on: Apr 12, 2021
 *      Author: macstein
 */

#ifndef LMP_MLIAP_SO3_MATH_H
#define LMP_MLIAP_SO3_MATH_H

#include "math_eigen_impl.h"
#include <cmath>

namespace SO3Math {

int jacobin(int n, double const *const *mat, double *eval, double **evec);
int invert_matrix(int n, double *A, double *Ainv);
int LUPdecompose(int n, double dtol, double *A, int *P);
void LUPSolve(int n, double *A, double *B, int *P);

}    // namespace SO3Math

using namespace MathEigen;

typedef Jacobi<double, double *, double **, double const *const *> Jacobi_v2;
inline int SO3Math::jacobin(int n, double const *const *mat, double *eval, double **evec)
{
  int *midx = new int[n];
  double **M = new double *[n];
  double **mat_cpy = new double *[n];

  for (int i = 0; i < n; i++) {
    mat_cpy[i] = new double[n];
    for (int j = 0; j < n; j++) mat_cpy[i][j] = mat[i][j];
    M[i] = &(mat_cpy[i][0]);
  }

  Jacobi_v2 ecalcn(n, M, midx);
  int ierror = ecalcn.Diagonalize(mat, eval, evec, Jacobi_v2::SORT_DECREASING_EVALS);

  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) std::swap(evec[i][j], evec[j][i]);
    delete[] mat_cpy[i];
  }
  delete[] mat_cpy;
  delete[] M;
  delete[] midx;

  return ierror;
}

inline int SO3Math::invert_matrix(int n, double *A, double *Ainv)
{

  int i, j;
  double dtol = 1.e-30;

  int *P;
  double *b, *Atemp;

  P = new int[n];
  b = new double[n];
  Atemp = new double[n * n];

  for (i = 0; i < n * n; i++) Atemp[i] = A[i];

  int rv = 0;
  if (LUPdecompose(n, dtol, Atemp, P) == 0) {

    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) b[j] = 0.0;

      b[i] = 1.0;
      LUPSolve(n, Atemp, b, P);

      for (j = 0; j < n; j++) Ainv[j * n + i] = b[j];
    }
  } else {
    rv = 1;
  }

  delete[] P;
  delete[] b;
  delete[] Atemp;

  return rv;
}

inline int SO3Math::LUPdecompose(int n, double dtol, double *A, int *P)
{
  int i, j, k, maxi;
  double maxA, Atemp;
  double *normi;

  maxi = 0;
  normi = new double[n];

  for (i = 0; i < n; i++) {
    maxA = 0.0;
    for (j = 0; j < n; j++) {
      Atemp = fabs(A[i * n + j]);
      if (Atemp > maxA) maxA = Atemp;
    }
    if (maxA < dtol) {
      delete[] normi;
      return 1;
    }

    normi[i] = 1.0 / maxA;
  }

  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++)
      for (k = 0; k < i; k++) A[i * n + j] -= A[i * n + k] * A[k * n + j];

    maxA = 0.0;
    for (i = j; i < n; i++) {
      for (k = 0; k < j; k++) A[i * n + j] -= A[i * n + k] * A[k * n + j];

      Atemp = fabs(A[i * n + j]) * normi[i];
      if (Atemp >= maxA) {
        maxA = Atemp;
        maxi = i;
      }
    }
    if (maxi != j) {
      if ((j == (n - 2)) && (A[j * n + j + 1] == 0.0))
        maxi = j;
      else {
        for (k = 0; k < n; k++) {
          Atemp = A[j * n + k];
          A[j * n + k] = A[maxi * n + k];
          A[maxi * n + k] = Atemp;
        }
        normi[maxi] = normi[j];
      }
    }

    P[j] = maxi;
    if (A[j * n + j] == 0.0) A[j * n + j] = dtol;

    if (j != (n - 1)) {
      Atemp = 1.0 / A[j * n + j];
      for (i = (j + 1); i < n; i++) A[i * n + j] *= Atemp;
    }
  }

  delete[] normi;
  return 0;
}

inline void SO3Math::LUPSolve(int n, double *A, double *B, int *P)
{
  int i, j;
  double dtemp;

  for (i = 0; i < n; i++) {

    dtemp = B[P[i]];
    B[P[i]] = B[i];
    for (j = (i - 1); j >= 0; j--) dtemp -= A[i * n + j] * B[j];

    B[i] = dtemp;
  }

  for (i = (n - 1); i >= 0; i--) {
    for (j = (i + 1); j < n; j++) B[i] -= A[i * n + j] * B[j];

    B[i] /= A[i * n + i];
  }
}

#endif /* LMP_MLIAP_SO3_MATH_H_ */
