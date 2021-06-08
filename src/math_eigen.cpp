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

/* ----------------------------------------------------------------------
   Contributing author: Andrew Jewett (Scripps Research)
------------------------------------------------------------------------- */

#include "math_eigen.h"
#include "math_eigen_impl.h"

#include <array>
#include <utility>
#include <vector>

using std::vector;
using std::array;
using namespace MathEigen;

// Special case: 3x3 matrices

typedef Jacobi<double,double*,double(*)[3],double const(*)[3]>  Jacobi_v1;
typedef Jacobi<double,double*,double**,double const*const*> Jacobi_v2;

int MathEigen::jacobi3(double const mat[3][3], double *eval, double evec[3][3])
{
  // make copy of const matrix

  double mat_cpy[3][3] = { {mat[0][0], mat[0][1], mat[0][2]},
                           {mat[1][0], mat[1][1], mat[1][2]},
                           {mat[2][0], mat[2][1], mat[2][2]} };
  double *M[3] = { &(mat_cpy[0][0]),  &(mat_cpy[1][0]), &(mat_cpy[2][0]) };
  int midx[3];

  // create instance of generic Jacobi class and get eigenvalues and -vectors

  Jacobi_v1 ecalc3(3, M, midx);
  int ierror = ecalc3.Diagonalize(mat, eval, evec, Jacobi_v1::SORT_DECREASING_EVALS);

  // transpose the evec matrix

  for (int i=0; i<3; i++)
    for (int j=i+1; j<3; j++)
      std::swap(evec[i][j], evec[j][i]);

  return ierror;
}

int MathEigen::jacobi3(double const* const* mat, double *eval, double **evec)
{
  // make copy of const matrix

  double mat_cpy[3][3] = { {mat[0][0], mat[0][1], mat[0][2]},
                           {mat[1][0], mat[1][1], mat[1][2]},
                           {mat[2][0], mat[2][1], mat[2][2]} };
  double *M[3] = { &(mat_cpy[0][0]),  &(mat_cpy[1][0]), &(mat_cpy[2][0]) };
  int midx[3];

  // create instance of generic Jacobi class and get eigenvalues and -vectors

  Jacobi_v2 ecalc3(3, M, midx);
  int ierror = ecalc3.Diagonalize(mat, eval, evec, Jacobi_v2::SORT_DECREASING_EVALS);

  // transpose the evec matrix

  for (int i=0; i<3; i++)
    for (int j=i+1; j<3; j++)
      std::swap(evec[i][j], evec[j][i]);

  return ierror;
}
