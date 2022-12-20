/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MATH_EIGEN_H
#define LMP_MATH_EIGEN_H

namespace MathEigen {

/** A specialized function which finds the eigenvalues and eigenvectors
 *  of a 3x3 matrix (in double ** format).
 *
 * \param  mat   the 3x3 matrix you wish to diagonalize
 * \param  eval  store the eigenvalues here
 * \param  evec  store the eigenvectors here...
 * \return       0 if eigenvalue calculation converged, 1 if it failed */

int jacobi3(double const *const *mat, double *eval, double **evec);

/** \overload */

int jacobi3(double const mat[3][3], double *eval, double evec[3][3]);

}    // namespace MathEigen

#endif    //#ifndef LMP_MATH_EIGEN_H
