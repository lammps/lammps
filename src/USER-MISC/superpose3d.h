/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.  (Some of the code in this file is also
   available using a more premissive license.  See below for details.)

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing author: Andrew Jewett (Scripps Research)
   Availability: https://github.com/jewettaij/superpose3d_cpp  (MIT license)
------------------------------------------------------------------------- */

/// @file     superpose3d.hpp
/// @brief    Calculate the optimal rotation, translation and scale needed to
///           optimally fit two different point clouds containing n points.
/// @author   Andrew Jewett
/// @license  MIT

#ifndef _SUPERPOSE3D_H
#define _SUPERPOSE3D_H

#include "math_eigen.h"   //functions to calculate eigenvalues and eigenvectors

namespace superpose3d {

using namespace math_eigen;

// -----------------------------------------------------------
// ------------------------ INTERFACE ------------------------
// -----------------------------------------------------------

/// @brief  Superpose3d is a class with only one important member function
///         Superpose().  It is useful for calculating the optimal
///         superposition (rotations, translations, and scale transformations)
///         between two point clouds of the same size.
template<typename Scalar,
         typename ConstArrayOfCoords,
         typename ConstArray=Scalar const*>
class Superpose3D {
private:
  size_t N;              //number of points in the point clouds
  Scalar *aWeights;      //weights applied to points when computing RMSD
  Jacobi<double, double*, double**> eigen_calc; // calculates eigenvectors
  Scalar **aaXf_shifted; //preallocated space for fixed point cloud (Nx3 array)
  Scalar **aaXm_shifted; //preallocated space for mobile point cloud (Nx3 array)

public:
  // The following data members store the rotation, translation and scale
  // after optimal superposition
  Scalar **R;  //!< store optimal rotation here (this is a 3x3 array).
  Scalar T[3]; //!< store optimal translation here
  Scalar c;  //!< store optimal scale (typically 1 unless requested by the user)
  Scalar q[4]; //!< quaternion corresponding to the rotation stored in R.
               //   The first entry of q is cos(θ/2).  The remaining 3 entries
               //   of q are the axis of rotation (with length sin(θ/2)).
               // (Note: This is not the same as "p" from Diamond's 1988 paper.)

  Superpose3D(size_t N = 0);  //!< N=number of points in both point clouds

  Superpose3D(size_t N,      //!< N = number of points in both point clouds
              ConstArray aWeights); //!< weight per point for computing RMSD

  ~Superpose3D();

  /// @brief specify he number of points in both point clouds
  void SetNumPoints(size_t N);
  /// @brief return the number of points in both point clouds
  size_t GetNumPoints() { return N; }
  /// @brief specify the weight applied to each point when computing RMSD
  void SetWeights(ConstArray aWeights);

  /// @brief Use rigid-body transformations (rotations, translations, and
  ///         optionally scale transformations) to superimpose two point clouds.
  ///
  /// @details
  /// This function takes two lists of xyz coordinates (of the same length) and 
  /// attempts to superimpose them using rotations, translations, and 
  /// (optionally) scale transformations.  These transformations are applied to
  /// to the coordinates in the "aaXm_orig" array (the "mobile" point cloud)
  /// in order to minimize the root-mean-squared-distance (RMSD) between the
  /// corresponding points in each cloud, where RMSD is defined as:
  ///
  /// @verbatim
  /// sqrt((Σ_n w[n]*Σ_i |X[n][i] - (Σ_j c*R[i][j]*x[n][j]+T[i])|^2)/(Σ_n w[n]))
  /// @endverbatim
  ///
  /// In this formula, the "X_i" and "x_i" are coordinates of the ith fixed and
  /// mobile point clouds (represented by "aaXf" and "aaXm" in the code below)
  /// and "w_i" are optional weights (represented by "aWeights" in the code).
  /// This function implements a more general variant of the method from:
  /// @verbatim
  /// R. Diamond, (1988) "A Note on the Rotational Superposition Problem", 
  /// Acta Cryst. A44, pp. 211-216
  /// @endverbatim
  ///
  /// @note:
  /// This code has been augmented with a new feature.  The version in the
  /// original paper only considers rotation and translation and does not allow
  /// coordinates of either cloud to be rescaled (multiplied by a scalar).
  /// To enable the ability to rescale the coordinates, set allow_rescale=true.
  /// (By default, this feature is disabled.)
  ///
  /// @returns
  /// The RMSD between the 2 pointclouds after optimal rotation, translation
  /// (and scaling if requested) was applied to the "mobile" point cloud.
  /// After this function is called, the optimal rotation, translation,
  /// and scale (if requested) will be stored in the "R", "T", and "c"
  /// public data members.
  Scalar Superpose(ConstArrayOfCoords aaXf, //!< coords for the "frozen" object
                   ConstArrayOfCoords aaXm, //!< coords for the "mobile" object
                   bool allow_rescale=false //!< rescale mobile object? (c≠1?)
                   );

  // C++ boilerplate: copy and move constructor, swap, and assignment operator
  Superpose3D(const Superpose3D<Scalar,ConstArrayOfCoords,ConstArray>& source);
  Superpose3D(Superpose3D<Scalar,ConstArrayOfCoords,ConstArray>&& other);
  void swap(Superpose3D<Scalar,ConstArrayOfCoords,ConstArray> &other);
  Superpose3D<Scalar,ConstArrayOfCoords,ConstArray>& operator = (Superpose3D<Scalar,ConstArrayOfCoords,ConstArray> source);

private:

  // memory management:
  void Alloc(size_t N);
  void Init();
  void Dealloc();

}; // class Superpose3D





// -------------- IMPLEMENTATION --------------


template<typename Scalar>
static inline Scalar SQR(Scalar x) {return x*x;}

template<typename Scalar, typename ConstArrayOfCoords, typename ConstArray>
Scalar Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>::
Superpose(ConstArrayOfCoords aaXf, // coords for the "frozen" object
          ConstArrayOfCoords aaXm, // coords for the "mobile" object
          bool allow_rescale)      // rescale mobile object? (c!=1?)
{
  assert(aaXf && aaXm);
  assert(aaXf_shifted && aaXm_shifted);
  assert(aWeights);
  assert(R && T);

  // Find the center of mass of each object:
  Scalar aCenter_f[3] = {0.0, 0.0, 0.0};
  Scalar aCenter_m[3] = {0.0, 0.0, 0.0};
  Scalar sum_weights = 0.0;
  for (size_t n=0; n < N; n++) {
    Scalar weight = aWeights[n];
    for (int d=0; d < 3; d++) {
      aCenter_f[d] += aaXf[n][d]*weight;
      aCenter_m[d] += aaXm[n][d]*weight;
    }
    sum_weights += weight;
  }
  assert(sum_weights != 0.0);
  for (int d=0; d < 3; d++) {
    aCenter_f[d] /= sum_weights;
    aCenter_m[d] /= sum_weights;
  }

  //Subtract the centers-of-mass from the original coordinates for each object
  for (size_t n=0; n < N; n++) {
    for (int d=0; d < 3; d++) {
      // shift the coordinates so that the new center of mass is at the origin
      aaXf_shifted[n][d] = aaXf[n][d] - aCenter_f[d];
      aaXm_shifted[n][d] = aaXm[n][d] - aCenter_m[d];
    }
  }

  // Calculate the "M" array from the Diamond paper (equation 16)
  Scalar M[3][3];
  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      M[i][j] = 0.0;

  for (size_t n=0; n < N; n++) {
    Scalar weight = aWeights[n];
    for (int i=0; i < 3; i++) {
      for (int j=0; j < 3; j++) {
        M[i][j] += weight * aaXm_shifted[n][i] * aaXf_shifted[n][j];
      }
    }
  }

  // Calculate Q (equation 17)
  Scalar traceM = 0.0;
  for (int i=0; i < 3; i++)
    traceM += M[i][i];
  Scalar Q[3][3];
  for (int i=0; i < 3; i++) {
    for (int j=0; j < 3; j++) {
      Q[i][j] = M[i][j] + M[j][i];
      if (i==j)
        Q[i][j] -= 2.0 * traceM;
    }
  }

  // Calculate V (equation 18)
  Scalar V[3];
  V[0] = M[1][2] - M[2][1];
  V[1] = M[2][0] - M[0][2];
  V[2] = M[0][1] - M[1][0];

  // Calculate "P" (equation 22)
  // First we must allocate space for the P matrix.  It's not safe to declare:
  // Scalar P[4][4];
  // ...because most matrix solvers expect arrays in pointer-to-pointer format.
  // (a different format).  Below I create a fixed size matrix P in this format.
  Scalar _P[4*4]; // Contiguous 1D array for storing contents of the 2D P array
  Scalar *P[4];   // This version of P has has ** (pointer-to-pointer) format.
  for (int i=0; i < 4; i++) // We must make sure that
    P[i] = &(_P[4*i]);      // P[i] points to the appropriate location in memory

  // Now fill the P array
  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      P[i][j] = Q[i][j];
  P[0][3] = V[0];
  P[3][0] = V[0];
  P[1][3] = V[1];
  P[3][1] = V[1];
  P[2][3] = V[2];
  P[3][2] = V[2];
  P[3][3] = 0.0;

  // The vector "p" contains the optimal rotation (backwards quaternion format)
  Scalar p[4] = {0.0, 0.0, 0.0, 1.0};  // default value
  Scalar pPp = 0.0;                    // = p^T * P * p  (zero by default)
  Scalar rmsd = 0.0;                   // default value

  bool singular = N<2;  // (it doesn't make sense to rotate a single point)

  if (! singular) {
    // Calculate the principal eigenvalue and eigenvector of matrix P.
    // Store the principal eigenvector in "p"
    // The vector "p" will contain the optimal rotation (in quaternion format)

    Scalar Evl[4];    // Store the eigenvalues of P here.
    Scalar *Evc[4];   // Store the eigevectors here. This version has ** format.
    Scalar _Evc[4*4]; // Contiguous 1D array for storing contents of "Evc" array
    for (int i=0; i < 4; i++) // We must make sure that
      Evc[i] = &(_Evc[4*i]);  // Evc[i] points to the correct location in memory

    eigen_calc.Diagonalize(P, Evl, Evc);

    // Note: The eigenvalues are sorted in decreasing order by default.
    pPp = Evl[0];  // = the maximum eigenvalue of P
    for (int i=0; i < 4; i++)
      p[i] = Evc[0][i]; //copy eigenvector corresponding to this eigenvalue to p
  } //if (! singular)

  // Now normalize p
  Scalar pnorm = 0.0;
  for (int i=0; i < 4; i++)
    pnorm += p[i]*p[i];
  pnorm = sqrt(pnorm);
  for (int i=0; i < 4; i++)
    p[i] /= pnorm;

  // Finally, calculate the rotation matrix corresponding to "p"
  // (convert a quaternion into a 3x3 rotation matrix)

  R[0][0] =  (p[0]*p[0])-(p[1]*p[1])-(p[2]*p[2])+(p[3]*p[3]);
  R[1][1] = -(p[0]*p[0])+(p[1]*p[1])-(p[2]*p[2])+(p[3]*p[3]);
  R[2][2] = -(p[0]*p[0])-(p[1]*p[1])+(p[2]*p[2])+(p[3]*p[3]);
  R[0][1] = 2*(p[0]*p[1] - p[2]*p[3]);
  R[1][0] = 2*(p[0]*p[1] + p[2]*p[3]);
  R[1][2] = 2*(p[1]*p[2] - p[0]*p[3]);
  R[2][1] = 2*(p[1]*p[2] + p[0]*p[3]);
  R[0][2] = 2*(p[0]*p[2] + p[1]*p[3]);
  R[2][0] = 2*(p[0]*p[2] - p[1]*p[3]);

  q[0] = p[3];  // Note: The "p" variable is not a quaternion in the
  q[1] = p[0];  //       conventional sense because its elements
  q[2] = p[1];  //       are in the wrong order.  I correct for that here.
  q[3] = p[2];  //       "q" is the quaternion correspond to rotation R.

  // Optional: Decide the scale factor, c
  c = 1.0;   // by default, don't rescale the coordinates

  if ((allow_rescale) && (! singular)) {
    Scalar Waxaixai = 0.0;
    Scalar WaxaiXai = 0.0;
    for (size_t a=0; a < N; a++) {
      Scalar weight = aWeights[a];
      for (int i=0; i < 3; i++) {
        Waxaixai += weight * aaXm_shifted[a][i] * aaXm_shifted[a][i];
        WaxaiXai += weight * aaXm_shifted[a][i] * aaXf_shifted[a][i];
      }
    }
    c = (WaxaiXai + pPp) / Waxaixai;

  } // if (allow_rescale)

  // Finally compute the RMSD between the two coordinate sets:
  // First compute E0 from equation 24 of the paper
  Scalar E0 = 0.0;
  for (size_t n=0; n < N; n++) {
    Scalar weight = aWeights[n];
    for (int d=0; d < 3; d++)
      // (remember to include the scale factor "c" that we inserted)
      E0 += weight * (SQR(aaXf_shifted[n][d] - c*aaXm_shifted[n][d]));
  }
  Scalar sum_sqr_dist = E0 - c*2.0*pPp;
  if (sum_sqr_dist < 0.0) //(edge case due to rounding error)
    sum_sqr_dist = 0.0;

  if (! singular)
    rmsd = sqrt(sum_sqr_dist/sum_weights);

  // Lastly, calculate the translational offset.
  // If c!=1, this is slightly more complicated than it seems.  Recall that:
  //RMSD=sqrt((Sum_i  w_i * |X_i - Sum_j(c*R_ij*x_j + T_i))|^2) / (Sum_j w_j))
  //    =sqrt((Sum_i  w_i * |X_i - x_i')|^2) / (Sum_j w_j))
  //  where
  // x_i' = Sum_j(c*R_ij*x_j) + T_i
  //      = Xcm_i + c*R_ij*(x_j - xcm_j)
  //  and Xcm and xcm = center_of_mass for the frozen and mobile point clouds
  //
  // Hence:
  //  T_i = Xcm_i - Sum_j c*R_ij*xcm_j
  // In the code, Xcm_i is represented by "aCenter_f[i]"
  //          and xcm_j is represented by "aCenter_m[j]"

  for (int i=0; i < 3; i++) {
    T[i] = aCenter_f[i];
    for (int j=0; j < 3; j++) {
      T[i] -= c*R[i][j]*aCenter_m[j];
    }
  }

  return rmsd;

} //Superpose3D::Superpose(aaXf, aaXm, allow_rescale)


template<typename Scalar, typename ConstArrayOfCoords, typename ConstArray>
void Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>::
SetNumPoints(size_t N) {
  Dealloc();
  Alloc(N);
}

template<typename Scalar, typename ConstArrayOfCoords, typename ConstArray>
void Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>::
SetWeights(ConstArray aWeights) {
  for (size_t i = 0; i < N; i++)
    this->aWeights[i] = aWeights[i];
}

template<typename Scalar, typename ConstArrayOfCoords, typename ConstArray>
Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>::Superpose3D(size_t N)
  :eigen_calc(4)
{
  Init();
  Alloc(N);
}

template<typename Scalar, typename ConstArrayOfCoords, typename ConstArray>
Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>::
Superpose3D(size_t N, ConstArray aWeights)
  :eigen_calc(4)
{
  Init();
  Alloc(N);
  SetWeights(aWeights);
}

template<typename Scalar, typename ConstArrayOfCoords, typename ConstArray>
Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>::~Superpose3D() {
  Dealloc();
}

template<typename Scalar, typename ConstArrayOfCoords, typename ConstArray>
void Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>::
Init() {
  R = nullptr;
  aWeights = nullptr;
  aaXf_shifted = nullptr;
  aaXm_shifted = nullptr;
}

// memory management:

template<typename Scalar, typename ConstArrayOfCoords, typename ConstArray>
void Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>::
Alloc(size_t N) {
  this->N = N;
  aWeights = new Scalar [N];
  for (size_t i = 0; i < N; i++)
    aWeights[i] = 1.0 / N;
  Alloc2D(3, 3, &R);
  Alloc2D(N, 3, &aaXf_shifted);
  Alloc2D(N, 3, &aaXm_shifted);
}

template<typename Scalar, typename ConstArrayOfCoords, typename ConstArray>
void Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>::
Dealloc() {
  if (R)
    Dealloc2D(&R);
  if (aWeights)
    delete [] aWeights;
  if (aaXf_shifted)
    Dealloc2D(&aaXf_shifted);
  if (aaXm_shifted)
    Dealloc2D(&aaXm_shifted);
}

// memory management: copy and move constructor, swap, and assignment operator:

template<typename Scalar, typename ConstArrayOfCoords, typename ConstArray>
Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>::
Superpose3D(const Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>& source)
  :eigen_calc(4)
{
  Init();
  Alloc(source.N);
  assert(N == source.N);
  for (int i = 0; i < N; i++) {
    std::copy(source.aaXf_shifted[i],
              source.aaXf_shifted[i] + 3,
              aaXf_shifted[i]);
    std::copy(source.aaXm_shifted[i],
              source.aaXm_shifted[i] + 3,
              aaXm_shifted[i]);
  }
}

template<typename Scalar, typename ConstArrayOfCoords, typename ConstArray>
void Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>::
swap(Superpose3D<Scalar, ConstArrayOfCoords, ConstArray> &other) {
  std::swap(N, other.N);
  std::swap(R, other.R);
  std::swap(aaXf_shifted, other.aaXf_shifted);
  std::swap(aaXm_shifted, other.aaXm_shifted);
}

// Move constructor (C++11)
template<typename Scalar, typename ConstArrayOfCoords, typename ConstArray>
Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>::
Superpose3D(Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>&& other) {
  Init();
  swap(*this, other);
}

// Using the "copy-swap" idiom for the assignment operator
template<typename Scalar, typename ConstArrayOfCoords, typename ConstArray>
Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>&
Superpose3D<Scalar, ConstArrayOfCoords, ConstArray>::
operator = (Superpose3D<Scalar, ConstArrayOfCoords, ConstArray> source) {
  this->swap(source);
  return *this;
}


} //namespace superposed3d



#endif //#ifndef _SUPERPOSE3D_H
