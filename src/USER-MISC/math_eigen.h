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
   Contributing authors: Yuya Kurebayashi (Tohoku University, Lanczos algorithm)
                         Andrew Jewett (Scripps Research, Jacobi algorith)
------------------------------------------------------------------------- */

#ifndef _MATH_EIGEN_H
#define _MATH_EIGEN_H

/// @file  This file contains a library of functions and classes which can
///        efficiently perform eigendecomposition for an extremely broad
///        range of matrix types: both real and complex, dense and sparse.
///        Matrices need not be of type "double **", for example.
///        In principle, almost any type of C++ container can be used.
///        Some general C++11 compatible functions for allocating matrices and
///        calculating norms of real and complex vectors are also provided.
/// @note
///        The "Jacobi" and "PEigenDense" classes are used for calculating
///        eigenvalues and eigenvectors of conventional dense square matrices.
/// @note
///        The "LambdaLanczos" class can calculate eigenalues and eigenvectors
///        of more general types of matrices, especially large, sparse matrices.
///        It uses C++ lambda expressions to simplify and generalize the way
///        matrices can be represented.  This allows it to be applied to
///        nearly any kind of sparse (or dense) matrix representation.
/// @note
///        The source code for Jacobi and LambdaLanczos is also available at:
///        https://github.com/jewettaij/jacobi_pd   (CC0-1.0 license)
///        https://github.com/mrcdr/lambda-lanczos  (MIT license)

#include <cassert>
#include <numeric>
#include <complex>
#include <limits>
#include <cmath>
#include <vector>
#include <random>
#include <functional>

namespace MathEigen {

// --- Memory allocation for matrices ---

/// @brief Allocate an arbitrary 2-dimensional array.  (Uses row-major order.)
/// @note  This function was intended for relatively small matrices (eg 4x4).
///        For large arrays, please use the 2d create() function from "memory.h"
template<typename Entry>
void Alloc2D(size_t nrows,          //!< size of the array (number of rows)
             size_t ncols,          //!< size of the array (number of columns)
             Entry ***paaX);        //!< pointer to a 2D C-style array

/// @brief Deallocate arrays that were created using Alloc2D().
template<typename Entry>
void Dealloc2D(Entry ***paaX);     //!< pointer to 2D multidimensional array

// --- Complex numbers ---

/// @brief  "realTypeMap" struct is used to the define "real_t<T>" type mapper
///         which returns the C++ type corresponding to the real component of T.
/// @details  Consider a function ("l2_norm()") that calculates the
///         (Euclidian) length of a vector of numbers (either real or complex):
/// @code
///    template <typename T>  real_t<T>  l2_norm(const std::vector<T>& vec);
/// @endcode
/// The l2_norm is always real by definition.
/// (See https://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm)
/// The return type of this function ("real_t<T>") indicates that
/// it returns a real number, even if the entries (of type T)
/// are complex numbers.  In other words, by default, real_t<T> returns T.
/// However real_t<std::complex<T>> returns T (not std::complex<T>).
/// We define "real_t<T>" below using C++ template specializations:
 
template <typename T>
struct realTypeMap {
  typedef T type;
};
template <typename T>
struct realTypeMap<std::complex<T>> {
  typedef T type;
};
template <typename T>
using real_t = typename realTypeMap<T>::type;


// --- Operations on vectors (of real and complex numbers) ---

/// @brief  Calculate the inner product of two vectors.
///         (For vectors of complex numbers, std::conj() is used.)
template <typename T>
T inner_prod(const std::vector<T>& v1, const std::vector<T>& v2);

/// @brief  Compute the sum of the absolute values of the entries in v
/// @returns a real number (of type real_t<T>).
template <typename T>
real_t<T> l1_norm(const std::vector<T>& v);

/// @brief  Calculate the l2_norm (Euclidian length) of vector v.
/// @returns a real number (of type real_t<T>).
template <typename T>
real_t<T> l2_norm(const std::vector<T>& v);

/// @brief  Multiply a vector (v) by a scalar (c).
template <typename T1, typename T2>
void scalar_mul(T1 c, std::vector<T2>& v);

/// @brief  Divide vector "v" in-place by it's length (l2_norm(v)).
template <typename T>
void normalize(std::vector<T>& v);


// ---- Eigendecomposition of small dense symmetric matrices ----

/// @class Jacobi
/// @brief Calculate the eigenvalues and eigevectors of a symmetric matrix
///        using the Jacobi eigenvalue algorithm.  Code for the Jacobi class
///        (along with tests and benchmarks) is available free of copyright at
///        https://github.com/jewettaij/jacobi_pd
/// @note  The "Vector" and "Matrix" type arguments can be any 
///        C or C++ object that support indexing, including pointers or vectors.
/// @details
/// -- Example: --
///
/// int n = 5;       // Matrix size
/// double **M;      // A symmetric n x n matrix you want to diagonalize
/// double *evals;   // Store the eigenvalues here.
/// double **evects; // Store the eigenvectors here.
/// // Allocate space for M, evals, and evects, and load contents of M (omitted)
///
/// // Now create an instance of Jacobi ("eigen_calc"). This will allocate space
/// // for storing intermediate calculations.  Once created, it can be reused
/// // multiple times without paying the cost of allocating memory on the heap.
/// 
/// Jacobi<double, double*, double**> eigen_calc(n);
///
/// // Note:
/// // If the matrix you plan to diagonalize (M) is read-only, use this instead:
/// // Jacobi<double, double*, double**, double const*const*> eigen_calc(n);
/// // If you prefer using vectors over C-style pointers, this works also:
/// // Jacobi<double, vector<double>&, vector<vector<double>>&> eigen_calc(n);
///
/// // Now, calculate the eigenvalues and eigenvectors of M
///
/// eigen_calc.Diagonalize(M, evals, evects);
///
/// --- end of example ---

template<typename Scalar,
         typename Vector,
         typename Matrix,
         typename ConstMatrix=Matrix>

class Jacobi
{
  int n;            //!< the size of the matrices you want to diagonalize
  Scalar **M;       //!< local copy of the current matrix being analyzed
  // Precomputed cosine, sine, and tangent of the most recent rotation angle:
  Scalar c;         //!< = cos(θ)
  Scalar s;         //!< = sin(θ)
  Scalar t;         //!< = tan(θ),  (note |t|<=1)
  int *max_idx_row; //!< = keep track of the the maximum element in row i (>i)

public:

  /// @brief  Specify the size of the matrices you want to diagonalize later.
  /// @param n  the size (ie. number of rows) of the (square) matrix.
  void SetSize(int n);

  Jacobi(int n = 0) { Init(); SetSize(n); }

  ~Jacobi() { Dealloc(); }

  // @typedef choose the criteria for sorting eigenvalues and eigenvectors
  typedef enum eSortCriteria {
    DO_NOT_SORT,
    SORT_DECREASING_EVALS,
    SORT_INCREASING_EVALS,
    SORT_DECREASING_ABS_EVALS,
    SORT_INCREASING_ABS_EVALS
  } SortCriteria;

  /// @brief Calculate the eigenvalues and eigevectors of a symmetric matrix
  ///        using the Jacobi eigenvalue algorithm.
  /// @returns The number_of_sweeps (= number_of_iterations / (n*(n-1)/2)).
  ///          If this equals max_num_sweeps, the algorithm failed to converge.
  /// @note  To reduce the computation time further, set calc_evecs=false.
  int
  Diagonalize(ConstMatrix mat, //!< the matrix you wish to diagonalize (size n)
              Vector eval,   //!< store the eigenvalues here
              Matrix evec,   //!< store the eigenvectors here (in rows)
              SortCriteria sort_criteria=SORT_DECREASING_EVALS,//!<sort results?
              bool calc_evecs=true,      //!< calculate the eigenvectors?
              int max_num_sweeps = 50);  //!< limit the number of iterations

private:
  // (Descriptions of private functions can be found in their implementation.)
  void CalcRot(Scalar const *const *M, int i, int j);
  void ApplyRot(Scalar **M, int i, int j);
  void ApplyRotLeft(Matrix E, int i, int j);
  int MaxEntryRow(Scalar const *const *M, int i) const;
  void MaxEntry(Scalar const *const *M, int& i_max, int& j_max) const;
  void SortRows(Vector v, Matrix M, int n, SortCriteria s=SORT_DECREASING_EVALS) const;
  void Init();
  void Alloc(int n);
  void Dealloc();
public:
  // C++ boilerplate: copy and move constructor, swap, and assignment operator
  Jacobi(const Jacobi<Scalar, Vector, Matrix, ConstMatrix>& source);
  Jacobi(Jacobi<Scalar, Vector, Matrix, ConstMatrix>&& other);
  void swap(Jacobi<Scalar, Vector, Matrix, ConstMatrix> &other);
  Jacobi<Scalar, Vector, Matrix, ConstMatrix>& operator = (Jacobi<Scalar, Vector, Matrix, ConstMatrix> source);

}; // class Jacobi





// ---- Eigendecomposition of large sparse (and dense) matrices ----

// The "LambdaLanczos" is a class useful for calculating eigenvalues
// and eigenvectors of large sparse matrices.  Unfortunately, before the
// LambdaLanczos class can be declared, several additional expressions,
// classes and functions that it depends on must be declared first.

// @brief  Create random vectors used at the beginning of the Lanczos algorithm.
// @note "Partially specialization of function" is not allowed, so
//       it is mimicked by wrapping the "init" function with a class template.
template <typename T>
struct VectorRandomInitializer {
public:
  static void init(std::vector<T>&);
};

template <typename T>
struct VectorRandomInitializer<std::complex<T>> {
public:
  static void init(std::vector<std::complex<T>>&);
};

/// @brief  Return the number of significant decimal digits of type T.
template <typename T>
inline constexpr int sig_decimal_digit() {
  return (int)(std::numeric_limits<T>::digits *
               std::log10(std::numeric_limits<T>::radix));
}

/// @brief  Return 10^-n where n=number of significant decimal digits of type T.
template <typename T>
inline constexpr T minimum_effective_decimal() {
  return std::pow(10, -sig_decimal_digit<T>());
}


/// @brief The LambdaLanczos class provides a general way to calculate
///        the smallest or largest eigenvalue and the corresponding eigenvector
///        of a symmetric (Hermitian) matrix using the Lanczos algorithm.
///        The characteristic feature of LambdaLanczos is that the matrix-vector
///        multiplication routine used in the Lanczos algorithm is adaptable.
/// @details
/// @code
///
/// //Example:
/// const int n = 3;
/// double M[n][n] = { {-1.0, -1.0, 1.0},
///                    {-1.0, 1.0, 1.0},
///                    { 1.0, 1.0, 1.0} };
/// // (Its eigenvalues are {-2, 1, 2})
///
/// // Specify the matrix-vector multiplication function
/// auto mv_mul = [&](const vector<double>& in, vector<double>& out) {
///   for(int i = 0;i < n;i++) {
///     for(int j = 0;j < n;j++) {
///       out[i] += M[i][j]*in[j];
///     }
///   } 
/// };
///
/// LambdaLanczos<double> engine(mv_mul, n, true);
///     // ("true" means to calculate the largest eigenvalue.)
/// engine.eigenvalue_offset = 3.0   # = max_i{Σ_j|Mij|}  (see below)
/// double eigenvalue;
/// vector<double> eigenvector(n);
/// int itern = engine.run(eigenvalue, eigenvector);
///
/// cout << "Iteration count: " << itern << endl;
/// cout << "Eigen value: " << setprecision(16) << eigenvalue << endl;
/// cout << "Eigen vector:";
/// for(int i = 0; i < n; i++) {
///   cout << eigenvector[i] << " ";
/// }
/// cout << endl;
///
/// @endcode
/// This feature allows you to use a matrix whose elements are partially given,
/// e.g. a sparse matrix whose non-zero elements are stored as a list of
/// {row-index, column-index, value} tuples.  You can also easily combine
/// LambdaLanczos with existing matrix libraries (e.g. Eigen)
///
/// @note
///   If the matrices you want to analyze are ordinary square matrices, (as in
///   the example) it might be easier to use "PEigenDense" instead.  (It is a
///   wrapper which takes care of all of the LambdaLanczos details for you.)
///
/// @note
///  IMPORTANT:
///  The Lanczos algorithm finds the largest magnitude eigenvalue, so you
///  MUST ensure that the eigenvalue you are seeking has the largest magnitude
///  (regardless of whether it is the maximum or minimum eigenvalue).
///  To insure that this is so, you can add or subtract a number to all
///  of the eigenvalues of the matrix by specifying the "eigenvalue_offset".
///  This number should exceed the largest magnitude eigenvalue of the matrix.
///  According to the Gershgorin theorem, you can estimate this number using
///     r = max_i{Σ_j|Mij|} = max_j{Σ_i|Mij|}
///  (where Mij are the elements of the matrix and Σ_j denotes the sum over j).
///  If find_maximum == true (if you are seeking the maximum eigenvalue), then
///    eigenvalue_offset = +r
///  If find_maximum == false, then
///    eigenvalue_offset = -r
///  The eigenvalue_offset MUST be specified by the user.  LambdaLanczos does
///  not have an efficient and general way to access the elements of the matrix.
///
///  (You can omit this step if you are seeking the maximum eigenvalue,
///   and the matrix is positive definite, or if you are seeking the minimum
///   eigenvalue and the matrix is negative definite.)
///
/// @note
///   LambdaLanczos is available under the MIT license and downloadable at:
///   https://github.com/mrcdr/lambda-lanczos

template <typename T>
class LambdaLanczos {
public:
  LambdaLanczos();
  LambdaLanczos(std::function<void(const std::vector<T>&, std::vector<T>&)> mv_mul, int matrix_size, bool find_maximum);
  LambdaLanczos(std::function<void(const std::vector<T>&, std::vector<T>&)> mv_mul, int matrix_size) : LambdaLanczos(mv_mul, matrix_size, true) {}

  /// @brief  Calculate the principal (largest or smallest) eigenvalue
  ///         of the matrix (and its corresponding eigenvector).
  int run(real_t<T>&, std::vector<T>&) const;

  // --- public data members ---

  /// @brief  Specify the size of the matrix you will analyze.
  ///         (This equals the size of the eigenvector which will be returned.)
  int matrix_size;

  /// @brief  Specify the function used for matrix*vector multiplication
  ///         used by the Lanczos algorithm.  For an ordinary dense matrix,
  ///         this function is the ordinary matrix*vector product. (See the
  ///         example above.  For a sparse matrix, it will be something else.)
  std::function<void(const std::vector<T>&, std::vector<T>&)> mv_mul;

  /// @brief  Are we searching for the maximum or minimum eigenvalue?
  /// @note   (Usually, you must also specify eigenvalue_offset.)
  bool find_maximum = false;

  /// @brief Shift all the eigenvalues by "eigenvalue_offset" during the Lanczos
  ///        iteration (ie. during LambdaLanczos::run()).  The goal is to insure
  ///        that the correct eigenvalue is selected (the one with the maximum
  ///        magnitude).
  /// @note  The eigevalue returned by LambdaLanczos::run() is not effected 
  ///        because after the iteration is finished, it will subtract this 
  ///        number from the eigenvalue before it is returned to the caller.
  /// @note  Unless your matrix is positive definite or negative definite,
  ///        you MUST specify eigenvalue_offset.  See comment above for details.
  real_t<T> eigenvalue_offset = 0.0;

  /// @brief This function sets "eigenvalue_offset" automatically.
  /// @note  Using this function is not recommended because it is very slow.
  ///        For efficiency, set the "eigenvalue_offset" yourself.
  void ChooseOffset();

  // The remaining data members usually can be left alone:
  int max_iteration;
  real_t<T> eps = minimum_effective_decimal<real_t<T>>() * 1e3;
  real_t<T> tridiag_eps_ratio = 1e-1;
  int initial_vector_size = 200;
  std::function<void(std::vector<T>&)> init_vector =
    VectorRandomInitializer<T>::init;

  // (for those who prefer "Set" functions...)
  int SetSize(int matrix_size);
  void SetMul(std::function<void(const std::vector<T>&,
                                 std::vector<T>&)> mv_mul);
  void SetInitVec(std::function<void(std::vector<T>&)> init_vector);
  void SetFindMax(bool find_maximum);
  void SetEvalOffset(T eigenvalue_offset);
  void SetEpsilon(T eps);
  void SetTriEpsRatio(T tridiag_eps_ratio);

private:
  static void schmidt_orth(std::vector<T>&, const std::vector<std::vector<T>>&);
  real_t<T> find_minimum_eigenvalue(const std::vector<real_t<T>>&,
                                    const std::vector<real_t<T>>&) const;
  real_t<T> find_maximum_eigenvalue(const std::vector<real_t<T>>&,
                                    const std::vector<real_t<T>>&) const;
  static real_t<T> tridiagonal_eigen_limit(const std::vector<real_t<T>>&,
                                           const std::vector<real_t<T>>&);
  static int num_of_eigs_smaller_than(real_t<T>,
                                      const std::vector<real_t<T>>&,
                                      const std::vector<real_t<T>>&);
  real_t<T> UpperBoundEvals() const;
};



/// @brief
///   PEigenDense is a class containing only one useful member function:
///   PrincipalEigen().  This function calculates the principal (largest
///   or smallest) eigenvalue and corresponding eigenvector of a square
///   n x n matrix.  This can be faster than diagionalizing the entire matrix.
///   (For example by using the Lanczos algorithm or something similar.)
/// @note
///   This code is a wrapper. Internally, it uses the "LambdaLanczos" class.
/// @note
///   For matrices larger than 13x13, PEigenDense::PrincipleEigen()
///   is usually faster than Jacobi::Diagonalize().)

template<typename Scalar, typename Vector, typename ConstMatrix>
class PEigenDense
{
  size_t n;                 // the size of the matrix
  std::vector<Scalar> evec; // preallocated vector

public:
  void SetSize(int matrix_size) {
    n = matrix_size;
    evec.resize(n);
  }

  PEigenDense(int matrix_size=0):evec(matrix_size) {
    SetSize(matrix_size);
  }

  /// @brief  Calculate the principal eigenvalue and eigenvector of a matrix.
  /// @return Return the principal eigenvalue of the matrix.
  ///         If you want the eigenvector, pass a non-null "evector" argument.
  Scalar
  PrincipalEigen(ConstMatrix matrix,   //!< the input patrix
                 Vector evector,       //!< the eigenvector is stored here
                 bool find_max=false); //!< want the max or min eigenvalue?

}; // class PEigenDense



// --------------------------------------
// ----------- IMPLEMENTATION -----------
// --------------------------------------




// --- Implementation: Memory allocation for matrices ---
template<typename Entry>
void Alloc2D(size_t nrows,          // size of the array (number of rows)
             size_t ncols,          // size of the array (number of columns)
             Entry ***paaX)         // pointer to a 2D C-style array
{
  assert(paaX);
  *paaX = new Entry* [nrows];  //conventional 2D C array (pointer-to-pointer)
  (*paaX)[0] = new Entry [nrows * ncols];  // 1D C array (contiguous memor)
  for(size_t iy=0; iy<nrows; iy++)
    (*paaX)[iy] = (*paaX)[0] + iy*ncols;
  // The caller can access the contents using (*paaX)[i][j]
}

template<typename Entry>
void Dealloc2D(Entry ***paaX)         // pointer to a 2D C-style array
{
  if (paaX && *paaX) {
    delete [] (*paaX)[0];
    delete [] (*paaX);
    *paaX = nullptr;
  }
}

/// @brief ConjugateProduct::prod(a,b) is a function which returns a*b by
/// default. If the arguments are complex numbers, it returns conj(a)*b instead.
template <typename T>
struct ConjugateProduct {
public:
  static T prod(T a, T b) { return a*b; }
};

template <typename T>
struct ConjugateProduct<std::complex<T>> {
public:
  static std::complex<T> prod(std::complex<T> a, std::complex<T> b) {
    return std::conj(a)*b;
  }
};


// --- Implementation: Operations on vectors (of real and complex numbers) ---

template <typename T>
inline T inner_prod(const std::vector<T>& v1, const std::vector<T>& v2) {
  return std::inner_product(std::begin(v1), std::end(v1),
                            std::begin(v2), T(),
                            [](T a, T b) -> T { return a+b; },
                            ConjugateProduct<T>::prod);
  // T() means zero value of type T
  // This spec is required because std::inner_product calculates
  // v1*v2 not conj(v1)*v2
}

template <typename T>
inline real_t<T> l2_norm(const std::vector<T>& vec) {
  return std::sqrt(std::real(inner_prod(vec, vec)));
  // The norm of any complex vector <v|v> is real by definition.
}

template <typename T1, typename T2>
inline void scalar_mul(T1 a, std::vector<T2>& vec) {
  int n = vec.size();
  for(int i = 0;i < n;i++)
    vec[i] *= a;
}

template <typename T>
inline void normalize(std::vector<T>& vec) {
  scalar_mul(1.0/l2_norm(vec), vec);
}


template <typename T>
inline real_t<T> l1_norm(const std::vector<T>& vec) {
  real_t<T> norm = real_t<T>(); // Zero initialization
  for(const T& element : vec)
    norm += std::abs(element);
  return norm;
}




// --- Implementation: Eigendecomposition of small dense matrices ---

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
int Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Diagonalize(ConstMatrix mat,    // the matrix you wish to diagonalize (size n)
            Vector eval,        // store the eigenvalues here
            Matrix evec,        // store the eigenvectors here (in rows)
            SortCriteria sort_criteria, // sort results?
            bool calc_evec,     // calculate the eigenvectors?
            int max_num_sweeps) // limit the number of iterations ("sweeps")
{
  // -- Initialization --
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)          //copy mat[][] into M[][]
      M[i][j] = mat[i][j];               //(M[][] is a local copy we can modify)

  if (calc_evec)
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        evec[i][j] = (i==j) ? 1.0 : 0.0; //Set evec equal to the identity matrix

  for (int i = 0; i < n-1; i++)          //Initialize the "max_idx_row[]" array 
    max_idx_row[i] = MaxEntryRow(M, i);  //(which is needed by MaxEntry())

  // -- Iteration --
  int n_iters;
  int max_num_iters = max_num_sweeps*n*(n-1)/2; //"sweep" = n*(n-1)/2 iters
  for (n_iters=0; n_iters < max_num_iters; n_iters++) {
    int i,j;
    MaxEntry(M, i, j); // Find the maximum entry in the matrix. Store in i,j

    // If M[i][j] is small compared to M[i][i] and M[j][j], set it to 0.
    if ((M[i][i] + M[i][j] == M[i][i]) && (M[j][j] + M[i][j] == M[j][j])) {
      M[i][j] = 0.0;
      max_idx_row[i] = MaxEntryRow(M,i); //must also update max_idx_row[i]
    }

    if (M[i][j] == 0.0)
      break;

    // Otherwise, apply a rotation to make M[i][j] = 0
    CalcRot(M, i, j);  // Calculate the parameters of the rotation matrix.
    ApplyRot(M, i, j); // Apply this rotation to the M matrix.
    if (calc_evec)     // Optional: If the caller wants the eigenvectors, then
      ApplyRotLeft(evec,i,j); // apply the rotation to the eigenvector matrix

  } //for (int n_iters=0; n_iters < max_num_iters; n_iters++)

  // -- Post-processing --
  for (int i = 0; i < n; i++)
    eval[i] = M[i][i];

  // Optional: Sort results by eigenvalue.
  SortRows(eval, evec, n, sort_criteria);

  return n_iters / (n*(n-1)/2); //returns the number of "sweeps" (converged?)
}


/// brief  Calculate the components of a rotation matrix which performs a
///        rotation in the i,j plane by an angle (θ) that (when multiplied on
///        both sides) will zero the ij'th element of M, so that afterwards
///        M[i][j] = 0.  The results will be stored in c, s, and t
///        (which store cos(θ), sin(θ), and tan(θ), respectively).

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
CalcRot(Scalar const *const *M,    //!< matrix
        int i,       //!< row index
        int j)       //!< column index
{
  t = 1.0; // = tan(θ)
  Scalar M_jj_ii = (M[j][j] - M[i][i]);
  if (M_jj_ii != 0.0) {
    // kappa = (M[j][j] - M[i][i]) / (2*M[i][j])
    Scalar kappa = M_jj_ii;
    t = 0.0;
    Scalar M_ij = M[i][j];
    if (M_ij != 0.0) {
      kappa /= (2.0*M_ij);
      // t satisfies: t^2 + 2*t*kappa - 1 = 0
      // (choose the root which has the smaller absolute value)
      t = 1.0 / (std::sqrt(1 + kappa*kappa) + std::abs(kappa));
      if (kappa < 0.0)
        t = -t;
    }
  }
  assert(std::abs(t) <= 1.0);
  c = 1.0 / std::sqrt(1 + t*t);
  s = c*t;
}


/// brief   Perform a similarity transformation by multiplying matrix M on both
///         sides by a rotation matrix (and its transpose) to eliminate M[i][j].
/// details This rotation matrix performs a rotation in the i,j plane by
///         angle θ.  This function assumes that c=cos(θ). s=som(θ), t=tan(θ)
///         have been calculated previously (using the CalcRot() function).
///         It also assumes that i<j.  The max_idx_row[] array is also updated.
///         To save time, since the matrix is symmetric, the elements
///         below the diagonal (ie. M[u][v] where u>v) are not computed.
///
/// verbatim
///
///   M' = R^T * M * R
/// where R the rotation in the i,j plane and ^T denotes the transpose.
///                 i         j
///       _                             _
///      |  1                            | 
///      |    .                          |
///      |      .                        |
///      |        1                      |
///      |          c   ...   s          |
///      |          .  .      .          |
/// R  = |          .    1    .          |
///      |          .      .  .          |
///      |          -s  ...   c          |
///      |                      1        |
///      |                        .      |
///      |                          .    |
///      |_                           1 _|
///
/// endverbatim
///
/// Let M' denote the matrix M after multiplication by R^T and R.
/// The components of M' are:
///   M'_uv =  Σ_w  Σ_z   R_wu * M_wz * R_zv
///
/// note
/// The rotation at location i,j will modify all of the matrix
/// elements containing at least one index which is either i or j
/// such as: M[w][i], M[i][w], M[w][j], M[j][w].
/// Check and see whether these modified matrix elements exceed the 
/// corresponding values in max_idx_row[] array for that row.
/// If so, then update max_idx_row for that row.
/// This is somewhat complicated by the fact that we must only consider
/// matrix elements in the upper-right triangle strictly above the diagonal.
/// (ie. matrix elements whose second index is > the first index).

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
ApplyRot(Scalar **M,  // matrix
         int i,       // row index
         int j)       // column index
{
  // Recall that c = cos(θ), s = sin(θ), t = tan(θ) (and t <= 1.0)

  // Compute the diagonal elements of M which have changed:
  M[i][i] -= t * M[i][j];
  M[j][j] += t * M[i][j];

  //Update the off-diagonal elements of M which will change (above the diagonal)
  assert(i < j);
  M[i][j] = 0.0;

  //compute M[w][i] and M[i][w] for all w!=i,considering above-diagonal elements
  for (int w=0; w < i; w++) {        // 0 <= w <  i  <  j < n
    M[i][w] = M[w][i]; //backup the previous value. store below diagonal (i>w)
    M[w][i] = c*M[w][i] - s*M[w][j]; //M[w][i], M[w][j] from previous iteration
    if (i == max_idx_row[w]) max_idx_row[w] = MaxEntryRow(M, w);
    else if (std::abs(M[w][i])>std::abs(M[w][max_idx_row[w]])) max_idx_row[w]=i;
    //assert(max_idx_row[w] == MaxEntryRow(M, w));
  }
  for (int w=i+1; w < j; w++) {      // 0 <= i <  w  <  j < n
    M[w][i] = M[i][w]; //backup the previous value. store below diagonal (w>i)
    M[i][w] = c*M[i][w] - s*M[w][j]; //M[i][w], M[w][j] from previous iteration
  }
  for (int w=j+1; w < n; w++) {      // 0 <= i < j+1 <= w < n
    M[w][i] = M[i][w]; //backup the previous value. store below diagonal (w>i)
    M[i][w] = c*M[i][w] - s*M[j][w]; //M[i][w], M[j][w] from previous iteration
  }

  // now that we're done modifying row i, we can update max_idx_row[i]
  max_idx_row[i] = MaxEntryRow(M, i);

  //compute M[w][j] and M[j][w] for all w!=j,considering above-diagonal elements
  for (int w=0; w < i; w++) {        // 0 <=  w  <  i <  j < n
    M[w][j] = s*M[i][w] + c*M[w][j]; //M[i][w], M[w][j] from previous iteration
    if (j == max_idx_row[w]) max_idx_row[w] = MaxEntryRow(M, w);
    else if (std::abs(M[w][j])>std::abs(M[w][max_idx_row[w]])) max_idx_row[w]=j;
    //assert(max_idx_row[w] == MaxEntryRow(M, w));
  }
  for (int w=i+1; w < j; w++) {      // 0 <= i+1 <= w <  j < n
    M[w][j] = s*M[w][i] + c*M[w][j]; //M[w][i], M[w][j] from previous iteration
    if (j == max_idx_row[w]) max_idx_row[w] = MaxEntryRow(M, w);
    else if (std::abs(M[w][j])>std::abs(M[w][max_idx_row[w]])) max_idx_row[w]=j;
    //assert(max_idx_row[w] == MaxEntryRow(M, w));
  }
  for (int w=j+1; w < n; w++) {      // 0 <=  i  <  j <  w < n
    M[j][w] = s*M[w][i] + c*M[j][w]; //M[w][i], M[j][w] from previous iteration
  }
  // now that we're done modifying row j, we can update max_idx_row[j]
  max_idx_row[j] = MaxEntryRow(M, j);
} //Jacobi::ApplyRot()


/// brief  Multiply matrix E on the left by the (previously calculated)
///        rotation matrix.
///
/// details
/// Multiply matrix M on the LEFT side by a transposed rotation matrix, R^T.
/// This matrix performs a rotation in the i,j plane by angle θ
/// (where the arguments "s" and "c" refer to cos(θ) and sin(θ), respectively).
///
/// verbatim
///   E'_uv = Σ_w  R_wu * E_wv
/// endverbatim

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
ApplyRotLeft(Matrix E,  // matrix
             int i,     // row index
             int j)     // column index
{
  // recall that c = cos(θ) and s = sin(θ)
  for (int v = 0; v < n; v++) {
    Scalar Eiv = E[i][v]; //backup E[i][v]
    E[i][v] = c*E[i][v] - s*E[j][v];
    E[j][v] = s*Eiv     + c*E[j][v];
  }
}

/// brief  Find the off-diagonal index in row i whose absolute value is largest
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
int Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
MaxEntryRow(Scalar const *const *M, int i) const {
  int j_max = i+1;
  for(int j = i+2; j < n; j++)
    if (std::abs(M[i][j]) > std::abs(M[i][j_max]))
      j_max = j;
  return j_max;
}

/// brief  Find the indices (i_max, j_max) marking the location of the
///        entry in the matrix with the largest absolute value.  This
///        uses the max_idx_row[] array to find the answer in O(n) time.
/// returns  This function does not return a avalue.  However after it is
///          invoked, the location of the largest matrix element will be
///          stored in the i_max and j_max arguments.
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
MaxEntry(Scalar const *const *M, int& i_max, int& j_max) const {
  // find the maximum entry in the matrix M in O(n) time
  i_max = 0;
  j_max = max_idx_row[i_max];
  Scalar max_entry = std::abs(M[i_max][j_max]);
  int nm1 = n-1;
  for (int i=1; i < nm1; i++) {
    int j = max_idx_row[i];
    if (std::abs(M[i][j]) > max_entry) {
      max_entry = std::abs(M[i][j]);
      i_max = i;
      j_max = j;
    }
  }
}

/// brief  Sort the rows in matrix "evec" according to the numbers in "eval".
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
SortRows(Vector eval,        // vector containing the keys used for sorting
         Matrix evec,        // matrix whose rows will be sorted according to v
         int n,              // size of the vector and matrix
         SortCriteria sort_criteria) const // sort eigenvalues?
{
  for (int i = 0; i < n-1; i++) {
    int i_max = i;
    for (int j = i+1; j < n; j++) {
      // find the "maximum" element in the array starting at position i+1
      switch (sort_criteria) {
      case SORT_DECREASING_EVALS:
        if (eval[j] > eval[i_max])
          i_max = j;
        break;
      case SORT_INCREASING_EVALS:
        if (eval[j] < eval[i_max])
          i_max = j;
        break;
      case SORT_DECREASING_ABS_EVALS:
        if (std::abs(eval[j]) > std::abs(eval[i_max]))
          i_max = j;
        break;
      case SORT_INCREASING_ABS_EVALS:
        if (std::abs(eval[j]) < std::abs(eval[i_max]))
          i_max = j;
        break;
      default:
        break;
      }
    }
    std::swap(eval[i], eval[i_max]); // sort "eval"
    for (int k = 0; k < n; k++)
      std::swap(evec[i][k], evec[i_max][k]); // sort "evec"
  }
}

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Init() {
  n = 0;
  M = nullptr;
  max_idx_row = nullptr;
}

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
SetSize(int n) {
  Dealloc();
  Alloc(n);
}

// Implementation: Jacobi memory management:

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Alloc(int n) {
  this->n = n;
  if (n > 0) {
    max_idx_row = new int[n];
    Alloc2D(n, n, &M);
  }
}

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Dealloc() {
  if (max_idx_row)
    delete [] max_idx_row;
  Dealloc2D(&M);
  Init();
}

// Jacobi copy and move constructor, swap, and assignment operator:

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Jacobi(const Jacobi<Scalar, Vector, Matrix, ConstMatrix>& source)
{
  Init();
  SetSize(source.n);
  assert(n == source.n);
  // The following lines aren't really necessary, because the contents
  // of source.M and source.max_idx_row are not needed (since they are
  // overwritten every time Jacobi::Diagonalize() is invoked).
  std::copy(source.max_idx_row,
            source.max_idx_row + n,
            max_idx_row);
  for (int i = 0; i < n; i++)
    std::copy(source.M[i],
              source.M[i] + n,
              M[i]);
}

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
swap(Jacobi<Scalar, Vector, Matrix, ConstMatrix> &other) {
  std::swap(n, other.n);
  std::swap(max_idx_row, other.max_idx_row);
  std::swap(M, other.M);
}

// Move constructor (C++11)
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Jacobi(Jacobi<Scalar, Vector, Matrix, ConstMatrix>&& other) {
  Init();
  swap(*this, other);
}

// Using the "copy-swap" idiom for the assignment operator
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>&
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
operator = (Jacobi<Scalar, Vector, Matrix, ConstMatrix> source) {
  this->swap(source);
  return *this;
}



// --- Implementation: Eigendecomposition of large matrices ----

template <typename T>
inline LambdaLanczos<T>::LambdaLanczos() {
  this->matrix_size = 0;
  this->max_iteration = 0;
  this->find_maximum = 0;
}


template <typename T>
inline LambdaLanczos<T>::
LambdaLanczos(std::function<void(const std::vector<T>&,
                                 std::vector<T>&)> mv_mul,
                                 int matrix_size,
              bool find_maximum)
{
  this->mv_mul = mv_mul;
  this->matrix_size = matrix_size;
  this->max_iteration = matrix_size;
  this->find_maximum = find_maximum;
}


template <typename T>
inline int LambdaLanczos<T>::
run(real_t<T>& eigvalue, std::vector<T>& eigvec) const
{ 
  assert(matrix_size > 0);
  assert(0 < this->tridiag_eps_ratio && this->tridiag_eps_ratio < 1);
  
  std::vector<std::vector<T>> u;     // Lanczos vectors
  std::vector<real_t<T>> alpha; // Diagonal elements of an approximated tridiagonal matrix
  std::vector<real_t<T>> beta;  // Subdiagonal elements of an approximated tridiagonal matrix

  const int n = this->matrix_size;

  u.reserve(this->initial_vector_size);
  alpha.reserve(this->initial_vector_size);
  beta.reserve(this->initial_vector_size);

  u.emplace_back(n, 0.0); // Same as u.push_back(std::vector<T>(n, 0.0))
  
  std::vector<T> vk(n, 0.0);
  
  real_t<T> alphak = 0.0;
  alpha.push_back(alphak);
  real_t<T> betak = 0.0;
  beta.push_back(betak);
  
  std::vector<T> uk(n);
  this->init_vector(uk);
  normalize(uk);
  u.push_back(uk);

  real_t<T> ev, pev; // Calculated eigenvalue and previous one
  pev = std::numeric_limits<real_t<T>>::max();

  int itern = this->max_iteration;
  for(int k = 1;k <= this->max_iteration;k++) {
    // vk = (A + offset*E)uk, here E is the identity matrix
    for(int i = 0;i < n;i++) {
      vk[i] = uk[i]*this->eigenvalue_offset;
    }
    this->mv_mul(uk, vk);
    
    alphak = std::real(inner_prod(u.back(), vk));
    
    // The inner product <uk|vk> is real.
    // Proof:
    //     <uk|vk> = <uk|A|uk>
    //   On the other hand its complex conjugate is
    //     <uk|vk>^* = <vk|uk> = <uk|A^*|uk> = <uk|A|uk>
    //   here the condition that matrix A is a symmetric (Hermitian) is used.
    //   Therefore
    //     <uk|vk> = <vk|uk>^*
    //   <uk|vk> is real.
    
    alpha.push_back(alphak);
    
    for(int i = 0;i < n; i++) {
      uk[i] = vk[i] - betak*u[k-1][i] - alphak*u[k][i];
    }
    
    schmidt_orth(uk, u);
    
    betak = l2_norm(uk);
    beta.push_back(betak);

    if(this->find_maximum) {
      ev = find_maximum_eigenvalue(alpha, beta);
    } else {
      ev = find_minimum_eigenvalue(alpha, beta);
    }

    const real_t<T> zero_threshold = minimum_effective_decimal<real_t<T>>()*1e-1;
    if(betak < zero_threshold) {
      u.push_back(uk);
      // This element will never be accessed,
      // but this "push" guarantees u to always have one more element than 
      // alpha and beta do.
      itern = k;
      break;
    }

    normalize(uk);
    u.push_back(uk);

    if(abs(ev-pev) < std::min(abs(ev), abs(pev))*this->eps) {
      itern = k;
      break;
    } else {
      pev = ev;
    }
  }

  eigvalue = ev - this->eigenvalue_offset;

  int m = alpha.size();
  std::vector<T> cv(m+1);
  cv[0] = 0.0;
  cv[m] = 0.0;
  cv[m-1] = 1.0;
  
  beta[m-1] = 0.0;

  if(eigvec.size() < n) {
    eigvec.resize(n);
  }
  
  for(int i = 0;i < n;i++) {
    eigvec[i] = cv[m-1]*u[m-1][i];
  }

  for(int k = m-2;k >= 1;k--) {
    cv[k] = ((ev - alpha[k+1])*cv[k+1] - beta[k+1]*cv[k+2])/beta[k];

    for(int i = 0;i < n;i++) {
      eigvec[i] += cv[k]*u[k][i];
    }
  }

  normalize(eigvec);

  return itern;

} //LambdaLancos::run()



template <typename T>
inline void LambdaLanczos<T>::
schmidt_orth(std::vector<T>& uorth, const std::vector<std::vector<T>>& u)
{
  // Vectors in u must be normalized, but uorth doesn't have to be.
  
  int n = uorth.size();
  
  for(int k = 0;k < u.size();k++) {
    T innprod = inner_prod(uorth, u[k]);
    for(int i = 0;i < n;i++)
      uorth[i] -= innprod * u[k][i];
  }
}


template <typename T>
inline real_t<T> LambdaLanczos<T>::
find_minimum_eigenvalue(const std::vector<real_t<T>>& alpha,
                        const std::vector<real_t<T>>& beta) const
{
  real_t<T> eps = this->eps * this->tridiag_eps_ratio;
  real_t<T> pmid = std::numeric_limits<real_t<T>>::max();
  real_t<T> r = tridiagonal_eigen_limit(alpha, beta);
  real_t<T> lower = -r;
  real_t<T> upper = r;
  real_t<T> mid;
  int nmid; // Number of eigenvalues smaller than the "mid"
  
  while(upper-lower > std::min(abs(lower), abs(upper))*eps) {
    mid = (lower+upper)/2.0;
    nmid = num_of_eigs_smaller_than(mid, alpha, beta);
    if(nmid >= 1) {
      upper = mid;
    } else {
      lower = mid;
    }
    
    if(mid == pmid) {
      break; // This avoids an infinite loop due to zero matrix
    }
    pmid = mid;
  }

  return lower; // The "lower" almost equals the "upper" here.
}


template <typename T>
inline real_t<T> LambdaLanczos<T>::
find_maximum_eigenvalue(const std::vector<real_t<T>>& alpha,
                        const std::vector<real_t<T>>& beta) const
{
  real_t<T> eps = this->eps * this->tridiag_eps_ratio;
  real_t<T> pmid = std::numeric_limits<real_t<T>>::max();
  real_t<T> r = tridiagonal_eigen_limit(alpha, beta);
  real_t<T> lower = -r;
  real_t<T> upper = r;
  real_t<T> mid;
  int nmid; // Number of eigenvalues smaller than the "mid"

  int m = alpha.size() - 1;  // Number of eigenvalues of the approximated 
                             // triangular matrix, which equals the rank of it
  
  
  while(upper-lower > std::min(abs(lower), abs(upper))*eps) {
    mid = (lower+upper)/2.0;
    nmid = num_of_eigs_smaller_than(mid, alpha, beta);
    
    if(nmid < m) {
      lower = mid;
    } else {
      upper = mid;
    }
    
    if(mid == pmid) {
      break; // This avoids an infinite loop due to zero matrix
    }
    pmid = mid;
  }

  return lower; // The "lower" almost equals the "upper" here.
}  


/// @brief
/// Compute the upper bound of the absolute value of eigenvalues
/// by Gerschgorin theorem. This routine gives a rough upper bound,
/// but it is sufficient because the bisection routine using
/// the upper bound converges exponentially.

template <typename T>
inline real_t<T> LambdaLanczos<T>::
tridiagonal_eigen_limit(const std::vector<real_t<T>>& alpha,
                        const std::vector<real_t<T>>& beta)
{
  real_t<T> r = l1_norm(alpha);
  r += 2*l1_norm(beta);
  
  return r;
}



// Algorithm from
//  Peter Arbenz et al.
// "High Performance Algorithms for Structured Matrix Problems"
//  Nova Science Publishers, Inc.

template <typename T>
inline int LambdaLanczos<T>::
num_of_eigs_smaller_than(real_t<T> c,
                         const std::vector<real_t<T>>& alpha,
                         const std::vector<real_t<T>>& beta)
{
  real_t<T> q_i = 1.0;
  int count = 0;
  int m = alpha.size();
  
  for(int i = 1;i < m;i++){
    q_i = alpha[i] - c - beta[i-1]*beta[i-1]/q_i;
    if(q_i < 0){
      count++;
    }
    if(q_i == 0){
      q_i = minimum_effective_decimal<real_t<T>>();
    }
  }

  return count;
}


template <typename T>
inline void LambdaLanczos<T>::ChooseOffset() {
  const auto n = this->matrix_size;
  std::vector<T> unit_vec_j(n);
  std::vector<T> matrix_column_j(n);
  real_t<T> eval_upper_bound = 0.0;
  /// According to Gershgorin theorem, the maximum (magnitude) eigenvalue should
  /// not exceed max_j{Σ_i|Mij|}.  We can infer the contents of each column in
  /// the matrix by multiplying it by different unit vectors.  This is slow.
  for (int j = 0; j < n; j++) {
    std::fill(unit_vec_j.begin(), unit_vec_j.end(), 0);  // fill with zeros
    unit_vec_j[j] = 1.0;  // = jth element is 1, all other elements are 0
    // Multiplying the matrix by a unit vector (a vector containing only one
    // non-zero element at position j) extracts the jth column of the matrix.
    this->mv_mul(unit_vec_j, matrix_column_j);
    real_t<T> sum_column = 0.0;  // compute Σ_i|Mij|
    for (int i = 0; i < n; i++)
      sum_column += std::abs(matrix_column_j[i]);
    if (eval_upper_bound < sum_column)
      eval_upper_bound = sum_column; // compute max_j{Σ_i|Mij|}
  }
  if (find_maximum)
    this->eigenvalue_offset = eval_upper_bound;
  else
    this->eigenvalue_offset = -eval_upper_bound;
}


template <typename T>
inline int LambdaLanczos<T>::SetSize(int matrix_size)
{
  this->matrix_size = matrix_size;
  this->max_iteration = matrix_size;
  return matrix_size;
}

template <typename T>
inline void LambdaLanczos<T>::SetMul(std::function<void(const std::vector<T>&, std::vector<T>&)> mv_mul)
{
  this->mv_mul = mv_mul;
}

template <typename T>
inline void LambdaLanczos<T>::SetInitVec(std::function<void(std::vector<T>&)> init_vector)
{
  this->init_vector = init_vector;
}

template <typename T>
inline void LambdaLanczos<T>::SetFindMax(bool find_maximum) {
  this->find_maximum = find_maximum;
}

template <typename T>
inline void LambdaLanczos<T>::SetEvalOffset(T offset)
{ 
  this->eigenvalue_offset = offset;
}

template <typename T>
inline void LambdaLanczos<T>::SetEpsilon(T epsilon)
{ 
  this->eps = epsilon;
}

template <typename T>
inline void LambdaLanczos<T>::SetTriEpsRatio(T tri_eps_ratio)
{ 
  this->tridiag_eps_ratio = tri_eps_ratio;
}





template <typename T>
inline void VectorRandomInitializer<T>::
init(std::vector<T>& v)
{
  std::random_device dev;
  std::mt19937 mt(dev());
  std::uniform_real_distribution<T> rand((T)(-1.0), (T)(1.0));

  int n = v.size();
  for(int i = 0;i < n;i++) {
    v[i] = rand(mt);
  }

  normalize(v);
}


template <typename T>
inline void VectorRandomInitializer<std::complex<T>>::
init(std::vector<std::complex<T>>& v)
{
  std::random_device dev;
  std::mt19937 mt(dev());
  std::uniform_real_distribution<T> rand((T)(-1.0), (T)(1.0));

  int n = v.size();
  for(int i = 0;i < n;i++) {
    v[i] = std::complex<T>(rand(mt), rand(mt));
  }

  normalize(v);
}


// --- Implementation of PEigenDense

template<typename Scalar, typename Vector, typename ConstMatrix>
Scalar PEigenDense<Scalar, Vector, ConstMatrix>::
PrincipalEigen(ConstMatrix matrix,
               Vector eigenvector,
               bool find_max)
{
  assert(n > 0);
  auto matmul = [&](const std::vector<Scalar>& in, std::vector<Scalar>& out) {
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        out[i] += matrix[i][j]*in[j];
      }
    } 
  };
  auto init_vec = [&](std::vector<Scalar>& vec) {
    for(int i = 0; i < n; i++)
      vec[i] = 0.0;
    vec[0] = 1.0;
  };

  // "ll_engine" calculates the eigenvalue and eigenvector.
  LambdaLanczos<Scalar> ll_engine(matmul, n, find_max);
  
  // The Lanczos algorithm selects the eigenvalue with the largest magnitude.
  // In order to insure that this is the one we want (maxima or minima), we can
  // add a constant to all of the eigenvalues by setting "eigenvalue_offset".
  Scalar eval_upper_bound = 0.0;
  for (int i = 0; i < n; i++) {
    Scalar sum_row = 0.0;
    for (int j = 0; j < n; i++)
      sum_row += std::abs(matrix[i][j]);
    if (eval_upper_bound < sum_row)
      eval_upper_bound = sum_row;
  }
  if (find_max)
    ll_engine.eigenvalue_offset = eval_upper_bound;
  else
    ll_engine.eigenvalue_offset = -eval_upper_bound;

  ll_engine.init_vector = init_vec;

  Scalar eval;

  // This line does all of the hard work:
  size_t itern = ll_engine.run(eval, evec);

  for (int i = 0; i < n; i++)
    eigenvector[i] = evec[i];

  return eval;
}


} //namespace MathEigen


#endif //#ifndef _MATH_EIGEN_H
