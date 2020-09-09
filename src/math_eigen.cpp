#include "math_eigen.h"
#include<array>

using std::vector;
using std::array;


using namespace MathEigen;

// --- Instantiate template class instances   ---
// --- to reduce bloat in the compiled binary ---

// When using one of these versions of Jacobi, you can reduce code bloat by
// using an "extern template class" declaration in your cpp file.  For example:
//
// #include"math_eigen.h"
// extern template class MathEigen::Jacobi<double, double*, double**, double const*const*>;
//
// ...This should (hopefully) use the version of Jacobi defined in this file
// instead of compiling a new version.


// template instantiations of Jacobi for pointer->pointer arrays
template class MathEigen::Jacobi<double, double*, double**, double const*const*>;

// template instantiations of Jacbi for fixed-size (3x3) arrays
template class MathEigen::Jacobi<double, double*, double (*)[3], double const (*)[3]>;

// template instantiations of Jacobi for vector-of-vectors:
template class MathEigen::Jacobi<double,
                                 vector<double>&,
                                 vector<vector<double> >&,
                                 const vector<vector<double> >& >;

// template instantiations of Jacobi for array-of-arrays:
template class MathEigen::Jacobi<double,
                                 array<double, 3>&,
                                 array<array<double, 3>, 3 >&,
                                 const array<array<double, 3>, 3 >& >;
// template instantiations of LambdaLanczos
template class MathEigen::LambdaLanczos<double>;

// template instantiations of PEidenDense for pointer->pointer arrays
template class MathEigen::PEigenDense<double, double*, double const*const*>;

// template instantiations of PEidenDense for vector<vector<>>
template class MathEigen::PEigenDense<double,
                                      vector<double>&,
                                      const vector<vector<double> >&>;

// If you plan to use other numeric types (eg floats), add them to this list.



// Special case: 3x3 matrices

int MathEigen::
jacobi3(double const mat[3][3], // the 3x3 matrix you wish to diagonalize
        double *eval,           // store the eigenvalues here
        double evec[3][3],      // store the eigenvectors here...
        bool evec_as_columns,   // ...as rows or columns?
        Jacobi<double,double*,double(*)[3],double const(*)[3]>::SortCriteria
        sort_criteria)
{
  // Create "ecalc3", an instance of the MathEigen::Jacobi class to calculate
  // the eigenvalues oand eigenvectors of matrix "mat".  It requires memory
  // to be allocated to store a copy of the matrix.  To avoid allocating
  // this memory on the heap, we create a fixed size 3x3 array on the stack
  // (and convert it to type double**).
  double mat_cpy[3][3] = { {mat[0][0], mat[0][1], mat[0][2]},
                           {mat[1][0], mat[1][1], mat[1][2]},
                           {mat[2][0], mat[2][1], mat[2][2]} };
  double *M[3] = { &(mat_cpy[0][0]),  &(mat_cpy[1][0]), &(mat_cpy[2][0]) };

  int midx[3];  // (another array which is preallocated to avoid using the heap)

  // variable "ecalc3" does all the work:
  Jacobi<double,double*,double(*)[3],double const(*)[3]> ecalc3(3, M, midx);
  int ierror = ecalc3.Diagonalize(mat, eval, evec, sort_criteria);

  // This will store the eigenvectors in the rows of "evec".
  // Do we want to return the eigenvectors in columns instead of rows?
  if (evec_as_columns)
    for (int i=0; i<3; i++)
      for (int j=i+1; j<3; j++)
        std::swap(evec[i][j], evec[j][i]);  // transpose the evec matrix

  return ierror;
}



int MathEigen::
jacobi3(double const* const* mat, // the 3x3 matrix you wish to diagonalize
        double *eval,             // store the eigenvalues here
        double **evec,            // store the eigenvectors here...
        bool evec_as_columns,     // ...as rows or columns?
        Jacobi<double,double*,double**,double const*const*>::SortCriteria
        sort_criteria)
{
  // Create "ecalc3", an instance of the MathEigen::Jacobi class to calculate
  // the eigenvalues oand eigenvectors of matrix "mat".  It requires memory
  // to be allocated to store a copy of the matrix.  To avoid allocating
  // this memory on the heap, we create a fixed size 3x3 array on the stack
  // (and convert it to type double**).
  double mat_cpy[3][3] = { {mat[0][0], mat[0][1], mat[0][2]},
                           {mat[1][0], mat[1][1], mat[1][2]},
                           {mat[2][0], mat[2][1], mat[2][2]} };
  double *M[3] = { &(mat_cpy[0][0]),  &(mat_cpy[1][0]), &(mat_cpy[2][0]) };

  int midx[3];  // (another array which is preallocated to avoid using the heap)

  // variable "ecalc3" does all the work:
  Jacobi<double,double*,double**,double const*const*> ecalc3(3, M, midx);
  int ierror = ecalc3.Diagonalize(mat, eval, evec, sort_criteria);

  // This will store the eigenvectors in the rows of "evec".
  // Do we want to return the eigenvectors in columns instead of rows?
  if (evec_as_columns)
    for (int i=0; i<3; i++)
      for (int j=i+1; j<3; j++)
        std::swap(evec[i][j], evec[j][i]);  // transpose the evec matrix

  return ierror;
}
