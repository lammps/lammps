
#ifndef SRC_SOT_MLIAP_SO3_MATH_H_
#define SRC_SOT_MLIAP_SO3_MATH_H_

#include<cassert>

#ifndef _MATRIX_ALLOC_JPD_H
#define _MATRIX_ALLOC_JPD_H

namespace matrix_alloc_jpd {

/// @brief  Allocate a 2-dimensional array.  (Uses row-major order.)
template<typename Entry>
void Alloc2D(size_t nrows, size_t ncols, Entry ***paaX);

/// @brief  Deallocate arrays that were created using Alloc2D().
template<typename Entry>
void Dealloc2D(Entry ***paaX);

// ---- implementation ----

template<typename Entry>
void Alloc2D(size_t nrows, size_t ncols, Entry ***paaX)
{
  assert(paaX);
  *paaX = new Entry* [nrows];
  (*paaX)[0] = new Entry [nrows * ncols];
  for(size_t iy=0; iy<nrows; iy++)
    (*paaX)[iy] = (*paaX)[0] + iy*ncols;
  // The caller can access the contents using (*paaX)[i][j]
}

template<typename Entry>
void Dealloc2D(Entry ***paaX)
{
  if (paaX && *paaX) {
    delete [] (*paaX)[0];
    delete [] (*paaX);
    *paaX = nullptr;
  }
}

} // namespace matrix_alloc_jpd

#endif //#ifndef _MATRIX_ALLOC_JPD_H


#ifndef _JACOBI_HPP
#define _JACOBI_HPP

#include <algorithm>
#include <cmath>


namespace jacobi_pd {

using namespace matrix_alloc_jpd;


template<typename Scalar,
         typename Vector,
         typename Matrix,
         typename ConstMatrix=Matrix>
class Jacobi
{
  int n;
  Scalar **M;
  Scalar c;
  Scalar s;
  Scalar t;
  int *max_idx_row;

 public:

  void SetSize(int n);

  Jacobi(int n = 0) {
    Init();
    SetSize(n);
  }

  ~Jacobi() {
    Dealloc();
  }

  typedef enum eSortCriteria {
    DO_NOT_SORT,
    SORT_DECREASING_EVALS,
    SORT_INCREASING_EVALS,
    SORT_DECREASING_ABS_EVALS,
    SORT_INCREASING_ABS_EVALS
  } SortCriteria;

  int Diagonalize(ConstMatrix mat, Vector eval, Matrix evec,
                  SortCriteria sort_criteria=SORT_DECREASING_EVALS,
                  bool calc_evecs=true, int max_num_sweeps = 50);

  int Diagonalize_1DArray(double *mat, Vector eval, double *evec,
                          SortCriteria sort_criteria=
                          SORT_DECREASING_EVALS,
                          bool calc_evecs=true,
                          int max_num_sweeps = 50);

 private:

  void CalcRot(Scalar const *const *M, int i, int j);

  void ApplyRot(Scalar **M, int i, int j);

  void ApplyRotLeft(Matrix E, int i, int j);

  void ApplyRotLeft_1DArray(double *E, int i, int j);

  int MaxEntryRow(Scalar const *const *M, int i) const;

  void MaxEntry(Scalar const *const *M, int& i_max, int& j_max) const;

  void SortRows(Vector v, Matrix M, int n,
                SortCriteria s=SORT_DECREASING_EVALS ) const;

  void Alloc(int N);
  void Init();
  void Dealloc();

  int LU_decompose( double *A, int *permute );
  void LU_substitution( double *A, double *B, int *permute );

 public:

  int invert_matrix( double *Am, double *Aminv );

  Jacobi(const Jacobi<Scalar, Vector, Matrix, ConstMatrix>& source);
  Jacobi(Jacobi<Scalar, Vector, Matrix, ConstMatrix>&& other);
  void swap(Jacobi<Scalar, Vector, Matrix, ConstMatrix> &other);
  Jacobi<Scalar, Vector, Matrix, ConstMatrix>& operator =
    (Jacobi<Scalar, Vector, Matrix, ConstMatrix> source);

}; // class Jacobi


template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
int Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
  invert_matrix( double *Am, double *Aminv )
{
  int i,j;

  int *permute;
  double *dxm, *Amtmp;
  permute=new int[n];
  dxm=new double[n];
  Amtmp=new double[n*n];

  for( i = 0 ; i < n*n ; i++ )
    Amtmp[i] = Am[i];


  if( LU_decompose( Amtmp,  permute ) != 0  )
    return(1);

  for( i = 0; i < n; i++ ) {

    for( j = 0 ; j < n ; j++ )
      dxm[j] = 0. ;

    dxm[i] = 1.;

    LU_substitution( Amtmp,  dxm, permute );

    for( j = 0 ; j < n ; j++ )
      Aminv[j*n+i] = dxm[j];

  }

  delete permute;
  delete dxm;
  delete Amtmp;

  return(0);
}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
int Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
  LU_decompose( double *A, int *permute )
{
  const  double absmin = 1.e-30;
  double *row_norm;
  double  absmax, maxtemp, mintemp;

  int i, j, k, max_row;

  row_norm=new double[n];

  max_row = 0;

  for( i = 0; i < n; i++ ) {
    absmax = 0.;

    for( j = 0; j < n ; j++ ) {

      maxtemp = fabs( A[i*n+j] );

      if( maxtemp > absmax )
        absmax = maxtemp;
    }

    if( absmax == 0. ) {
      fprintf(stderr, "LU_decompose(): row-wise singular matrix!\n");
      return(1);
    }

    row_norm[i] = 1. / absmax ;
  }


  for( j = 0; j < n; j++ ) {

    for( i = 0; i < j; i++ ) {
      for( k = 0; k < i; k++ )
        A[i*n+j] -= A[i*n+k] * A[k*n+j];
    }

    absmax = 0.0;

    for( i = j; i < n; i++ ) {

      for (k = 0; k < j; k++)
        A[i*n+j] -= A[i*n+k] * A[k*n+j];

      maxtemp = fabs(A[i*n+j]) * row_norm[i] ;

      if( maxtemp >= absmax ) {
        absmax = maxtemp;
        max_row = i;
      }

    }

    if( max_row != j ) {

      if( (j == (n-2)) && (A[j*n+j+1] == 0.) )
        max_row = j;

      else {

        for( k = 0; k < n; k++ ) {

          maxtemp       = A[   j*n+k] ;
          A[   j*n+k] = A[max_row*n+k] ;
          A[max_row*n+k] = maxtemp;

        }

        row_norm[max_row] = row_norm[j] ;

      }
    }

    permute[j] = max_row;

    if( A[j*n+j] == 0. )
      A[j*n+j] = absmin;

    if( j != (n-1) ) {
      maxtemp = 1. / A[j*n+j]  ;

      for( i = (j+1) ; i < n; i++ )
        A[i*n+j] *= maxtemp;
    }

  }

  delete row_norm;
  return(0);

}


template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
  LU_substitution( double *A, double *B, int *permute)
{
  int i, j ;
  double tmpvar,tmpvar2;

  for(i = 0; i < n; i++) {

    tmpvar        = B[permute[i]];
    B[permute[i]] = B[    i     ];
    for( j = (i-1); j >= 0 ; j-- )
      tmpvar -=  A[i*n+j] * B[j];

    B[i] = tmpvar;
  }


  for( i = (n-1); i >= 0; i-- ) {
    for( j = (i+1); j < n ; j++ )
      B[i] -=  A[i*n+j] * B[j];

    B[i] /= A[i*n+i] ;
  }

}


// -------------- implementation --------------

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
int Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Diagonalize_1DArray(double *mat, Vector eval, double *evec,
                    SortCriteria sort_criteria,
                    bool calc_evec, int max_num_sweeps)
{
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      M[i][j] = mat[i*n+j];

  if (calc_evec)
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        evec[i*n+j] = (i==j) ? 1.0 : 0.0;

  for (int i = 0; i < n-1; i++)
    max_idx_row[i] = MaxEntryRow(M, i);

  int n_iters;
  int max_num_iters = max_num_sweeps*n*(n-1)/2;
  for (n_iters=1; n_iters <= max_num_iters; n_iters++) {
    int i,j;
    MaxEntry(M, i, j);

    if ((M[i][i] + M[i][j] == M[i][i]) &&
      (M[j][j] + M[i][j] == M[j][j])) {
        M[i][j] = 0.0;
        max_idx_row[i] = MaxEntryRow(M,i);
    }

    if (M[i][j] == 0.0)
      break;

    CalcRot(M, i, j);
    ApplyRot(M, i, j);
    if (calc_evec)
      ApplyRotLeft_1DArray(evec,i,j);

  }

  for (int i = 0; i < n; i++)
    eval[i] = M[i][i];

  if ((n_iters > max_num_iters) && (n>1))
    return 0;

  return n_iters;
}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
int Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
  Diagonalize(ConstMatrix mat,
            Vector eval,
            Matrix evec,
            SortCriteria sort_criteria,
            bool calc_evec,
            int max_num_sweeps)
{
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      M[i][j] = mat[i][j];

  if (calc_evec)
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        evec[i][j] = (i==j) ? 1.0 : 0.0;

  for (int i = 0; i < n-1; i++)
    max_idx_row[i] = MaxEntryRow(M, i);

  // -- Iteration --
  int n_iters;
  int max_num_iters = max_num_sweeps*n*(n-1)/2;
  for (n_iters=1; n_iters <= max_num_iters; n_iters++) {
    int i,j;
    MaxEntry(M, i, j);


    if ((M[i][i] + M[i][j] == M[i][i]) && (M[j][j] + M[i][j] == M[j][j])) {
      M[i][j] = 0.0;
      max_idx_row[i] = MaxEntryRow(M,i);
    }

    if (M[i][j] == 0.0)
      break;

    CalcRot(M, i, j);
    ApplyRot(M, i, j);
    if (calc_evec)
      ApplyRotLeft(evec,i,j);

  }

  for (int i = 0; i < n; i++)
    eval[i] = M[i][i];

  SortRows(eval, evec, n, sort_criteria);

  if ((n_iters > max_num_iters) && (n>1))
    return 0;

  return n_iters;
}


template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
  CalcRot(Scalar const *const *M, int i, int j)
{
  t = 1.0;
  Scalar M_jj_ii = (M[j][j] - M[i][i]);
  if (M_jj_ii != 0.0) {

    Scalar kappa = M_jj_ii;
    t = 0.0;
    Scalar M_ij = M[i][j];
    if (M_ij != 0.0) {
      kappa /= (2.0*M_ij);
      t = 1.0 / (std::sqrt(1 + kappa*kappa) + std::abs(kappa));
      if (kappa < 0.0)
        t = -t;
    }
  }

  c = 1.0 / std::sqrt(1 + t*t);
  s = c*t;

}

template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
ApplyRot(Scalar **M,
         int i,
         int j)
{
  M[i][i] -= t * M[i][j];
  M[j][j] += t * M[i][j];

  M[i][j] = 0.0;

  for (int w=0; w < i; w++) {
    M[i][w] = M[w][i];
    M[w][i] = c*M[w][i] - s*M[w][j];
    if (i == max_idx_row[w])
      max_idx_row[w] = MaxEntryRow(M, w);
    else if (std::abs(M[w][i])>std::abs(M[w][max_idx_row[w]]))
      max_idx_row[w]=i;

  }

  for (int w=i+1; w < j; w++) {
    M[w][i] = M[i][w];
    M[i][w] = c*M[i][w] - s*M[w][j];
  }

  for (int w=j+1; w < n; w++) {
    M[w][i] = M[i][w];
    M[i][w] = c*M[i][w] - s*M[j][w];
  }

  max_idx_row[i] = MaxEntryRow(M, i);

  for (int w=0; w < i; w++) {
    M[w][j] = s*M[i][w] + c*M[w][j];
    if (j == max_idx_row[w])
      max_idx_row[w] = MaxEntryRow(M, w);
    else if (std::abs(M[w][j])>std::abs(M[w][max_idx_row[w]]))
      max_idx_row[w]=j;

  }

  for (int w=i+1; w < j; w++) {
    M[w][j] = s*M[w][i] + c*M[w][j];
    if (j == max_idx_row[w])
      max_idx_row[w] = MaxEntryRow(M, w);
    else if (std::abs(M[w][j])>std::abs(M[w][max_idx_row[w]]))
      max_idx_row[w]=j;
  }

  for (int w=j+1; w < n; w++)
    M[j][w] = s*M[w][i] + c*M[j][w];

  max_idx_row[j] = MaxEntryRow(M, j);

}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
  ApplyRotLeft(Matrix E, int i, int j)
{
  for (int v = 0; v < n; v++) {
    Scalar Eiv = E[i][v];
    E[i][v] = c*E[i][v] - s*E[j][v];
    E[j][v] = s*Eiv + c*E[j][v];
  }
}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
  ApplyRotLeft_1DArray(double *E, int i, int j)
{
  for (int v = 0; v < n; v++) {
    Scalar Eiv = E[i*n+v];
    E[i*n+v] = c*E[i*n+v] - s*E[j*n+v];
    E[j*n+v] = s*Eiv + c*E[j*n+v];
  }
}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
int Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
  MaxEntryRow(Scalar const *const *M, int i) const
{
  int j_max = i+1;
  for(int j = i+2; j < n; j++)
    if (std::abs(M[i][j]) > std::abs(M[i][j_max]))
      j_max = j;
  return j_max;
}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
  MaxEntry(Scalar const *const *M,
           int& i_max, int& j_max) const {

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

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
  SortRows(Vector eval, Matrix evec, int n,
           SortCriteria sort_criteria) const
{
  for (int i = 0; i < n-1; i++) {
    int i_max = i;
    for (int j = i+1; j < n; j++) {
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
    std::swap(eval[i], eval[i_max]);
    for (int k = 0; k < n; k++)
      std::swap(evec[i][k], evec[i_max][k]);
  }
}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::Init()
{
  n = 0;
  M = nullptr;
  max_idx_row = nullptr;
}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::SetSize(int n)
{
  Dealloc();
  Alloc(n);
}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::Alloc(int N)
{
  n = N;
  if (n > 0) {
    max_idx_row = new int[n];
    Alloc2D(n, n, &M);
  }
}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Dealloc() {
  if (max_idx_row) {
    delete [] max_idx_row;
    max_idx_row = nullptr;
  }
  Dealloc2D(&M);
  Init();
}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
  Jacobi(const Jacobi<Scalar, Vector, Matrix, ConstMatrix>& source)
{
  Init();
  SetSize(source.n);

  std::copy(source.max_idx_row,
    source.max_idx_row + n,max_idx_row);
  for (int i = 0; i < n; i++)
    std::copy(source.M[i],source.M[i] + n,M[i]);

}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
  swap(Jacobi<Scalar, Vector, Matrix, ConstMatrix> &other)
{
  std::swap(n, other.n);
  std::swap(max_idx_row, other.max_idx_row);
  std::swap(M, other.M);
}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
  Jacobi(Jacobi<Scalar, Vector, Matrix, ConstMatrix>&& other)
{
  Init();
  this->swap(other);
}

template<typename Scalar,typename Vector,typename Matrix,
  typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>&
  Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
    operator = (Jacobi<Scalar, Vector, Matrix, ConstMatrix>
      source)
{
  this->swap(source);
  return *this;
}

} // namespace jacobi


#endif //#ifndef _JACOBI_HPP


#endif /* SRC_SOT_MLIAP_SO3_MATH_H_ */
