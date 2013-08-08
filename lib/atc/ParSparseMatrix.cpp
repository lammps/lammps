#include "ParSparseMatrix.h"
#include <fstream>
#ifdef TIMING_ON
double time_diff(timespec &start, timespec &end)
{
  return (double)(1e9 * (end.tv_sec - start.tv_sec) +
                  end.tv_nsec - start.tv_nsec) / 1e9;
}
#endif
using namespace MPI_Wrappers;
using namespace std;
namespace ATC_matrix {

// All the same constructors as for SparseMatrix
ParSparseMatrix<double>::ParSparseMatrix(MPI_Comm comm, INDEX rows, INDEX cols)
  : SparseMatrix<double>(rows, cols), _comm(comm) {}

ParSparseMatrix<double>::ParSparseMatrix(MPI_Comm comm,
    const SparseMatrix<double> &c) :
    SparseMatrix<double>(c), _comm(comm) {}

ParSparseMatrix<double>::ParSparseMatrix(MPI_Comm comm,
    INDEX* rows, INDEX* cols, double* vals, INDEX size,
    INDEX nRows, INDEX nCols, INDEX nRowsCRS)
  : SparseMatrix<double>(rows, cols, vals, size, nRows,
    nCols, nRowsCRS), _comm(comm) {}

//============================================================
void ParSparseMatrix<double>::MultMv(const Vector<double>& v,
    DenseVector<double>& c) const
{
  int numProcs = MPI_Wrappers::size(_comm);
#ifdef DISABLE_PAR_HEURISTICS
  // Use much more lenient heuristics to exercise parallel code
  if (numProcs == 1 ||  _size < 300) {
#else  
  // These are simple heuristics to perform multiplication in serial if
  //   parallel will be slower. They were determined experimentally.
  if ( numProcs == 1 ||
      (_size < 50000  || _size > 10000000) ||
     ((_size < 150000 || _size > 5000000) && numProcs > 8) ||
     ((_size < 500000 || _size > 2500000) && numProcs > 16 ) ||
     (numProcs > 32)) {
#endif
    SparseMatrix<double>::MultMv(v, c);
    return;
  }
 

  SparseMatrix<double>::compress(*this);
  GCK(*this, v, this->nCols() != v.size(), "ParSparseMatrix * Vector")

  SparseMatrix<double> A_local;

  // Split the sparse matrix. partition() takes a ParSparMat, so we cast.
  partition(*static_cast<ParSparseMatrix<double>*>(&A_local));

  // actually do multiplication - end up with partial result vector
  // on each processor
#ifdef TIMING_ON
  timespec before, after;
//  barrier(MPI_COMM_WORLD);
  clock_gettime(CLOCK_MONOTONIC, &before);
#endif
  DenseVector<double> c_local = A_local * v;
#ifdef TIMING_ON
//  barrier(MPI_COMM_WORLD);
  clock_gettime(CLOCK_MONOTONIC, &after);
  cout << "P" << MPI_Wrappers::rank(MPI_COMM_WORLD) << " " << time_diff(before,after) << " mat.vec time\n";
  //LammpsInterface::instance()->all_print((after-before),"mat.vec time");
  barrier(MPI_COMM_WORLD);
#endif

  // destroy A_local intelligently
  static_cast<ParSparseMatrix<double>*>(&A_local)->finalize();

  // Add all the result vectors together on each processor.
#ifdef TIMING_ON
  barrier(MPI_COMM_WORLD);
//barrier(MPI_COMM_WORLD);
  clock_gettime(CLOCK_MONOTONIC, &before);
#endif
  allsum(_comm, c_local.ptr(), c.ptr(), c_local.size());
#ifdef TIMING_ON
//barrier(MPI_COMM_WORLD);
  clock_gettime(CLOCK_MONOTONIC, &after);
  cout << "P" << MPI_Wrappers::rank(MPI_COMM_WORLD) << " " << time_diff(before,after) << " allsum time\n";
  //LammpsInterface::instance()->print_msg_once((after-before),"allsum time");
#endif
}

DenseVector<double> ParSparseMatrix<double>::transMat(
    const Vector<double>& v) const {
  SparseMatrix<double>::compress(*this);
  GCK(*this, v, this->nRows() != v.size(), "ParSparseMatrix transpose * Vector")

  DenseVector<double> c(nCols(), true);

  SparseMatrix<double> A_local;
  partition(*static_cast<ParSparseMatrix<double>*>(&A_local));

  // actually do multiplication - end up with partial result vector
  // on each processor
  DenseVector<double> c_local = A_local.transMat(v);

  static_cast<ParSparseMatrix<double>*>(&A_local)->finalize();

  // Add all the result vectors together on each processor.
  allsum(_comm, c_local.ptr(), c.ptr(), c_local.size());

  return c;
}

void ParSparseMatrix<double>::MultAB(const Matrix<double>& B,
    DenseMatrix<double>& C) const {
  SparseMatrix<double>::compress(*this);
  GCK(*this, B, this->nCols() != B.nRows(), "ParSparseMatrix * Matrix")

  SparseMatrix<double> A_local;
  partition(*static_cast<ParSparseMatrix<double>*>(&A_local));

  // actually do multiplication - end up with partial result matrix
  // on each processor

#ifdef TIMING_ON
  timespec before, after;
  barrier(MPI_COMM_WORLD);
  clock_gettime(CLOCK_MONOTONIC, &before);
#endif
  DenseMatrix<double> C_local = A_local * B;
#ifdef TIMING_ON
  barrier(MPI_COMM_WORLD);
  clock_gettime(CLOCK_MONOTONIC, &after);
  cout << "P" << MPI_Wrappers::rank(MPI_COMM_WORLD) << " " << time_diff(after,before) << " mat.vec time\n";
  //LammpsInterface::instance()->all_print((after-before),"mat.vec time");
#endif

  static_cast<ParSparseMatrix<double>*>(&A_local)->finalize();

  // Add all the result vectors together on each processor.
#ifdef TIMING_ON
  barrier(MPI_COMM_WORLD);
  clock_gettime(CLOCK_MONOTONIC, &before);
#endif
  allsum(_comm, C_local.ptr(), C.ptr(), C_local.size());
#ifdef TIMING_ON
  barrier(MPI_COMM_WORLD);
  clock_gettime(CLOCK_MONOTONIC, &after);
  cout << "P" << MPI_Wrappers::rank(MPI_COMM_WORLD) << " " << time_diff(after,before) << " allsum time\n";
  //LammpsInterface::instance()->print_msg_once((after-before),"allsum time");
#endif
}

DenseMatrix<double> ParSparseMatrix<double>::transMat(
    const DenseMatrix<double>& B) const {
  SparseMatrix<double>::compress(*this);
  GCK(*this, B, this->nRows() != B.nRows(), "ParSparseMatrix transpose * Matrix")

  DenseMatrix<double> C(nCols(), B.nCols(), true);

  SparseMatrix<double> A_local;
  partition(*static_cast<ParSparseMatrix<double>*>(&A_local));

  // actually do multiplication - end up with partial result matrix
  // on each processor
  DenseMatrix<double> C_local = A_local.transMat(B);

  static_cast<ParSparseMatrix<double>*>(&A_local)->finalize();

  // Add all the result vectors together on each processor.
  allsum(_comm, C_local.ptr(), C.ptr(), C_local.size());

  return C;
}

/*
 The two commented-out functions both need to return SparseMatrices. It's hard
 to combine sparse matrices between processors, so this has not yet been completed.

 void ParMultAB(const SparseMatrix<double> &B, SparseMatrix<double> &C) const
 {
 //SparseMatrix<T>::compress(*this);
 GCK(*this, B, this->nCols()!=B.nRows(), "ParSparseMatrix * SparseMatrix")

 ParSparseMatrix<double> A_local(this->_comm);
 this->partition(A_local);

 // actually do multiplication - end up with partial result matrix
 // on each processor

 SparseMatrix<double> C_local = ((SparseMatrix<double>)A_local) * B;

 // destroy newA intelligently
 static_cast<ParSparseMatrix<double>*>(&A_local)->finalize();

 // Add all the result vectors together on each processor.
 sumSparse(C_local, C);
 }*/

DenseMatrix<double> ParSparseMatrix<double>::transMat(
    const SparseMatrix<double>& B) const {
  SparseMatrix<double>::compress(*this);
  GCK(*this, B, this->nRows() != B.nRows(), "ParSparseMatrix transpose * SparseMatrix")

  DenseMatrix<double> C(nCols(), B.nCols(), true);

  SparseMatrix<double> A_local;
  partition(*static_cast<ParSparseMatrix<double>*>(&A_local));

  // actually do multiplication - end up with partial result matrix
  // on each processor
  DenseMatrix<double> C_local = A_local.transMat(B);

  static_cast<ParSparseMatrix<double>*>(&A_local)->finalize();

  // Add all the result vectors together on each processor.
  allsum(_comm, C_local.ptr(), C.ptr(), C_local.size());

  return C;
}

/*void ParMultAB(const DiagonalMatrix<double> &B, SparseMatrix<double> &C) const
 {
 //SparseMatrix<T>::compress(*this);
 GCK(*this, B, this->nCols()!=B.nRows(), "ParSparseMatrix * DiagonalMatrix")

 ParSparseMatrix<double> A_local(this->_comm);
 this->partition(A_local);

 // actually do multiplication - end up with partial result matrix
 // on each processor

 SparseMatrix<double> C_local = ((SparseMatrix<double>)A_local) * B;

 // destroy newA intelligently
 A_local._val = NULL;
 A_local._ja = NULL;

 // Add all the result vectors together on each processor.
 sumSparse(C_local, C);
 }*/

void ParSparseMatrix<double>::partition(
    ParSparseMatrix<double>& A_local) const {
  // create new sparse matrix on each processor, with same size and
  // a disjoint subset of A's elements.
  //
  // Ex: on two processors,
  //
  // |0 1 0|             |0 1 0|               |0 0 0|
  // |2 6 0| splits into |2 0 0| on proc 1 and |0 6 0| on proc 2
  // |0 0 3|             |0 0 0|               |0 0 3|
  //
  // We compute the subproducts individually on each processor, then
  // sum up all the vectors to get our final result.
  //

  // decide which elements will be in each submatrix
  INDEX startIndex = (MPI_Wrappers::rank(_comm) * size()) / MPI_Wrappers::size(_comm);
  INDEX endIndex = ((MPI_Wrappers::rank(_comm) + 1) * size()) / MPI_Wrappers::size(_comm);

  // update number of elements
  A_local._nRows = _nRows;
  A_local._nCols = _nCols;
  A_local._size = endIndex - startIndex;
  A_local._nRowsCRS = _nRowsCRS;
  // use pointer arithmetic to:
  // set newA's _val (to inside A's _val)
  A_local._val = _val + startIndex;
  // set newA's _ja (to inside A's _ja)
  A_local._ja = _ja + startIndex;
  // set newA's _ia (from scratch)
  A_local._ia = new INDEX[nRowsCRS() + 1];
  INDEX numRows = nRowsCRS();
  if (A_local._size > 0) {
    for (INDEX i = 0; i < numRows + 1; i++) {
      A_local._ia[i] = std::min(std::max((_ia[i] - startIndex), 0),
          endIndex - startIndex);
    }
  } else {
    A_local._nRowsCRS = 0;
  }
}

// Prepare an A_local matrix for deletion after it has been loaded with
//   data members from another matrix.
void ParSparseMatrix<double>::finalize() {
  _val = NULL;
  _ja = NULL;
}

void ParSparseMatrix<double>::operator=(const SparseMatrix<double> &source)
{
  copy(source);
}


/*void sumSparse(SparseMatrix<double> &C_local, SparseMatrix<double> &C)
 {
 }*/

}
