// Header file for this class
#include "LinearSolver.h"
#include <sstream>
using std::stringstream;
using std::set;







namespace ATC {
const double  kPenalty = 1.0e4;
const double  kTol     = 1.0e-8;
const int kMaxDirect   = 1000;

// ====================================================================
//  LinearSolver
// ====================================================================
LinearSolver::LinearSolver(
    const SPAR_MAT & A,
    const BC_SET & bcs,
    const int solverType,
    const int constraintHandlerType,
    bool parallel
)
  : solverType_(solverType),
    constraintHandlerType_(constraintHandlerType),
    nVariables_(0),
    initialized_(false),
    initializedMatrix_(false),
    initializedInverse_(false),
    matrixModified_(false),
    allowReinitialization_(false),
    homogeneousBCs_(false),
    bcs_(&bcs),
    rhs_(nullptr),
    rhsDense_(),
    b_(nullptr),
    matrix_(A),
    matrixDense_(),
    matrixFreeFree_(), matrixFreeFixed_(),matrixInverse_(),
    penalty_(1),
    maxIterations_(0), maxRestarts_(0), tol_(0),
    parallel_(parallel)
{
  // deep copy
  matrixCopy_ = A;
  matrixSparse_ = &matrixCopy_;
  setup();
}

LinearSolver::LinearSolver(
    const SPAR_MAT & A,
    const int solverType,
    bool parallel
)
  : solverType_(solverType),
    constraintHandlerType_(NO_CONSTRAINTS),
    nVariables_(0),
    initialized_(false),
    initializedMatrix_(true),
    initializedInverse_(false),
    matrixModified_(false),
    allowReinitialization_(false),
    homogeneousBCs_(false),
    bcs_(nullptr), // null implies no constraints will be added later
    rhs_(nullptr),
    rhsDense_(), b_(nullptr),
    matrix_(A),
    matrixDense_(),
    matrixFreeFree_(), matrixFreeFixed_(),matrixInverse_(),
    penalty_(1),
    maxIterations_(0), maxRestarts_(0), tol_(0),
    parallel_(parallel)
{
  // shallow copy
  matrixSparse_ = &A;
  setup();
}


// --------------------------------------------------------------------
//  Setup
// --------------------------------------------------------------------
void LinearSolver::setup()
{
  tol_     = kTol;
  nVariables_ = matrix_.nRows();
  maxIterations_=2*nVariables_;
  maxRestarts_=nVariables_;


  // switch method based on size
  if (solverType_ < 0) {
    if (nVariables_ > kMaxDirect ) {
      solverType_ = ITERATIVE_SOLVE_SYMMETRIC;
      constraintHandlerType_ = PENALIZE_CONSTRAINTS;
    }
    else {
      solverType_ = DIRECT_SOLVE;
    }
  }
  if (constraintHandlerType_ < 0) {
    constraintHandlerType_ = PENALIZE_CONSTRAINTS;
    if (solverType_ == DIRECT_SOLVE) constraintHandlerType_ = CONDENSE_CONSTRAINTS;
  }
  if ( solverType_ == DIRECT_SOLVE && constraintHandlerType_ == CONDENSE_CONSTRAINTS ) allowReinitialization_ = true;
  if ( solverType_ == ITERATIVE_SOLVE_SYMMETRIC && constraintHandlerType_ == CONDENSE_CONSTRAINTS )  { throw ATC_Error("LinearSolver::unimplemented method"); }
}

// --------------------------------------------------------------------
//  Initialize
// --------------------------------------------------------------------
void LinearSolver::allow_reinitialization()
{
  if (constraintHandlerType_ == PENALIZE_CONSTRAINTS) {
    if (matrixModified_ ) throw ATC_Error("LinearSolver: can't allow reinitialization after matrix has been modified");
    matrixOriginal_ = *matrixSparse_;
  }
  allowReinitialization_  = true;
}

void LinearSolver::initialize(const BC_SET * bcs)
{
  if (bcs) {
    if (! allowReinitialization_ ) throw ATC_Error("LinearSolver: reinitialization not allowed");
   //if (! bcs_ ) throw ATC_Error("LinearSolver: adding constraints after constructing without constraints is not allowed");
    // shallow --> deep copy
    if (! bcs_ ) { // constraintHandlerType_ == NO_CONSTRAINTS
      if (matrixModified_) {
        throw ATC_Error("LinearSolver: adding constraints after constructing without constraints is not allowed if matrix has been modified");
      }
      else {
        matrixCopy_ = *matrixSparse_;
        matrixSparse_ = &matrixCopy_;
        constraintHandlerType_ = -1;
        setup();
      }
    }
    bcs_ = bcs;
    initializedMatrix_ = false;
    initializedInverse_   = false;
    if (matrixModified_) {
      matrixCopy_ = matrixOriginal_;
      matrixSparse_ = &matrixCopy_;
    }
  }
  initialize_matrix();
  initialize_inverse();
  initialize_rhs();

  initialized_ = true;
}

// --------------------------------------------------------------------
//  initialize_matrix
// --------------------------------------------------------------------
void LinearSolver::initialize_matrix()
{
  if ( initializedMatrix_ ) return;
  if       (constraintHandlerType_ == PENALIZE_CONSTRAINTS) {
    add_matrix_penalty();
  }
  else if  (constraintHandlerType_ == CONDENSE_CONSTRAINTS) {
    partition_matrix();
  }
  initializedMatrix_ = true;
}

// --------------------------------------------------------------------
//  initialize_inverse
// --------------------------------------------------------------------
void LinearSolver::initialize_inverse()
{
  if ( initializedInverse_ ) return;
  if (solverType_ == ITERATIVE_SOLVE_SYMMETRIC
   || solverType_ == ITERATIVE_SOLVE ) {
    matrixDiagonal_ = matrixSparse_->diag(); // preconditioner
  }
  else { // DIRECT_SOLVE
    if (constraintHandlerType_ == CONDENSE_CONSTRAINTS) {
      if( num_unknowns() > 0 ) {
        matrixInverse_ = inv(matrixFreeFree_);
      }
    }
    else { // NO_CONSTRAINTS || PENALIZE_CONSTRAINTS
      matrixDense_ = matrixSparse_->dense_copy(); // need dense for lapack
      matrixInverse_ = inv(matrixDense_);
    }
  }
  initializedInverse_ = true;
}

// --------------------------------------------------------------------
//  initialize_rhs
// --------------------------------------------------------------------
void LinearSolver::initialize_rhs()
{
  if (! rhs_ ) return;
  if (! bcs_ ) {
    b_ = rhs_;
    return;
  }
  if      (constraintHandlerType_ == PENALIZE_CONSTRAINTS) {
    add_rhs_penalty();
  }
  else if (constraintHandlerType_ == CONDENSE_CONSTRAINTS) {
    add_rhs_influence();
  }
}

// --------------------------------------------------------------------
// add matrix penalty
// - change matrix for Dirichlet conditions: add penalty
// --------------------------------------------------------------------
void LinearSolver::add_matrix_penalty()
{
  penalty_ = kPenalty; // relative to matrix diagonal
  SPAR_MAT & A = matrixCopy_;
  penalty_ *= (A.diag()).maxabs();
  BC_SET::const_iterator itr;
  for (itr = bcs_->begin(); itr != bcs_->end(); itr++) {
    int i = itr->first;
    A.add(i,i,penalty_); // modifies matrix
  }
  A.compress();
  matrixModified_ = true;
}

// --------------------------------------------------------------------
// partition matrix
// - partition matrix based on Dirichlet constraints
// --------------------------------------------------------------------
void LinearSolver::partition_matrix()
{
  fixedSet_.clear();
  BC_SET::const_iterator itr;
  for (itr = bcs_->begin(); itr != bcs_->end(); itr++) {
    int i = itr->first;
    fixedSet_.insert(i);
  }
  freeSet_.clear();
  freeGlobalToCondensedMap_.clear();
  int j = 0; // local index
  for (int i = 0; i < nVariables_;  i++) {
    if (fixedSet_.find(i) == fixedSet_.end() ) {
       freeSet_.insert(i);
       freeGlobalToCondensedMap_[i] = j++;
    }
  }

  if (matrixDense_.nRows() == 0) matrixDense_ =matrixSparse_->dense_copy();
  DENS_MAT & K = matrixDense_;
  K.row_partition(freeSet_,matrixFreeFree_,matrixFreeFixed_);
}

// --------------------------------------------------------------------
//  add_rhs_penalty
// --------------------------------------------------------------------
void LinearSolver::add_rhs_penalty()
{

  // deep copy
  VECTOR & b = rhsDense_;

  const VECTOR & r = *rhs_;
  int size = r.nRows();
  b.reset(size);
  for (int i = 0; i < size;  i++) {
    b(i) = r(i);
  }

  if ( ! homogeneousBCs_ ){
    BC_SET::const_iterator itr;
    for (itr = bcs_->begin(); itr != bcs_->end(); itr++) {
      int i = itr->first;
      double v = itr->second;
      b(i) += penalty_ * v;
    }
  }
  b_ = &rhsDense_;
}

// --------------------------------------------------------------------
//  add_rhs_influence
// --------------------------------------------------------------------
void LinearSolver::add_rhs_influence()
{
  if (! initializedMatrix_ ) partition_matrix();

  // rhs = rhs + K_free,fixed * x_fixed
  int nbcs = bcs_->size();
  if (nbcs == 0) { // no bcs to handle
    b_ = rhs_;
  }
  else {
    DENS_VEC & b = rhsDense_;
    if ( ! homogeneousBCs_ ){
      DENS_VEC xFixed(nbcs);
      BC_SET::const_iterator itr;
      int i = 0;
      for (itr = bcs_->begin(); itr != bcs_->end(); itr++,i++) {
        double v = itr->second;
        xFixed(i,0) = -v;
      }
      b = matrixFreeFixed_*xFixed; // matrix and bcs have same ordering
    }
    else {
      b.reset(matrixFreeFixed_.nRows());
    }
    const VECTOR & r = *rhs_;
    set<int>::const_iterator iter;
    int i = 0;
    for (iter = freeSet_.begin(); iter != freeSet_.end(); iter++,i++) {
      b(i) += r(*iter);
    }
    b_ = &rhsDense_;
  }
}

// --------------------------------------------------------------------
//  set fixed values
//  - {x_i = y_i}
// --------------------------------------------------------------------
void LinearSolver::set_fixed_values(VECTOR & X)
{
  BC_SET::const_iterator itr;
  for (itr = bcs_->begin(); itr != bcs_->end(); itr++) {
    int i = itr->first;
    double v = 0;
    if ( ! homogeneousBCs_ ) v = itr->second;
    X(i) = v;
  }
}

// --------------------------------------------------------------------
//  Eigensystem
// --------------------------------------------------------------------
// calls lapack

void LinearSolver::eigen_system( DENS_MAT & eigenvalues, DENS_MAT & eigenvectors, const DENS_MAT * M) /* const */
{
  initialize_matrix(); // no inverse needed
  const DENS_MAT * Kp = nullptr;
  const DENS_MAT * Mp =M;
  DENS_MAT MM;
  DENS_MAT KM;
  if (constraintHandlerType_ == CONDENSE_CONSTRAINTS) {
    Kp = &matrixFreeFree_;
    if (M) {
      DENS_MAT MfreeFixed; // not used
      M->row_partition(freeSet_,MM,MfreeFixed);
      Mp = &MM;
    }
  }
  else {
    if (matrixDense_.nRows() == 0) matrixDense_ =matrixSparse_->dense_copy();
    Kp = &matrixDense_;
  }
  if (!M) {
    MM.identity(Kp->nRows());
    Mp = &MM;
  }

  DENS_MAT eVecs, eVals;
  eVecs = eigensystem(*Kp,*Mp,eVals);
  eigenvalues.reset(nVariables_,1);
  eigenvectors.reset(nVariables_,nVariables_);
  set<int>::const_iterator itr;

  for (int i = 0; i < Kp->nRows(); i++) { // ordering is by energy not node
    eigenvalues(i,0) = eVals(i,0);
    int j = 0;
    for (itr = freeSet_.begin(); itr != freeSet_.end(); itr++,j++) {
      int jj = *itr;
      eigenvectors(jj,i) = eVecs(j,i); // transpose
    }
  }
}

// --------------------------------------------------------------------
// solve
// - solves A x = b
// - if a "b" is provided it is used as the new rhs
// --------------------------------------------------------------------

bool LinearSolver::solve(VECTOR & x, const VECTOR & b)
{
  SPAR_MAT * A = nullptr;

  rhs_ = &b;
  initialized_ = false;
  initialize();
  if (num_unknowns() == 0) {
    set_fixed_values(x);
    return true;
  }
  const VECTOR & r = *b_;
  if (solverType_ == ITERATIVE_SOLVE_SYMMETRIC) {



    if (parallel_) {
      A = new PAR_SPAR_MAT(LammpsInterface::instance()->world(), *matrixSparse_);
    }
    else {
      A = new SPAR_MAT(*matrixSparse_);
    }
    DIAG_MAT & PC = matrixDiagonal_;
    int niter = maxIterations_;
    double tol = tol_;
    int convergence = CG(*A, x, r, PC, niter, tol);// CG changes niter, tol
    if (convergence>0) {
       stringstream ss;
       ss << "CG solve did not converge,";
       ss << " iterations: " << niter;
       ss << " residual: " << tol;
       throw ATC_Error(ss.str());
    }
  }
  else if (solverType_ == ITERATIVE_SOLVE) {
    if (parallel_) {
      A = new PAR_SPAR_MAT(LammpsInterface::instance()->world(), *matrixSparse_);
    }
    else {
      A = new SPAR_MAT(*matrixSparse_);
    }
    const DIAG_MAT & PC = matrixDiagonal_;
    int iterations = maxIterations_;
    int restarts = maxRestarts_;
    double tol = tol_;
    DENS_MAT H(maxRestarts_+1, maxRestarts_);
    DENS_VEC xx(nVariables_);
    DENS_VEC bb;
    bb = b;
    int convergence = GMRES(*A, xx, bb, PC, H, restarts, iterations, tol);
    if (convergence>0) {
       stringstream ss;
       ss << "GMRES greens_function solve did not converge,";
       ss << " iterations: " << iterations;
       ss << " residual: " << tol;
       throw ATC_Error(ss.str());
    }
    x.copy(xx.ptr(),xx.nRows());
  }
  else { // DIRECT_SOLVE
    const DENS_MAT & invA = matrixInverse_;
    if (constraintHandlerType_ == CONDENSE_CONSTRAINTS) {
      DENS_MAT xx = invA*r;
      int i = 0;
      set<int>::const_iterator itr;
      for (itr = freeSet_.begin(); itr != freeSet_.end(); itr++,i++) {
        int ii = *itr;
        x(ii) = xx(i,0);
      }
      set_fixed_values(x);
    }
    else {

      DENS_VEC xx = invA*r;
      for (int i = 0; i < xx.nRows(); i++) {
        x(i) = xx(i);
      }
    }
  }
  delete A;


  return true;
}

// --------------------------------------------------------------------
// greens function
// - returns the solution to a Kronecker delta rhs b = {0 0 .. 1 .. 0 0}
//   and with homogeneous constraints {x_i = 0}
// --------------------------------------------------------------------

void LinearSolver::greens_function(int I, VECTOR & G_I)
{
  SPAR_MAT * A = nullptr;



  initialize_matrix();
  initialize_inverse();
  G_I.reset(nVariables_);
  VECTOR & x = G_I;

  if (solverType_ == ITERATIVE_SOLVE_SYMMETRIC) {
    DENS_VEC b(nVariables_); b = 0.0; b(I) = 1.0;
    if (parallel_) {
      A = new PAR_SPAR_MAT(LammpsInterface::instance()->world(), *matrixSparse_);
    }
    else {
      A = new SPAR_MAT(*matrixSparse_);
    }
    const DIAG_MAT & PC = matrixDiagonal_;
    int niter = maxIterations_;
    double tol = tol_;
    int convergence = CG(*A, x, b, PC, niter, tol);
    if (convergence>0) {
       stringstream ss;
       ss << "CG greens_function solve did not converge,";
       ss << " iterations: " << niter;
       ss << " residual: " << tol;
       throw ATC_Error(ss.str());
    }
  }
  else if (solverType_ == ITERATIVE_SOLVE) {
    DENS_VEC b(nVariables_); b = 0.0; b(I) = 1.0;
    //  VECTOR & bb = b;
    if (parallel_) {
      A = new PAR_SPAR_MAT(LammpsInterface::instance()->world(), *matrixSparse_);
    }
    else {
      A = new SPAR_MAT(*matrixSparse_);
    }
    //  const DENS_MAT A = matrixSparse_->dense_copy();
    const DIAG_MAT & PC = matrixDiagonal_;
    int iterations = maxIterations_;
    int restarts = maxRestarts_;
    double tol = tol_;
    DENS_MAT H(maxRestarts_+1, maxRestarts_);
    DENS_VEC xx(nVariables_);
    int convergence = GMRES(*A, xx, b, PC, H, restarts, iterations, tol);
    if (convergence>0) {
       stringstream ss;
       ss << "GMRES greens_function solve did not converge,";
       ss << " iterations: " << iterations;
       ss << " residual: " << tol;
       throw ATC_Error(ss.str());
    }
    x.copy(xx.ptr(),xx.nRows());
  }
  else {
    const DENS_MAT & invA = matrixInverse_;
    if (constraintHandlerType_ == CONDENSE_CONSTRAINTS) {
      set<int>::const_iterator itr;
      for (itr = fixedSet_.begin(); itr != fixedSet_.end(); itr++) {
        int ii = *itr;
        x(ii) = 0;
      }
      itr = freeSet_.find(I);
      if (itr !=freeSet_.end() ) {
        int j = freeGlobalToCondensedMap_[I];
        int i = 0;
        for (itr = freeSet_.begin(); itr != freeSet_.end(); itr++,i++) {
          int ii = *itr;
          x(ii) = invA(j,i);
        }
      }
    }
    else {
      for (int i = 0; i < nVariables_; ++i) x(i) = invA(I,i);
    }
  }

  delete A;
}


} // namespace ATC
