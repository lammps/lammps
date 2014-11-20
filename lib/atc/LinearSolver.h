#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

// ATC includes
#include "ATC_TypeDefs.h"
#include "MatrixLibrary.h"
#include "ATC_Error.h"

// other includes
#include <set>
#include <map>

#include "LammpsInterface.h"
#include "CG.h"
#include "GMRES.h"


namespace ATC {

/**
 *  @class LinearSolver
 *  @brief a class to solve a system of linear equations
 *  A x = b subject to a set of constraints { x_i = y_i }
 */
  
class LinearSolver {

 public:
  enum LinearSolveType {
    AUTO_SOLVE=-1, 
    DIRECT_SOLVE=0,
    ITERATIVE_SOLVE,
    ITERATIVE_SOLVE_SYMMETRIC
  };

  enum LinearSolveConstraintHandlingType {
    AUTO_HANDLE_CONSTRAINTS=-1, 
    NO_CONSTRAINTS=0,
    CONDENSE_CONSTRAINTS,
    PENALIZE_CONSTRAINTS
  };

  /** Constructor */
  LinearSolver( // does not assume that A is persistent
    const SPAR_MAT & A,  // lhs matrix "deep" copy
    const BC_SET & bcs,     // constraints
    const int solverType = AUTO_SOLVE, 
    const int bcHandlerType = -1,
    bool parallel = false
  );
  LinearSolver( // assumes A is persistent
    const SPAR_MAT & A,  // lhs matrix "shallow" copy
    const int solverType = AUTO_SOLVE,
    bool parallel = false
  );

  /** Destructor */
  virtual ~LinearSolver() {};

  /** (re)initialize 
      - if bcs are provided the lhs matrix is re-configured 
        for the new constraints 
      - if the class is to be reused with new constraints 
         allow_reinitialization must be called before first solve, etc */
  void allow_reinitialization(void); // depending on method save a copy of A
  void set_homogeneous_bcs(void) { homogeneousBCs_ = true;} // for nonlinear solver, solve for increment
  void initialize(const BC_SET * bcs = NULL);

  /** solve
      - solves A x = b
      - if a "b" is provided it is used as the new rhs */
  bool solve(VECTOR & x, const VECTOR & b);

  /** greens function 
      - returns the solution to a Kronecker delta rhs b = {0 0 .. 1 .. 0 0}
        and with homogeneous constraints {x_i = 0} */
  void greens_function(int I, VECTOR & G_I);

  /** eigensystem
      - returns the e-values & e-vectors for constrained system Ax + v x = 0 
      - if M is provided the eval problem : ( A + v M ) x = 0 is solved*/
  void eigen_system(DENS_MAT & eigenvalues, DENS_MAT & eigenvectors, 
   const DENS_MAT * M = NULL);

  /** access to penalty coefficient
      - if a penalty method is not being used this returns zero */
  double penalty_coefficient(void) const {return penalty_;};

  /** change iterative solver parameters */
  void set_max_iterations(const int maxIter) { 
    if (solverType_ != ITERATIVE_SOLVE && solverType_ != ITERATIVE_SOLVE_SYMMETRIC ) throw ATC_Error("inappropriate parameter set in LinearSolver");
    maxIterations_=maxIter;
  }
  void set_tolerance(const double tol) { tol_=tol;}

  
  /* access to number of unknowns */
  int num_unknowns(void) const 
  { 
    int nUnknowns = nVariables_;
    if (bcs_) { nUnknowns -= bcs_->size(); }
    return nUnknowns;
  }


 protected:
  /** flavors */
  int solverType_;
  int constraintHandlerType_ ;

  /** number of variables = number of rows of matrix */
  int nVariables_;

  /** initialize methods */
  bool initialized_,initializedMatrix_,initializedInverse_;
  bool matrixModified_,allowReinitialization_;
  bool homogeneousBCs_;
  void setup(void);
  void initialize_matrix(void);
  void initialize_inverse(void);
  void initialize_rhs(void);

  /** constraint handling methods to modify the RHS */
  void add_matrix_penalty(void); /** penalty */
  void partition_matrix(void); /** condense */

  /** constraint handling methods to modify the RHS */
  void add_rhs_penalty(void); /** penalty */
  void add_rhs_influence(void); /** condense */

  /** set fixed values */
  void set_fixed_values(VECTOR & x);

  /** constraints container */
  const BC_SET * bcs_; 

  /** rhs vector/s */
  const VECTOR * rhs_; 
  DENS_VEC rhsDense_; // modified
  const VECTOR * b_; // points to appropriate rhs

  /** lhs matrix */
  const SPAR_MAT & matrix_; 
  SPAR_MAT matrixCopy_; // a copy that will be modified by penalty methods

  SPAR_MAT matrixOriginal_; // a copy that is used for re-initialization
  const SPAR_MAT * matrixSparse_; // points to matrix_ or matrixCopy_

  DENS_MAT matrixDense_; // a dense copy for lapack 

  /** partitioned matrix - condense constraints */
  DENS_MAT matrixFreeFree_, matrixFreeFixed_;

  /** maps for free and fixed variables for partitioned matrix - condense */
  std::set<int> freeSet_, fixedSet_; 
  std::map<int,int> freeGlobalToCondensedMap_;

  /** inverse matrix matrix - direct solve */

  DENS_MAT matrixInverse_;

  /** pre-conditioner diagonal of the matrix - iterative solve */
  DIAG_MAT matrixDiagonal_;

  /** penalty coefficient - penalty constraints */
  double penalty_;

  /** max iterations - iterative solve */
  int maxIterations_;

  /** max restarts - GMRES solve */
  int maxRestarts_;

  /** tolerance - iterative solve */
  double tol_;

  /** run solve in parallel */
  bool parallel_;
};

} // namespace ATC

#endif
