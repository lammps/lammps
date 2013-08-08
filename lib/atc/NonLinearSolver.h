#ifndef NON_LINEAR_SOLVER_H
#define NON_LINEAR_SOLVER_H

// ATC includes
#include "ATC_TypeDefs.h"
#include "MatrixLibrary.h"
#include "ATC_Error.h"

// other includes
#include <set>
#include <map>
using std::pair;

namespace ATC {

/**
 *  @class TangentOperator
 *  @brief an adaptor to allow NonLinearSolver to work with a generic function
 */
class TangentOperator {
 public:
  TangentOperator(){};
  virtual ~TangentOperator(){};
  virtual void function(const VECTOR & x, DENS_VEC & f) {}; // =0;
  virtual void tangent(const VECTOR & x, DENS_VEC & f, MATRIX & dfdx) =0;
  //virtual void function(const VECTOR & x, VECTOR & f) {}; // =0;
  //virtual void tangent(const VECTOR & x, VECTOR & f, MATRIX & dfdx) {}; // =0;
};

/**
 *  @class NonLinearSolver
 *  @brief a class to solve a system of non-linear equations
 *  f(x) = 0 
 */


class NonLinearSolver {

 public:
  enum NonLinearSolveType {
    NEWTON_RAPHSON=0,
  };

  /** Constructor */
  NonLinearSolver( 
    TangentOperator * f, // provides f and f' at x, pointer for polymorphism
    const BC_SET * bcs = NULL,
    const int dof = 0,
    bool parallel = false
  );

  /** Destructor */
  virtual ~NonLinearSolver() {};

  /** residual norm */
  double residual_norm(VECTOR & x); 

  /** solve */
  bool solve(VECTOR & x); // incoming: initial guess, outgoing: solution

  /** line search */
  bool line_search(VECTOR & x);

  
  /** access to current state */
  DENS_MAT & tangent()  { return A_;}
  DENS_VEC & residual() { return r_;}


  /** change solver parameters */
  void set_max_iterations(const int maxIter) { maxIterations_=maxIter; }
  void set_residual_tolerance(const double tol) { tol_=tol;}
  void set_solution_tolerance(const double tol) { tolx_=tol;}

 protected:
  /** function & tangent */
  TangentOperator * f_;
  DENS_VEC r_;
  DENS_MAT A_;
  DENS_VEC dx_;

  /** equality constraints */
  const BC_SET * bcs_;

  /** degree of freedom */
  int dof_;

  /** flavors */
  int solverType_;

  /** state */ 
  double rNorm0_, rNorm0P_, rNorm_, rNormP_;

  /** parameters & tolerances */
  double tol_; // tolerance on f
  double tolx_; // tolerance on dx
  double tol0_; // tolerance on initial f
  int maxIterations_;

  /** run solve in parallel */
  bool parallel_;
};

} // namespace ATC

#endif
