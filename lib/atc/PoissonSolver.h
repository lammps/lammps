#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

// ATC includes
#include "Array2D.h"
#include "LinearSolver.h"
#include "PhysicsModel.h"
#include "NonLinearSolver.h"

// other includes
#include <vector>
#include <map>

namespace ATC {

// Forward class declarations
class ATC_Coupling;
class FE_Engine;
class PrescribedDataManager;
class PhysicsModel;

/**
 *  @class PoissonSolver
 *  @brief a class to solve the Poisson equation of electro-statics
 */
class PoissonSolver {

 public:

  /** Constructor */
  PoissonSolver(
    const FieldName fieldName,
    const PhysicsModel * physicsModel,
    const FE_Engine * feEngine,
    const PrescribedDataManager * prescribedDataMgr,
    /*const*/ ATC_Coupling * atc,
    const Array2D<bool> & rhsMask,
    const int solverType = LinearSolver::DIRECT_SOLVE,
    bool parallel = false
  );

  /** Destructor */
  ~PoissonSolver();

  /** parser */
  bool modify(int narg, char **arg);

  /** initialize */
  void initialize(void);

  /** compute_rhs */
  void compute_rhs(const FIELDS & fields, FIELDS & rhs);

  /** solve */
  bool solve(FIELDS & fields, FIELDS & rhs); // rhs is a return value
  bool solve(DENS_MAT & field, const DENS_MAT & rhs); // rhs is input

  /** greens function */
  void greens_function(const int I, VECTOR & inv_stiffness_I) const
  {
    solver_->greens_function(I,inv_stiffness_I);
  }

  /** access to penalty coefficient */
  double penalty_coefficient() const
  {
    return solver_->penalty_coefficient();
  }

  /** set tolerance for underlying solver */
  void set_tolerance(double tol)
  {
    solverTol_ = tol;
  }

  /** set max iterations for underlying solver */
  void set_max_iterations(int maxIter)
  {
    solverMaxIter_ = maxIter;
  }

 protected:

  /** set atomic charges from electron model */
  void set_charges(FIELDS & fields);

  /** Pointer to ATC_Tranfer */
  ATC_Coupling * atc_;

  /** Pointer to FE_Engine */
  const FE_Engine * feEngine_;

  /** Pointer to PrescribedDataManager */
  const PrescribedDataManager * prescribedDataMgr_;

  /** Pointer to FE_Engine */
  const PhysicsModel * physicsModel_;

  /** field to solve for */
  FieldName fieldName_;

  /** number of nodes */
  int nNodes_;

  /** degree-of-freedom */
  int dof_;

  /** righthand side vector */
  Array2D<bool> rhsMask_;

  /** linear Poisson equation */
  bool linear_;

  /** solver */
  LinearSolver * solver_;
  NonLinearSolver *solverNL_;
  PhysicsModelTangentOperator * tangent_;
  int solverType_;
  double solverTol_;
  int solverMaxIter_;

  /** source quadrature */
  IntegrationDomainType integrationType_;

  /** stiffness matrix */

  SPAR_MAT stiffness_;

  /** use owned grid */
  bool useOwnGrid_;

  /** run solve in parallel */
  bool parallel_;
};

} // namespace ATC
#endif
