#ifndef IMPLICIT_SOLVE_OPERATOR_H
#define IMPLICIT_SOLVE_OPERATOR_H

// ATC includes
#include "Array2D.h"
#include "MatrixLibrary.h"
#include "PhysicsModel.h"

// other includes
#include <vector>
#include <map>

namespace ATC {

// Forward class declarations
class ATC_Coupling;
class FE_Engine;

/**
 *  @class ImplicitSolveOperator
 *  @brief Helper class to compute matrix-free product for use with IML++ solvers
 */
class ImplicitSolveOperator {
 public:
  /** Constructor */
  ImplicitSolveOperator(double alpha, double dt);
  /** Destructor */
  virtual ~ImplicitSolveOperator() {};
  /** pure virtual operator to compute Ax, for equation Ax=b */
  virtual DENS_VEC operator * (const DENS_VEC &x) const;
  /** pure virtual method to return the rhs vector b */
  virtual void R(const DENS_VEC &f, DENS_VEC &v) const = 0;
  /** pure virtual method to return the rhs vector b */
  virtual DENS_VEC r() const;
  /** pure virtual method to return preconditioner */
  virtual DiagonalMatrix<double> preconditioner() const;
  /** finalize */
  virtual void solution(const DENS_MAT & dx, DENS_MAT &x) const = 0;
  virtual void rhs(const DENS_MAT & r, DENS_MAT &rhs) const = 0;
 protected:
  int n_,dof_;
  DENS_VEC x0_;             // condensed previous
  mutable DENS_VEC x_;      // condensed current
  DENS_VEC R0_;        // condensed previous
  mutable DENS_VEC R_; // condensed current
  double dt_; // timestep
  double alpha_; // implicit/explicit parameter (0 -> explicit, else implicit)
  double epsilon0_; // small parameter to compute increment
};

/**
 *  @class FieldImplicitSolveOperator
 *  @brief Class to perform A*x operation for electron temperature solution
 */
class FieldImplicitSolveOperator : public ImplicitSolveOperator {
 public:
  /** Constructor */
  FieldImplicitSolveOperator(ATC_Coupling * atc,
                             FIELDS & fields,
                             const FieldName f,
                             const Array2D< bool > & rhsMask,
                             const PhysicsModel * physicsModel,
                             double simTime, double dt, double alpha = 0.5);
  /** Destructor */
  virtual ~FieldImplicitSolveOperator() {};
  virtual void R(const DENS_VEC &f, DENS_VEC &v) const;
  virtual void solution(const DENS_MAT & dx, DENS_MAT &x) const;
  virtual void rhs(const DENS_MAT & r, DENS_MAT &rhs) const;
 protected:
  void to_all(const VECTOR& free, MATRIX& all) const;
  void to_free(const MATRIX& all, VECTOR& free) const;
  FieldName fieldName_; // field name of ODE to solve
  ATC_Coupling * atc_; /** Pointer to atc */
  const PhysicsModel * physicsModel_; /** Pointer to PhysicsModel */
  FIELDS & fields0_; // Reference to current fields (passed in constructor)
  // Local fields
  mutable FIELDS fields_; // full
  // set of fixed nodes
  std::map<int,int> freeNodes_;
  // masks
  Array2D<bool> rhsMask_;
  Array<FieldName> massMask_;
  // right hand side of dot x = R(x)
  mutable FIELDS rhs_; // full
  // simulation time
  double time_;
};

} // namespace ATC

#endif
