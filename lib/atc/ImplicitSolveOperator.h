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
  ImplicitSolveOperator(ATC_Coupling * atc,
                        /*const*/ FE_Engine * feEngine,
                        const PhysicsModel * physicsModel);

  /** Destructor */
  virtual ~ImplicitSolveOperator() {};

  /** pure virtual operator to compute Ax, for equation Ax=b */
  virtual DENS_VEC operator * (DENS_VEC x) const = 0;

  /** pure virtual method to return the rhs vector b */
  virtual DENS_VEC rhs() = 0;

  /** pure virtual method to return preconditioner */
  virtual DiagonalMatrix<double> preconditioner(FIELDS & fields) = 0;

 protected:

  /** Pointer to atc */
  ATC_Coupling * atc_;

  /** Pointer to FE_Engine */
  /*const*/ FE_Engine * feEngine_;

  /** Pointer to PhysicsModel */
  const PhysicsModel * physicsModel_;

};

/**
 *  @class FieldImplicitSolveOperator
 *  @brief Class to perform A*x operation for electron temperature solution
 */
class FieldImplicitSolveOperator : public ImplicitSolveOperator {

 public:

  /** Constructor */
  FieldImplicitSolveOperator(ATC_Coupling * atc,
                             /*const*/ FE_Engine * fe_Engine,
                             FIELDS & fields,
                             const FieldName electronField,
                             const Array2D< bool > & rhsMask,
                             const PhysicsModel * physicsModel,
                             double simTime,
                             double dt,
                             double alpha);

  /** Destructor */
  virtual ~FieldImplicitSolveOperator() {};

  /** operator to compute A*x for the electron temperature equation */
  virtual DENS_VEC operator * (DENS_VEC x) const;

  /** method to return the rhs vector b */
  virtual DENS_VEC rhs();

  /** method to return preconditioner (identity matrix) */
  virtual DIAG_MAT preconditioner(FIELDS & fields);

 protected:
  // field name of ODE to solve
  FieldName fieldName_;

  // Reference to current fields (passed in constructor)
  FIELDS & fields_;

  // Local fields
  mutable FIELDS fieldsNp1_;

  // Vector to hold current temperature
  DENS_VEC TnVect_;

  // Old and new RHS maps (not including inverse mass)
  FIELDS RnMap_;
  mutable FIELDS RnpMap_;

  // Matrices/vectors to hold electron temperature components of RHS
  // vectors (including inverse mass)
  DENS_VEC RnVect_;
  mutable DENS_VEC RnpVect_;

  Array2D<bool> rhsMask_;
  Array<FieldName> massMask_;

  // simulation time
  double time_;

  // timestep
  double dt_;

  // implicit/explicit parameter (0 -> explicit, else implicit)
  double alpha_;

  // small parameter to compute increment 
  double epsilon0_;

};

} // namespace ATC

#endif
