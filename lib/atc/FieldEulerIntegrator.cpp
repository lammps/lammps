#include "FieldEulerIntegrator.h"
#include "ATC_Coupling.h"
#include "FE_Engine.h"
#include "PhysicsModel.h"
#include "GMRES.h"
#include "CG.h"
#include "ImplicitSolveOperator.h"

namespace ATC {

// ====================================================================
// FieldEulerIntegrator
// ====================================================================
FieldEulerIntegrator::FieldEulerIntegrator(
  const FieldName fieldName,
  const PhysicsModel * physicsModel, 
  /*const*/ FE_Engine * feEngine,
  /*const*/ ATC_Coupling * atc,
  const Array2D< bool > & rhsMask  // copy 
)
  : atc_(atc),
    feEngine_(feEngine),
    physicsModel_(physicsModel),
    fieldName_(fieldName),
    rhsMask_(rhsMask)
{
  nNodes_ = feEngine->num_nodes();
}

// ====================================================================
// FieldImplicitIntegrator
// ====================================================================
FieldExplicitEulerIntegrator::FieldExplicitEulerIntegrator(
  const FieldName fieldName,
  const PhysicsModel * physicsModel, 
  /*const*/ FE_Engine * feEngine,
  /*const*/ ATC_Coupling * atc,
  const Array2D< bool > & rhsMask  // copy 
) : FieldEulerIntegrator(fieldName,physicsModel,feEngine,atc,rhsMask)
{
}

// --------------------------------------------------------------------
// update 
// --------------------------------------------------------------------
void FieldExplicitEulerIntegrator::update(const double dt, double time,
  FIELDS & fields, FIELDS & rhs)
{ 
  // write and add update mass matrix to handled time variation
  // update mass matrix to be consistent/lumped, and handle this in apply_inverse_mass_matrix
  atc_->compute_rhs_vector(rhsMask_, fields, rhs,
    FULL_DOMAIN, physicsModel_);
  DENS_MAT & myRhs(rhs[fieldName_].set_quantity());
  atc_->apply_inverse_mass_matrix(myRhs,fieldName_);
  fields[fieldName_] += dt*myRhs;
}

// ====================================================================
// FieldImplicitEulerIntegrator
// ====================================================================
FieldImplicitEulerIntegrator::FieldImplicitEulerIntegrator(
  const FieldName fieldName,
  const PhysicsModel * physicsModel, 
  /*const*/ FE_Engine * feEngine,
  /*const*/ ATC_Coupling * atc,
  const Array2D< bool > & rhsMask,  // copy 
  const double alpha
) : FieldEulerIntegrator(fieldName,physicsModel,feEngine,atc,rhsMask),
  alpha_(alpha),
  dT_(1.0e-6), 
  maxRestarts_(50),
  maxIterations_(1000),
  tol_(1.0e-8)
{
}

// --------------------------------------------------------------------
// update 
// --------------------------------------------------------------------
void FieldImplicitEulerIntegrator::update(const double dt, double time,
  FIELDS & fields, FIELDS & rhs) 
{ // solver handles bcs
  FieldImplicitSolveOperator solver(atc_, feEngine_, 
    fields, fieldName_, rhsMask_, physicsModel_,
    time, dt, alpha_);
  DiagonalMatrix<double> preconditioner = solver.preconditioner(fields);
  DENS_VEC myRhs = solver.rhs();
  DENS_VEC dT(nNodes_); dT = dT_; 
  DENS_MAT H(maxRestarts_+1, maxRestarts_);
  double tol = tol_; // tol returns the residual
  int iterations = maxIterations_; // iterations returns number of iterations
  int restarts = maxRestarts_;
  int convergence = GMRES(solver,
    dT, myRhs, preconditioner, H, restarts, iterations, tol);
  if (convergence != 0) {
    throw ATC_Error(field_to_string(fieldName_) + " evolution did not converge");
  }
  fields[fieldName_] += dT;
  rhs[fieldName_] = myRhs;
}

} // namespace ATC
