// Header file for this class
#include "FieldEulerIntegrator.h"

// Other ATC includes
#include "ATC_Transfer.h"
#include "FE_Engine.h"
#include "PhysicsModel.h"
#include "GMRES.h"
#include "ImplicitSolveOperator.h"

namespace ATC {

// ====================================================================
// FieldEulerIntegrator
// ====================================================================
FieldEulerIntegrator::FieldEulerIntegrator(
  const FieldName fieldName,
  const PhysicsModel * physicsModel, 
  /*const*/ FE_Engine * feEngine,
  /*const*/ ATC_Transfer * atcTransfer,
  const Array2D< bool > & rhsMask  // copy 
)
  : atc_(atcTransfer),
    feEngine_(feEngine),
    physicsModel_(physicsModel),
    fieldName_(fieldName),
    rhsMask_(rhsMask)
{
  nNodes_ = feEngine->get_nNodes();
}

// ====================================================================
// FieldImplicitIntegrator
// ====================================================================
FieldExplicitEulerIntegrator::FieldExplicitEulerIntegrator(
  const FieldName fieldName,
  const PhysicsModel * physicsModel, 
  /*const*/ FE_Engine * feEngine,
  /*const*/ ATC_Transfer * atcTransfer,
  const Array2D< bool > & rhsMask  // copy 
) : FieldEulerIntegrator(fieldName,physicsModel,feEngine,atcTransfer,rhsMask)
{
}

// --------------------------------------------------------------------
// update 
// --------------------------------------------------------------------
void FieldExplicitEulerIntegrator::update(const double dt, double time,
  FIELDS & fields, FIELDS & rhs)
{ // NOTE time is not used
  atc_->compute_rhs_vector(rhsMask_, fields, rhs,
    atc_->FULL_DOMAIN, physicsModel_);
  atc_->apply_inverse_mass_matrix(rhs[fieldName_],fieldName_);
  fields[fieldName_] += dt*rhs[fieldName_];
}

// ====================================================================
// FieldImplicitEulerIntegrator
// ====================================================================
FieldImplicitEulerIntegrator::FieldImplicitEulerIntegrator(
  const FieldName fieldName,
  const PhysicsModel * physicsModel, 
  /*const*/ FE_Engine * feEngine,
  /*const*/ ATC_Transfer * atcTransfer,
  const Array2D< bool > & rhsMask,  // copy 
  const double alpha
) : FieldEulerIntegrator(fieldName,physicsModel,feEngine,atcTransfer,rhsMask),
  alpha_(alpha),
  dT_(1.0e-6), 
  maxRestarts_(50),
  maxIterations_(200),
  tol_(1.0e-8)
{
}

// --------------------------------------------------------------------
// update 
// --------------------------------------------------------------------
void FieldImplicitEulerIntegrator::update(const double dt, double time,
  FIELDS & fields, FIELDS & rhs) 
{ // solver handles bcs
  FieldImplicitSolveOperator solver(atc_, feEngine_, // NOTE make persistent
    fields, fieldName_, rhsMask_, physicsModel_,
    time, dt, alpha_);
  DiagonalMatrix<double> preconditioner = solver.get_preconditioner(fields);
  DENS_VEC myRhs = solver.get_rhs();
  DENS_VEC dT(nNodes_); dT = dT_; 
  DENS_MAT H(maxRestarts_+1, maxRestarts_);
  double tol = tol_; // tol returns the residual
  int iterations = maxIterations_; // iterations returns number of iterations
  int restarts = maxRestarts_;
  int convergence = GMRES(solver,
    dT, myRhs, preconditioner, H, restarts, iterations, tol);
  if (convergence != 0) {
    throw ATC_Error(0,field_to_string(fieldName_) + " evolution did not converge");
  }
  fields[fieldName_] += dT;
  rhs[fieldName_] = myRhs;
}

} // namespace ATC
