#include "FieldEulerIntegrator.h"
#include "ATC_Coupling.h"
#include "FE_Engine.h"
#include "PhysicsModel.h"
#include "PrescribedDataManager.h"
//#include "GMRES.h"
//#include "CG.h"
#include "ImplicitSolveOperator.h"
#include "MatrixDef.h"
#include "LinearSolver.h"

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
  FieldImplicitSolveOperator solver(atc_, 
    fields, fieldName_, rhsMask_, physicsModel_,
    time, dt, alpha_);
  DiagonalMatrix<double> preconditioner = solver.preconditioner();
  DENS_VEC rT = solver.r();
  DENS_VEC dT(nNodes_); dT = dT_; 
  DENS_MAT H(maxRestarts_+1, maxRestarts_);
  double tol = tol_; // tol returns the residual
  int iterations = maxIterations_; // iterations returns number of iterations
  int restarts = maxRestarts_;
  int convergence = GMRES(solver,
    dT, rT, preconditioner, H, restarts, iterations, tol);
  if (convergence != 0) {
    throw ATC_Error(field_to_string(fieldName_) + " evolution did not converge");
  }
  solver.solution(dT,fields[fieldName_].set_quantity());
}

// ====================================================================
// FieldImplicitDirectEulerIntegrator
// ====================================================================
FieldImplicitDirectEulerIntegrator::FieldImplicitDirectEulerIntegrator(
  const FieldName fieldName,
  const PhysicsModel * physicsModel, 
  /*const*/ FE_Engine * feEngine,
  /*const*/ ATC_Coupling * atc,
  const Array2D< bool > & rhsMask,  // copy 
  const double alpha
) : FieldEulerIntegrator(fieldName,physicsModel,feEngine,atc,rhsMask),
  alpha_(alpha),solver_(NULL)
{
   rhsMask_(fieldName_,FLUX) = false; // handle laplacian term with stiffness
   const BC_SET & bcs = (atc_->prescribed_data_manager()->bcs(fieldName_))[0];
   solver_ = new LinearSolver(_lhsMK_,bcs);
   solver_->allow_reinitialization();
}
FieldImplicitDirectEulerIntegrator::~FieldImplicitDirectEulerIntegrator()
{
   if (solver_) delete solver_;
}

// --------------------------------------------------------------------
// initialize 
// --------------------------------------------------------------------
void FieldImplicitDirectEulerIntegrator::initialize(const double dt, double time,
  FIELDS & fields) 
{ 
   std::pair<FieldName,FieldName> p(fieldName_,fieldName_);
   Array2D <bool>  rmask = atc_->rhs_mask();
   rmask(fieldName_,FLUX) = true; 
   atc_->tangent_matrix(p,rmask,physicsModel_,_K_);
   _lhsMK_ = (1./dt)*(_M_)-     alpha_*(_K_);
   _rhsMK_ = (1./dt)*(_M_)+(1.+alpha_)*(_K_);
}
// --------------------------------------------------------------------
// update 
// --------------------------------------------------------------------
void FieldImplicitDirectEulerIntegrator::update(const double dt, double time,
  FIELDS & fields, FIELDS & rhs) 
{ 
  atc_->compute_rhs_vector(rhsMask_, fields, rhs,
    FULL_DOMAIN, physicsModel_);
  CLON_VEC myRhs   = column(   rhs[fieldName_].set_quantity(),0);
  CLON_VEC myField = column(fields[fieldName_].set_quantity(),0);
  myRhs += _rhsMK_*myField; // f = (1/dt M + (1+alpha) K) T + f
  solver_->solve(myField,myRhs); // (1/dt M -alpha K)^-1 f
}
} // namespace ATC
