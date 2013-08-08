// Header file for this class
#include "ImplicitSolveOperator.h"

// Other ATC includes
#include "ATC_Coupling.h"
#include "FE_Engine.h"
#include "PhysicsModel.h"
#include "PrescribedDataManager.h"

namespace ATC {

// --------------------------------------------------------------------
// --------------------------------------------------------------------
//  ImplicitSolveOperator
// --------------------------------------------------------------------
// --------------------------------------------------------------------
ImplicitSolveOperator::
ImplicitSolveOperator(ATC_Coupling * atc,
                      /*const*/ FE_Engine * feEngine,
                      const PhysicsModel * physicsModel)
  : atc_(atc),
    feEngine_(feEngine),
    physicsModel_(physicsModel)
{
  // Nothing else to do here
}

// --------------------------------------------------------------------
// --------------------------------------------------------------------
//  FieldImplicitSolveOperator
// --------------------------------------------------------------------
// --------------------------------------------------------------------
FieldImplicitSolveOperator::
FieldImplicitSolveOperator(ATC_Coupling * atc,
                           /*const*/ FE_Engine * feEngine,
                           FIELDS & fields,
                           const FieldName fieldName,
                           const Array2D< bool > & rhsMask,
                           const PhysicsModel * physicsModel,
                           double simTime,
                           double dt,
                           double alpha)
  : ImplicitSolveOperator(atc, feEngine, physicsModel),
    fieldName_(fieldName),
    fields_(fields), // ref to fields
    time_(simTime),
    dt_(dt),
    alpha_(alpha),
    epsilon0_(1.0e-8)
{
  // find field associated with ODE
  rhsMask_.reset(NUM_FIELDS,NUM_FLUX);
  rhsMask_ = false;
  for (int i = 0; i < rhsMask.nCols(); i++) {
    rhsMask_(fieldName_,i) =  rhsMask(fieldName_,i);
  }
  massMask_.reset(1);
  massMask_(0) = fieldName_;

  // Save off current field
  TnVect_ = column(fields_[fieldName_].quantity(),0); 

  // Allocate vectors for fields and rhs
  int nNodes = atc_->num_nodes();
  // copy fields
  fieldsNp1_ = fields_;
  // size rhs
  int dof = fields_[fieldName_].nCols();
  RnMap_ [fieldName_].reset(nNodes,dof);
  RnpMap_[fieldName_].reset(nNodes,dof);
 
  // Compute the RHS vector R(T^n) 
  // Set BCs on Rn, multiply by inverse mass and then extract its vector
  atc_->compute_rhs_vector(rhsMask_, fields_, RnMap_,
                                   FULL_DOMAIN, physicsModel_);
  DENS_MAT & Rn = RnMap_[fieldName_].set_quantity();
  atc_->prescribed_data_manager()
    ->set_fixed_dfield(time_, fieldName_, Rn);
  atc_->apply_inverse_mass_matrix(Rn,fieldName_);
  RnVect_ = column(Rn,0);
}


// --------------------------------------------------------------------
//  operator *(Vector)
// --------------------------------------------------------------------
DENS_VEC
FieldImplicitSolveOperator::operator * (DENS_VEC x) const
{
  // This method uses a matrix-free approach to approximate the
  // multiplication by matrix A in the matrix equation Ax=b, where the
  // matrix equation results from an implicit treatment of the
  // fast field solve for the Two Temperature Model.  In
  // brief, if the ODE for the fast field can be written:
  //
  //  dT/dt = R(T)
  // 
  // A generalized discretization can be written:
  //
  //  1/dt * (T^n+1 - T^n) = alpha * R(T^n+1) + (1-alpha) * R(T^n)
  //
  // Taylor expanding the R(T^n+1) term and rearranging gives the
  // equation to be solved for dT at each timestep:
  //
  //  [1 - dt * alpha * dR/dT] * dT = dt * R(T^n)
  //
  // The operator defined in this method computes the left-hand side,
  // given a vector dT.  It uses a finite difference, matrix-free
  // approximation of dR/dT * dT, giving:
  //
  //  [1 - dt * alpha * dR/dT] * dT = dt * R(T^n)
  //      ~=  dT - dt*alpha/epsilon * ( R(T^n + epsilon*dT) - R(T^n) )
  //
  
  
  // Compute epsilon
  double epsilon = (x.norm() > 0.0) ? epsilon0_ * TnVect_.norm()/x.norm() : epsilon0_;

  // Compute incremented vector = T + epsilon*dT
  fieldsNp1_[fieldName_] = TnVect_ + epsilon * x;

  // Evaluate R(b)
  atc_->compute_rhs_vector(rhsMask_, fieldsNp1_, RnpMap_,
                                   FULL_DOMAIN, physicsModel_);
  DENS_MAT & Rnp = RnpMap_[fieldName_].set_quantity();
  atc_->prescribed_data_manager()
    ->set_fixed_dfield(time_, fieldName_, Rnp);
  atc_->apply_inverse_mass_matrix(Rnp,fieldName_);
  RnpVect_ = column(Rnp,0);

  // Compute full left hand side and return it
  DENS_VEC Ax = x - dt_ * alpha_ / epsilon * (RnpVect_ - RnVect_);
  return Ax;
}

// --------------------------------------------------------------------
//  rhs
// --------------------------------------------------------------------
DENS_VEC
FieldImplicitSolveOperator::rhs()
{
  // Return dt * R(T^n)
  return dt_ * RnVect_;
}

// --------------------------------------------------------------------
//  preconditioner
// --------------------------------------------------------------------
DIAG_MAT
FieldImplicitSolveOperator::preconditioner(FIELDS & fields)
{
  // Just create and return identity matrix
  int nNodes = atc_->num_nodes();
  DENS_VEC ones(nNodes);
  ones = 1.0;
  DIAG_MAT identity(ones);
  return identity;
}

} // namespace ATC
