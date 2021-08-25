// Header file for this class
#include "ImplicitSolveOperator.h"

// Other ATC includes
#include "ATC_Coupling.h"
#include "FE_Engine.h"
#include "PhysicsModel.h"
#include "PrescribedDataManager.h"


using std::set;
using std::vector;

namespace ATC {

// --------------------------------------------------------------------
// --------------------------------------------------------------------
//  ImplicitSolveOperator
// --------------------------------------------------------------------
// --------------------------------------------------------------------
ImplicitSolveOperator::
ImplicitSolveOperator(double alpha, double dt)
  : n_(0),
    dof_(0),
    dt_(dt),
    alpha_(alpha),
    epsilon0_(1.0e-8)
{
  // Nothing else to do here
}

// --------------------------------------------------------------------
//  operator *
// --------------------------------------------------------------------
DENS_VEC
ImplicitSolveOperator::operator * (const DENS_VEC & x) const
{
  // This method uses a matrix-free approach to approximate the
  // multiplication by matrix A in the matrix equation Ax=b, where the
  // matrix equation results from an implicit treatment of the
  // fast field. In brief, if the ODE for the fast field can be written:
  //
  //  dx/dt = R(x)
  //
  // A generalized discretization can be written:
  //
  //  1/dt * (x^n+1 - x^n) = alpha * R(x^n+1) + (1-alpha) * R(x^n)
  //
  // Taylor expanding the R(x^n+1) term and rearranging gives the
  // equation to be solved for dx at each timestep:
  //
  //  [1 - dt * alpha * dR/dx] * dx = dt * R(x^n)
  //
  // The operator defined in this method computes the left-hand side,
  // given a vector dx.  It uses a finite difference, matrix-free
  // approximation of dR/dx * dx, giving:
  //
  //  [1 - dt * alpha * dR/dx] * dx = dt * R(x^n)
  //      ~=  dx - dt*alpha/epsilon * ( R(x^n + epsilon*dx) - R(x^n) )
  //
  // Compute epsilon
  double epsilon = (x.norm()>0.0) ? epsilon0_*x0_.norm()/x.norm():epsilon0_;
  // Compute incremented vector x^n+1 = x^n + epsilon*dx
  x_ = x0_ + epsilon * x;
  // Evaluate R(x)
  this->R(x_,R_);
  // Compute full left hand side and return it
  DENS_VEC Ax = x - dt_ * alpha_ / epsilon * (R_ - R0_);
  return Ax;
}

// --------------------------------------------------------------------
//  rhs of Ax = r
// --------------------------------------------------------------------
DENS_VEC
ImplicitSolveOperator::r() const
{
  return dt_ * R0_; // dt * R(T^n)
}

// --------------------------------------------------------------------
//  preconditioner
// --------------------------------------------------------------------
DIAG_MAT
ImplicitSolveOperator::preconditioner() const
{
  DENS_VEC diag(n_);
  diag = 1.0;
  DIAG_MAT preconditioner(diag);
  return preconditioner;
}

// --------------------------------------------------------------------
// --------------------------------------------------------------------
//  FieldImplicitSolveOperator
// --------------------------------------------------------------------
// --------------------------------------------------------------------
FieldImplicitSolveOperator::
FieldImplicitSolveOperator(ATC_Coupling * atc,
                           FIELDS & fields,
                           const FieldName fieldName,
                           const Array2D< bool > & rhsMask,
                           const PhysicsModel * physicsModel,
                           double simTime,
                           double dt,
                           double alpha)
  : ImplicitSolveOperator(alpha, dt),
    fieldName_(fieldName),
    atc_(atc),
    physicsModel_(physicsModel),
    fields0_(fields), // ref to fields
    fields_ (fields), // copy of fields
    rhsMask_(rhsMask),
    time_(simTime)
{
  const DENS_MAT & f = fields0_[fieldName_].quantity();
  dof_ = f.nCols();
  if (dof_ > 1) throw ATC_Error("Implicit solver operator can only handle scalar fields");
  // create all to free map
  int nNodes = f.nRows();
  set<int> fixedNodes_ = atc_->prescribed_data_manager()->fixed_nodes(fieldName_);
  n_ = nNodes;
  vector<bool> tag(nNodes);
  set<int>::const_iterator it;  int i = 0;
  for (i = 0; i < nNodes; ++i) { tag[i] = true; }
  for (it=fixedNodes_.begin();it!=fixedNodes_.end();++it) {tag[*it]=false;}
  int m = 0;
  for (i = 0; i < nNodes; ++i) { if (tag[i]) freeNodes_[i]= m++; }
//std::cout << " nodes " << n_ << " " << nNodes << "\n";

  // Save current field
  x0_.reset(n_);
  to_free(f,x0_);
  x_  = x0_; // initialize

  // righthand side/forcing vector
  rhsMask_.reset(NUM_FIELDS,NUM_FLUX);
  rhsMask_ = false;
  for (int i = 0; i < rhsMask.nCols(); i++) {
    rhsMask_(fieldName_,i) =  rhsMask(fieldName_,i);
  }
//std::cout << print_mask(rhsMask_) << "\n";
  massMask_.reset(1);
  massMask_(0) = fieldName_;
  rhs_[fieldName_].reset(nNodes,dof_);
  // Compute the RHS vector R(T^n)
  R0_.reset(n_);
  R_ .reset(n_);
  R(x0_, R0_);
}

void FieldImplicitSolveOperator::to_all(const VECTOR &x, MATRIX &f) const
{
  f.reset(x.nRows(),1);
  for (int i = 0; i < x.nRows(); ++i) {
    f(i,0) = x(i);
  }
}
void FieldImplicitSolveOperator::to_free(const MATRIX &r, VECTOR &v) const
{
  v.reset(r.nRows());
  for (int i = 0; i < r.nRows(); ++i) {
    v(i) = r(i,0);
  }
}
void
FieldImplicitSolveOperator::R(const DENS_VEC &x, DENS_VEC &v ) const
{
  DENS_MAT & f = fields_[fieldName_].set_quantity();
  atc_->prescribed_data_manager()->set_fixed_field(time_, fieldName_, f);
  to_all(x,f);
  atc_->compute_rhs_vector(rhsMask_,fields_,rhs_,FULL_DOMAIN,physicsModel_);
  DENS_MAT & r = rhs_[fieldName_].set_quantity();
  atc_->prescribed_data_manager()->set_fixed_dfield(time_, fieldName_, r);
  atc_->apply_inverse_mass_matrix(r,fieldName_);
  to_free(r,v);
#if 0
int n = 6;
//std::cout << "# x "; for (int i = 0; i < n_; ++i)  std::cout << x(i) << " "; std::cout << "\n";
//std::cout << "# f "; for (int i = 0; i < n; ++i)  std::cout << f(i,0) << " "; std::cout << "\n";
std::cout << "# r "; for (int i = 0; i < n; ++i)  std::cout << r(i,0) << " "; std::cout << "\n";
//std::cout << "# v "; for (int i = 0; i < n; ++i)  std::cout << v(i) << " "; std::cout << "\n";
#endif
}

void FieldImplicitSolveOperator::solution(const DENS_MAT & dx, DENS_MAT &f) const
{
  DENS_MAT & df = fields_[fieldName_].set_quantity();
  to_all(column(dx,0),df);
  atc_->prescribed_data_manager()->set_fixed_dfield(time_, fieldName_, df);
  f += df;
}
void FieldImplicitSolveOperator::rhs(const DENS_MAT & r, DENS_MAT &rhs) const
{

  to_all(column(r,0),rhs);
}


} // namespace ATC
