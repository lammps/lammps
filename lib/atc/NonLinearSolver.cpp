// ATC headers 
#include "NonLinearSolver.h"
#include "LinearSolver.h"
#include "ATC_Error.h"
#include "LammpsInterface.h"


using std::stringstream;

namespace ATC {
  //===================================================================
  //  TangentOperator
  //===================================================================

  //===================================================================
  //  NonLinearSolver
  //===================================================================
  NonLinearSolver::NonLinearSolver(TangentOperator * f,
                                   const BC_SET * bcs, const int dof,
                                   bool parallel):
    f_(f),
    bcs_(bcs),
    dof_(dof),
    rNorm0P_(1.0),
    tol_(1.e-10),
    tolx_(1.e-8),
    tol0_(1.e-6),
    maxIterations_(20),
    parallel_(parallel)
  {
  }
  //--------------------------------------------------------------------
  double NonLinearSolver::residual_norm(VECTOR & r)
  {
    
    if (bcs_) {
      DENS_VEC R = r;
      BC_SET::const_iterator itr;
      for (itr = bcs_->begin(); itr != bcs_->end(); itr++) {
        int i = itr->first;
        R(i) = 0;
      }
      return R.norm();
    }
    else { return r.norm(); }
  }
  //--------------------------------------------------------------------
  bool NonLinearSolver::solve(VECTOR & x)
  {
    f_->function(x, r_);
    rNorm0_ = residual_norm(r_);
    if (rNorm0_ < tol_*rNorm0P_) { // if a "solution" does pass here rNorm0_ will be too small to allow for convergence
      return true; // note abs vs rel tol
    }
    if (rNorm0_ == 0.0) rNorm0_ = 1.0;
    if (rNorm0_ < tol0_ ) rNorm0_ = rNorm0P_;
    if (rNorm0P_ == 1.0) rNorm0P_ = rNorm0_; 
    rNormP_ = rNorm0_;
    dx_.reset(r_.nRows()); // needs to be sized for linear solver
    // newton's method  
    for (int iter = 0; iter < maxIterations_ ; iter++ ) {
      // compute tangent
      f_->tangent(x, r_, A_);
      rNorm_ = residual_norm(r_);
      rNorm_ /= rNorm0_;
      if (rNorm_ < tol_) {
        return true;
      }
      SPAR_MAT Asparse(A_); 
      LinearSolver linearSolver(Asparse, LinearSolver::AUTO_SOLVE, parallel_);
      if (bcs_) {
        linearSolver.allow_reinitialization();
        linearSolver.initialize(bcs_);
        if (iter > 0) linearSolver.set_homogeneous_bcs();
        else { x.zero(); } // linear solve w/ bcs will replace guess
      }
      r_ *= -1; 
      linearSolver.solve(dx_,r_);

      if (iter > 0 && rNorm_ > rNormP_) {
        bool descent = line_search(x);
        if (! descent ) {
//        return false;
        }
      }
      rNormP_ = rNorm_;
      x += dx_;
    }
    stringstream ss;
    ss << "WARNING NonLinearSolver: did not converge, iterations="<< maxIterations_ <<" error= " << rNorm_;
    ATC::LammpsInterface::instance()->print_msg_once(ss.str());

    return false;
  }

  //--------------------------------------------------------------------



  bool NonLinearSolver::line_search(VECTOR & x) 
  {
    double rNormP = rNormP_; 
    double dxnorm = dx_.norm();
    while ( dxnorm > tolx_) {
      dx_ *= 0.5; // bisection
      dxnorm = dx_.norm();
      f_->function(x+dx_,r_);
      rNorm_ = residual_norm(r_)/rNorm0_;
      if (rNorm_ < rNormP) return true;
    }
    return false; // no descent
  }

} // end namespace ATC
