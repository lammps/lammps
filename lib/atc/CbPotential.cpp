#include "CbPotential.h"
#include <cmath>
namespace ATC
{
  static const double EPS    = 1.0e-8;

  // Approximates the derivative of phi
  double CbPotential::phi_r(const double &r) const
  {
    const double dr = r*EPS;
    return (phi(r+dr)-phi(r)) / dr;
  }
  // Approximates the second derivative of phi
  double CbPotential::phi_rr(const double &r) const
  {
    const double dr = r*EPS;
    return (phi_r(r+dr)-phi_r(r)) / dr;
  }
  // Approximates the third derivative of phi
  double CbPotential::phi_rrr(const double &r) const 
  {
    const double dr = r*EPS;
    return (phi_rr(r+dr)-phi_rr(r)) / dr;
  }
  // Approximates the derivative of rho
  double CbPotential::rho_r(const double &r) const
  {
    const double dr = r*EPS;
    return (rho(r+dr)-rho(r)) / dr;
  }
  // Approximates the second derivative of rho
  double CbPotential::rho_rr(const double &r) const
  {
    const double dr = r*EPS;
    return (rho_r(r+dr)-rho_r(r)) / dr;
  }
  // Approximates the third derivative of rho 
  double CbPotential::rho_rrr(const double &r) const 
  {
    const double dr = r*EPS;
    return (rho_rr(r+dr)-rho_rr(r)) / dr;
  }
  // Approximates the derivative of the embedding function
  double CbPotential::F_p(const double &p) const
  {
    const double dp = p*EPS;
    return (F(p+dp)-F(p)) / dp; 
  }
  // Approximates the second derivative of the embedding function
  double CbPotential::F_pp(const double &p) const
  {
    const double dp = p*EPS;
    return (F_p(p+dp)-F_p(p)) / dp;
  }
  // Approximates the third derivative of the embedding function
  double CbPotential::F_ppp(const double &p) const
  {
    const double dp = p*EPS;
    return (F_pp(p+dp)-F_pp(p)) / dp;
  }
  // Approximates the derivative of phi3.
  double CbPotential::phi3_q (const double &q) const
  {
    const double dq = q*EPS;
    return (phi3(q+dq)-phi3(q)) / dq;
  }
  // Approximates the second derivative of phi3.
  double CbPotential::phi3_qq(const double &q) const
  {
    const double dq = q*EPS;
    return (phi3_q(q+dq)-phi3_q(q)) / dq;
  }
  // Compute bond angle jik from the squared length of vectors ij,ik,kj.
  double calculate_theta(double ij2, double ik2, double jk2)
  {
    return acos( 0.5*(ik2+ij2-jk2)/sqrt(ij2*ik2) );
  }
  // Initializes atomic interactions for up to three different terms.
  Interactions::Interactions(int a, int b, int c)
  {
    // bitwise OR combines the terms that are listed
    const int abc = a|b|c;  
    pairwise      = (abc&PAIRWISE)>0;
    embedding     = (abc&EAM)>0;
    three_body    = (abc&THREE_BDY)>0;
    angle_bending = (abc&ANGLE_BND)>0;
  }
}
