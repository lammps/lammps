#ifndef CBLJSMOOTHLINEAR_H
#define CBLJSMOOTHLINEAR_H

#include <iostream>
#include <cstdlib>
#include "CbPotential.h"
namespace ATC
{

  /**
   *  @class  CbLjSmoothLinear
   *  @brief  Class for computing Cauchy-Born quantities for a Lennard-Jones/smooth/linear  material
   *          (A factor of one-half is already included to split the
   *           bond energy between atoms)
   */

  class CbLjSmoothLinear : public CbPotential
  {
  public:
    //! Constructor - initializes coefficients and enables PAIRWISE term.
    CbLjSmoothLinear(double eps, double sig, double cutoff_radius)
      : CbPotential(Interactions(PAIRWISE)),
        A  (4.0*eps*pow(sig, 12)),
        B  (4.0*eps*pow(sig,  6)),
        rcut(cutoff_radius)
    { 
      ricut = 1.0/rcut;
      r6icut = (ricut*ricut*ricut)*(ricut*ricut*ricut); 
      phicut = r6icut*(A*r6icut - B);
      dphicut = r6icut*ricut*(6.0*B - 12.0*A*r6icut); 
    }
    
    //! Returns the cutoff readius of the LJ potential.
    double cutoff_radius() const { return rcut; }
    //! Returns the LJ energy between two interacting atoms (6,12).
    double phi(const double &r) const 
    {
      const double r6i = 1.0/((r*r*r)*(r*r*r));
      if (r < rcut) { 
        return (r6i*(A*r6i - B) - phicut - (r-rcut)*dphicut);
      }
      else {
        return 0.0;
      }
    }
    //! Returns the first derivative of the LJ energy (7,13).
    double phi_r(const double &r) const 
    {
      const double ri  = 1.0/r;
      const double r6i = (ri*ri*ri)*(ri*ri*ri);
      if (r < rcut) { 
        return (r6i*ri*(6.0*B - 12.0*A*r6i) - dphicut);
      }
      else {
        return 0.0;
      }
    }
    //! Returns the second derivative of the LJ energy (8,14).
    double phi_rr(const double &r) const 
    {
      const double r2i = 1.0/(r*r);
      const double r6i = r2i*r2i*r2i;
      if (r < rcut) {
        return r6i*r2i*(13.0*12.0*A*r6i - 7.0*6.0*B);
      }
      else {
        return 0.0;
      }
    }
    //! Returns the third derivative of the LJ bond energy (9,15).
    double phi_rrr(const double &r) const
    {
      const double r3i = 1.0/(r*r*r);
      const double r9i = r3i*r3i*r3i;
      if (r < rcut) {
        return r9i*(8.0*7.0*6.0*B - 14.0*13.0*12.0*A*r3i*r3i);
      }
      else {
        return 0.0;
      }
    }

    const double A, B;  //!< phi  =  Ar^-12 +    Br^-6
    const double rcut;  //!< cutoff radius
    double ricut, r6icut, phicut, dphicut;
  };
}
#endif
