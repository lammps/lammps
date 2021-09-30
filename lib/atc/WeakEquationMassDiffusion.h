#ifndef WEAK_EQUATION_MASS_DIFFUSION_H
#define WEAK_EQUATION_MASS_DIFFUSION_H

#include <set>
#include <string>

#include "WeakEquation.h"

namespace ATC{

 /**
   *  @class  WeakEquationMassDiffusion
   *  @brief  Mass diffusion
   *   c rho,t = div rho  -->
   *   int M c rho,t = int B rho
   */

class WeakEquationMassDiffusion : public WeakEquation {

  public:

  // constructor
  WeakEquationMassDiffusion();

  // destructor
  virtual ~WeakEquationMassDiffusion();

  /** density that used to form the mass matrix */
  virtual bool has_M_integrand(void) const {return true;}
  virtual void M_integrand(const FIELD_MATS &fields,
                           const Material * material,
                           DENS_MAT &density ) const ;

  /** flux that is integrated with B = Grad N as its weight */
  virtual bool has_B_integrand(void) const {return true;}
  virtual void B_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const ;

  /** necessary interfaces */
  virtual std::set<std::string> needs_material_functions(void) const
  {
    std::string list[1] = {"mass_density"};
    std::set<std::string> needs(list,list+1);
    return needs;
  }

};

};
#endif
