#ifndef WEAK_EQUATION_DIFFUSION_H
#define WEAK_EQUATION_DIFFUSION_H

#include <set>
#include <string>

#include "WeakEquation.h"

namespace ATC{

 /**
   *  @class  WeakEquationDiffusion
   *  @brief  species diffusion
   *  c q,t = div q  -->
   *  int M c q,t = int B q 
   */

class WeakEquationDiffusion : public WeakEquation {

  public:
  
  // constructor 
  WeakEquationDiffusion();

  // destructor
  virtual ~WeakEquationDiffusion();
  
  /** density that used to form the mass matrix */
  virtual bool has_M_integrand(void) const {return true;}
  virtual void M_integrand(const FIELD_MATS &fields,
                           const Material * material,
                           DENS_MAT &density ) const ;

  /** flux that is integrated with B = Grad N as its weight */
  virtual bool has_B_integrand(void) const {return false;} // note true
  virtual void B_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const ;

  /** necessary interfaces */
  virtual std::set<std::string> needs_material_functions(void)  const
  {
    std::set<std::string> needs;
    return needs;
  } 
};

};
#endif
