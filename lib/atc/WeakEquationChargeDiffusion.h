#ifndef WEAK_EQUATION_CHARGE_DIFFUSION_H
#define WEAK_EQUATION_CHARGE_DIFFUSION_H

#include <set>
#include <string>

#include "WeakEquation.h"

namespace ATC{

 /**
   *  @class  WeakEquationChargeDiffusion
   *  @brief  Charge density computation
   *  int M q = sum_s N q_s
   */

class WeakEquationChargeDiffusion : public WeakEquation {

  public:

  // constructor
  WeakEquationChargeDiffusion();

  // destructor
  virtual ~WeakEquationChargeDiffusion();

  /** density that used to form the mass matrix */
  virtual bool has_M_integrand(void) const {return true;}
  virtual void M_integrand(const FIELD_MATS &fields,
                           const Material * material,
                           DENS_MAT &density ) const ;
  /** necessary interfaces */
  virtual std::set<std::string> needs_material_functions(void)  const
  {
    std::set<std::string> needs;
    return needs;
  }
};

};
#endif
