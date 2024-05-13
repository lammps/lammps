#ifndef WEAK_EQUATION_POISSON_H
#define WEAK_EQUATION_POISSON_H

#include <set>
#include <string>

#include "WeakEquation.h"

namespace ATC{

 /**
   *  @class  WeakEquationPoisson
   *  @brief  Poisson equation
   *  0 = perm grad^2 \phi + \rho -->
   *  \int B^T perm B dv \phi  = \int N^T \rho dv
   */

class WeakEquationPoisson : public WeakEquation {

  public:

  // constructor
  WeakEquationPoisson();

  // destructor
  virtual ~WeakEquationPoisson(){};

  /** flux that is integrated with B = Grad N as its weight */
  virtual bool has_B_integrand(void) const {return true;}
  virtual void B_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const ;

  /** flux that is integrated with N as its weight */
  //virtual bool has_N_integrand(void) const {return false;}
  virtual bool has_N_integrand(void) const {return true;}
  virtual bool N_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT &flux) const ;


  /** linear RHS stiffness matrix */
  virtual bool has_BB_tangent_coefficients(void) const {return true;}
  virtual void BB_tangent_coefficients(const FieldName field,
                                       const FIELD_MATS &fields,
                                       const Material * material,
                                       DENS_MAT &coefs) const;

  /** linear RHS stiffness matrix */
  virtual bool has_NN_tangent_coefficients(void) const {return true;}
  virtual void NN_tangent_coefficients(const FieldName field,
                                       const FIELD_MATS &fields,
                                       const Material * material,
                                       DENS_MAT &coefs) const;

  /** necessary interfaces */
  virtual std::set<std::string> needs_material_functions(void) const
  {
    std::set<std::string> needs;
    needs.insert("electric_displacement");
    needs.insert("electron_charge_density");
    return needs;
  }
};

 /**
   *  @class  WeakEquationPoissonConstantRHS
   *  @brief  Poisson equation with constant RHS
   */


class WeakEquationPoissonConstantRHS : public WeakEquationPoisson {

  public:

  // constructor
  WeakEquationPoissonConstantRHS();

  // destructor
  virtual ~WeakEquationPoissonConstantRHS(){};

  /** flux that is integrated with B = Grad N as its weight */
  virtual bool has_B_integrand(void) const {return true;}
  virtual void B_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const
  { WeakEquationPoisson::B_integrand(fields, grad_fields, material, flux) ; }

  /** flux that is integrated with N as its weight */
  virtual bool has_N_integrand(void) const {return true;}
  virtual bool N_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT &flux) const ;

  /** rhs is constant */
  virtual bool has_NN_tangent_coefficients(void) const {return false;}
  virtual void NN_tangent_coefficients(const FieldName /* field */,
                                       const FIELD_MATS & /* fields */,
                                       const Material * /* material */,
                                       DENS_MAT & /* coefs */) const {};

  virtual std::set<std::string> needs_material_functions(void) const
  {
    std::set<std::string> needs;
    needs.insert("electric_displacement");
    return needs;
  }

};

}; // namespace
#endif
