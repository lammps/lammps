#ifndef WEAK_EQUATION_MOMENTUM_H
#define WEAK_EQUATION_MOMENTUM_H

#include "WeakEquation.h"

namespace ATC{

class Material;

 /**
   *  @class  WeakEquationMomentum
   *  @brief  Momentum 
   *   rho v,t = div P  -->
   *   int M rho v,t = int B P 
   */

class WeakEquationMomentum : public WeakEquation {

  public:
  // constructor
  WeakEquationMomentum();

  // destructor
  virtual ~WeakEquationMomentum();


  /** integrand  that used to form the energy */
  virtual bool has_E_integrand(void) const {return true;}
  virtual void E_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT &energy ) const;
  
  /** density that used to form the mass matrix */
  virtual bool has_M_integrand(void) const {return true;}
  virtual void M_integrand(const FIELD_MATS &fields,
                           const Material * material,
                           DENS_MAT &density ) const ;

  /** flux that is integrated with N as its weight */
  virtual bool has_N_integrand(void) const {return true;}
  virtual bool N_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT &flux) const ;

  /** flux that is integrated with B = Grad N as its weight */
  virtual bool has_B_integrand(void) const {return true;}
  virtual void B_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const ;

  /** necessary interfaces */
  virtual set<string> needs_material_functions(void) const
  {
    string list[4] = {"mass_density","stress","elastic_energy","body_force"};
    set<string> needs(list,list+4);
    return needs;
  }


};

 /**
   *  @class  WeakEquationMomentumElectrostatic
   *  @brief  Momentum with electrostatic forces
   *   rho v,t = div P  + b -->
   *   int M rho v,t = int B P + N b
   */

class WeakEquationMomentumElectrostatic : public WeakEquationMomentum {

  public:
  // constructor
  WeakEquationMomentumElectrostatic();

  // destructor
  virtual ~WeakEquationMomentumElectrostatic(){};
  
  /** density that used to form the mass matrix */
  virtual bool has_M_integrand(void) const {return true;}
  virtual void M_integrand(const FIELD_MATS &fields,
                           const Material * material,
                           DENS_MAT &density ) const
  { WeakEquationMomentum::M_integrand(fields, material, density ); }

  /** flux that is integrated with B = Grad N as its weight */
  virtual bool has_B_integrand(void) const {return true;}
  virtual void B_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const 
  { WeakEquationMomentum::B_integrand(fields, grad_fields, material, flux); }

  /** flux that is integrated with N as its weight */
  virtual bool has_N_integrand(void) const {return true;}
  virtual bool N_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT &flux) const ; 

  /** necessary interfaces */
  virtual set<string> needs_material_functions(void) const
  {
    set<string> needs
      = WeakEquationMomentum::needs_material_functions();
    needs.insert("electric_field");
    return needs;
  }


  protected:

  /** conversion factor charge * electric field --> force */
  double qE2f_;

  mutable DENS_MAT_VEC _E_;
};

class WeakEquationMomentumDiffusion : public WeakEquation {

  public:
  // constructor
  WeakEquationMomentumDiffusion();

  // destructor
  virtual ~WeakEquationMomentumDiffusion();


  /** integrand  that used to form the energy */
  virtual bool has_E_integrand(void) const {return false;}
  
  /** density that used to form the mass matrix */
  virtual bool has_M_integrand(void) const {return true;}
  virtual void M_integrand(const FIELD_MATS &fields,
                           const Material * material,
                           DENS_MAT &density ) const ;

  /** flux that is integrated with N as its weight */
  virtual bool has_N_integrand(void) const {return true;}
  virtual bool N_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT &flux) const ;

  /** flux that is integrated with B = Grad N as its weight */
  virtual bool has_B_integrand(void) const {return true;}
  virtual void B_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const ;

  /** linear RHS stiffness matrix */
  virtual bool has_BB_tangent_coefficients(void) const {return true;}
  virtual void BB_tangent_coefficients(const FieldName field,
                                       const FIELD_MATS &fields,
                                       const Material * material,
                                       DENS_MAT &coefs) const;

  /** necessary interfaces */
  virtual set<string> needs_material_functions(void) const
  {
    string list[4] = {"mass_density","viscous_stress","body_force"};
    set<string> needs(list,list+3);
    return needs;
  }

};
} // NAMESPACE
#endif
