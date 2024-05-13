#ifndef WEAK_EQUATION_PHONON_TEMPERATURE_H
#define WEAK_EQUATION_PHONON_TEMPERATURE_H

#include <set>
#include <string>

#include "WeakEquation.h"

namespace ATC{

 /**
   *  @class  WeakEquationPhononTemperature
   *  @brief  Phonon temperature
   *   c T_p,t = div q_p  -->
   *   int M c T_p,t = int B q_p
   */

class WeakEquationPhononTemperature : public WeakEquation {

  public:

  // constructor
  WeakEquationPhononTemperature();

  // destructor
  virtual ~WeakEquationPhononTemperature();

  /** integrand  that used to form the energy */
  virtual bool has_E_integrand(void) const {return true;}
  virtual void E_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &gradFields,
                           const Material * material,
                           DENS_MAT &energy ) const ;

  /** density that used to form the mass matrix */
  virtual bool has_M_integrand(void) const {return true;}
  virtual void M_integrand(const FIELD_MATS &fields,
                           const Material * material,
                           DENS_MAT &density ) const ;

  /** flux that is integrated with B = Grad N as its weight */
  virtual bool has_B_integrand(void) const {return true;}
  virtual void B_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &gradFields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const ;

   /** necessary interfaces */
  virtual std::set<std::string> needs_material_functions(void) const
  {
    std::string list[3] = {"thermal_energy","heat_capacity","heat_flux"};
    std::set<std::string> needs(list,list+3);
    return needs;

  }

};

 /**
   *  @class  WeakEquationPhononTemperatureExchange
   *  @brief  Phonon temperature with exchange to electrons
   *   c T_p,t = div q_p  + g(T_p-T_e) -->
   *   int M c T_p,t = int B q_p + int N g
   */

class WeakEquationPhononTemperatureExchange :
  public WeakEquationPhononTemperature {

  public:

  // constructor
  WeakEquationPhononTemperatureExchange();

  // destructor
  virtual ~WeakEquationPhononTemperatureExchange();

  /** integrand  that used to form the energy */
  virtual bool has_E_integrand(void) const {return true;}
  virtual void E_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &gradFields,
                           const Material * material,
                           DENS_MAT &energy ) const;

  /** density that used to form the mass matrix */
  virtual bool has_M_integrand(void) const {return true;}
  virtual void M_integrand(const FIELD_MATS &fields,
                           const Material * material,
                           DENS_MAT &density ) const;

  /** flux that is integrated with B = Grad N as its weight */
  virtual bool has_B_integrand(void) const {return true;}
  virtual void B_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &gradFields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const;

  /** flux that is integrated with N as its weight */
  virtual bool has_N_integrand(void) const {return true;}
  virtual bool N_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &gradFields,
                           const Material * material,
                           DENS_MAT &flux) const;

  /** necessary interfaces */
  virtual std::set<std::string> needs_material_functions(void) const
  {
    std::set<std::string> needs
      = WeakEquationPhononTemperature::needs_material_functions();
    needs.insert("electron_phonon_exchange");
    return needs;
  }

};

}; // namespace
#endif
