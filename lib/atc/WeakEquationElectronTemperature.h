#ifndef WEAK_EQUATION_ELECTRON_TEMPERATURE_H
#define WEAK_EQUATION_ELECTRON_TEMPERATURE_H

#include <set>
#include <string>

#include "WeakEquation.h"

namespace ATC{

 /**
   *  @class  WeakEquationElectronTemperature
   *  @brief  Electron temperature 
   *   c T_e,t = div q_e + g(T_e-T_p) --> 
   *   int M c T_e,t = int B q_e + int N g
   */

class WeakEquationElectronTemperature : public WeakEquation {

  public:
  
  // constructor 
  WeakEquationElectronTemperature();

  // destructor
  virtual ~WeakEquationElectronTemperature();
  
  /** integrand  that used to form the energy */
  virtual bool has_E_integrand(void) const {return true;}
  virtual void E_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
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
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const ;

  /** flux that is integrated with N as its weight */
  virtual bool has_N_integrand(void) const {return true;}
  virtual bool N_integrand(const FIELD_MATS &fields, 
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT &flux) const ;

  /** necessary interfaces */
  virtual std::set<std::string> needs_material_functions(void) const
  {
    std::string list[4] = {"electron_thermal_energy",
                      "electron_heat_capacity",
                      "electron_phonon_exchange",
                      "electron_heat_flux"};
    std::set<std::string> needs(list,list+4);
    return needs;
  }

};

 /**
   *  @class  WeakEquationElectronTemperatureJouleHeating
   *  @brief  Electron temperature with Joule heating
   *   c T_e,t = div q_e + g(T_e-T_p) + J.E --> 
   *   int M c T_e,t = int B q_e + int N ( g + J.E)
   */

class WeakEquationElectronTemperatureJouleHeating : 
  public WeakEquationElectronTemperature  
{
  public:
  
  // constructor 
  WeakEquationElectronTemperatureJouleHeating();

  // destructor
  virtual ~WeakEquationElectronTemperatureJouleHeating();
  
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
                           DENS_MAT &density ) const;

  /** flux that is integrated with B = Grad N as its weight */
  virtual bool has_B_integrand(void) const {return true;}
  virtual void B_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const;

  /** flux that is integrated with N as its weight */
  virtual bool has_N_integrand(void) const {return true;}
  virtual bool N_integrand(const FIELD_MATS &fields, 
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT &flux) const;

  /** necessary interfaces */
  virtual std::set<std::string> needs_material_functions(void) const
  {
    std::set<std::string> needs 
      = WeakEquationElectronTemperature::needs_material_functions();
    needs.insert("electric_field");
    return needs;
  }


  protected:
  
  /** conversion factor for charge * voltage --> mass length^2 / time^2 */
  double eV2E_;

  mutable DENS_MAT_VEC _J_, _E_;
};

 /**
   *  @class  WeakEquationElectronTemperatureConvection
   *  @brief  Electron temperature with convection
   *   c ( T_e,t + grad m_e J.E T_e / e ) = div q_e + g(T_e-T_p) + c n T_e grad J.E/(n e) --> 
   *   int M c T_e,t = int B ( q_e - c m_e J.E T_e / e ) + int N ( g + c n T_e grad J.E/(n e) )
   */

class WeakEquationElectronTemperatureConvection : 
  public WeakEquationElectronTemperatureJouleHeating
{
  public:
  
  // constructor 
  WeakEquationElectronTemperatureConvection();

  // destructor
  virtual ~WeakEquationElectronTemperatureConvection();

  /** flux that is integrated with B = Grad N as its weight */
  virtual bool has_B_integrand(void) const {return true;}
  virtual void B_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const;

  /** flux that is integrated with N as its weight */
  virtual bool has_N_integrand(void) const {return true;}
  virtual bool N_integrand(const FIELD_MATS &fields, 
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT &flux) const;

  /** necessary interfaces */
  virtual std::set<std::string> needs_material_functions(void) const
  {
    std::set<std::string> needs 
      = WeakEquationElectronTemperature::needs_material_functions();
    needs.insert("electron_drag_power");
    return needs;
  }

  protected:

  /** workspace variable for the convective flux */
  mutable DENS_MAT_VEC _convectiveFlux_;
  
};

}; // end namespace
#endif
