#ifndef WEAK_EQUATION_ELECTRON_MOMENTUM_H
#define WEAK_EQUATION_ELECTRON_MOMENTUM_H

#include <set>
#include <string>

#include "WeakEquation.h"

namespace ATC{

  class Material;

 /**
   *  @class  WeakEquationElectronMomentum
   *  @brief  Electron momentum
   *  rho v,t = div P  -->
   *  int M rho v,t = int B P
   */

  class WeakEquationElectronMomentum : public WeakEquation {

  public:
    // constructor
    WeakEquationElectronMomentum();

    // destructor
    virtual ~WeakEquationElectronMomentum();

    /** density that used to form the mass matrix */
    virtual bool has_M_integrand(void) const {return true;}
    virtual void M_integrand(const FIELD_MATS &fields,
                             const Material * material,
                             DENS_MAT &density ) const ;

    /** convection for flux */
    virtual bool has_B_integrand(void) const {return true;}
    virtual void B_integrand(const FIELD_MATS &fields,
                             const GRAD_FIELD_MATS &grad_fields,
                             const Material * material,
                             DENS_MAT_VEC &flux) const ;

    /** necessary interfaces */
    virtual std::set<std::string> needs_material_functions(void) const
    {
      std::string list[2] = {"inv_effective_mass","electron_heat_capacity"};
      std::set<std::string> needs(list,list+2);
      return needs;
    }

  protected:
    /** computes standard velocity convective fluxes */
    virtual void convection(const FIELD_MATS &fields,
                            const Material * material,
                            DENS_MAT_VEC & flux) const ;

  };

 /**
   *  @class  WeakEquationElectronMomentumDDM
   *  @brief  Electron momentum - drift diffusion
   *  rho v,t = div P  -->
   *  int M rho v,t = int B P
   */

  class WeakEquationElectronMomentumDDM : public WeakEquationElectronMomentum {

  public:
    // constructor
    WeakEquationElectronMomentumDDM();

    // destructor
    virtual ~WeakEquationElectronMomentumDDM();

    /** density that used to form the mass matrix */
    virtual void M_integrand(const FIELD_MATS &fields,
                             const Material * material,
                             DENS_MAT &density ) const ;

    /** flux that is integrated with grad N as its weight */
    virtual bool has_B_integrand(void) const {return false;}
    virtual void B_integrand(const FIELD_MATS & /* fields */,
                             const GRAD_FIELD_MATS & /* grad_fields */,
                             const Material * /* material */,
                             DENS_MAT_VEC &/* flux */) const {};

    /** flux that is integrated with N as its weight */
    virtual bool has_N_integrand(void) const {return true;}
    virtual bool N_integrand(const FIELD_MATS &fields,
                             const GRAD_FIELD_MATS &grad_fields,
                             const Material * material,
                             DENS_MAT &flux) const ;

    /** necessary interfaces */
    virtual std::set<std::string> needs_material_functions(void) const
    {
      std::set<std::string> needs
        = WeakEquationElectronMomentum::needs_material_functions();
      needs.insert("electron_drag_coefficient");
      needs.insert("electron_heat_capacity");
      needs.insert("electric_displacement");
      return needs;
    }

  protected:
    /** computes thermal stresses arising from convection */
    virtual void thermal_stress(const FIELD_MATS &fields,
                                const GRAD_FIELD_MATS &gradFields,
                                const Material * material,
                                DENS_MAT &flux) const ;

    /** workspace variables */
    mutable DENS_MAT_VEC _dnCp_;
    mutable DENS_MAT_VEC _electricForce_;
  };

};
#endif
