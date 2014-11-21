#ifndef WEAK_EQUATION_ELECTRON_CONTINUITY_H
#define WEAK_EQUATION_ELECTRON_CONTINUITY_H

#include <set>
#include <string>

#include "WeakEquation.h"

namespace ATC{

/**
  *  @class  WeakEquationElectronContinuity
  *  @brief  Electron continuity 
  *  n,t = div J + (G-R)  -->
  *  M n,t = int B J + int N (G-R)
  */

class WeakEquationElectronContinuity : public WeakEquation {

  public:
  
  // constructor 
  WeakEquationElectronContinuity();

  // destructor
  virtual ~WeakEquationElectronContinuity();
  
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
    std::string list[2] = {"electron_flux","electron_recombination"};
    std::set<std::string> needs(list,list+2);
    return needs;
  }

};

/**
  *  @class  WeakEquationElectronEquilibrium
  *  @brief  Electron continuity from equilibrium 
  *  n = n(\phi)
  *  M n = int N n(\phi)
 */

class WeakEquationElectronEquilibrium : public WeakEquation {

  public:
  
  // constructor 
  WeakEquationElectronEquilibrium();

  // destructor
  virtual ~WeakEquationElectronEquilibrium();
  
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
                           DENS_MAT &flux) const;
 
  
  /** flux that is integrated with B = Grad N as its weight */
  virtual bool has_B_integrand(void) const {return true;}
  virtual void B_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const;


  /** necessary interfaces */
  virtual std::set<std::string> needs_material_functions(void) const
  {
    std::set<std::string> needs;
    needs.insert("electron_charge_density");
    return needs;
  }
};
};
#endif

