#ifndef WEAK_EQUATION_H
#define WEAK_EQUATION_H

#include "ATC_TypeDefs.h"

namespace ATC{

class Material;

 /**
   *  @class  WeakEquation
   *  @brief  WeakEquation is a template for adapting FE weak equations
   *  for the FE_Engine. It wrappers the correct material functions for
   *  integration over elements
   */

class WeakEquation {

  public:
  /** type of equation */
  enum PDE_Type {DYNAMIC_PDE=0, STATIC_PDE, EIGENVALUE_PDE, PROJECTION_PDE};
  
  // constructor 
  WeakEquation(PDE_Type type, FieldName fieldName, int fieldSize) : 
  type_(type), fieldName_(fieldName), fieldSize_(fieldSize), nsd_(3) {};

  // destructor
  virtual ~WeakEquation(){};
  
  /** integrand  that used to form the energy */
  virtual bool has_E_integrand(void) const {return false;}
  virtual void E_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT &energy ) const {};

  /** density that used to form the mass matrix */
  virtual bool has_M_integrand(void) const {return false;}
  virtual void M_integrand(const FIELD_MATS &fields,
                           const Material * material,
                           DENS_MAT &density ) const {}; 

  /** flux that is integrated with B = Grad N as its weight */
  virtual bool has_B_integrand(void) const {return false;}
  virtual void B_integrand(const FIELD_MATS &fields,
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT_VEC &flux) const {};

  /** flux that is integrated with N as its weight */
  virtual bool has_N_integrand(void) const {return false;}
  // N_integrand bool is for masking in FE_Engine
  virtual bool N_integrand(const FIELD_MATS &fields, 
                           const GRAD_FIELD_MATS &grad_fields,
                           const Material * material,
                           DENS_MAT &flux) const {return false;};

  /** stiffness matrix */ 

  // linear
  virtual void BB_tangent_coefficients(const FieldName field,
                                       const Material * material,
                                       DENS_MAT &coefs) const {};
  virtual void NN_tangent_coefficients(const FieldName field,
                                       const Material * material,
                                       DENS_MAT &coefs) const {};
  // non-linear
  virtual bool has_BB_tangent_coefficients(void) const {return false;}
  virtual void BB_tangent_coefficients(const FieldName field,
                                       const FIELD_MATS &fields, 
                                       const Material * material,
                                       DENS_MAT &coefs) const {};
  virtual bool has_NN_tangent_coefficients(void) const {return false;}
  virtual void NN_tangent_coefficients(const FieldName field,
                                       const FIELD_MATS &fields, 
                                       const Material * material,
                                       DENS_MAT &coefs) const {};
  
  /** type of equation */
  PDE_Type type(void) const {return type_;}

  /** primary field */
  FieldName field_name(void) const {return fieldName_;}
  int field_size(void) const {return fieldSize_;}

  /** list of require interfaces */
  virtual set<string> needs_material_functions(void) const = 0;

  protected:
  /** type of equation */
  PDE_Type type_;

  /** field */
  FieldName fieldName_;

  /** field size */
  int fieldSize_;

  /** number of spatial dimensions */
  int nsd_;

};

};
#endif
