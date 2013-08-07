#ifndef WEAK_EQUATION_SCHRODINGER_H
#define WEAK_EQUATION_SCHRODINGER_H 

#include "WeakEquation.h"

namespace ATC{

/**
  *  @class  WeakEquationSchrodinger
  *  @brief  Schrodinger equation
  *  -hbar^2 /2 1/m^* psi,xx + phi psi = E psi             
  *  (-c*B^T B + phi N^T N) psi = E N^T N psi
  */

class WeakEquationSchrodinger : public WeakEquation {

  public:
  
  // constructor 
  WeakEquationSchrodinger();

  // destructor
  virtual ~WeakEquationSchrodinger();
  
  /** RHS stiffness matrix  */
  virtual bool has_BB_tangent_coefficients(void) const {return true;}
  virtual void BB_tangent_coefficients(const FieldName field,
                                       const FIELD_MATS & fields,
                                       const Material * material,
                                       DENS_MAT &coefs) const;
  virtual bool has_NN_tangent_coefficients(void) const {return true;}
  virtual void NN_tangent_coefficients(const FieldName field,
                                       const FIELD_MATS & fields,
                                       const Material * material,
                                       DENS_MAT &coefs) const;

  /** necessary interfaces */
  virtual set<string> needs_material_functions(void) const
  {
    set<string> needs;
    needs.insert("inv_effective_mass");
    needs.insert("band_edge_potential");
    return needs;
  }

};

};
#endif
