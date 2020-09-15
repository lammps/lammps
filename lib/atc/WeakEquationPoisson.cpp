#include "WeakEquationPoisson.h"
#include "Material.h"
#include <iostream>
#include <fstream>

namespace ATC {

//==============================================================
//  Class WeakEquationPoisson
//==============================================================

// R = B^T flux       + N source = 0
// 0 = B^T perm B phi - N rho_+ 

  
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationPoisson::WeakEquationPoisson()
  : WeakEquation(STATIC_PDE,ELECTRIC_POTENTIAL,1)
{}

//---------------------------------------------------------------------
void WeakEquationPoisson::B_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT_VEC &flux) const
{
  // flux = - perm D_x phi
  material->electric_displacement(fields, grad_fields, flux);
  flux[0] *= -1;
  flux[1] *= -1;
  flux[2] *= -1;
}

//---------------------------------------------------------------------
void WeakEquationPoisson::BB_tangent_coefficients(
  const FieldName /* field */,
  const FIELD_MATS &fields,
  const Material* material,
  DENS_MAT &coefs) const
{
  material->permittivity(fields,coefs); // convention evalues (B^T perm B) > 0
}

//---------------------------------------------------------------------
bool WeakEquationPoisson::N_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS & /* grad_fields */,
  const Material * material,
  DENS_MAT &flux) const
{
  return material->electron_charge_density(fields,flux); // - N rho_+ = N n
}

//---------------------------------------------------------------------
void WeakEquationPoisson::NN_tangent_coefficients(
  const FieldName /* field */,
  const FIELD_MATS &fields,
  const Material* material,
  DENS_MAT &coefs) const
{
  material->D_electron_charge_density(ELECTRIC_POTENTIAL,fields,coefs);
}

//==============================================================
//  Class WeakEquationPoissonConstantRHS
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationPoissonConstantRHS::WeakEquationPoissonConstantRHS()
  : WeakEquationPoisson()
{}

//---------------------------------------------------------------------
bool WeakEquationPoissonConstantRHS::N_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS & /* grad_fields */,
  const Material * /* material */,
  DENS_MAT &flux) const

{
  FIELD_MATS::const_iterator   nIter = fields.find(ELECTRON_DENSITY);
  if (nIter != fields.end()) {
    const DENS_MAT &  n = nIter->second;
    flux = -1.0*n; 
    return true;
  }
  return false;
}
};
