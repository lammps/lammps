#include "WeakEquationElectronContinuity.h"
#include "Material.h"
#include <iostream>
#include <fstream>

namespace ATC {

//==============================================================
//  Class WeakEquationElectronContinuity
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationElectronContinuity::WeakEquationElectronContinuity()
  : WeakEquation(DYNAMIC_PDE,ELECTRON_DENSITY,1)
{}
//--------------------------------------------------------------
//  Destructor
//---------------------------------------------------------------------
WeakEquationElectronContinuity::~WeakEquationElectronContinuity(void)
{}

//---------------------------------------------------------------------
void WeakEquationElectronContinuity::M_integrand(
  const FIELD_MATS &fields,
  const Material * material,
  DENS_MAT & density ) const
{
  FIELD_MATS::const_iterator nField = fields.find(ELECTRON_DENSITY);
  const DENS_MAT &  n = nField->second;
  density.resize(n.nRows(),n.nCols()); 
  density = 1;
}

//---------------------------------------------------------------------
void WeakEquationElectronContinuity::B_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT_VEC &flux) const
{
  material->electron_flux(fields, grad_fields, flux);
}

//---------------------------------------------------------------------
bool WeakEquationElectronContinuity::N_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT &recombination) const
{
  return material->electron_recombination(fields, grad_fields, recombination);
}

//==============================================================
//  Class WeakEquationElectronEquilbrium
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationElectronEquilibrium::WeakEquationElectronEquilibrium()
  : WeakEquation(PROJECTION_PDE,ELECTRON_DENSITY,1)
{}
//--------------------------------------------------------------
//  Destructor
//---------------------------------------------------------------------
WeakEquationElectronEquilibrium::~WeakEquationElectronEquilibrium(void)
{}

//---------------------------------------------------------------------
void WeakEquationElectronEquilibrium::M_integrand(
  const FIELD_MATS &fields,
  const Material * material,
  DENS_MAT & density ) const
{
  FIELD_MATS::const_iterator nField = fields.find(ELECTRON_DENSITY);
  const DENS_MAT &  n = nField->second;
  density.reset(n.nRows(),n.nCols()); 
  density = 1;
}

//---------------------------------------------------------------------
bool WeakEquationElectronEquilibrium::N_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT &flux) const
{
  bool flag = material->electron_charge_density(fields, flux);
  flux *= -1.; // transform from charge density to number density
  return flag;
}

//---------------------------------------------------------------------
void WeakEquationElectronEquilibrium::B_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT_VEC &flux) const
{
  material->electron_flux(fields, grad_fields, flux);
}
};
