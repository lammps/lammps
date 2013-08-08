#include "ATC_Transfer.h"
#include "WeakEquationMassDiffusion.h"
#include "Material.h"
#include <string>
#include <iostream>
#include <fstream>

namespace ATC {

//==============================================================
//  Class WeakEquationMassDiffusion
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationMassDiffusion::WeakEquationMassDiffusion()
  : WeakEquation(DYNAMIC_PDE,MASS_DENSITY,1)
{}
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
WeakEquationMassDiffusion::~WeakEquationMassDiffusion(void)
{}
//---------------------------------------------------------------------
//   compute capacity
//---------------------------------------------------------------------
void WeakEquationMassDiffusion::M_integrand(
  const FIELD_MATS &fields,
  const Material * material,
  DENS_MAT & capacity ) const
{
  FIELD_MATS::const_iterator dField = fields.find(MASS_DENSITY);
  const DENS_MAT &  rho = dField->second;
  capacity.reset(rho.nRows(),rho.nCols());
  capacity = 1.;
}
//--------------------------------------------------------------
//  compute flux
//--------------------------------------------------------------
void WeakEquationMassDiffusion::B_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT_VEC &flux) const
{
//  material->mass_flux(fields, grad_fields, flux[MASS_DENSITY]);
}
}; // end namespace

  
