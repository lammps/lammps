#include "ATC_Transfer.h"
#include "WeakEquationDiffusion.h"
#include "Material.h"
#include <string>
#include <iostream>
#include <fstream>

namespace ATC {

//==============================================================
//  Class WeakEquationDiffusion
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationDiffusion::WeakEquationDiffusion()
  : WeakEquation(DYNAMIC_PDE,SPECIES_CONCENTRATION,1) 
{}
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
WeakEquationDiffusion::~WeakEquationDiffusion(void)
{}
//---------------------------------------------------------------------
//   compute capacity
//---------------------------------------------------------------------
void WeakEquationDiffusion::M_integrand(
  const FIELD_MATS &fields,
  const Material * material,
  DENS_MAT & capacity ) const
{
  FIELD_MATS::const_iterator rhoField = fields.find(SPECIES_CONCENTRATION);
  const DENS_MAT &  rho = rhoField->second;
  capacity.reset(rho.nRows(),rho.nCols());
  capacity = 1.;
}
//--------------------------------------------------------------
//  compute flux
//--------------------------------------------------------------
void WeakEquationDiffusion::B_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT_VEC &flux) const
{
//  material->diffusion_flux(fields, grad_fields, flux[SPECIES_CONCENTRATION]);
}
}; // end namespace

  
