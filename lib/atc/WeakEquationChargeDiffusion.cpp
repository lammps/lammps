#include "ATC_Transfer.h"
#include "WeakEquationChargeDiffusion.h"
#include "Material.h"
#include <iostream>
#include <fstream>

namespace ATC {

//==============================================================
//  Class WeakEquationChargeDiffusion
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationChargeDiffusion::WeakEquationChargeDiffusion()
  : WeakEquation(PROJECTION_PDE,CHARGE_DENSITY,1)
{}
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
WeakEquationChargeDiffusion::~WeakEquationChargeDiffusion(void)
{}
//---------------------------------------------------------------------
//   compute capacity
//---------------------------------------------------------------------
void WeakEquationChargeDiffusion::M_integrand(
  const FIELD_MATS &fields,
  const Material * /* material */,
  DENS_MAT & capacity ) const
{
  FIELD_MATS::const_iterator rhoField = fields.find(CHARGE_DENSITY);
  const DENS_MAT &  rho = rhoField->second;
  capacity.reset(rho.nRows(),rho.nCols());
  capacity = 1.;
}
}; // end namespace


