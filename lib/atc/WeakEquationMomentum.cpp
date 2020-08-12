#include "WeakEquationMomentum.h"
#include "Material.h"
#include "LammpsInterface.h"

namespace ATC {

//==============================================================
//  Class WeakEquationMomentum
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationMomentum::WeakEquationMomentum()
  : WeakEquation(DYNAMIC_PDE,VELOCITY,3)
{}

//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
WeakEquationMomentum::~WeakEquationMomentum()
{}

//---------------------------------------------------------------------
//   compute mass density
//---------------------------------------------------------------------
void WeakEquationMomentum::M_integrand(
  const FIELD_MATS &fields,
  const Material * material,
  DENS_MAT & density ) const
{
  material->mass_density(fields, density);
}

//--------------------------------------------------------------
void WeakEquationMomentum::B_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT_VEC &flux) const
{
  material->stress(fields, grad_fields, flux);
}

//--------------------------------------------------------------
bool WeakEquationMomentum::N_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS & /* grad_fields */,
  const Material * material,
  DENS_MAT &flux) const
{
  return material->body_force(fields, flux); 
}

//--------------------------------------------------------------
void WeakEquationMomentum::E_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT &energy) const
{
  material->elastic_energy(fields, grad_fields, energy);
}
//==============================================================
//  Class WeakEquationMomentumElectrostatic
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationMomentumElectrostatic::WeakEquationMomentumElectrostatic()
  : WeakEquationMomentum()
{
  // convert charge * electric field --> force
  qE2f_ = (ATC::LammpsInterface::instance()->qe2f());
  _E_.assign(3, DENS_MAT()); 
}

//--------------------------------------------------------------
bool WeakEquationMomentumElectrostatic::N_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT &flux) const 
{
  material->electric_field(fields, grad_fields, _E_);
  // "conversion" of grad scalar to vector field
  int nsd = _E_.size();
  flux.resize(_E_[0].nRows(),nsd); 

  FIELD_MATS::const_iterator nField = fields.find(ELECTRON_DENSITY);
  const DENS_MAT &  n = nField->second;
  for (int i=0; i < nsd; i++) {
    CLON_VEC fi = column(flux,i);
    fi = _E_[i];
    fi *= -qE2f_*n; 
  }
  return true; 
}
//==============================================================
//  Class WeakEquationMomentumDiffusion
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationMomentumDiffusion::WeakEquationMomentumDiffusion()
  : WeakEquation(DYNAMIC_PDE,VELOCITY,3)
{}

//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
WeakEquationMomentumDiffusion::~WeakEquationMomentumDiffusion()
{}

//---------------------------------------------------------------------
//   compute mass density
//---------------------------------------------------------------------
void WeakEquationMomentumDiffusion::M_integrand(
  const FIELD_MATS &fields,
  const Material * material,
  DENS_MAT & density ) const
{
  material->mass_density(fields, density);
}

//--------------------------------------------------------------
void WeakEquationMomentumDiffusion::B_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT_VEC &flux) const
{
  material->viscous_stress(fields, grad_fields, flux);
}

//---------------------------------------------------------------------
void WeakEquationMomentumDiffusion::BB_tangent_coefficients(
  const FieldName /* field */,
  const FIELD_MATS &fields,
  const Material* material,
  DENS_MAT &coefs) const
{
  material->viscosity(fields,coefs);
}

//--------------------------------------------------------------
bool WeakEquationMomentumDiffusion::N_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS & /* grad_fields */,
  const Material * material,
  DENS_MAT &flux) const
{
  return material->body_force(fields, flux); 
}
}; // END namespace
