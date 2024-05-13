#include "WeakEquationElectronTemperature.h"
#include "Material.h"
#include <iostream>
#include <fstream>

namespace ATC {
//==============================================================
//  Class WeakEquationElectronTemperature
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationElectronTemperature::WeakEquationElectronTemperature()
  : WeakEquation(DYNAMIC_PDE,ELECTRON_TEMPERATURE,1)
{
}
//--------------------------------------------------------------
//  Destructor
//---------------------------------------------------------------------
WeakEquationElectronTemperature::~WeakEquationElectronTemperature()
{}

//---------------------------------------------------------------------
//   compute energy
//---------------------------------------------------------------------
void WeakEquationElectronTemperature::E_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS & /* grad_fields */,
  const Material * material,
  DENS_MAT & energy ) const
{
  material->electron_thermal_energy(fields, energy);
}

//---------------------------------------------------------------------
//   compute heat capacities
//---------------------------------------------------------------------
void WeakEquationElectronTemperature::M_integrand(
  const FIELD_MATS &fields,
  const Material * material,
  DENS_MAT & capacity ) const
{
  material->electron_heat_capacity(fields, capacity);
}

//---------------------------------------------------------------------
//   compute heat fluxes
//---------------------------------------------------------------------
void WeakEquationElectronTemperature::B_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT_VEC &flux) const
{
  material->electron_heat_flux(fields, grad_fields, flux);
}

//---------------------------------------------------------------------
//   compute exchange fluxes
//---------------------------------------------------------------------
bool WeakEquationElectronTemperature::N_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS & /* grad_fields */,
  const Material * material,
  DENS_MAT &flux) const
{

  DENS_MAT exchange_flux;
  bool has = material->electron_phonon_exchange(fields, exchange_flux);
  if (has) flux = -1.*exchange_flux;
  return has;
}

//==============================================================
//  Class WeakEquationElectronJouleHeating
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationElectronTemperatureJouleHeating::WeakEquationElectronTemperatureJouleHeating()
  : WeakEquationElectronTemperature()
{
  // convert charge * voltage --> mass length^2 / time^2
  //eV2E_ = (ATC::LammpsInterface::instance()->qe2f())
  //      * (ATC::LammpsInterface::instance()->ftm2v());
  eV2E_ = ATC::LammpsInterface::instance()->qv2e();
  int nSD = 3;
  _J_.assign(nSD, DENS_MAT());
  _E_.assign(nSD, DENS_MAT());
}
//--------------------------------------------------------------
//  Destructor
//---------------------------------------------------------------------
WeakEquationElectronTemperatureJouleHeating::~WeakEquationElectronTemperatureJouleHeating()
{}

//---------------------------------------------------------------------
void WeakEquationElectronTemperatureJouleHeating::E_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT &energy) const
{
  WeakEquationElectronTemperature::E_integrand(fields, grad_fields, material, energy);
}
//---------------------------------------------------------------------
void WeakEquationElectronTemperatureJouleHeating::M_integrand(
  const FIELD_MATS &fields,
  const Material * material,
  DENS_MAT &capacity) const
{
  WeakEquationElectronTemperature::M_integrand(fields, material, capacity);
}
//---------------------------------------------------------------------
void WeakEquationElectronTemperatureJouleHeating::B_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT_VEC &flux) const
{
  WeakEquationElectronTemperature::B_integrand(fields, grad_fields, material, flux);
}

//---------------------------------------------------------------------
bool WeakEquationElectronTemperatureJouleHeating::N_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT &flux) const
{

  // call base class to get electron_temperature terms
  WeakEquationElectronTemperature::N_integrand(fields, grad_fields, material, flux);
  // Joule heating = -I.grad Psi = J.grad Psi \approx J.E
  DENS_MAT jouleHeating;
  material->electron_flux (fields, grad_fields, _J_);
  material->electric_field(fields, grad_fields, _E_);
  jouleHeating = _J_[0].mult_by_element(_E_[0]);
  for (DENS_MAT_VEC::size_type i=1; i < _J_.size(); i++)
    jouleHeating += _J_[i].mult_by_element(_E_[i]);
  jouleHeating *= eV2E_;
  flux -= jouleHeating;
  return true;
}

//==============================================================
//  Class WeakEquationElectronConvection
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationElectronTemperatureConvection::WeakEquationElectronTemperatureConvection()
  : WeakEquationElectronTemperatureJouleHeating()
{
  int nSD = 3;
  _convectiveFlux_.assign(nSD, DENS_MAT());
}

//--------------------------------------------------------------
//  Destructor
//---------------------------------------------------------------------
WeakEquationElectronTemperatureConvection::~WeakEquationElectronTemperatureConvection()
{
  // do nothing
}

//---------------------------------------------------------------------
void WeakEquationElectronTemperatureConvection::B_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT_VEC &flux) const
{
  // add diffusion term

  WeakEquationElectronTemperatureJouleHeating::B_integrand(fields, grad_fields, material, flux);
  //flux[0] = 0.;
  //flux[1] = 0.;
  //flux[2] = 0.;

  // add convection term
  DENS_MAT_VEC convectiveFlux;
  material->electron_heat_convection(fields,_convectiveFlux_);
  flux[0] += _convectiveFlux_[0];
  flux[1] += _convectiveFlux_[1];
  flux[2] += _convectiveFlux_[2];
}

//---------------------------------------------------------------------
bool WeakEquationElectronTemperatureConvection::N_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT &flux) const
{
  // call base class to get electron_temperature terms
  WeakEquationElectronTemperatureJouleHeating::N_integrand(fields, grad_fields, material, flux);
#ifdef TEST
  // add exchange with kinetic energy
  DENS_MAT keExchange;
  DENS_MAT capacity;
  material->electron_heat_capacity(fields, capacity);
  capacity *= 2./3.; // correction in DDM equations


  //FIELD_MATS::const_iterator dField = fields.find(ELECTRON_DENSITY);
  FIELD_MATS::const_iterator tField = fields.find(ELECTRON_TEMPERATURE);
  //const DENS_MAT & density = dField->second;
  const DENS_MAT & temperature = tField->second;

  GRAD_FIELD_MATS::const_iterator velocityGradients = grad_fields.find(ELECTRON_VELOCITY);
  const DENS_MAT_VEC & dv = velocityGradients->second;
  CLON_VEC vxx(dv[0],CLONE_COL,0);
  CLON_VEC vyy(dv[1],CLONE_COL,1);
  CLON_VEC vzz(dv[2],CLONE_COL,2);

  keExchange = vxx + vyy + vzz;
  //keExchange *= density;
  keExchange *= temperature;
  keExchange *= capacity;

  flux -= keExchange;
#endif
  return true;
}

}; // end namespace
