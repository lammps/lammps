#include "ATC_Transfer.h"
#include "WeakEquationPhononTemperature.h"
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
WeakEquationPhononTemperature::WeakEquationPhononTemperature()
  : WeakEquation(DYNAMIC_PDE,TEMPERATURE,1)
{}
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
WeakEquationPhononTemperature::~WeakEquationPhononTemperature(void)
{}
//---------------------------------------------------------------------
//   compute total energy
//---------------------------------------------------------------------
void WeakEquationPhononTemperature::E_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS & /* gradFields */,
  const Material * material,
  DENS_MAT &energy) const
{
  material->thermal_energy(fields, energy);
}
//---------------------------------------------------------------------
//   compute heat capacity
//---------------------------------------------------------------------
void WeakEquationPhononTemperature::M_integrand(
  const FIELD_MATS &fields,
  const Material * material,
  DENS_MAT & capacity ) const
{
  material->heat_capacity(fields, capacity);
}
//--------------------------------------------------------------
//  compute heat flux
//--------------------------------------------------------------
void WeakEquationPhononTemperature::B_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &gradFields,
  const Material * material,
  DENS_MAT_VEC &flux) const
{
  material->heat_flux(fields, gradFields, flux);
}

//==============================================================
//  Class WeakEquationPhononTemperatureExchange
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationPhononTemperatureExchange::WeakEquationPhononTemperatureExchange()
  : WeakEquationPhononTemperature()
{}
//--------------------------------------------------------------
//  Destructor
//---------------------------------------------------------------------
WeakEquationPhononTemperatureExchange::~WeakEquationPhononTemperatureExchange(void)
{}

//---------------------------------------------------------------------
//   compute energy
//---------------------------------------------------------------------
void WeakEquationPhononTemperatureExchange::E_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &gradFields,
  const Material * material,
  DENS_MAT &energy ) const
{
  WeakEquationPhononTemperature::E_integrand(fields,gradFields,material,energy);
}

//---------------------------------------------------------------------
//   compute heat capacities
//---------------------------------------------------------------------
void WeakEquationPhononTemperatureExchange::M_integrand(
  const FIELD_MATS &fields,
  const Material * material,
  DENS_MAT & density ) const
{
  WeakEquationPhononTemperature::M_integrand(fields,material,density);
}

//---------------------------------------------------------------------
//   compute heat fluxes
//---------------------------------------------------------------------
void WeakEquationPhononTemperatureExchange::B_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &gradFields,
  const Material * material,
  DENS_MAT_VEC &flux) const
{
  WeakEquationPhononTemperature::B_integrand(fields,gradFields,material,flux);
}

//---------------------------------------------------------------------
//   compute exchange fluxes
//---------------------------------------------------------------------
bool WeakEquationPhononTemperatureExchange::N_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &gradFields,
  const Material * material,
  DENS_MAT &flux) const
{
  bool has =   material->electron_phonon_exchange(fields, flux);
  has = has || material->electron_drag_power(fields, gradFields, flux);
  return has;
}

}; // end namespace


