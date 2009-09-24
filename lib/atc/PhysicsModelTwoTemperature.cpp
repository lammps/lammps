#include "PhysicsModelTwoTemperature.h"
#include "Material.h"
#include <string>
#include <iostream>
#include <fstream>

namespace ATC {
//==============================================================
//  Class PhysicsModelTwoTemperature
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
PhysicsModelTwoTemperature::PhysicsModelTwoTemperature(string matFileName,
  ATC_Transfer * atcTransfer)
  : PhysicsModel(matFileName,atcTransfer)
{
  initialize();
}
//--------------------------------------------------------------
//  Destructor
//---------------------------------------------------------------------
PhysicsModelTwoTemperature::~PhysicsModelTwoTemperature(void)
{}

//---------------------------------------------------------------------
//  intialization 
//---------------------------------------------------------------------
void PhysicsModelTwoTemperature::initialize(void)
{
  string list[5] = {"heat_capacity",
           "electron_heat_capacity",
                    "heat_flux",
           "electron_heat_flux",
           "electron_phonon_exchange"};
  set<string> needs(list,list+5);
  vector< Material* >::iterator iter;
  for (iter = materials_.begin(); iter != materials_.end(); iter++) {
    Material * mat = *iter;
    if (! (mat->check_registry(needs)) ) {
      throw ATC_Error(0,"material " + mat->label() + " does not provide all interfaces for physics");
    }
  }
}

//---------------------------------------------------------------------
//   compute energy
//---------------------------------------------------------------------
void PhysicsModelTwoTemperature::E_integrand(const Array<FieldName> &mask, 
                                             const FIELDS &fields, 
                                             const GRAD_FIELDS &grad_fields,
                                             FIELDS &energy,
                                             const int matIndex) const
{
  Material* material = materials_[matIndex]; 
  for (int n = 0; n < mask.get_length(); n++) {
    if (mask(n) == TEMPERATURE) {
      material->thermal_energy(fields, energy[TEMPERATURE]);
    }
    else if (mask(n) == ELECTRON_TEMPERATURE) {
      material->electron_thermal_energy(fields, energy[ELECTRON_TEMPERATURE]);
    }
  }
}

//---------------------------------------------------------------------
//   compute heat capacities
//---------------------------------------------------------------------
void PhysicsModelTwoTemperature::M_integrand(const Array<FieldName> &mask, 
                                             const FIELDS &fields, 
                                             FIELDS &capacity,
                                             const int matIndex) const
{
  Material* material = materials_[matIndex]; 
  for (int n = 0; n < mask.get_length(); n++) {
    if (mask(n) == TEMPERATURE) {
      material->heat_capacity(fields, capacity[TEMPERATURE]);
    }
    else if (mask(n) == ELECTRON_TEMPERATURE) {
      material->electron_heat_capacity(fields, capacity[ELECTRON_TEMPERATURE]);
    }
  }
}

//---------------------------------------------------------------------
//   compute heat fluxes
//---------------------------------------------------------------------
void PhysicsModelTwoTemperature::B_integrand(const Array2D<bool> & mask, 
                                             const FIELDS &fields,
                                             const GRAD_FIELDS &grad_fields,
                                             GRAD_FIELDS &flux,
                                             const int matIndex) const
{
  Material * material = materials_[matIndex];
  if (mask(TEMPERATURE,FLUX)) {
    material->heat_flux(fields, grad_fields, flux[TEMPERATURE]);
  }
  if (mask(ELECTRON_TEMPERATURE,FLUX)) {
    material->electron_heat_flux(fields, grad_fields, flux[ELECTRON_TEMPERATURE]);
  }
}

//---------------------------------------------------------------------
//   compute exchange fluxes
//---------------------------------------------------------------------
void PhysicsModelTwoTemperature::N_integrand(const Array2D<bool> &mask, 
                                             const FIELDS &fields, 
                                             const GRAD_FIELDS &grad_fields,
                                             FIELDS &flux,
                                             const int matIndex) const
{
  Material * material = materials_[matIndex]; 
  FIELD exchange_flux;
  material->electron_phonon_exchange(fields, exchange_flux);
  if (mask(TEMPERATURE,SOURCE)) {
    flux[TEMPERATURE] = exchange_flux;
  }
  if (mask(ELECTRON_TEMPERATURE,SOURCE)) {
    flux[ELECTRON_TEMPERATURE] = -1.*exchange_flux;
  }
}
}// end namespace
