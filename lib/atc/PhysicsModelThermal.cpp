#include "ATC_Transfer.h"
#include "PhysicsModelThermal.h"
#include "Material.h"
#include <string>
#include <iostream>
#include <fstream>

namespace ATC {

  PhysicsModelThermal::PhysicsModelThermal(string matFileName,
                                           ATC_Transfer * atcTransfer)
    : PhysicsModel(matFileName,atcTransfer)
  {
    initialize();
  }

  PhysicsModelThermal::~PhysicsModelThermal(void)
  {
  }

  //--------------------------------------------------------------
  //  initialize
  //--------------------------------------------------------------
  void PhysicsModelThermal::initialize(void)
  {
    string list[2] = {"heat_capacity",
                      "heat_flux"};
    set<string> needs(list,list+2);
    vector< Material* >::iterator iter;
    for (iter = materials_.begin(); iter != materials_.end(); iter++) {
      Material * mat = *iter;
      if (! (mat->check_registry(needs)) ) {
        throw ATC_Error(0,"material does not provide all interfaces for physics");
      }
    }
  }
  //---------------------------------------------------------------------
  //   compute heat capacity
  //---------------------------------------------------------------------
  void PhysicsModelThermal::M_integrand(const Array<FieldName> &mask, 
                                        const FIELDS &fields, 
                                        FIELDS &capacity,
                                        const int matIndex) const
  {
    Material* material = materials_[matIndex]; 
    for (int n = 0; n < mask.get_length(); n++) {
      if (mask(n) == TEMPERATURE) {
        material->heat_capacity(fields, capacity[TEMPERATURE]);
      }
    }
  }

  //---------------------------------------------------------------------
  //   compute total energy
  //---------------------------------------------------------------------
  void PhysicsModelThermal::E_integrand(const Array<FieldName> &mask, 
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
    }
  }

  //--------------------------------------------------------------
  //  compute heat flux
  //--------------------------------------------------------------
  void PhysicsModelThermal::B_integrand(const Array2D<bool> & mask, 
					const FIELDS & fields,
					const GRAD_FIELDS & gradFields, 
					GRAD_FIELDS & flux,
          const int matIndex) const
  {
    Material *material = materials_[matIndex];
    material->heat_flux(fields, gradFields, flux[TEMPERATURE]);
  }
};

  
