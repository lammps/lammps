#include "Material.h"
#include "ATC_Transfer.h"
#include "LammpsInterface.h"
#include "ElectronHeatCapacity.h"
#include "ElectronHeatFlux.h"
#include "ElectronPhononExchange.h"

namespace ATC {
using namespace ATC_STRING;

  Material::Material()
    : rhoCp_(0), 
      rho_(0),
      heatCapacity_(0),
      massDensity_(0),
      chargeDensity_(0),
      electronHeatCapacity_(NULL),
      heatConductivity_(0),
      electronHeatFlux_(NULL),
      stress_(NULL),
      electronPhononExchange_(NULL),
      permittivity_(1.),
      electronEquilibriumDensity_(0),
      electronRecombinationInvTau_(0),
      donorConcentration_(0)
  {
  }
  //--------------------------------------------------------------
  //  Constructor (parser)
  //--------------------------------------------------------------
  // Example:
  // material Cu
  //   heat_capacity  constant
  //     capacity  1.0
  //   end
  //   heat_flux linear
  //     conductivity    1.0
  //   end
  //   electron_heat_flux linear
  //     conductivity  1.0
  //   end
  //   electron_heat_capacity linear
  //     capacity 1.0
  //   end
  //   electron_phonon_exchange linear
  //     coefficient   0.1
  //   end
  // end
  Material::Material(string & tag, fstream &fileId)
    : tag_(tag),
      rhoCp_(0), 
      rho_(0),
      heatCapacity_(0),
      massDensity_(0),
      chargeDensity_(0),
      electronHeatCapacity_(NULL),
      heatConductivity_(0),
      electronHeatFlux_(NULL),
      stress_(NULL),
      electronPhononExchange_(NULL),
      permittivity_(1.),
      electronEquilibriumDensity_(0),
      electronRecombinationInvTau_(0),
      donorConcentration_(0)
  {
    linearFlux_.reset(NUM_FIELDS);
    linearFlux_ = false;
    linearSource_.reset(NUM_FIELDS);
    linearSource_ = true;
    constantDensity_.reset(NUM_FIELDS);
    constantDensity_ = false;

    // NOTE these next two rely on many assumptions
    rhoCp_ = ATC::LammpsInterface::instance()->heat_capacity();
    parameters_["heat_capacity"] = rhoCp_;
    heatCapacity_ = rhoCp_;
    registry_.insert("heat_capacity");
    constantDensity_(TEMPERATURE) = true;

    rho_   = ATC::LammpsInterface::instance()->mass_density();
    parameters_["mass_density"] = rho_;
    massDensity_ = rho_;
    registry_.insert("mass_density");
    constantDensity_(DISPLACEMENT) = true;
    constantDensity_(VELOCITY) = true;


    vector<string> line;
    while(fileId.good()) {
      get_command_line(fileId, line);
      if (line.size() == 0) continue;
      if (line.size() == 1) {
        if (line[0] == "end") {
          return;
        }
      }
      if      (line[0] == "heat_capacity") { // over-ride default
        registry_. insert("heat_capacity");
        registry_. insert("thermal_energy");
        if (line[1] == "constant") {
          while(fileId.good()) {
            get_command_line(fileId, line);
            if (line.size() == 0) continue;
            if (line[0] == "end") break;
            double value = str2dbl(line[1]);
            if (line[0] == "capacity") {
              heatCapacity_ = value;
              parameters_["heat_capacity"] = heatCapacity_;
            }
          }
        }
      }
      else if (line[0] == "heat_flux") {
        registry_. insert("heat_flux");
        if (line[1] == "linear") {
          linearFlux_(TEMPERATURE) = true;
          while(fileId.good()) {
            get_command_line(fileId, line);
            if (line.size() == 0) continue;
            if (line[0] == "end") break;
            double value = str2dbl(line[1]);
            if (line[0] == "conductivity") {
              heatConductivity_ = value;
            }
          }
        }
      }
      else if (line[0] == "electron_heat_flux") {
        registry_. insert("electron_heat_flux");
        if      (line[1] == "linear") {
          linearFlux_(ELECTRON_TEMPERATURE) = true;
          electronHeatFlux_ = new ElectronHeatFluxLinear(fileId, parameters_);
        }
        else if (line[1] == "power_law") {
          electronHeatFlux_ = new ElectronHeatFluxPowerLaw(fileId, parameters_);
        }
      }
      else if (line[0] == "electron_heat_capacity") {
        registry_. insert("electron_heat_capacity");
        registry_. insert("electron_thermal_energy");
        if      (line[1] == "constant") {
          constantDensity_(ELECTRON_TEMPERATURE) = true;
          electronHeatCapacity_ = new ElectronHeatCapacityConstant(fileId,
                                                                parameters_);
        }
        else if (line[1] == "linear") {
          electronHeatCapacity_ = new ElectronHeatCapacityLinear(fileId,
                                                                parameters_);
        }
      }
      else if (line[0] == "electron_phonon_exchange") {
        registry_. insert("electron_phonon_exchange");
        if      (line[1] == "linear") {
          electronPhononExchange_ = new ElectronPhononExchangeLinear(fileId,  
                                                                parameters_);
        }
        else if (line[1] == "power_law") {
          linearSource_(TEMPERATURE) = false;
          linearSource_(ELECTRON_TEMPERATURE) = false;
          electronPhononExchange_ = new ElectronPhononExchangePowerLaw(fileId, 
                                                                parameters_);
        }
      }
      else if (line[0] == "mass_density") { // over-ride default
        registry_. insert("mass_density");
        registry_. insert("kinetic_energy");
        if (line[1] == "constant") {
          while(fileId.good()) {
            get_command_line(fileId, line);
            if (line.size() == 0) continue;
            if (line[0] == "end") break;
            double value = str2dbl(line[1]);
            if (line[0] == "density") {
              massDensity_ = value;
              parameters_["mass_density"] = massDensity_;
            }
          }
        }
      }
      else if (line[0] == "electron_recombination") {
        registry_. insert("electron_recombination");
        if (line[1] == "linear") {
          while(fileId.good()) {
            get_command_line(fileId, line);
            if (line.size() == 0) continue;
            if (line[0] == "end") break;
            double value = str2dbl(line[1]);
            if (line[0] == "inv_relaxation_time") {
              electronRecombinationInvTau_ = value;
              parameters_["inv_relaxation_time"] = electronRecombinationInvTau_;
            }
            else if (line[0] == "equilibrium_carrier_density") {
              electronEquilibriumDensity_ = value;
              parameters_["equilibrium_carrier_density"] 
                = electronEquilibriumDensity_;
            }
          }
        }
      }
      else if (line[0] == "charge_density") {
        registry_. insert("charge_density");
        if (line[1] == "constant") {
          while(fileId.good()) {
            get_command_line(fileId, line);
            if (line.size() == 0) continue;
            if (line[0] == "end") break;
            double value = str2dbl(line[1]);
            if (line[0] == "donor_concentration") {
              donorConcentration_ = value;
              parameters_["donor_concentration"] = donorConcentration_;
            }
            else if (line[0] == "density") {
              chargeDensity_ = value;
              parameters_["charge_density"] = chargeDensity_;
            }
          }
        }
      }
      else {
        throw ATC_Error(0, "unrecognized material function: "+line[0]);
      }
    }
  }
//---------------------------------------------------------------------
void Material::heat_capacity(
  const FIELDS & fields,
  FIELD & capacity)
{
  const FIELD & T = (fields.find(TEMPERATURE))->second;
  int nNodes  = T.nRows();
  capacity.reset(nNodes,1);
  capacity = heatCapacity_;
};
//---------------------------------------------------------------------
void Material::thermal_energy(
  const FIELDS &fields,
  FIELD &energy)
{
  const FIELD & T = (fields.find(TEMPERATURE))->second;
  energy = heatCapacity_ * T;
};
//---------------------------------------------------------------------
void Material::electron_heat_capacity(
  const FIELDS & fields,
  FIELD & capacity)
{
  electronHeatCapacity_->electron_heat_capacity(fields,capacity);
};
//---------------------------------------------------------------------
void Material::electron_thermal_energy(
  const FIELDS &fields,
  FIELD &energy)
{
  electronHeatCapacity_->electron_thermal_energy(fields,energy);
};
//---------------------------------------------------------------------
void Material::mass_density(
  const FIELDS &fields,
  FIELD &density)
{
  int nNodes = 0;
  FIELDS::const_iterator field = fields.find(MASS_DENSITY);
  if (field != fields.end()) {
    const FIELD & d = field->second;
    nNodes  = d.nRows();
  }
  else {
    FIELDS::const_iterator field = fields.find(VELOCITY);
    if (field != fields.end()) {
      const FIELD & v = field->second;
      nNodes  = v.nRows();
    }
  }
  density.reset(nNodes,1);
  density = massDensity_;
};
//---------------------------------------------------------------------
void Material::kinetic_energy(
  const FIELDS &fields,
  FIELD &energy)
{
  FIELDS::const_iterator field = fields.find(VELOCITY);
  if (field != fields.end()) {
    const FIELD & v = field->second;
    energy = 0.5*massDensity_*v;
    energy *= v;
  }
  else {
    energy = 0.;
  }
};
//---------------------------------------------------------------------
void Material::heat_flux(
  const FIELDS & fields,
  const GRAD_FIELDS & gradFields,
  GRAD_FIELD & flux)
{
  const GRAD_FIELD & dT = (gradFields.find(TEMPERATURE))->second;
  flux.push_back(-heatConductivity_* dT[0]);
  flux.push_back(-heatConductivity_* dT[1]);
  flux.push_back(-heatConductivity_* dT[2]);
}
//---------------------------------------------------------------------
void Material::electron_heat_flux(
  const FIELDS & fields,
  const GRAD_FIELDS & gradFields,
  GRAD_FIELD & flux)
{
  electronHeatFlux_->electron_heat_flux(fields,gradFields,flux);
}
//---------------------------------------------------------------------
void Material::electron_phonon_exchange(
  const FIELDS & fields,
  FIELD & flux)
{
  electronPhononExchange_->electron_phonon_exchange(fields,flux);
}
//---------------------------------------------------------------------
void Material::electron_recombination(
  const FIELDS &fields,
  const GRAD_FIELDS &gradFields,
  FIELD & flux)
{
  // 1/tau (n - n0)
  const FIELD & n   = (fields.find(ELECTRON_DENSITY))->second;
  flux  = n;
  flux -= electronEquilibriumDensity_;
  flux *= -electronRecombinationInvTau_;
}

//---------------------------------------------------------------------
void Material::charge_density(
  const FIELDS &fields,
  const GRAD_FIELDS &gradFields,
  FIELD & flux)
{
  // Extrinsic/electron charge density 
  FIELDS::const_iterator field = fields.find(ELECTRON_DENSITY);
  if (field != fields.end()) {
    // (n - nD) , where n > 0
    const FIELD & n   = field->second;
    flux  = n;
    flux -= donorConcentration_;// NOTE uniform
    flux  *= -1.0; // account for negative charge
  }
};

} // end namespace

