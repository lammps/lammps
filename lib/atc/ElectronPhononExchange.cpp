#include "Material.h"
#include "ElectronPhononExchange.h"
#include "ATC_Error.h"
#include "LammpsInterface.h"

#include <iostream>
#include <fstream>
#include <math.h>

using ATC_Utility::command_line;
using ATC_Utility::str2dbl;
using ATC_Utility::str2int;

namespace ATC {

ElectronPhononExchangeLinear::ElectronPhononExchangeLinear(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronPhononExchange(),
  exchangeCoef_(0)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    else if (line[0] == "coefficient") {
      exchangeCoef_ = str2dbl(line[1]);
      parameters["electron_phonon_exchange_coefficient"] = exchangeCoef_;
    }
    else {
      throw ATC_Error( "unrecognized material function "+line[0]);
    }
  }
}

ElectronPhononExchangePowerLaw::ElectronPhononExchangePowerLaw(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronPhononExchange(),
  exchangeCoef_(0),
  exponent_(1)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    else if (line[0] == "coefficient") {
      exchangeCoef_ = str2dbl(line[1]);
      parameters["electron_phonon_exchange_coefficient"] = exchangeCoef_;
    }
    else if (line[0] == "exponent") {
      exponent_ = str2int(line[1]);
    }
    else {
      throw ATC_Error( "unrecognized material function "+line[0]);
    }
  }
}

ElectronPhononExchangeHertel::ElectronPhononExchangeHertel(fstream &fileId,
                                                           map<string,double> & parameters,
                                                           Material * material) 
  : ElectronPhononExchange(),
    exchangeCoef_(0),
    debeyeTemperature_(1),
    massEnhancement_(0),
    material_(material)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") break;
    else if (line[0] == "debeye_temperature") {
      debeyeTemperature_ = str2dbl(line[1]);
      parameters["debeye_temperature"] = debeyeTemperature_;
    }
    else if (line[0] == "mass_enhancement") {
      massEnhancement_ = str2dbl(line[1]);
      parameters["mass_enhancement"] = massEnhancement_;
    }
    else {
      throw ATC_Error( "unrecognized material function "+line[0]);
    }
  }
  // coupling coefficient, eqn. 15 of Hertel 2002
  double kb = LammpsInterface::instance()->kBoltzmann();
  double hbar = LammpsInterface::instance()->hbar();
  double PI = 3.141592653589793238; 
  exchangeCoef_ = 144.*1.0369*kb/(PI*hbar);
  exchangeCoef_ *= massEnhancement_/pow(debeyeTemperature_,2);
}

bool ElectronPhononExchangeHertel::electron_phonon_exchange(const FIELD_MATS &fields,
                                                            DENS_MAT &flux)
{
  FIELD_MATS::const_iterator tField = fields.find(TEMPERATURE);
  FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
  const DENS_MAT & T  = tField->second;
  const DENS_MAT & Te = etField->second;

  // flux = g C_e (T_e - T_p)^5 / T_e
  flux = (Te - T).pow(5);
  flux /= Te;
  flux *= exchangeCoef_;
  material_->electron_heat_capacity(fields,capacityWorkspace_);
  flux*= capacityWorkspace_;
  return true;
}

}

