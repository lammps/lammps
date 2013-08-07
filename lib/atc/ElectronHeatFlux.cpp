#include "ElectronHeatFlux.h"
#include "ATC_Error.h"

#include <iostream>
#include <fstream>

using ATC_Utility::command_line;
using ATC_Utility::str2dbl;

namespace ATC {

ElectronHeatFlux::ElectronHeatFlux(ElectronHeatCapacity * electronHeatCapacity)
  :
  electronHeatCapacity_(electronHeatCapacity)
{
  // do nothing
}

ElectronHeatFluxLinear::ElectronHeatFluxLinear(fstream &fileId, map<string,double> & parameters,
                                               ElectronHeatCapacity * electronHeatCapacity) 
  : ElectronHeatFlux(electronHeatCapacity),
  conductivity_(0)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    else if (line[0] == "conductivity") {
      conductivity_ = str2dbl(line[1]);
      parameters["electron_thermal_conductivity"] = conductivity_;
    }
    else {
      throw ATC_Error( "unrecognized material function "+line[0]);
    }
  }
}

ElectronHeatFluxPowerLaw::ElectronHeatFluxPowerLaw(fstream &fileId, map<string,double> & parameters,
                                                   ElectronHeatCapacity * electronHeatCapacity) 
  : ElectronHeatFlux(electronHeatCapacity),
  conductivity_(0)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    else if (line[0] == "conductivity") {
      conductivity_ = str2dbl(line[1]);
      parameters["electron_thermal_conductivity"] = conductivity_;
    }
    else {
      throw ATC_Error( "unrecognized material function "+line[0]);
    }
  }
}

ElectronHeatFluxThermopower::ElectronHeatFluxThermopower(
  fstream &fileId, map<string,double> & parameters,
  /*const*/ ElectronFlux * electronFlux,
  ElectronHeatCapacity * electronHeatCapacity) 
  : ElectronHeatFlux(electronHeatCapacity),
    conductivity_(0),
    seebeckCoef_(0),
    electronFlux_(electronFlux)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    double value = str2dbl(line[1]);
    if (line[0] == "conductivity") {
      conductivity_ = value;
      parameters["electron_thermal_conductivity"] = conductivity_;
    }
    else {
      throw ATC_Error( "unrecognized material function "+line[0]);
    }
    
    seebeckCoef_ = parameters["seebeck_coefficient"];
  }
}

}

