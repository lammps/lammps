#include "ElectronHeatFlux.h"
#include "StringManip.h"
#include "ATC_Error.h"

#include <iostream>
#include <fstream>

namespace ATC {
using namespace ATC_STRING;

ElectronHeatFluxLinear::ElectronHeatFluxLinear(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronHeatFlux(),
  conductivity_(0)
{
  if (!fileId.is_open()) throw ATC_Error(0,"cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    get_command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    else if (line[0] == "conductivity") {
      conductivity_ = str2dbl(line[1]);
      parameters["electron_thermal_conductivity"] = conductivity_;
    }
    else {
      throw ATC_Error(0, "unrecognized material function "+line[0]);
    }
  }
}

ElectronHeatFluxPowerLaw::ElectronHeatFluxPowerLaw(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronHeatFlux(),
  conductivity_(0)
{
  if (!fileId.is_open()) throw ATC_Error(0,"cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    get_command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    else if (line[0] == "conductivity") {
      conductivity_ = str2dbl(line[1]);
      parameters["electron_thermal_conductivity"] = conductivity_;
    }
    else {
      throw ATC_Error(0, "unrecognized material function "+line[0]);
    }
  }
}

ElectronHeatFluxThermopower::ElectronHeatFluxThermopower(
  fstream &fileId, map<string,double> & parameters,
  /*const*/ ElectronFlux * electronFlux) 
  : ElectronHeatFlux(),
  electronFlux_(electronFlux),
  conductivity_(0),
  seebeckCoef_(0)
{
  if (!fileId.is_open()) throw ATC_Error(0,"cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    get_command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    double value = str2dbl(line[1]);
    if (line[0] == "conductivity") {
      conductivity_ = value;
      parameters["electron_thermal_conductivity"] = conductivity_;
    }
    else {
      throw ATC_Error(0, "unrecognized material function "+line[0]);
    }
    seebeckCoef_ = parameters["seebeck_coefficient"];
  }
}

}

