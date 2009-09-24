#include "ElectronFlux.h"
#include "StringManip.h"
#include "ATC_Error.h"

#include <iostream>
#include <fstream>

namespace ATC {
using namespace ATC_STRING;

ElectronFluxLinear::ElectronFluxLinear(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronFlux(), 
  electronMobility_(0),
  electronDiffusivity_(0)
{
  if (!fileId.is_open()) throw ATC_Error(0,"cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    get_command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    double value = str2dbl(line[1]);
    if (line[0] == "mobility") {
      electronMobility_ = value;
      parameters["electron_mobility"] = electronMobility_;
    }
    else if (line[0] == "diffusivity") {
      electronDiffusivity_ = value;
      parameters["electron_diffusivity"] = electronDiffusivity_;
    }
    else {
      throw ATC_Error(0, "unrecognized material function "+line[0]);
    }
  }
}

ElectronFluxThermopower::ElectronFluxThermopower(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronFlux(),
  electronMobility_(0),
  seebeckCoef_(0)
{
  if (!fileId.is_open()) throw ATC_Error(0,"cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    get_command_line(fileId, line);
    if (line.size() == 0) continue; 
    if (line[0] == "end") return;
    double value = str2dbl(line[1]);
    if (line[0] == "mobility") {
      electronMobility_ = value;
      parameters["electron_mobility"] = electronMobility_;
    }
    else if (line[0] == "seebeck") {
      seebeckCoef_ = value;
      parameters["seebeck_coefficient"] = seebeckCoef_;
    }
    else {
      throw ATC_Error(0, "unrecognized material function "+line[0]);
    }
  }
}

}

