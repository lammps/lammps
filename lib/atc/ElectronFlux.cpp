#include "ElectronFlux.h"
#include "ATC_Error.h"

#include <iostream>
#include <fstream>

using ATC_Utility::command_line;
using ATC_Utility::str2dbl;

namespace ATC {

ElectronFlux::ElectronFlux() :
  maskX_(false),maskY_(false),maskZ_(false)
{}


ElectronFluxLinear::ElectronFluxLinear(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronFlux(), 
  electronMobility_(0),
  electronDiffusivity_(0)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
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
    else if (line[0] == "mask_x") { maskX_ = true; }
    else if (line[0] == "mask_y") { maskY_ = true; }
    else if (line[0] == "mask_z") { maskZ_ = true; }
    else {
      throw ATC_Error( "unrecognized material function "+line[0]);
    }
  }
}

ElectronFluxThermopower::ElectronFluxThermopower(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronFlux(),
  electronMobility_(0),
  seebeckCoef_(0)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
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
      throw ATC_Error( "unrecognized material function "+line[0]);
    }
  }
}

ElectronFluxConvection::ElectronFluxConvection(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronFlux()
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    else if (line[0] == "mask_x") { maskX_ = true; }
    else if (line[0] == "mask_y") { maskY_ = true; }
    else if (line[0] == "mask_z") { maskZ_ = true; }
    else {
      throw ATC_Error( "unrecognized material function "+line[0]);
    }
  }
}

}

