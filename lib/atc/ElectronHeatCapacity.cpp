#include "ElectronHeatCapacity.h"
#include "StringManip.h"
#include "ATC_Error.h"

#include <iostream>
#include <fstream>

namespace ATC {
using namespace ATC_STRING;

ElectronHeatCapacityConstant::ElectronHeatCapacityConstant(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronHeatCapacity(),
  electronHeatCapacity_(0)
{
  if (!fileId.is_open()) throw ATC_Error(0,"cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    get_command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    else if (line[0] == "capacity") {
      electronHeatCapacity_ = str2dbl(line[1]);
      parameters["electron_heat_capacity"] = electronHeatCapacity_;
    }
    else {
      throw ATC_Error(0, "unrecognized material function:" + line[0]);
    }
  }
}

ElectronHeatCapacityLinear::ElectronHeatCapacityLinear(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronHeatCapacity(),
  electronHeatCapacity_(0)
{
  if (!fileId.is_open()) throw ATC_Error(0,"cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    get_command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    else if (line[0] == "capacity") {
      electronHeatCapacity_ = str2dbl(line[1]);
      parameters["electron_heat_capacity"] = electronHeatCapacity_;
    }
    else {
      throw ATC_Error(0, "unrecognized material function: " + line[0]);
    }
  }
}

}

