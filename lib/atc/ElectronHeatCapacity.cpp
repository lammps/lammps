#include "ElectronHeatCapacity.h"
#include "ATC_Error.h"

#include <iostream>
#include <fstream>
#include <vector>

using ATC_Utility::command_line;
using ATC_Utility::str2dbl;
using std::fstream;
using std::map;
using std::string;
using std::vector;

namespace ATC {

ElectronHeatCapacityConstant::ElectronHeatCapacityConstant(
  fstream &fileId, map<string,double> & parameters)
  : ElectronHeatCapacity(),
  electronHeatCapacity_(0)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    else if (line[0] == "capacity") {
      electronHeatCapacity_ = str2dbl(line[1]);
      parameters["electron_heat_capacity"] = electronHeatCapacity_;
    }
    else {
      throw ATC_Error( "unrecognized material function:" + line[0]);
    }
  }
}

ElectronHeatCapacityLinear::ElectronHeatCapacityLinear(
  fstream &fileId, map<string,double> & parameters)
  : ElectronHeatCapacity(),
  electronHeatCapacity_(0)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    else if (line[0] == "capacity") {
      electronHeatCapacity_ = str2dbl(line[1]);
      parameters["electron_heat_capacity"] = electronHeatCapacity_;
    }
    else {
      throw ATC_Error( "unrecognized material function: " + line[0]);
    }
  }
}

ElectronHeatCapacityConstantAddDensity::ElectronHeatCapacityConstantAddDensity(fstream &fileId,
                                                                               map<string,double> & parameters,
                                                                               Material * material)
  : ElectronHeatCapacityConstant(fileId, parameters),
    material_(material)
{
  // do nothing
}

ElectronHeatCapacityLinearAddDensity::ElectronHeatCapacityLinearAddDensity(fstream &fileId,
                                                                           map<string,double> & parameters,
                                                                           Material * material)
  : ElectronHeatCapacityLinear(fileId, parameters),
    material_(material)
{
  // do nothing
}

}

