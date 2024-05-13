#include "BodyForce.h"
#include "ATC_Error.h"

#include <iostream>
#include <fstream>
#include <vector>

using ATC_Utility::command_line;
using ATC_Utility::str2dbl;
using ATC_Utility::str2int;
using std::fstream;
using std::string;
using std::map;
using std::vector;

namespace ATC {

BodyForceViscous::BodyForceViscous(
  fstream &fileId, std::map<std::string,double> & parameters)
  : BodyForce(), gamma_(0)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  std::vector<std::string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    double value = str2dbl(line[1]);
    if (line[0] == "gamma") {
      gamma_ = value;
      parameters["gamma"] = gamma_;
    }

  }
}

}
