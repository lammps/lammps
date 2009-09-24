#include "ElectronPhononExchange.h"
#include "StringManip.h"
#include "ATC_Error.h"

#include <iostream>
#include <fstream>

namespace ATC {
using namespace ATC_STRING;

ElectronPhononExchangeLinear::ElectronPhononExchangeLinear(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronPhononExchange(),
  exchangeCoef_(0)
{
  if (!fileId.is_open()) throw ATC_Error(0,"cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    get_command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") return;
    else if (line[0] == "coefficient") {
      exchangeCoef_ = str2dbl(line[1]);
      parameters["electron_phonon_exchange_coefficient"] = exchangeCoef_;
    }
    else {
      throw ATC_Error(0, "unrecognized material function "+line[0]);
    }
  }
}

ElectronPhononExchangePowerLaw::ElectronPhononExchangePowerLaw(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronPhononExchange(),
  exchangeCoef_(0),
  exponent_(1)
{
  if (!fileId.is_open()) throw ATC_Error(0,"cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    get_command_line(fileId, line);
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
      throw ATC_Error(0, "unrecognized material function "+line[0]);
    }
  }
}

}

