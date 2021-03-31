#include "ElectronChargeDensity.h"
#include "ATC_Error.h"

#include <iostream>
#include <vector>

using ATC_Utility::command_line;
using ATC_Utility::str2dbl;
using ATC_Utility::str2int;
using std::fstream;
using std::map;
using std::string;
using std::vector;

namespace ATC {
ElectronChargeDensityInterpolation::ElectronChargeDensityInterpolation(
  fstream &fileId, map<string,double> & /* parameters */) 
  : ElectronChargeDensity(), n_()
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  int npts = 0;
  double coef = 1.;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue; 
    if (line[0] == "end") return;
    else if (line[0] == "scale") coef = str2dbl(line[1]);
    else if (line[0] == "number_of_points") {
      npts = str2int(line[1]);
      n_.initialize(npts,fileId,coef);
    }
  }
}

ElectronChargeDensityLinear::ElectronChargeDensityLinear(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronChargeDensity()
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue; 
    if (line[0] == "end") return;
    double value = str2dbl(line[1]);
    if (line[0] == "coefficient") {
      C_ = value;
      parameters["coefficient"] = C_;
    }
  }
}

ElectronChargeDensityExponential::ElectronChargeDensityExponential(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronChargeDensity(),
  intrinsicConcentration_(0),
  intrinsicEnergy_(0),
  referenceTemperature_(0)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue; 
    if (line[0] == "end") return;
    double value = str2dbl(line[1]);
    if (line[0] == "intrinsic_concentration") {
      intrinsicConcentration_ = value;
      parameters["intrinsic_concentration"] = intrinsicConcentration_;
    }
    else if (line[0] == "intrinsic_energy") {
      intrinsicEnergy_ = value;
      parameters["intrinsic_energy"] = intrinsicEnergy_;
    }
    else if (line[0] == "reference_temperature") {
      referenceTemperature_ = value;
      parameters["reference_temperature"] = referenceTemperature_;
    }
    else {
      throw ATC_Error( "unrecognized material function "+line[0]);
    }
  }
}

ElectronChargeDensityFermiDirac::ElectronChargeDensityFermiDirac(
  fstream &fileId, map<string,double> & parameters) 
  : ElectronChargeDensity(),
  Ef_(0),
    referenceTemperature_(0),
    Ed_(0), Nd_(0), Eb_(0),
    hasReferenceTemperature_(false),
    donorIonization_(false)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue; 
    if (line[0] == "end") return;
    double value = str2dbl(line[1]);
    if (line[0] == "fermi_energy") {
      Ef_ = value;
      parameters["fermi_energy"] = Ef_;
    }
    else if (line[0] == "reference_temperature") {
      hasReferenceTemperature_ = true;
      referenceTemperature_ = value;
      parameters["reference_temperature"] = referenceTemperature_;
    }
    else if (line[0] == "band_edge") {
      Eb_ = value;
      parameters["band_edge_potential"] = Eb_;
    }
    else if (line[0] == "donor_ionization_energy") {
      donorIonization_ = true;
      Ed_ = value;
      parameters["donor_ionization_energy"] = Ed_;
    }
    else if (line[0] == "donor_concentration") {
      donorIonization_ = true;
      Nd_ = value;
      parameters["donor_concentration"] = Nd_; 
    }
    else {
      throw ATC_Error( "unrecognized material function "+line[0]);
    }
  }
}

}
