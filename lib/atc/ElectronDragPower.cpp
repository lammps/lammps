#include "ElectronDragPower.h"
#include "Material.h"
#include "ATC_Error.h"

#include <iostream>
#include <vector>

using ATC_Utility::command_line;
using ATC_Utility::str2dbl;
using std::map;
using std::string;
using std::fstream;
using std::vector;

namespace ATC {

ElectronDragPowerLinear::ElectronDragPowerLinear(fstream &fileId,
                                                 map<string,double> & parameters,
                                                 Material * material)
  : ElectronDragPower(),
    electronDragInvTau_(0),
    material_(material)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line.size() == 0) continue;
    if (line[0] == "end") break;
    double value = str2dbl(line[1]);
    if (line[0] == "inv_momentum_relaxation_time") {
      electronDragInvTau_ = value;
      parameters["inv_momentum_relaxation_time"] = electronDragInvTau_;
    }
    else {
      throw ATC_Error( "unrecognized material function "+line[0]);
    }
  }
}

  bool ElectronDragPowerLinear::electron_drag_power(const FIELD_MATS &fields,
                                                    const GRAD_FIELD_MATS & /* gradFields */,
                                                    DENS_MAT & flux)
{

  FIELD_MATS::const_iterator evField = fields.find(ELECTRON_VELOCITY);
  const DENS_MAT & v = evField->second;

  // -1/tau (m_e * n * v)
  electron_drag_velocity_coefficient(fields,dragCoefWorkspace_);
  for (int i = 0; i < v.nRows(); i++) {
    double velocityMagnitude = 0.;
    for (int j = 0; j < v.nCols(); j++)
      velocityMagnitude -= v(i,j)*v(i,j);
    flux(i,0) += velocityMagnitude*dragCoefWorkspace_(i,0); // adds flux to phonon temperature
  }

  return true;
}

void ElectronDragPowerLinear::electron_drag_velocity_coefficient(const FIELD_MATS &fields,
                                                                 DENS_MAT & dragCoef)
{
  FIELD_MATS::const_iterator enField = fields.find(ELECTRON_DENSITY);
  const DENS_MAT & n = enField->second;
  //n.print("DENS");
  // -1/tau (m_e * n)
  material_->inv_effective_mass(fields,invEffMassWorkspace_);
  //invEffMassWorkspace_.print("INV MASS");
  dragCoef = n;
  dragCoef /= invEffMassWorkspace_;
  dragCoef *= -electronDragInvTau_;
  //dragCoef.print("DRAG COEF");
}

}

