#include "ViscousStress.h"
#include "ATC_Error.h"
#include "LammpsInterface.h"
#include <iostream>

using ATC_Utility::command_line;
using ATC_Utility::str2dbl;
using std::string;
using std::vector;
using std::fstream;

namespace ATC {

//=============================================================================
// isotropic constant viscosity
//=============================================================================
ViscousStressConstant::ViscousStressConstant(fstream &fileId) 
  : ViscousStress(), viscosity_(0)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line[0] == "end") {
      if (viscosity_ < 0.0) 
        throw ATC_Error("ViscousStressConstant:: bad constant viscosity");
      return;
    }
    else if (line[0]=="viscosity")        viscosity_ = str2dbl(line[1]);
    else throw ATC_Error( "viscosity constant - unrecognized material function");
  }
}
//=============================================================================
// compute the stress at N integration points from the velocity gradients
// T_{ij} = viscosity * du_i/dx_j
//=============================================================================
void ViscousStressConstant::viscous_stress(const FIELD_MATS      &fields,
                                           const GRAD_FIELD_MATS &gradFields,
                                           DENS_MAT_VEC &sigma)
{
  GRAD_FIELD_MATS::const_iterator du_itr = gradFields.find(VELOCITY);
  const DENS_MAT_VEC &du = du_itr->second;

  CLON_VEC DuxDx(du[0],CLONE_COL,0);
  CLON_VEC DuxDy(du[1],CLONE_COL,0);
  CLON_VEC DuxDz(du[2],CLONE_COL,0);
  CLON_VEC DuyDx(du[0],CLONE_COL,1);
  CLON_VEC DuyDy(du[1],CLONE_COL,1);
  CLON_VEC DuyDz(du[2],CLONE_COL,1);
  CLON_VEC DuzDx(du[0],CLONE_COL,2);
  CLON_VEC DuzDy(du[1],CLONE_COL,2);
  CLON_VEC DuzDz(du[2],CLONE_COL,2);

  const INDEX N = DuxDx.size();          // # of integration pts
  sigma.assign(3, DENS_MAT(N,3));

  // precompute the pressure and copy to the diagonal
  column(sigma[0],0) = DuxDx;
  column(sigma[0],1) = DuyDx;
  column(sigma[0],2) = DuzDx;
  column(sigma[1],0) = DuxDy;
  column(sigma[1],1) = DuyDy;
  column(sigma[1],2) = DuzDy;
  column(sigma[2],0) = DuxDz;
  column(sigma[2],1) = DuyDz;
  column(sigma[2],2) = DuzDz;

  sigma[0] *= -1.*viscosity_;
  sigma[1] *= -1.*viscosity_;
  sigma[2] *= -1.*viscosity_;
}

void ViscousStressConstant::viscosity(const FIELD_MATS &fields,
                                      DENS_MAT &coefs) const
{
  const DENS_MAT & v = (fields.find(VELOCITY))->second;
  
  coefs.resize(v.nRows(),v.nCols());
  coefs = -1.*viscosity_;
}

}// end atc namespace
