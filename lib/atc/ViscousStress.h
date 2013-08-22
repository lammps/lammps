#ifndef VISCOUS_STRESS_H
#define VISCOUS_STRESS_H

#include <map>
#include <string>
#include <vector>
#include <fstream>
#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"

namespace ATC {

  /**
   * @class ViscousStress
   * @brief Base class that defines interface for a constitutive law 
   * @brief that computes viscous stresses given all field and gradient information.
   */
  class ViscousStress
  {
    public:
      ViscousStress()  {};
      virtual ~ViscousStress() {};
      virtual void initialize(void){};
      //* Returns parameter values, (Nothing uses this).
      virtual void parameters(std::map<std::string,double> &parameters) {}
      //* Computes viscous stress given a strain rate tensor.
      //* Units: mvv/L^3 (i.e. for units Real: g/(mol ps^2 A^2) )
      virtual void viscous_stress(const FIELD_MATS &fields,
                                  const GRAD_FIELD_MATS &gradFields,
                                  DENS_MAT_VEC &stress)=0; 
      virtual void viscosity(const FIELD_MATS & fields,
                             DENS_MAT & coefs) const
        {throw ATC_Error("ViscousStress::viscosity: unimplemented function");}
      //* Returns the derivative of the stress tensor for a given strain-rate tensor.
      virtual void tangent(const MATRIX &F, MATRIX &C) const 
        {throw ATC_Error("ViscousStress::tangent: unimplemented function");}
  };


  /**
   *  @class  ViscousStressConstant
   *  @brief  Class for computing stress for a constant viscosity material
   *          assuming divergence-free flow
   */ 

  class ViscousStressConstant : public ViscousStress
  {
    public:
      ViscousStressConstant():viscosity_(0.){};
      ViscousStressConstant(std::fstream &matfile);
      ViscousStressConstant(double viscosity)
        : viscosity_(viscosity) {};
      void viscous_stress(const FIELD_MATS &fields,
                          const GRAD_FIELD_MATS &gradFields,
                          DENS_MAT_VEC &flux);
      void viscosity(const FIELD_MATS & fields,
                     DENS_MAT & coefs) const;
    protected:
      double viscosity_;
      void set_tangent();
  };

}
#endif 
