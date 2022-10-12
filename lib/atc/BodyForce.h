#ifndef BODY_FORCE_H
#define BODY_FORCE_H

#include <map>
#include <string>
#include "ATC_TypeDefs.h"
#include "Function.h"

namespace ATC {

  /**
   *  @class  BodyForce
   *  @brief  Base class for models of body forces in the momentum eqn
   */

  class BodyForce
  {
    public:
      BodyForce()   {};
      virtual ~BodyForce() {};
      virtual bool body_force(const FIELD_MATS & /* fields */,
                              DENS_MAT & /* flux */) const { return false; };
  };

  /**
   *  @class  BodyForceViscous
   *  @brief  viscous body forces
   */
  class BodyForceViscous : public BodyForce
  {
    public:
      BodyForceViscous(std::fstream &matfile,std::map<std::string,double> & parameters);
      virtual ~BodyForceViscous() {};
      virtual bool body_force(const FIELD_MATS &fields,
                                    DENS_MAT &flux) const
      {
        FIELD_MATS::const_iterator v_field = fields.find(VELOCITY);
        const DENS_MAT & v = v_field->second;
        flux = -gamma_*v;
        return true;
       }
     protected:
       double gamma_;
  };
  /**
   *  @class  BodyForceElectricField
   *  @brief  electric field body forces
   */
  class BodyForceElectricField : public BodyForce
  {
    public:
    BodyForceElectricField(std::fstream & /* matfile */,std::map<std::string,double> & /* parameters */)
        { throw ATC_Error("unimplemented due to issues with accessing electric field"); }
      virtual ~BodyForceElectricField() {};
      virtual bool body_force(const FIELD_MATS &fields,
                                    DENS_MAT &flux) const
      {
        FIELD_MATS::const_iterator v_field = fields.find(VELOCITY);
        const DENS_MAT & v = v_field->second;
        int nNodes  = v.nRows();
        flux.reset(nNodes,1);
        return true;
       }
  };
}
#endif


