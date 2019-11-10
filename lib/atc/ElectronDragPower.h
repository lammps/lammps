#ifndef ELECTRON_DRAG_POWER_H
#define ELECTRON_DRAG_POWER_H

#include <map>
#include <string>
#include <fstream>

#include "Material.h"
#include "ATC_TypeDefs.h"

namespace ATC {

  /**
   *  @class  ElectronDragPower 
   *  @brief  Base class for defining the lattice drag power from electrons
   */

  class ElectronDragPower
  {
    public:
    ElectronDragPower() {};
    virtual ~ElectronDragPower() {};
    /** computes drag power */
    virtual bool electron_drag_power(const FIELD_MATS & /* fields */,
                                     const GRAD_FIELD_MATS & /* gradFields */,
                                     DENS_MAT & /* flux */)
    {
      return false;
    };
    virtual void electron_drag_velocity_coefficient(const FIELD_MATS &fields,
                                                          DENS_MAT & dragCoef)
    {
      FIELD_MATS::const_iterator t_field = fields.find(TEMPERATURE);
      dragCoef.reset((t_field->second).nRows(),1); // zero out matrix, resize if necessary
    };
  };
  //-------------------------------------------------------------------
  
  /**
   *  @class  ElectronDragPowerLinear
   *  @brief  Class for electron drag that linearly depends on the difference between the electron and lattice velocities
   */  
  
  class ElectronDragPowerLinear : public ElectronDragPower
  {
    public:
      ElectronDragPowerLinear(std::fstream &matfile,
                              std::map<std::string,double> & parameters,
                              Material * material_);
      virtual ~ElectronDragPowerLinear() {};
      virtual bool electron_drag_power(const FIELD_MATS &fields,
                                       const GRAD_FIELD_MATS &gradFields,
                                             DENS_MAT & flux);
      virtual void electron_drag_velocity_coefficient(const FIELD_MATS &fields,
                                                            DENS_MAT & dragCoef);

  protected:

      double electronDragInvTau_;
      Material * material_;

      // used to avoid unnecessary resizing
      DENS_MAT dragCoefWorkspace_;
      DENS_MAT invEffMassWorkspace_;

  };
}

#endif 


