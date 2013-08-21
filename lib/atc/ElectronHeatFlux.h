#ifndef ELECTRON_HEAT_FLUX_H
#define ELECTRON_HEAT_FLUX_H

#include <map>
#include <string>
#include "ATC_TypeDefs.h"
#include "ElectronFlux.h"
#include "ElectronHeatCapacity.h"

namespace ATC {

  /**
   *  @class  ElectronHeatFlux
   *  @brief  Base class for the electron heat flux
   */

  class ElectronHeatFlux
  {
    public:
      ElectronHeatFlux(/*const*/ ElectronHeatCapacity * electronHeatCapacity = NULL);
      virtual ~ElectronHeatFlux() {};
      /** computes heat flux */
      virtual void electron_heat_flux(const FIELD_MATS &fields,
                                      const GRAD_FIELD_MATS &gradFields,
                                            DENS_MAT_VEC &flux)
      {
         
        FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
        const DENS_MAT & Te = etField->second;
        zeroWorkspace_.reset(Te.nRows(),Te.nCols());
        flux[0] = zeroWorkspace_;
        flux[1] = zeroWorkspace_;
        flux[2] = zeroWorkspace_;
      };
      void electron_heat_convection(const FIELD_MATS &fields,
                                          DENS_MAT_VEC & flux)
      {
        FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
        FIELD_MATS::const_iterator evField = fields.find(ELECTRON_VELOCITY);
        const DENS_MAT & Te = etField->second;
        const DENS_MAT & v = evField->second;
        electronHeatCapacity_->electron_heat_capacity(fields,cpTeWorkspace_);
        cpTeWorkspace_ *= Te;
        const CLON_VEC vx(v,CLONE_COL,0);
        const CLON_VEC vy(v,CLONE_COL,1);
        const CLON_VEC vz(v,CLONE_COL,2);
        flux[0] = vx;
        flux[1] = vy;
        flux[2] = vz;
        // scale by thermal energy 
        flux[0] *= cpTeWorkspace_;
        flux[1] *= cpTeWorkspace_;
        flux[2] *= cpTeWorkspace_;
      };
  protected:
      ElectronHeatCapacity * electronHeatCapacity_;
      DENS_MAT zeroWorkspace_;
      DENS_MAT cpTeWorkspace_; // hopefully avoid resizing
  };
  //-----------------------------------------------------------------------
  
  /**
   *  @class  ElectronHeatFluxLinear
   *  @brief  Class for an electron heat flux proportional to the temperature gradient with constant conductivity
   */  
  
  class ElectronHeatFluxLinear : public ElectronHeatFlux
  {
    public:
    ElectronHeatFluxLinear(std::fstream &matfile,std::map<std::string,double> & parameters,
                           /*const*/ ElectronHeatCapacity * electronHeatCapacity = NULL);
      virtual ~ElectronHeatFluxLinear() {};
      virtual void electron_heat_flux(const FIELD_MATS &fields,
                                      const GRAD_FIELD_MATS &gradFields,
                                            DENS_MAT_VEC &flux)
      {
        GRAD_FIELD_MATS::const_iterator dEtField = gradFields.find(ELECTRON_TEMPERATURE);
         // flux = -ke dTe/dx
         const DENS_MAT_VEC & dT = dEtField->second;
         flux[0] = -conductivity_ * dT[0];
         flux[1] = -conductivity_ * dT[1];
         flux[2] = -conductivity_ * dT[2];
      }; 
    protected:
      double conductivity_;
  };
  //-----------------------------------------------------------------------
  
  /**
   *  @class  ElectronHeatFluxPowerLaw
   *  @brief  Class for an electron heat flux proportional to the temperature gradient but with a conductivity proportional to the ratio of the electron and phonon temperatures
   */  
  
  class ElectronHeatFluxPowerLaw : public ElectronHeatFlux
  {
    public:
    ElectronHeatFluxPowerLaw(std::fstream &matfile,std::map<std::string,double> &parameters,
                             /*const*/ ElectronHeatCapacity * electronHeatCapacity = NULL);
      virtual ~ElectronHeatFluxPowerLaw() {};
      virtual void electron_heat_flux(const FIELD_MATS &fields,
                                      const GRAD_FIELD_MATS &gradFields,
                                            DENS_MAT_VEC &flux)
      {
        FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
        FIELD_MATS::const_iterator tField  = fields.find(TEMPERATURE);
        GRAD_FIELD_MATS::const_iterator dEtField = gradFields.find(ELECTRON_TEMPERATURE);
        const DENS_MAT_VEC & dT = dEtField->second;
        const DENS_MAT & T = tField->second;
        const DENS_MAT & Te = etField->second;

        // flux = -ke * ( Te / T ) dT;
        flux[0] = dT[0];
        flux[1] = dT[1];
        flux[2] = dT[2];
        electronConductivity_ = (-conductivity_* Te) / T;
        flux[0] *= electronConductivity_;
        flux[1] *= electronConductivity_;
        flux[2] *= electronConductivity_;
      }; 
    protected:
      double conductivity_;
      DENS_MAT electronConductivity_; // hopefully avoid resizing
  };
  //-----------------------------------------------------------------------
  
  /**
   *  @class  ElectronHeatFluxThermopower
   *  @brief  Class for an electron heat flux proportional to the temperature gradient but with a condu
ctivity proportional to the ratio of the electron and phonon temperatures with the thermopower from teh electric current included
   */  
  
  class ElectronHeatFluxThermopower : public ElectronHeatFlux
  {
    public:
      ElectronHeatFluxThermopower(std::fstream &matfile, 
                                  std::map<std::string,double> & parameters,
                                  /*const*/ ElectronFlux * electronFlux = NULL,
                                  /*const*/ ElectronHeatCapacity * electronHeatCapacity = NULL);
      virtual ~ElectronHeatFluxThermopower() {};
      virtual void electron_heat_flux(const FIELD_MATS &fields,
                                      const GRAD_FIELD_MATS &gradFields,
                                            DENS_MAT_VEC &flux)
      {
        FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
        FIELD_MATS::const_iterator tField  = fields.find(TEMPERATURE);
        GRAD_FIELD_MATS::const_iterator dEtField = gradFields.find(ELECTRON_TEMPERATURE);
        const DENS_MAT_VEC & dT = dEtField->second;
        const DENS_MAT & T = tField->second;
        const DENS_MAT & Te = etField->second;
        
        // flux = -ke * ( Te / T ) dT + pi J_e;
        flux[0] = dT[0];
        flux[1] = dT[1];
        flux[2] = dT[2];
        elecCondWorkspace_ = (-conductivity_* Te) / T;
        flux[0] *= elecCondWorkspace_;
        flux[1] *= elecCondWorkspace_;
        flux[2] *= elecCondWorkspace_;

        electronFlux_->electron_flux(fields, gradFields, tmp_);
        tmp_[0] *=  Te;
        tmp_[1] *=  Te;
        tmp_[2] *=  Te;
        flux[0] += seebeckCoef_*tmp_[0];
        flux[1] += seebeckCoef_*tmp_[1];
        flux[2] += seebeckCoef_*tmp_[2]; 
      }; 
    protected:
      double conductivity_,seebeckCoef_;
      ElectronFlux * electronFlux_;
      DENS_MAT elecCondWorkspace_; // hopefully avoid resizing
      DENS_MAT_VEC tmp_; 
  };
}

#endif 


