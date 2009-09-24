#ifndef ELECTRON_HEAT_FLUX_H
#define ELECTRON_HEAT_FLUX_H

#include <map>
#include <string>

using std::map;
using std::string;

#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"
#include "ElectronFlux.h"

namespace ATC {
  class ElectronHeatFlux
  {
    public:
      ElectronHeatFlux() 	{};
      virtual ~ElectronHeatFlux() {};
      /** computes heat flux */
      virtual void electron_heat_flux(const FIELDS &fields,
                                      const GRAD_FIELDS &gradFields,
                                      GRAD_FIELD &flux)=0; 
  };
  //-----------------------------------------------------------------------
  class ElectronHeatFluxLinear : public ElectronHeatFlux
  {
    public:
      ElectronHeatFluxLinear(fstream &matfile,map<string,double> & parameters);
      virtual ~ElectronHeatFluxLinear() {};
      virtual void electron_heat_flux(const FIELDS &fields,
                                      const GRAD_FIELDS &gradFields,
                                      GRAD_FIELD &flux)
      {
         // flux = -ke dTe/dx
         const GRAD_FIELD& dT = (gradFields.find(ELECTRON_TEMPERATURE))->second;
         flux.push_back(-conductivity_ * dT[0]);
         flux.push_back(-conductivity_ * dT[1]);
         flux.push_back(-conductivity_ * dT[2]);
      }; 
    protected:
      double conductivity_;
  };
  //-----------------------------------------------------------------------
  class ElectronHeatFluxPowerLaw : public ElectronHeatFlux
  {
    public:
      ElectronHeatFluxPowerLaw(fstream &matfile,map<string,double> &parameters);
      virtual ~ElectronHeatFluxPowerLaw() {};
      virtual void electron_heat_flux(const FIELDS &fields,
                                      const GRAD_FIELDS &gradFields,
                                      GRAD_FIELD &flux)
      {
         // flux = -ke * ( Te / T ) dT;
         const GRAD_FIELD dT  = (gradFields.find(ELECTRON_TEMPERATURE))->second;
         const FIELD & T  = (fields.find(TEMPERATURE))->second;
         const FIELD & Te = (fields.find(ELECTRON_TEMPERATURE))->second;
         flux.push_back(dT[0]);
         flux.push_back(dT[1]);
         flux.push_back(dT[2]);
         FIELD k_e;
         k_e = (-conductivity_* Te) / T;
         flux[0] *= k_e;
         flux[1] *= k_e;
         flux[2] *= k_e;
      }; 
    protected:
      double conductivity_;
  };
  //-----------------------------------------------------------------------
  class ElectronHeatFluxThermopower : public ElectronHeatFlux
  {
    public:
      ElectronHeatFluxThermopower(fstream &matfile, 
                                  map<string,double> & parameters,
                                  /*const*/ ElectronFlux * electronFlux = NULL);
      virtual ~ElectronHeatFluxThermopower() {};
      virtual void electron_heat_flux(const FIELDS &fields,
                                      const GRAD_FIELDS &gradFields,
                                      GRAD_FIELD &flux)
      {
         // flux = -ke * ( Te / T ) dT + pi J_e;
         const GRAD_FIELD dT  = (gradFields.find(ELECTRON_TEMPERATURE))->second;
         const FIELD & T  = (fields.find(TEMPERATURE))->second;
         const FIELD & Te = (fields.find(ELECTRON_TEMPERATURE))->second;
         flux.push_back(dT[0]);
         flux.push_back(dT[1]);
         flux.push_back(dT[2]);
         FIELD k_e;
         k_e = (-conductivity_* Te) / T;
         flux[0] *= k_e;
         flux[1] *= k_e;
         flux[2] *= k_e;
         GRAD_FIELD tmp;
         electronFlux_->electron_flux(fields, gradFields, tmp);
         tmp[0] *=  Te;
         tmp[1] *=  Te;
         tmp[2] *=  Te;
         flux[0] += seebeckCoef_*tmp[0];
         flux[1] += seebeckCoef_*tmp[1];
         flux[2] += seebeckCoef_*tmp[2]; 

      }; 
    protected:
      double conductivity_,seebeckCoef_;
      ElectronFlux * electronFlux_;
  };
}

#endif 


