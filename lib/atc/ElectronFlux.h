#ifndef ELECTRON_FLUX_H
#define ELECTRON_FLUX_H

#include <map>
#include <string>

using std::map;
using std::string;

#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"

namespace ATC {
  class ElectronFlux
  {
    public:
      ElectronFlux() 	{};
      virtual ~ElectronFlux() {};
      /** computes flux */
      virtual void electron_flux(const FIELDS &fields,
                                 const GRAD_FIELDS &gradFields,
                                 GRAD_FIELD &flux)=0; 
  };
  //-----------------------------------------------------------------------
  class ElectronFluxLinear : public ElectronFlux
  {
    public:
      ElectronFluxLinear(fstream &matfile, map<string,double> & parameters);
      virtual ~ElectronFluxLinear() {};
      virtual void electron_flux(const FIELDS &fields,
                                      const GRAD_FIELDS &gradFields,
                                      GRAD_FIELD &flux)
      {
         // J_n = - \mu n grad \phi  - D grad n
         const FIELD & n   = (fields.find(ELECTRON_DENSITY))->second;
         const GRAD_FIELD & dn   =(gradFields.find(ELECTRON_DENSITY))->second;
         const GRAD_FIELD & dphi =(gradFields.find(ELECTRIC_POTENTIAL))->second;
         // NOTE use electron velocity instead
         flux.push_back(-electronMobility_*dphi[0]);
         flux.push_back(-electronMobility_*dphi[1]);
         flux.push_back(-electronMobility_*dphi[2]);   
         flux[0] *= n; // scale by n to get : -n \mu grad(\phi)
         flux[1] *= n;
         flux[2] *= n;
         flux[0] += -electronDiffusivity_* dn[0];
         flux[1] += -electronDiffusivity_* dn[1];
         flux[2] += -electronDiffusivity_* dn[2];
      }; 
    protected:
      double electronMobility_, electronDiffusivity_;
  };
  //-----------------------------------------------------------------------
  class ElectronFluxThermopower : public ElectronFlux
  {
    public:
      ElectronFluxThermopower(fstream &matfile,map<string,double> & parameters);
      virtual ~ElectronFluxThermopower() {};
      virtual void electron_flux(const FIELDS &fields,
                                      const GRAD_FIELDS &gradFields,
                                      GRAD_FIELD &flux)
      {
         static const double kB_ = 8.617343e-5;// [eV/K]
         // J_n = - \mu n grad \phi  - \mu kB/e T_e grad n 
         //       - \mu S n grad T_e - \mu kB/e n grad T_e
         const FIELD & n   = (fields.find(ELECTRON_DENSITY))->second;
         const GRAD_FIELD & dn   =(gradFields.find(ELECTRON_DENSITY))->second;
         const GRAD_FIELD & dphi =(gradFields.find(ELECTRIC_POTENTIAL))->second;
         const GRAD_FIELD & dT =(gradFields.find(ELECTRON_TEMPERATURE))->second;
         flux.push_back(-electronMobility_*dphi[0]);
         flux.push_back(-electronMobility_*dphi[1]);
         flux.push_back(-electronMobility_*dphi[2]);   
         double coef = -electronMobility_*(seebeckCoef_ + kB_);
         flux[0] += coef* dT[0];
         flux[1] += coef* dT[1];
         flux[2] += coef* dT[2];
         flux[0] *= n; // scale by n 
         flux[1] *= n;
         flux[2] *= n;
         GRAD_FIELD tmp = dn;
         const FIELD & Te   = (fields.find(ELECTRON_TEMPERATURE))->second;
         tmp[0] *= Te;
         tmp[1] *= Te;
         tmp[2] *= Te;
         coef = -electronMobility_*kB_; // kB: eV/K/e -> V/K
         flux[0] += tmp[0];
         flux[1] += tmp[1];
         flux[2] += tmp[2];
      }; 
    protected:
      double electronMobility_, seebeckCoef_;
  };
}

#endif 


