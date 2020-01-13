#ifndef ELECTRON_FLUX_H
#define ELECTRON_FLUX_H

#include <map>
#include <string>
#include "ATC_TypeDefs.h"

namespace ATC {

  /**
   *  @class  ElectronFlux
   *  @brief  Base class for the flux appearing in the electron density transport equation
   */

  class ElectronFlux
  {
    public:
      ElectronFlux();
      virtual ~ElectronFlux() {};
      /** computes flux */
      virtual void electron_flux(const FIELD_MATS &fields,
                                 const GRAD_FIELD_MATS & /* gradFields */,
                                       DENS_MAT_VEC &flux) 
      {
         
        
        
        
        FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
        const DENS_MAT & etMat = etField->second;
        zeroWorkspace_.reset(etMat.nRows(),etMat.nCols());
        flux[0] = zeroWorkspace_;
        flux[1] = zeroWorkspace_;
        flux[2] = zeroWorkspace_;
      };
      void electron_convection(const FIELD_MATS &fields,
                                     DENS_MAT_VEC &flux)
      {
        FIELD_MATS::const_iterator edField = fields.find(ELECTRON_DENSITY);
        FIELD_MATS::const_iterator evField = fields.find(ELECTRON_VELOCITY);
        const DENS_MAT & n = edField->second;
        const DENS_MAT & v = evField->second;
        const CLON_VEC vx(v,CLONE_COL,0);
        const CLON_VEC vy(v,CLONE_COL,1);
        const CLON_VEC vz(v,CLONE_COL,2);
        zeroWorkspace_.reset(v.nRows(),1);

        if (maskX_) {
          flux[0] = zeroWorkspace_;
        }
        else {
          flux[0] = vx;
          flux[0] *= n;  // scale by n
        }

        if (maskY_) {
          flux[1] = zeroWorkspace_;
        }
        else {
          flux[1] = vy;
          flux[1] *= n;
        }

        if (maskZ_) {
          flux[2] = zeroWorkspace_;
        }
        else {
          flux[2] = vz;
          flux[2] *= n;
        }

      };
    protected:
      bool maskX_, maskY_, maskZ_;
      DENS_MAT zeroWorkspace_;
  };
  //-----------------------------------------------------------------------
  
  /**
   *  @class  ElectronFluxLinear
   *  @brief  Class for drift-diffusion electron flux with linear dependency on the electron density gradient
   */  
  
  class ElectronFluxLinear : public ElectronFlux
  {
    public:
      ElectronFluxLinear(std::fstream &matfile, std::map<std::string,double> & parameters);
      virtual ~ElectronFluxLinear() {};
      virtual void electron_flux(const FIELD_MATS &fields,
                                 const GRAD_FIELD_MATS &gradFields,
                                       DENS_MAT_VEC &flux)
      {
        FIELD_MATS::const_iterator edField = fields.find(ELECTRON_DENSITY);
        GRAD_FIELD_MATS::const_iterator dEdField = gradFields.find(ELECTRON_DENSITY);
        GRAD_FIELD_MATS::const_iterator dPhiField = gradFields.find(ELECTRIC_POTENTIAL);
         // J_n = - \mu n E  - D grad n
        // note electrons move counter to electric field grad E = - grad \phi
        const DENS_MAT & n = edField->second;
        const DENS_MAT_VEC & dn   = dEdField->second;
        const DENS_MAT_VEC & dphi = dPhiField->second;
        
        //n.print("DENSITY");
        //for (int i = 0; i < 3; i++) {
        //  dn[i].print("GRAD N");
        //  dphi[i].print("GRAD PHI");
        //}
        //cout << "------------------------------------------------====\n";

        flux[0] = n;
        flux[1] = n;
        flux[2] = n;

        if (maskX_)
          flux[0] = 0.;
        else {
          flux[0] *= electronMobility_*dphi[0]; // scale by n to get : n \mu grad(\phi)
          flux[0] += -electronDiffusivity_* dn[0];
        }

        if (maskY_)
          flux[1] = 0.;
        else {
          flux[1] *= electronMobility_* dphi[1] ;
          flux[1] += -electronDiffusivity_* dn[1];
        }

        if (maskZ_)
          flux[2] = 0.;
        else {
          flux[2] *= electronMobility_*dphi[2];
          flux[2] += -electronDiffusivity_* dn[2];
        }

      }; 
    protected:
      double electronMobility_, electronDiffusivity_;
  };
  //-----------------------------------------------------------------------
  
  /**
   *  @class  ElectronFluxThermopower
   *  @brief  Class for defining the electron flux (i.e., current) to include the elctron velocity or have a electron temperature-dependent mobility
   */  
  
  class ElectronFluxThermopower : public ElectronFlux
  {
    public:
      ElectronFluxThermopower(std::fstream &matfile,std::map<std::string,double> & parameters);
      virtual ~ElectronFluxThermopower() {};
      virtual void electron_flux(const FIELD_MATS &fields,
                                 const GRAD_FIELD_MATS &gradFields,
                                 DENS_MAT_VEC &flux)
      {
        if (fields.find(ELECTRON_VELOCITY)!=fields.end()) {
          // J_n = - e n v, note the electron charge e is unity
          electron_convection(fields,flux);
          flux[0] *= -1;
          flux[1] *= -1;
          flux[2] *= -1;
        }
        else {
          FIELD_MATS::const_iterator edField = fields.find(ELECTRON_DENSITY);
          FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
          GRAD_FIELD_MATS::const_iterator dEdField = gradFields.find(ELECTRON_VELOCITY);
          GRAD_FIELD_MATS::const_iterator dPhiField = gradFields.find(ELECTRIC_POTENTIAL);
          GRAD_FIELD_MATS::const_iterator dEtField = gradFields.find(ELECTRON_TEMPERATURE);

          // J_n = - \mu n grad \phi  - \mu kB/e T_e grad n 
          //       - \mu S n grad T_e - \mu kB/e n grad T_e
          const DENS_MAT & n   = edField->second;
          const DENS_MAT_VEC & dn   = dEdField->second;
          const DENS_MAT_VEC & dphi = dPhiField->second;
          const DENS_MAT_VEC & dT = dEtField->second;
          
          flux[0] = -electronMobility_*dphi[0];
          flux[1] = -electronMobility_*dphi[1];
          flux[2] = -electronMobility_*dphi[2];
          double coef = -electronMobility_*(seebeckCoef_ + kBeV_);
          flux[0] += coef* dT[0];
          flux[1] += coef* dT[1];
          flux[2] += coef* dT[2];
          flux[0] *= n; // scale by n 
          flux[1] *= n;
          flux[2] *= n;

          //GRAD_FIELD tmp = dn;
          const DENS_MAT & Te   = etField->second;
          //tmp[0] *= Te;
          //tmp[1] *= Te;
          //tmp[2] *= Te;
          coef = -electronMobility_*kBeV_; 
          //flux[0] += tmp[0];
          flux[0] += dn[0].mult_by_element(Te);
          flux[1] += dn[1].mult_by_element(Te);
          flux[2] += dn[2].mult_by_element(Te);
        }
      }; 
  protected:
      double electronMobility_, seebeckCoef_;
  };
  //-----------------------------------------------------------------------
  
  /**
   *  @class  ElectronFluxConvection
   *  @brief  Class for electron flux based on the standard convection term
   */  
  
  class ElectronFluxConvection : public ElectronFlux
  {
    public:
      ElectronFluxConvection(std::fstream &matfile,std::map<std::string,double> & parameters);
      virtual ~ElectronFluxConvection() {};
      virtual void electron_flux(const FIELD_MATS &fields,
                                 const GRAD_FIELD_MATS & /* gradFields */,
                                       DENS_MAT_VEC &flux)
      {
        // flux = n v
        electron_convection(fields,flux);
      };
  };
}

#endif 


