#ifndef ELECTRON_PHONON_EXCHANGE_H
#define ELECTRON_PHONON_EXCHANGE_H

#include <map>
#include <string>
#include <fstream>
#include "ATC_TypeDefs.h"

namespace ATC {
  
  class Material;

  /**
   *  @class  ElectronPhononExchange
   *  @brief  Base class for energy exchange between the electron and phonon temperatures
   */

  class ElectronPhononExchange
  {
    public:
      ElectronPhononExchange()  {};
      virtual ~ElectronPhononExchange() {};
      /** computes heat capacity */
      virtual bool electron_phonon_exchange(const FIELD_MATS & /* fields */,
                                            DENS_MAT & /* flux */) { return false; }
  };
  //-------------------------------------------------------------------

  /**
   *  @class  ElectronPhononExchangeLinear
   *  @brief  Class for electron-phonon energy exchange proportional to the difference between the two temperatures
   */

  class ElectronPhononExchangeLinear : public ElectronPhononExchange
  {
    public:
      ElectronPhononExchangeLinear(std::fstream &matfile,
                                   std::map<std::string,double> & parameters);
      virtual ~ElectronPhononExchangeLinear() {};
      virtual bool electron_phonon_exchange(const FIELD_MATS &fields,
                                                  DENS_MAT &flux)
      {
        FIELD_MATS::const_iterator tField = fields.find(TEMPERATURE);
        FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
        const DENS_MAT & T  = tField->second;
        const DENS_MAT & Te = etField->second;

        // flux = g * ( T- Te)
        flux = Te - T;
        flux *= exchangeCoef_;
        return true;
      };
    protected:
      double exchangeCoef_;
  };
  //-------------------------------------------------------------------
  
  /**
   *  @class  ElectronPhononExchangePowerLaw
   *  @brief  Class for electron-phonon exchange proportional to the temperature difference raised to a constant power
   */  
  
  class ElectronPhononExchangePowerLaw : public ElectronPhononExchange
  {
    public:
      ElectronPhononExchangePowerLaw(std::fstream &matfile,  
                                     std::map<std::string,double> & parameters);
      virtual ~ElectronPhononExchangePowerLaw() {};
      virtual bool electron_phonon_exchange(const FIELD_MATS &fields,
                                                  DENS_MAT &flux)
      {
        FIELD_MATS::const_iterator tField = fields.find(TEMPERATURE);
        FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
        const DENS_MAT & T  = tField->second;
        const DENS_MAT & Te = etField->second;

        // flux = g c_e T_e (T_e - T_p)^5 / T_e
        flux = (Te - T).pow(exponent_);
        flux *= exchangeCoef_;
        return true;
      }; 
    protected:
      double exchangeCoef_;
      int exponent_;
  };
  //-------------------------------------------------------------------
  
  /**
   *  @class  ElectronPhononExchangeHertel
   *  @brief  Class for electron-phonon exchange based on the formulation of Hertel for Cu
   */  
  
  class ElectronPhononExchangeHertel : public ElectronPhononExchange
  {
    public:
      ElectronPhononExchangeHertel(std::fstream &matfile,  
                                   std::map<std::string,double> & parameters,
                                   Material * material);
      virtual ~ElectronPhononExchangeHertel() {};
      virtual bool electron_phonon_exchange(const FIELD_MATS &fields,
                                                  DENS_MAT &flux);
    protected:
      double exchangeCoef_;
      double debeyeTemperature_;
      double massEnhancement_;
      Material * material_;
      
  private:
      ElectronPhononExchangeHertel() {};
      DENS_MAT capacityWorkspace_;
  };
}

#endif 


