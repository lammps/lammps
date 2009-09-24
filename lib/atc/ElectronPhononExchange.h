#ifndef ELECTRON_PHONON_EXCHANGE_H
#define ELECTRON_PHONON_EXCHANGE_H

#include <map>
#include <string>

using std::map;
using std::string;

#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"

namespace ATC {
  class ElectronPhononExchange
  {
    public:
      ElectronPhononExchange() 	{};
      virtual ~ElectronPhononExchange() {};
      /** computes heat capacity */
      virtual void electron_phonon_exchange(const FIELDS &fields,
                                     DENS_MAT &flux)=0; 
  };
  //-------------------------------------------------------------------
  class ElectronPhononExchangeLinear : public ElectronPhononExchange
  {
    public:
      ElectronPhononExchangeLinear(fstream &matfile,
                                   map<string,double> & parameters);
      virtual ~ElectronPhononExchangeLinear() {};
      virtual void electron_phonon_exchange(const FIELDS &fields,
                                     DENS_MAT &flux)
      {
         // flux = g * ( T- Te)
         const FIELD & T  = (fields.find(TEMPERATURE))->second;
         const FIELD & Te = (fields.find(ELECTRON_TEMPERATURE))->second;
         flux = Te;
         flux -= T;
         flux *= exchangeCoef_;
      }; 
    protected:
      double exchangeCoef_;
  };
  //-------------------------------------------------------------------
  class ElectronPhononExchangePowerLaw : public ElectronPhononExchange
  {
    public:
      ElectronPhononExchangePowerLaw(fstream &matfile,  
                                     map<string,double> & parameters);
      virtual ~ElectronPhononExchangePowerLaw() {};
      virtual void electron_phonon_exchange(const FIELDS &fields,
                                     DENS_MAT &flux)
      {
         // flux = g c_e T_e (T_e - T_p)^5 / T_e
         const FIELD & T  = (fields.find(TEMPERATURE))->second;
         const FIELD & Te = (fields.find(ELECTRON_TEMPERATURE))->second;
         flux = Te;
         flux -= T;
         flux = flux.pow(exponent_);
         flux *= exchangeCoef_;
      }; 
    protected:
      double exchangeCoef_;
      int exponent_;
  };
}

#endif 


