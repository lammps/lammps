#ifndef ELECTRON_HEAT_CAPACITY_H
#define ELECTRON_HEAT_CAPACITY_H

#include <map>
#include <string>

using std::map;
using std::string;

#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"

namespace ATC {
  class ElectronHeatCapacity
  {
    public:
      ElectronHeatCapacity() 	{};
      virtual ~ElectronHeatCapacity() {};
      /** computes heat capacity */
      virtual void electron_heat_capacity(const FIELDS &fields,
                                     FIELD &capacity)=0; 
      /** computes thermal energy */
      virtual void electron_thermal_energy(const FIELDS &fields,
                                     FIELD &energy)=0; 
  };
  //-------------------------------------------------------------------
  class ElectronHeatCapacityConstant : public ElectronHeatCapacity
  {
    public:
      ElectronHeatCapacityConstant(fstream &matfile, 
                                   map<string,double> & parameters);
      virtual ~ElectronHeatCapacityConstant() {};
      virtual void electron_heat_capacity(const FIELDS &fields,
                                     FIELD &capacity)
      {
        const FIELD & T = (fields.find(ELECTRON_TEMPERATURE))->second;
        capacity.reset(T.nRows(),T.nCols());
        capacity = electronHeatCapacity_;
      }; 
      virtual void electron_thermal_energy(const FIELDS &fields,
                                     FIELD &energy)
      {
        const FIELD & T = (fields.find(ELECTRON_TEMPERATURE))->second;
        energy = electronHeatCapacity_ * T;
      }; 
    protected:
      double electronHeatCapacity_;
  };
  //-------------------------------------------------------------------
  class ElectronHeatCapacityLinear : public ElectronHeatCapacity
  {
    public:
      ElectronHeatCapacityLinear(fstream &matfile,
                                 map<string,double> & parameters);
      virtual ~ElectronHeatCapacityLinear() {};
      virtual void electron_heat_capacity(const FIELDS &fields,
                                     FIELD &capacity)
      {
        const FIELD & T = (fields.find(ELECTRON_TEMPERATURE))->second;
        capacity = electronHeatCapacity_*T;
      }; 
      virtual void electron_thermal_energy(const FIELDS &fields,
                                     FIELD &energy)
      {
        const FIELD & T = (fields.find(ELECTRON_TEMPERATURE))->second;
        energy = 0.5 * electronHeatCapacity_ * T;
        energy *= T;
      }; 
    protected:
      double electronHeatCapacity_;
  };
}

#endif 


