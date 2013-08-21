#ifndef ELECTRON_HEAT_CAPACITY_H
#define ELECTRON_HEAT_CAPACITY_H

#include <map>
#include <string>
#include "ATC_TypeDefs.h"
#include "Material.h"

namespace ATC {

  /**
   *  @class  ElectronHeatCapacity
   *  @brief  Base class for defining the heat capcity of the electron gas
   */

  class ElectronHeatCapacity
  {
    public:
      ElectronHeatCapacity()    {};
      virtual ~ElectronHeatCapacity() {};
      /** computes heat capacity */
      virtual void electron_heat_capacity(const FIELD_MATS &fields,
                                                DENS_MAT &capacity)=0;
      /** derivative of electron heat capacity */
      virtual void D_electron_heat_capacity(const FIELD_MATS &fields,
                                            const GRAD_FIELD_MATS &gradFields,
                                                  DENS_MAT_VEC & Dcapacity)=0;
      /** computes thermal energy */
      virtual void electron_thermal_energy(const FIELD_MATS &fields,
                                                 DENS_MAT &energy)=0; 
  };
  //-------------------------------------------------------------------
  
  /**
   *  @class  ElectronHeatCapacityConstant
   *  @brief  Class for a constant electron heat capacity
   */  
  
  class ElectronHeatCapacityConstant : public ElectronHeatCapacity
  {
    public:
      ElectronHeatCapacityConstant(std::fstream &matfile, 
                                   std::map<std::string,double> & parameters);
      virtual ~ElectronHeatCapacityConstant() {};
      virtual void electron_heat_capacity(const FIELD_MATS &fields,
                                                DENS_MAT &capacity)
      {
        FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
        const DENS_MAT & T = etField->second;
        capacity.resize(T.nRows(),T.nCols());
        capacity = electronHeatCapacity_;
      };
      virtual void D_electron_heat_capacity(const FIELD_MATS &fields,
                                            const GRAD_FIELD_MATS &gradFields,
                                                  DENS_MAT_VEC & Dcapacity)
      {
        FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
        zeroWorkspace_.reset((etField->second).nRows(),(etField->second).nCols());
        Dcapacity[0] = zeroWorkspace_;
        Dcapacity[1] = zeroWorkspace_;
        Dcapacity[2] = zeroWorkspace_;
      }                                     
      virtual void electron_thermal_energy(const FIELD_MATS &fields,
                                           DENS_MAT &energy)
      {
        FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
        const DENS_MAT & T = etField->second;
        energy = electronHeatCapacity_ * T;
      }; 
    protected:
      double electronHeatCapacity_;
      DENS_MAT zeroWorkspace_;
  };
  //-------------------------------------------------------------------
  
  /**
   *  @class  ElectronHeatCapacityLinear
   *  @brief  Class for an electron capacity that is directly proportional to the electron temperature
   */  
  
  class ElectronHeatCapacityLinear : public ElectronHeatCapacity
  {
    public:
      ElectronHeatCapacityLinear(std::fstream &matfile,
                                 std::map<std::string,double> & parameters);
      virtual ~ElectronHeatCapacityLinear() {};
      virtual void electron_heat_capacity(const FIELD_MATS &fields,
                                                DENS_MAT &capacity)
      {
        FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
        const DENS_MAT & T = etField->second;
        capacity = electronHeatCapacity_*T;
      };
      virtual void D_electron_heat_capacity(const FIELD_MATS &fields,
                                            const GRAD_FIELD_MATS &gradFields,
                                                  DENS_MAT_VEC &Dcapacity)
      {
        GRAD_FIELD_MATS::const_iterator dEtField = gradFields.find(ELECTRON_TEMPERATURE);
        const DENS_MAT_VEC & dT = dEtField->second;
        Dcapacity[0] = electronHeatCapacity_ * dT[0];
        Dcapacity[1] = electronHeatCapacity_ * dT[1];
        Dcapacity[2] = electronHeatCapacity_ * dT[2];
      }
      virtual void electron_thermal_energy(const FIELD_MATS &fields,
                                           DENS_MAT &energy)
      {
        FIELD_MATS::const_iterator etField = fields.find(ELECTRON_TEMPERATURE);
        const DENS_MAT & T = etField->second;
        energy = electronHeatCapacity_ * T;
        energy *= T;
      }; 
    protected:
      double electronHeatCapacity_;
  };
  //-------------------------------------------------------------------
  
  /**
   *  @class  ElectronHeatCapacityConstantAddDensity
   *  @brief  Class for a constant electron specific heat capacity (i.e, does not include the electron density)
   */  
  
  class ElectronHeatCapacityConstantAddDensity : public ElectronHeatCapacityConstant
  {
    public:
      ElectronHeatCapacityConstantAddDensity(std::fstream &matfile, 
                                             std::map<std::string,double> & parameters,
                                             Material * material);
      virtual ~ElectronHeatCapacityConstantAddDensity() {};
      virtual void electron_heat_capacity(const FIELD_MATS &fields,
                                                DENS_MAT &capacity)
      {
        ElectronHeatCapacityConstant::electron_heat_capacity(fields,capacity);
        FIELD_MATS::const_iterator edField = fields.find(ELECTRON_DENSITY);
        const DENS_MAT & density = edField->second;
        capacity = capacity.mult_by_element(density);
      };
      virtual void D_electron_heat_capacity(const FIELD_MATS &fields,
                                            const GRAD_FIELD_MATS &gradFields,
                                                  DENS_MAT_VEC &Dcapacity)
      {
        ElectronHeatCapacityConstant::D_electron_heat_capacity(fields,gradFields,Dcapacity);
        FIELD_MATS::const_iterator edField = fields.find(ELECTRON_DENSITY);
        const DENS_MAT & density = edField->second;
        Dcapacity[0] *= density;
        Dcapacity[1] *= density;
        Dcapacity[2] *= density;

        GRAD_FIELD_MATS::const_iterator dEdField = gradFields.find(ELECTRON_DENSITY);
        const DENS_MAT_VEC & Ddensity = dEdField->second;
        ElectronHeatCapacityConstant::electron_heat_capacity(fields,capacityMat_);
        Dcapacity[0] += Ddensity[0].mult_by_element(capacityMat_);
        Dcapacity[1] += Ddensity[1].mult_by_element(capacityMat_);
        Dcapacity[2] += Ddensity[2].mult_by_element(capacityMat_);
      }
      virtual void electron_thermal_energy(const FIELD_MATS &fields,
                                           DENS_MAT &energy)
      {
        ElectronHeatCapacityConstant::electron_thermal_energy(fields,energy);
        FIELD_MATS::const_iterator edField = fields.find(ELECTRON_DENSITY);
        const DENS_MAT & density = edField->second;
        energy *= density;
      };
    protected:
      Material * material_;
      DENS_MAT capacityMat_;  // avoid resizing if possible
  };
  //-------------------------------------------------------------------
  
  /**
   *  @class  ElectronHeatCapacityLinearAddDensity
   *  @brief  Class for a electron specific heat capacity that is proportional to the temperature (i.e., does not include density)
   */  
  
  class ElectronHeatCapacityLinearAddDensity : public ElectronHeatCapacityLinear
  {
    public:
      ElectronHeatCapacityLinearAddDensity(std::fstream &matfile,
                                           std::map<std::string,double> & parameters,
                                           Material * material);
      virtual ~ElectronHeatCapacityLinearAddDensity() {};
      virtual void electron_heat_capacity(const FIELD_MATS &fields,
                                                DENS_MAT &capacity)
      {
        ElectronHeatCapacityLinear::electron_heat_capacity(fields,capacity);
        FIELD_MATS::const_iterator edField = fields.find(ELECTRON_DENSITY);
        const DENS_MAT & density = edField->second;
        capacity *= density;
      };
       virtual void D_electron_heat_capacity(const FIELD_MATS &fields,
                                             const GRAD_FIELD_MATS &gradFields,
                                                   DENS_MAT_VEC &Dcapacity)
      {
        ElectronHeatCapacityLinear::D_electron_heat_capacity(fields,gradFields,Dcapacity);
        FIELD_MATS::const_iterator edField = fields.find(ELECTRON_DENSITY);
        const DENS_MAT & density = edField->second;
        Dcapacity[0] *= density;
        Dcapacity[1] *= density;
        Dcapacity[2] *= density;

        GRAD_FIELD_MATS::const_iterator dEdField = gradFields.find(ELECTRON_DENSITY);
        const DENS_MAT_VEC & Ddensity = dEdField->second;
        ElectronHeatCapacityLinear::electron_heat_capacity(fields,capacityWorkspace_);  
        Dcapacity[0] += Ddensity[0].mult_by_element(capacityWorkspace_);
        Dcapacity[1] += Ddensity[1].mult_by_element(capacityWorkspace_);
        Dcapacity[2] += Ddensity[2].mult_by_element(capacityWorkspace_);
      }
      virtual void electron_thermal_energy(const FIELD_MATS &fields,
                                           DENS_MAT &energy)
      {
        ElectronHeatCapacityLinear::electron_thermal_energy(fields,energy);
     
        FIELD_MATS::const_iterator edField = fields.find(ELECTRON_DENSITY);
        const DENS_MAT & density = edField->second;

        energy *= density;
      };
    protected:
      Material * material_;
      DENS_MAT capacityWorkspace_;  // avoid resizing if possible
  };
}

#endif 


