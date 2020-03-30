#ifndef MATERIAL_H
#define MATERIAL_H

#include <map>
#include <set>
#include <string>
#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"
#include "ATC_Error.h"
#include "LammpsInterface.h"

namespace ATC 
{
  class BodyForce;
  class Stress;
  class ViscousStress;
  class ElectronChargeDensity;
  class ElectronHeatCapacity;
  class ElectronHeatFlux;
  class ElectronFlux;
  class ElectronPhononExchange;
  class ElectronDragPower;

/**
 *  @class Material 
 *  @brief Base class for computing and storing properties and fields for a material
 */

  class Material
  {
  public:
    Material();
    virtual ~Material();
    /** this constructor parses material file */
    Material(std::string & tag, std::fstream & fileId);

    /** initialize */
    virtual void initialize();

    /** return label */
    std::string label(void) const {return tag_;}

    /** check material has required interfaces */
    bool check_registry(const std::set<std::string> functionList) const
    { 
      std::set<std::string>::const_iterator itr;
      for (itr=functionList.begin(); itr!=functionList.end(); itr++) {
        if (registry_.find(*itr) == registry_.end()) {
          std::stringstream ss;
          ss << "WARNING: material: [" << tag_ << "] cannot find " << *itr ;
          ATC::LammpsInterface::instance()->print_msg_once(ss.str());
        } 
        if (registry_.find(*itr) == registry_.end()) return false;
      }
      return true;
    }

    /** access to material parameters */
    bool parameter(const std::string name, double & value) const
    { 
      std::map<std::string,double>::const_iterator iter = parameters_.find(name);
      if ( iter == parameters_.end()) {
        value = 0.0;
        return false;
      }
      value = iter->second;
      return true;
    } 
    /** true if rhs flux is linear (per field) */
    bool linear_flux(FieldName name) const { 
     return linearFlux_(name);
    };

    /** true if rhs source is linear (per field) */
    bool linear_source(FieldName name) const { 
     return linearSource_(name);
    };

    /** true if lhs density is constant (per field) */
    bool constant_density(FieldName name) const { 
      return constantDensity_(name);
    };

    /** each of these is a field function computed at a set of points */
    /** if there is only one function it is in the base class 
     ** otherwise, a subsidiary class is setup */
    /* -----------------------------------------------------------------*/
    /** densities */
    /* -----------------------------------------------------------------*/
    /** thermal energy */
    void thermal_energy(const FIELD_MATS & fields,
                        DENS_MAT & energy) const; 
    /** heat capacity */
    void heat_capacity(const FIELD_MATS & fields, 
                       DENS_MAT & capacity) const;
    /** thermal energy */
    void electron_thermal_energy(const FIELD_MATS & fields, 
                                 DENS_MAT & energy) const; 
    /** electron capacities */
    void electron_mass_density(const FIELD_MATS &fields,
                               DENS_MAT &density) const;
    void electron_heat_capacity(const FIELD_MATS &fields,
                                DENS_MAT &capacity) const;
    /** derivative of electron heat capacity */
    void D_electron_heat_capacity(const FIELD_MATS &fields,
                                  const GRAD_FIELD_MATS &gradFields,
                                  DENS_MAT_VEC &Dcapacity) const;
    /** kinetic energy */
    void kinetic_energy(const FIELD_MATS & fields, 
                        DENS_MAT & energy) const; 
    /** mass density */
    void mass_density(const FIELD_MATS &fields,
                      DENS_MAT &density) const;
    /** elastic energy */
    void elastic_energy(const FIELD_MATS & fields, 
                        const GRAD_FIELD_MATS & gradFields,
                        DENS_MAT & energy) const; 
    /** permitivity */
    void permittivity(const FIELD_MATS & fields, 
                      DENS_MAT & energy) const; 
    /** inverse effective mass */
    void inv_effective_mass(const FIELD_MATS & fields, 
                            DENS_MAT & energy) const; 
    /** band-edge potential */
    void band_edge_potential(const FIELD_MATS & fields, 
                             DENS_MAT & energy) const; 
    /** viscosity */
    void viscosity(const FIELD_MATS & fields, 
                   DENS_MAT & energy) const; 
    /* -----------------------------------------------------------------*/
    /** fluxes */
    /* -----------------------------------------------------------------*/
    /** heat_flux */
    void heat_flux(const FIELD_MATS & fields,
                   const GRAD_FIELD_MATS & gradFields,
                   DENS_MAT_VEC  & heatFlux) const;
    /** electron conduction flux */
    void electron_heat_flux(const FIELD_MATS &fields,
                            const GRAD_FIELD_MATS &gradFields,
                            DENS_MAT_VEC &flux) const;
    /** electron heat convection */
    void electron_heat_convection(const FIELD_MATS &fields,
                                  DENS_MAT_VEC &flux) const;
    /** electron momentum convection */
    void electron_momentum_convection(const FIELD_MATS &fields,
                                      DENS_MAT_VEC &flux) const;
    /** stress */
    void stress(const FIELD_MATS &fields,
                const GRAD_FIELD_MATS &gradFields,
                DENS_MAT_VEC &stress) const;

    /** viscous stress */
    void viscous_stress(const FIELD_MATS &fields,
                        const GRAD_FIELD_MATS &gradFields,
                        DENS_MAT_VEC &viscousStress) const;

    /** computes electron flux */
    void electron_flux(const FIELD_MATS &fields,
                       const GRAD_FIELD_MATS &gradFields,
                       DENS_MAT_VEC &flux) const;
    void electron_thermal_stress(const FIELD_MATS &fields,
                                 const GRAD_FIELD_MATS &gradFields,
                                 DENS_MAT_VEC &stress) const;
    /** computes electric displacement */
    void electric_displacement(const FIELD_MATS &fields,
                               const GRAD_FIELD_MATS &gradFields,
                               DENS_MAT_VEC &flux) const;
    /** computes electric field */
    void electric_field(const FIELD_MATS &fields,
                        const GRAD_FIELD_MATS &gradFields,
                        DENS_MAT_VEC &flux) const;
    /* -----------------------------------------------------------------*/
    /** sources */
    /* -----------------------------------------------------------------*/
    /** electron-phonon exchange flux */
    bool electron_phonon_exchange(const FIELD_MATS &fields,
                                  DENS_MAT &flux) const;
    bool electron_drag_power(const FIELD_MATS &fields,
                             const GRAD_FIELD_MATS &gradFields,
                             DENS_MAT &power) const;
    void electron_drag_velocity_coefficient(const FIELD_MATS &fields,
                                            DENS_MAT &dragCoef) const;
    /** computes net generation */
    virtual bool electron_recombination(const FIELD_MATS &fields,
                                        const GRAD_FIELD_MATS &gradFields,
                                        DENS_MAT &recombination) const;
    /** computes drift diffusion charge density */ 
    virtual bool electron_charge_density(const FIELD_MATS &fields,
                                         DENS_MAT &density) const;
    virtual void D_electron_charge_density(const FieldName fieldName, 
                                           const FIELD_MATS &fields,
                                           DENS_MAT &D_density) const;
    /** computes momentum source */ 
    virtual bool body_force(const FIELD_MATS &fields,
                                  DENS_MAT &density) const;


  protected:
    /** material's label  */
    std::string tag_;
    /** dictionary of material parameters */
    std::map<std::string,double> parameters_;
    /** dictionary of instantiated functions */
    std::set<std::string> registry_;
    /** per eqn flag of constant density */
    Array<bool> constantDensity_;
    /** per eqn flag of linearity/nonlinearity */
    Array<bool> linearFlux_, linearSource_;
    /** default heat capacity */
    double rhoCp_;
    /** heat capacity */
    double heatCapacity_;
    /** electron heat capacity */
    ElectronHeatCapacity * electronHeatCapacity_;
    /** mass density */
    double massDensity_;
    /** thermal conductivity */
    double heatConductivity_;
    /** electron heat flux */
    ElectronHeatFlux * electronHeatFlux_;
    /** stress */
    Stress * stress_;
    /** viscous stress */
    ViscousStress * viscousStress_;
    /** body force */
    BodyForce * bodyForce_;
    /** electron-phonon exchange */
    ElectronPhononExchange * electronPhononExchange_;
    /** electron drag power */
    ElectronDragPower * electronDragPower_;
    /** electron flux */
    ElectronFlux * electronFlux_;
    /** electric permittivity */
    double permittivity_;
    /** inverse effective mass */
    double invEffectiveMass_;
    /** equilibrium carrier density */
    double electronEquilibriumDensity_;
    /** relaxation time */
    double electronRecombinationInvTau_;
    /** electron charge density model */
    ElectronChargeDensity * electronChargeDensity_;
  };

}
#endif // Material.h


