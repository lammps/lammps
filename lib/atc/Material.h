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
  class ATC_Transfer;
  class Stress;
  class ElectronHeatCapacity;
  class ElectronHeatFlux;
  class ElectronFlux;
  class ElectronPhononExchange;

  class Material
  {
  public:
    Material();
    virtual ~Material() {};
    /** this constructor parses material file */
    Material(string & tag, fstream & fileId);

    /** return label */
    string label(void) {return tag_;}

    /** check material has required interfaces */
    bool check_registry(const set<string> functionList) const
    { 
      set<string>::const_iterator itr;
      for (itr=functionList.begin(); itr!=functionList.end(); itr++) {
        if (registry_.find(*itr) == registry_.end()) return false;
      }
      return true;
    }

    /** access to material parameters */
    bool get_parameter(const string name, double & value) const
    { 
      map<string,double>::const_iterator iter = parameters_.find(name);
      if ( iter == parameters_.end()) {
        value = 0.0;
        return false;
      }
      value = iter->second;
      return true;
    } 
    /** true if rhs flux is linear (per field) */
    bool linear_flux(FieldName name) { 
     return linearFlux_(name);
    };

    /** true if rhs source is linear (per field) */
    bool linear_source(FieldName name) { 
     return linearSource_(name);
    };

    /** true if lhs density is constant (per field) */
    bool constant_density(FieldName name) { 
      return constantDensity_(name);
    };

    /** each of these is a field function computed at a set of points */
    /** if there is only one function it is in the base class 
     ** otherwise, a subsidary class is setup */
    /* -----------------------------------------------------------------*/
    /** densitities */
    /* -----------------------------------------------------------------*/
    /** thermal energy */
    void thermal_energy(const FIELDS & fields, 
                FIELD & energy); 
    /** heat capacity */
    void heat_capacity(const FIELDS & fields, 
                       FIELD & capacity);
    /** thermal energy */
    void electron_thermal_energy(const FIELDS & fields, 
                FIELD & energy); 
    /** electron heat capacity */
    void electron_heat_capacity(const FIELDS &fields,
                                FIELD &capacity);
    /** kinetic energy */
    void kinetic_energy(const FIELDS & fields, 
                FIELD & energy); 
    /** mass density */
    void mass_density(const FIELDS &fields,
                      FIELD &density);
    /** elastic energy */
    void elastic_energy(const FIELDS & fields, 
                const GRAD_FIELDS & gradFields,
                FIELD & energy); 
    /* -----------------------------------------------------------------*/
    /** fluxes */
    /* -----------------------------------------------------------------*/
    /** heat_flux */
    void heat_flux(const FIELDS & fields,
                   const GRAD_FIELDS & gradFields,
                   GRAD_FIELD  & heatFlux);
    /** electron conduction flux */
    void electron_heat_flux(const FIELDS &fields,
                            const GRAD_FIELDS &gradFields,
                            GRAD_FIELD &flux);
    /** stress */
    void stress(const FIELDS &fields,
                const GRAD_FIELDS &gradFields,
                GRAD_FIELD &stress);
    /** computes electron flux */
    void electron_flux(const FIELDS &fields,
                       const GRAD_FIELDS &gradFields,
                       GRAD_FIELD &flux);
    /** computes electric displacement */
    void electric_displacement(const FIELDS &fields,
                        const GRAD_FIELDS &gradFields,
                        GRAD_FIELD &flux);
    /** computes electric field */
    void electric_field(const FIELDS &fields,
                        const GRAD_FIELDS &gradFields,
                        GRAD_FIELD &flux);
    /* -----------------------------------------------------------------*/
    /** sources */
    /* -----------------------------------------------------------------*/
    /** electron-phonon exchange flux */
    void electron_phonon_exchange(const FIELDS &fields,
                                  FIELD &flux);
    /** computes net generation */
    virtual void electron_recombination(const FIELDS &fields,
                                  const GRAD_FIELDS &gradFields,
                                  FIELD &flux);
    /** computes drift diffusion charge density */
    virtual void charge_density(const FIELDS &fields,
                                const GRAD_FIELDS &gradFields,
                                FIELD &flux);


  protected:
    /** material's label  */
    string tag_;
    /** dictionary of material parameters */
    map<string,double> parameters_;
    /** dictionary of instantiated functions */
    set<string> registry_;
    /** per eqn flag of constant density */
    Array<bool> constantDensity_;
    /** per eqn flag of linearity/nonlinearity */
    Array<bool> linearFlux_, linearSource_;
    /** default heat capacity */
    double rhoCp_;
    /** default mass density */
    double rho_;
    /** heat capacity */
    double heatCapacity_;
    /** electron heat capacity */
    ElectronHeatCapacity * electronHeatCapacity_;
    /** mass density */
    double massDensity_;
    /** charge density */
    double chargeDensity_;
    /** thermal conductivity */
    double heatConductivity_;
    /** electron heat flux */
    ElectronHeatFlux * electronHeatFlux_;
    /** stress */
    Stress * stress_;
    /** electron-phonon exchange */
    ElectronPhononExchange * electronPhononExchange_;
    /** electron heat flux */
    ElectronFlux * electronFlux_;
    /** electric permittivity */
    double permittivity_;
    /** equilibrium carrier density */
    double electronEquilibriumDensity_;
    /** relaxation time */
    double electronRecombinationInvTau_;
    /** uniform donor concentration */
    double donorConcentration_; // NOTE only for uniform
  };

}
#endif // Material.h


