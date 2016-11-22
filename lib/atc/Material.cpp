#include "Material.h"
#include "ATC_Transfer.h"
#include "LammpsInterface.h"
#include "ElectronChargeDensity.h"
#include "ElectronHeatCapacity.h"
#include "ElectronHeatFlux.h"
#include "ElectronPhononExchange.h"
#include "ElectronDragPower.h"
#include "Stress.h"
#include "ViscousStress.h"
#include "BodyForce.h"
#include "ElectronFlux.h"
#include <sstream>
#include <fstream>
#include <vector>

using ATC_Utility::command_line;
using ATC_Utility::str2dbl;
using ATC_Utility::str2int;
using std::stringstream;
using std::set;
using std::fstream;
using std::string;
using std::vector;

namespace ATC {

  Material::Material()
    : rhoCp_(0),
      heatCapacity_(0),
      electronHeatCapacity_(NULL),
      massDensity_(0),
      heatConductivity_(0),
      electronHeatFlux_(NULL),
      stress_(NULL),
      viscousStress_(NULL),
      bodyForce_(NULL),
      electronPhononExchange_(NULL),
      electronDragPower_(NULL),
      electronFlux_(NULL),
      permittivity_(1.),
      invEffectiveMass_(1.),
      electronEquilibriumDensity_(0),
      electronRecombinationInvTau_(0),
      electronChargeDensity_(NULL)
  {
  }
  //--------------------------------------------------------------
  //  Constructor (parser)
  //--------------------------------------------------------------
  // Example:
  // material Cu
  //   heat_capacity  constant
  //     capacity  1.0
  //   end
  //   heat_flux linear
  //     conductivity    1.0
  //   end
  //   electron_heat_flux linear
  //     conductivity  1.0
  //   end
  //   electron_heat_capacity linear
  //     capacity 1.0
  //   end
  //   electron_phonon_exchange linear
  //     coefficient   0.1
  //   end
  // end
  Material::Material(string & tag, fstream &fileId)
    : tag_(tag),
      rhoCp_(0),
      heatCapacity_(0),
      electronHeatCapacity_(NULL),
      massDensity_(0),
      heatConductivity_(0),
      electronHeatFlux_(NULL),
      stress_(NULL),
      viscousStress_(NULL),
      bodyForce_(NULL),
      electronPhononExchange_(NULL),
      electronDragPower_(NULL),
      electronFlux_(NULL),
      permittivity_(1.),
      invEffectiveMass_(1.),
      electronEquilibriumDensity_(0),
      electronRecombinationInvTau_(0),
      electronChargeDensity_(NULL)
  {
    /*! \page man_material material
      \section syntax
        material <tag> <units> \n
           <commands> \n
        end \n
        tag - a unique identifier for the material type which can be referenced in input decks.  Multiple materials are specified using different tag regions, terminated with an 'end', in the material file.
        units - the LAMMPS units system the material is based on, used as a check against the actual LAMMPS units.  AtC units are consistent units using the LAMMPS length, mass, time, charge, and volts.  The only units conversion occuring within AtC are LAMMPS to AtC units and charge to volts units. 
        \section examples
        material Argon real
           -------
        end
        \section description
        Starts a section in which material properties can be specified.  Materials are organized by material, identified by a tag, and all associated material models are specified within its scope.  Unspecified material properties use defaults as indicated or are considered as null.  Null material properties contribute no value to integrals using them.  Material properties defined which are not part of the physics model are ignored.  Functions which are specified correspond to those implemented in the code and there is no mechanism for user-specified material models unless they are added to the main code.\n 
        \section restrictions
        Material models are only used for evaluating finite element integrals with for physics models they are associated with.
        \section related
        \section default
        Default for all material properties is null.  The null material using the tag 'null' is the only material defined by default. \n
      */
    linearFlux_.reset(NUM_FIELDS);
    linearFlux_ = false;
    linearSource_.reset(NUM_FIELDS);
    linearSource_ = true;
    constantDensity_.reset(NUM_FIELDS);
    constantDensity_ = false;

    
    rhoCp_ = ATC::LammpsInterface::instance()->heat_capacity();
    parameters_["heat_capacity"] = rhoCp_;
    heatCapacity_ = rhoCp_;
    registry_.insert("heat_capacity");
    registry_.insert("thermal_energy");
    constantDensity_(TEMPERATURE) = true;

    constantDensity_(DISPLACEMENT) = true;
    constantDensity_(VELOCITY) = true;
    electronDragPower_ = new ElectronDragPower(); 

    vector<string> line;
    while(fileId.good()) {
      command_line(fileId, line);
      if (line.size() == 0) continue;
      if (line.size() == 1) {
        if (line[0] == "end") {
          return;
        }
      }
      /*! \page man_mat_heat_capacity material heat_capcity 
        \section syntax
        heat_capacity constant\n
          capacity <value> \n
        end \n 
        \section description
        Overrides use of lattice heat capacity using Dulong-Petit law for continuum regions. \n 
        \section restrictions
        Only valid with AtC models incorporating a phonon temperature:  thermal, two-temperature, drift-diffusion
        \section related
        material
        \section default
        If no value is given, the Dulong-Petit value for the lattice is used. \n
      */
      if      (line[0] == "heat_capacity") { // over-ride default
        registry_. insert("heat_capacity");
        registry_. insert("thermal_energy");
        if (line[1] == "constant") {
          while(fileId.good()) {
            command_line(fileId, line);
            if (line.size() == 0) continue;
            if (line[0] == "end") break;
            double value = str2dbl(line[1]);
            if (line[0] == "capacity") {
              heatCapacity_ = value;
              parameters_["heat_capacity"] = heatCapacity_;
            }
          }
        }
      }
      /*! \page man_mat_heat_flux material heat_flux 
        \section syntax
        heat_flux linear\n
          conductivity <value> \n
        end \n 
        \section description
        Specifies a heat flux proportional to the temperature gradient. \n 
        \section restrictions
        Only valid with AtC models incorporating a phonon temperature:  thermal, two-temperature, drift-diffusion
        \section related
        material
        \section default
        Null. \n
      */
      else if (line[0] == "heat_flux") {
        registry_. insert("heat_flux");
        if (line[1] == "linear") {
          linearFlux_(TEMPERATURE) = true;
          while(fileId.good()) {
            command_line(fileId, line);
            if (line.size() == 0) continue;
            if (line[0] == "end") break;
            double value = str2dbl(line[1]);
            if (line[0] == "conductivity") {
              heatConductivity_ = value;
            }
          }
        }
      }
      /*! \page man_mat_electron_heat_flux material electron_heat_flux 
        \section syntax
        electron_heat_flux <null|linear|power_law|thermopower>\n
          <parameter> <value> \n
        end \n 
        null - no electron heat flux contributions \n
        linear - a heat flux proportional to the temperature gradient, parameter is 'conductivity'\n
        power_law - a heat flux proportional to the temperature gradient and ratio of electron to phonon temperatures, parameter is 'conductivity'\n
        thermopower - same as power_law but with an addition proportional to the electron current, parameters are 'conductivity' but it also uses the Seebeck coefficient defined elsewhere
        \section description
        Specifies the form for the electron heat flux. \n 
        \section restrictions
        Only valid with AtC models incorporating an electron temperature:  two-temperature, drift-diffusion
        \section related
        material
        \section default
        Null. \n
      */
      else if (line[0] == "electron_heat_flux") {
        registry_. insert("electron_heat_flux");
        if      (line[1] == "null") {
          linearFlux_(ELECTRON_TEMPERATURE) = true;
          if (electronHeatFlux_) delete electronHeatFlux_;
          electronHeatFlux_ = new ElectronHeatFlux();
        }
        else if      (line[1] == "linear") {
          linearFlux_(ELECTRON_TEMPERATURE) = true;
          if (electronHeatCapacity_) {
            if (electronHeatFlux_) delete electronHeatFlux_;
            electronHeatFlux_ = new ElectronHeatFluxLinear(fileId, parameters_,electronHeatCapacity_);
          }
          else {
            if (electronHeatFlux_) delete electronHeatFlux_;
            electronHeatFlux_ = new ElectronHeatFluxLinear(fileId, parameters_);
          }
        }
        else if (line[1] == "power_law") {
          if (electronHeatCapacity_) {
            if (electronHeatFlux_) delete electronHeatFlux_;
            electronHeatFlux_ = new ElectronHeatFluxPowerLaw(fileId, parameters_,electronHeatCapacity_);
          }
          else {
            if (electronHeatFlux_) delete electronHeatFlux_;
            electronHeatFlux_ = new ElectronHeatFluxPowerLaw(fileId, parameters_);
          }
        }
        else if (line[1] == "thermopower") {
          
          if (! electronFlux_) {
            throw ATC_Error( "for thermopower please define electron_flux before electron_heat_flux");
          }
          if (electronHeatFlux_) delete electronHeatFlux_;
          electronHeatFlux_ = new ElectronHeatFluxThermopower(fileId, 
            parameters_, electronFlux_);
        }
      }
      /*! \page man_mat_electron_heat_capacity material electron_heat_capacity 
        \section syntax
        electron_heat_capacity <constant|linear> <no_density>\n
          capacity <value> \n
        end \n 
        no_density - if this keyword is present, the electron density does not multiply the capacity\n
        constant - a constant electron heat flux \n
        linear - a heat flux proportional to the electron temperature\n
        \section description
        Specifies the form for the electron heat capacity. \n 
        \section restrictions
        Only valid with AtC models incorporating an electron temperature:  two-temperature, drift-diffusion
        \section related
        material
        \section default
        Null. \n
      */
      else if (line[0] == "electron_heat_capacity") {
        registry_. insert("electron_heat_capacity");
        registry_. insert("electron_thermal_energy");
        if (line[1] == "constant") {
          if ((line.size() == 3) && (line[2] == "no_density")) {
            if (electronHeatCapacity_) delete electronHeatCapacity_;
            electronHeatCapacity_ = new ElectronHeatCapacityConstantAddDensity(fileId,
                                                                               parameters_,
                                                                               this);
          }
          else {
            constantDensity_(ELECTRON_TEMPERATURE) = true;
            if (electronHeatCapacity_) delete electronHeatCapacity_;
            electronHeatCapacity_ = new ElectronHeatCapacityConstant(fileId,
                                                                     parameters_);
          }
        }
        else if (line[1] == "linear") {
          if ((line.size() == 3) && line[2] == "no_density") {
            if (electronHeatCapacity_) delete electronHeatCapacity_;
            electronHeatCapacity_ = new ElectronHeatCapacityLinearAddDensity(fileId,
                                                                             parameters_,
                                                                             this);
          }
          else {
            if (electronHeatCapacity_) delete electronHeatCapacity_;
            electronHeatCapacity_ = new ElectronHeatCapacityLinear(fileId,
                                                                   parameters_);
          }
        }
      }
      /*! \page man_mat_electron_phonon_exchange material electron_phonon_exchange 
        \section syntax
        electron_phonon_exchange <null|linear|power_law|hertel>\n
          <parameter> <value> \n
        end \n 
        null - no electron heat flux contributions \n
        linear - an energy exchange proportional to the temperature difference between the electron and phonon temperatures, parameter is 'coefficient'\n
        power_law - same as linear, but the temperature difference is raised to a specified power, parameters are 'coefficient' and 'exponent'\n
        hertel - exchange proportional to temperature difference to the 5th divided by the electron temperature, the coefficient is a function of the mass enhancement and Debeye temperature, parameters are 'debeye_temperature' and 'mass_enhancement'
        \section description
        Specifies the form for the electron/phonon heat exchange. \n 
        \section restrictions
        Only valid with AtC models incorporating an electron temperature:  two-temperature, drift-diffusion
        \section related
        material
        \section default
        Null. \n
      */
      else if (line[0] == "electron_phonon_exchange") {
        registry_. insert("electron_phonon_exchange");
        if      (line[1] == "null") {
          if (electronPhononExchange_) delete electronPhononExchange_;
          electronPhononExchange_ = new ElectronPhononExchange();
        }
        else if      (line[1] == "linear") {
          if (electronPhononExchange_) delete electronPhononExchange_;
          electronPhononExchange_ = new ElectronPhononExchangeLinear(fileId,  
                                                                     parameters_);
        }
        else if (line[1] == "power_law") {
          linearSource_(TEMPERATURE) = false;
          linearSource_(ELECTRON_TEMPERATURE) = false;
          if (electronPhononExchange_) delete electronPhononExchange_;
          electronPhononExchange_ = new ElectronPhononExchangePowerLaw(fileId, 
                                                                       parameters_);
        }
        else if (line[1] == "hertel") {
          linearSource_(TEMPERATURE) = false;
          linearSource_(ELECTRON_TEMPERATURE) = false;
          if (electronPhononExchange_) delete electronPhononExchange_;
          electronPhononExchange_ = new ElectronPhononExchangeHertel(fileId,parameters_,this);
        }
      }
      /*! \page man_mass_density material mass_density 
        \section syntax
        mass_density <no entry|basis|constant>\n
          <keyword> <values> \n
        end \n 
        no entry - compute mass density from the lattice using the mass of the first type, no keyword or values\n
        basis - compute mass density for the given number of atoms of each type in the lattice, no keyword, values are one integer per type specifying the number of atoms of that type in the lattice\n
        constant - prescribed mass density, keyword = density, value = desired mass density
        \section description
        Specifies the mass density of the system. \n 
        \section restrictions
        Valid for all AtC physics models.
        \section related
        material
        \section default
        Compute from the basis. \n
      */
      else if (line[0] == "mass_density") { // over-ride default
        registry_. insert("mass_density");
        registry_. insert("kinetic_energy");
        if (line.size() == 1 ) { // try to calculate from lattice
          massDensity_ = LammpsInterface::instance()->mass_density();
          parameters_["mass_density"] = massDensity_;
          stringstream ss;
          ss << "computed mass density : " << massDensity_ ;
          ATC::LammpsInterface::instance()->print_msg_once(ss.str());
        }
        else if (line[1] == "basis") {
          while(fileId.good()) {
            command_line(fileId, line);
            if (line.size() == 0) continue;
            if (line[0] == "end") break;
            int n = line.size();
            int* numPerType = new int[n];
            for (int i = 0; i < n; i++) {
              numPerType[i] = str2int(line[i]);
            }
            massDensity_ = LammpsInterface::instance()->mass_density(numPerType);
            delete [] numPerType;
            parameters_["mass_density"] = massDensity_;
            stringstream ss;
            ss << "computed mass density (from given basis) : " << massDensity_ ;
            ATC::LammpsInterface::instance()->print_msg_once(ss.str());
          }
        }
        else if (line[1] == "constant") {
          while(fileId.good()) {
            command_line(fileId, line);
            if (line.size() == 0) continue;
            if (line[0] == "end") break;
            double value = str2dbl(line[1]);
            if (line[0] == "density") {
              massDensity_ = value;
              parameters_["mass_density"] = massDensity_;
            }
          }
        }
      }
      /*! \page man_mat_stress material stress
        \section syntax
        stress <linear|cubic|damped_cubic|cauchy_born>\n
          <keyword> <values> \n
        end \n 
        null - no electron heat flux contributions \n
        linear - a stress tensor proportional to the displacements, keywords are 'modulus' and 'poissons_ratio'\n
        cubic - an anisotropic linear stress tensor, keywords are 'c11', 'c12', and 'c44'\n
        damped_cubic - same as cubic, with a damping term proportional to the velocity vector, keywords are 'c11', 'c12', 'c44', and the damping parameter 'gamma'\n
        cauchy_born - stress tensor is computed using the Cauchy-Born formalism from the lattice and given potential, keywords are 'pairstyle', 'linear' (linearizes the Cauchy-Born relationship), or 'temperature' (the temperature used to determine the Cauchy-Born stress).  The 'pairstyle' lines are followed by values of 'lj/cut', 'lj/smooth/linear', and 'eam', the latter two of which are followed on the line by the value for the cut-off radius.  The 'lj/cut' and 'lj/smooth/linear' pairstyles are followed on the next line using the keyword 'pair_coeff' followed by value of the pair-coefficients \sigma and \epsilon.
        \section description
        Specifies the form for the mechanical stress tensor. \n 
        \section restrictions
        Only valid with AtC models incorporating a mechanical stress:  elastic
        \section related
        material
        \section default
        Null. \n
      */
      else if (line[0] == "stress") {
        registry_. insert("stress");
        registry_. insert("elastic_energy");
        if      (line[1] == "linear") {
          linearFlux_(VELOCITY) = true;
          linearFlux_(DISPLACEMENT) = true;
          if (stress_) delete stress_;
          stress_ = new StressLinearElastic(fileId);
        }
        else if (line[1] == "cubic") {
          linearFlux_(VELOCITY) = true;
          linearFlux_(DISPLACEMENT) = true;
          if (stress_) delete stress_;
          stress_ = new StressCubicElastic(fileId);
        }
        else if (line[1] == "damped_cubic") {
          linearFlux_(VELOCITY) = true;
          linearFlux_(DISPLACEMENT) = true;
          if (stress_) delete stress_;
          stress_ = new StressCubicElasticDamped(fileId);
        }
        else if (line[1] == "cauchy-born") {
          CbData cb;
          LammpsInterface *lmp = LammpsInterface::instance();
          lmp->lattice(cb.cell_vectors, cb.basis_vectors);
          cb.inv_atom_volume = 1.0 / lmp->volume_per_atom();
          cb.e2mvv           = 1.0 / lmp->mvv2e();
          cb.boltzmann       = lmp->boltz();
          cb.atom_mass       = lmp->atom_mass(1); 
          if (stress_) delete stress_;
          stress_ = new StressCauchyBorn(fileId, cb);
        }
      }
      else if (line[0] == "viscous_stress") {
        registry_.insert("viscous_stress");
        if (line[1] == "constant") {
          linearFlux_(VELOCITY) = true;
          if (viscousStress_) delete viscousStress_;
          viscousStress_ = new ViscousStressConstant(fileId);
        }
      }
      /*! \page man_body_force material body_force 
        \section syntax
        body_force <electric_field|viscous>\n
          <keyword> <values> \n
        end \n 
        electric_field - adds body force proportional to the electric field and charge density, no keywords or values\n
        viscous - adds a body force proportional to the velocity vector, keyword = gamma (damping parameter) followed by its value\n
        \section description
        Specifies body forces acting on the system. \n 
        \section restrictions
        Valid for all AtC mechanical models:  elastic
        \section related
        material
        \section default
        Null. \n
      */
      else if (line[0] == "body_force") {
        registry_. insert("body_force");
        if (line.size() > 1) {
          if      (line[1] == "electric_field") {
            if (bodyForce_) delete bodyForce_;
            bodyForce_ = new BodyForceElectricField(fileId,parameters_);
          }
          else if (line[1] == "viscous") {
            if (bodyForce_) delete bodyForce_;
            bodyForce_ = new BodyForceViscous(fileId,parameters_);
          }
          else {
            if (bodyForce_) delete bodyForce_;
            bodyForce_ = new BodyForce();
          }
        }
        else {
          if (bodyForce_) delete bodyForce_;
          bodyForce_ = new BodyForce();
        }
      }
      else if (line[0] == "electron_flux") {
        registry_. insert("electron_flux");
        if      (line[1] == "null") {
          linearFlux_(ELECTRON_DENSITY) = true;
          if (electronFlux_) delete electronFlux_;
          electronFlux_ = new ElectronFlux();
        }
        else if      (line[1] == "linear") {
          linearFlux_(ELECTRON_DENSITY) = true;
          if (electronFlux_) delete electronFlux_;
          electronFlux_ = new ElectronFluxLinear(fileId, parameters_);
        }
        else if (line[1] == "thermopower") {
          if (electronFlux_) delete electronFlux_;
          electronFlux_ = new ElectronFluxThermopower(fileId, parameters_);
        }
        else if (line[1] == "convection") {
          if (electronFlux_) delete electronFlux_;
          electronFlux_ = new ElectronFluxConvection(fileId, parameters_);
        }
      }
      /*! \page man_electric_field material electric_field
        \section syntax
        electric_field linear\n
          permittivity <value> \n
        end \n 
        Provide a value for the permittivity or use LAMMPS' value if no value is given.\n
        \section description
        Specifies the electric displacement vector to be proportional to the electric field. \n 
        \section restrictions
        Valid for AtC physics models using electric fields:  fem_efield, drift-diffusion
        \section related
        material
        \section default
        Use LAMMPS' permittivity. \n
      */
      else if (line[0] == "electric_field") {
        registry_. insert("electric_field");
        registry_. insert("electric_displacement");
        if (line[1] == "linear") {
          linearFlux_(ELECTRIC_POTENTIAL) = true;
          while(fileId.good()) {
            command_line(fileId, line);
            if (line.size() == 0) continue;
            if (line[0] == "end") break;
            if (line[0] == "permittivity") {
              // if no value given use lammps dielectric constant
              if (line.size() == 1 ) {
                permittivity_ = LammpsInterface::instance()->dielectric();
              }
              else {
                double value = str2dbl(line[1]);
                permittivity_ = value;
              }
              // convert relative permitivity (dielectric) to abs internal units
              stringstream ss;
              ss << "permittivity: relative= " << permittivity_ ;
              permittivity_ *= LammpsInterface::instance()->epsilon0();
              ss << ", absolute= " << permittivity_ ;
              ss << ", lattice constant= " << LammpsInterface::instance()->max_lattice_constant()  ;
              ATC::LammpsInterface::instance()->print_msg_once(ss.str());
              LammpsInterface::UnitsType utype = LammpsInterface::instance()->units_style();
              if ( utype != LammpsInterface::REAL 
                && utype != LammpsInterface::METAL) {
                ATC::LammpsInterface::instance()->print_msg_once("WARNING: must use a unit system where: [Energy/force] = [Length] and [charge] = e");
              // note this is so that: perm0/e = 1 / (Epotential_units * Length)
              // our charge densities are multiples of the elemental charge
              }
              parameters_["permittivity"]    = permittivity_;
            }
          }
        }
      }
      else if (line[0] == "effective_mass") {
        registry_. insert("inv_effective_mass");
        if (line[1] == "constant") {
          linearFlux_(ELECTRON_WAVEFUNCTION) = true;
          while(fileId.good()) {
            command_line(fileId, line);
            if (line.size() == 0) continue;
            if (line[0] == "end") break;
            if (line[0] == "inv_effective_mass") {
              double value = str2dbl(line[1]);
              invEffectiveMass_ = value; // 1/m* = inv_eff_mass/m_e
              // convert to hbar^2 / 2 / m* / e
              //double hbar = LammpsInterface::instance()->hbar();
              //invEffectiveMass_ *= 0.5*hbar*hbar;
              // m_e in units [eV-ps^2/A^2] : 5.68562958414706e-32
              double scale = 3.80998192145007; // units [V A^2]
              invEffectiveMass_ *= scale;
              parameters_["inverse_effective_mass"] = invEffectiveMass_;
            }
          }
        }
      }
      else if (line[0] == "electron_drag") {
        registry_.insert("electron_drag_power");
        registry_.insert("electron_drag_coefficient");
        if      (line[1] == "null") {
          if (electronDragPower_) delete electronDragPower_;
          electronDragPower_ = new ElectronDragPower();
        }
        else if      (line[1] == "linear") {
          if (electronDragPower_) delete electronDragPower_;
          electronDragPower_ = new ElectronDragPowerLinear(fileId,  
                                                           parameters_,
                                                           this);
        }        
      }
      else if (line[0] == "electron_recombination") {
        registry_. insert("electron_recombination");
        if (line[1] == "linear") {
          while(fileId.good()) {
            command_line(fileId, line);
            if (line.size() == 0) continue;
            if (line[0] == "end") break;
            double value = str2dbl(line[1]);
            if (line[0] == "inv_relaxation_time") {
              electronRecombinationInvTau_ = value;
              parameters_["inv_relaxation_time"] = electronRecombinationInvTau_;
            }
            else if (line[0] == "equilibrium_carrier_density") {
              electronEquilibriumDensity_ = value;
              parameters_["equilibrium_carrier_density"] 
                = electronEquilibriumDensity_;
            }
          }
        }
      }
      /*! \page man_mat_electron_density material electron_density 
        \section syntax
        electron_density <null|linear|interpolation|exponential|fermi_dirac>\n
          <keyword> <values> \n
        end \n 
        null - no electron density constitutive model, uses the state variable \n
        linear - density is proportional to the electric field, keyword is 'coefficient' followed by its value\n
        interpolation - interpolates in a look-up table contained in a file provided after the 'interpolation' word, keywords are 'scale' and 'number_of_points', followed by their values \n
        exponential - density is based on Boltzmann statistics for the electric potential above an activation energy, keywords are 'intrinsic_concentration', 'intrinsic_energy', and reference_temperature', followed by their values\n
        fermi_dirac - density is based on Fermi-Dirac statistics for the electric potential relative to an activation energy, keywords are 'fermi_energy', 'reference_temperature', 'band_edge', 'donor_ionization_energy', and 'donor_concentration'
        \section description
        Specifies the form for the electron density. \n 
        \section restrictions
        Only valid with AtC models incorporating an electrons:  electrostatic, two-temperature, drift-diffusion
        \section related
        material
        \section default
        Null. \n
      */
      else if (line[0] == "electron_density") { // density is converted to charge
        registry_. insert("electron_charge_density");
        if (line[1] == "null") {
          if (electronChargeDensity_) delete electronChargeDensity_;
          electronChargeDensity_ = new ElectronChargeDensity();
        }
        else if (line[1] == "linear") {
          linearSource_(ELECTRIC_POTENTIAL) = false;
          if (electronChargeDensity_) delete electronChargeDensity_;
          electronChargeDensity_ = new ElectronChargeDensityLinear(fileId, parameters_);
        }
        else if (line[1] == "interpolation") {
          linearSource_(ELECTRIC_POTENTIAL) = false;
          if (electronChargeDensity_) delete electronChargeDensity_;
          electronChargeDensity_ = new ElectronChargeDensityInterpolation(fileId, parameters_);
        }
        else if (line[1] == "exponential") {
          linearSource_(ELECTRIC_POTENTIAL) = false;
          if (electronChargeDensity_) delete electronChargeDensity_;
          electronChargeDensity_ = new ElectronChargeDensityExponential(fileId, parameters_);
        }
        else if (line[1] == "fermi_dirac") {
          registry_. insert("band_edge_potential");  
          //linearSource_(ELECTRIC_POTENTIAL) = false;  // treated as constant
          if (electronChargeDensity_) delete electronChargeDensity_;
          electronChargeDensity_ = new ElectronChargeDensityFermiDirac(fileId, parameters_);
        }
        else {
          throw ATC_Error("unrecognized material function type: "+line[0]+" - "+line[1]);
        }
      }
      else {
        throw ATC_Error( "unrecognized material function: "+line[0]);
      }
    }
  }
  //--------------------------------------------------------------------
  Material::~Material()
  {
    if (electronDragPower_) delete electronDragPower_;
    if (electronChargeDensity_) delete electronChargeDensity_;
    if (electronHeatCapacity_) delete electronHeatCapacity_;
    if (electronHeatFlux_) delete electronHeatFlux_;
    if (electronFlux_) delete electronFlux_;
    if (stress_) delete stress_;
    if (viscousStress_) delete viscousStress_;
    if (bodyForce_) delete bodyForce_;
    if (electronPhononExchange_) delete electronPhononExchange_;
  }
//---------------------------------------------------------------------
void Material::initialize(){if (stress_) stress_->initialize();}
void Material::heat_capacity(
  const FIELD_MATS & fields,
  DENS_MAT & capacity) const
{
  const DENS_MAT & T = (fields.find(TEMPERATURE))->second;
  int nNodes  = T.nRows();
  capacity.reset(nNodes,1);
  capacity = heatCapacity_;
};
//---------------------------------------------------------------------
void Material::thermal_energy(
  const FIELD_MATS &fields,
  DENS_MAT &energy) const
{
  
  const DENS_MAT & T = (fields.find(TEMPERATURE))->second;
  energy = heatCapacity_ * T;
};
//---------------------------------------------------------------------
void Material::electron_heat_capacity(
  const FIELD_MATS & fields,
  DENS_MAT & capacity) const
{
  electronHeatCapacity_->electron_heat_capacity(fields,capacity);
};
//---------------------------------------------------------------------
void Material::D_electron_heat_capacity(
  const FIELD_MATS & fields,
  const GRAD_FIELD_MATS & gradFields,
  DENS_MAT_VEC & Dcapacity) const
{
  electronHeatCapacity_->D_electron_heat_capacity(fields,gradFields,Dcapacity);
};
//---------------------------------------------------------------------
void Material::electron_thermal_energy(
  const FIELD_MATS &fields,
  DENS_MAT &energy) const
{
  electronHeatCapacity_->electron_thermal_energy(fields,energy);
};
//---------------------------------------------------------------------
void Material::mass_density(
  const FIELD_MATS &fields,
  DENS_MAT &density) const
{
  
  int nNodes = 0;
  FIELD_MATS::const_iterator field = fields.find(MASS_DENSITY);
  if (field != fields.end()) {
    const DENS_MAT & d = field->second;
    nNodes  = d.nRows();
  }
  else {
    FIELD_MATS::const_iterator field = fields.find(VELOCITY);
    if (field != fields.end()) {
      const DENS_MAT & v = field->second;
      nNodes  = v.nRows();
    }
  }
  density.reset(nNodes,1);
  density = massDensity_;
};
//---------------------------------------------------------------------
void Material::electron_mass_density(
  const FIELD_MATS &fields,
  DENS_MAT &density) const
{
  
  
  int nNodes = 0;
  FIELD_MATS::const_iterator field = fields.find(ELECTRON_DENSITY);
  //if (field != fields.end()) {
  const DENS_MAT & d = field->second;
    nNodes  = d.nRows();
    //}
  density.reset(nNodes,1);
  inv_effective_mass(fields,density);
  density = d.div_by_element(density);
};
//---------------------------------------------------------------------
void Material::kinetic_energy(
  const FIELD_MATS &fields,
  DENS_MAT &energy) const
{
  FIELD_MATS::const_iterator field = fields.find(VELOCITY);
  if (field != fields.end()) {
    const DENS_MAT & v = field->second;
    energy = 0.5*massDensity_*v;
    energy *= v;
  }
  else {
    energy = 0.;
  }
};
//---------------------------------------------------------------------
void Material::permittivity(
  const FIELD_MATS &fields,
  DENS_MAT &density) const
{
  const DENS_MAT & phi = (fields.find(ELECTRIC_POTENTIAL))->second;
  int nNodes  = phi.nRows();
  density.reset(nNodes,1);
  density = permittivity_;
};
//---------------------------------------------------------------------
void Material::band_edge_potential(
  const FIELD_MATS &fields,
  DENS_MAT &density) const
{
  electronChargeDensity_->band_edge_potential(fields,density);
};
//---------------------------------------------------------------------
void Material::inv_effective_mass(
  const FIELD_MATS &fields,
  DENS_MAT &density) const
{
  const DENS_MAT & phi = (fields.find(ELECTRON_DENSITY))->second;
  int nNodes  = phi.nRows();
  density.reset(nNodes,1);
  density = invEffectiveMass_;
};
//---------------------------------------------------------------------
void Material::heat_flux(
  const FIELD_MATS & fields,
  const GRAD_FIELD_MATS & gradFields,
  DENS_MAT_VEC & flux) const
{
  const DENS_MAT_VEC & dT = (gradFields.find(TEMPERATURE))->second;
  flux[0] = -heatConductivity_* dT[0];
  flux[1] = -heatConductivity_* dT[1];
  flux[2] = -heatConductivity_* dT[2];
}
//---------------------------------------------------------------------
void Material::electron_heat_flux(
  const FIELD_MATS & fields,
  const GRAD_FIELD_MATS & gradFields,
  DENS_MAT_VEC & flux) const
{
  electronHeatFlux_->electron_heat_flux(fields,gradFields,flux);
}
//---------------------------------------------------------------------
void Material::electron_heat_convection(
  const FIELD_MATS & fields,
  DENS_MAT_VEC & flux) const
{
  electronHeatFlux_->electron_heat_convection(fields,flux);
}
//---------------------------------------------------------------------
bool Material::electron_phonon_exchange(
  const FIELD_MATS & fields,
  DENS_MAT & flux) const
{
  return electronPhononExchange_->electron_phonon_exchange(fields,flux);
}
//---------------------------------------------------------------------
void Material::electron_drag_velocity_coefficient(
  const FIELD_MATS &fields,
  DENS_MAT & dragCoef) const
{
  electronDragPower_->electron_drag_velocity_coefficient(fields,dragCoef);
}
//---------------------------------------------------------------------
bool  Material::electron_drag_power(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &gradFields,
  DENS_MAT & power) const
{
  return electronDragPower_->electron_drag_power(fields,gradFields,power);
}
//---------------------------------------------------------------------
bool Material::electron_recombination(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &gradFields,
  DENS_MAT & recombination) const
{
  // 1/tau (n - n0)
  const DENS_MAT & n   = (fields.find(ELECTRON_DENSITY))->second;
  recombination  = n;
  recombination -= electronEquilibriumDensity_;
  recombination *= -electronRecombinationInvTau_;
  return true; 
}
//---------------------------------------------------------------------
bool Material::electron_charge_density(
  const FIELD_MATS &fields,
  DENS_MAT & density) const
{
  return electronChargeDensity_->electron_charge_density(fields,density);
};
//---------------------------------------------------------------------
void Material::D_electron_charge_density(const FieldName fieldName,
                                         const FIELD_MATS &fields,
                                         DENS_MAT & D_density) const
{
  electronChargeDensity_->D_electron_charge_density(fieldName,fields,D_density);
};
//---------------------------------------------------------------------
bool Material::body_force(
  const FIELD_MATS &fields,
  DENS_MAT & density) const
{
  return bodyForce_->body_force(fields,density);
};
//---------------------------------------------------------------------
void Material::stress(
  const FIELD_MATS & fields,
  const GRAD_FIELD_MATS & gradFields,
  DENS_MAT_VEC & stress) const
{
  stress_->stress(fields,gradFields,stress);
}
//---------------------------------------------------------------------
void Material::elastic_energy(
  const FIELD_MATS & fields,
  const GRAD_FIELD_MATS & gradFields,
  DENS_MAT & energy) const
{
  stress_->elastic_energy(fields,gradFields,energy);
}
//---------------------------------------------------------------------
void Material::viscous_stress(
  const FIELD_MATS & fields,
  const GRAD_FIELD_MATS & gradFields,
  DENS_MAT_VEC & stress) const
{
  viscousStress_->viscous_stress(fields,gradFields,stress);
}
//---------------------------------------------------------------------
void Material::viscosity(
  const FIELD_MATS &fields,
  DENS_MAT &coefs) const
{
  viscousStress_->viscosity(fields,coefs);
}
//---------------------------------------------------------------------
void Material::electron_flux(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &gradFields,
  DENS_MAT_VEC &flux) const
{
  electronFlux_->electron_flux(fields,gradFields,flux);
}
//---------------------------------------------------------------------
void Material::electric_field(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &gradFields,
  DENS_MAT_VEC &flux) const
{
  // E = - grad \phi
  const DENS_MAT_VEC & dphi = (gradFields.find(ELECTRIC_POTENTIAL))->second;

  flux[0] = -1.0* dphi[0];
  flux[1] = -1.0* dphi[1];
  flux[2] = -1.0* dphi[2];
}
//---------------------------------------------------------------------
void Material::electric_displacement(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &gradFields,
  DENS_MAT_VEC &flux) const
{
  // D = - permitivity grad \phi
  const DENS_MAT_VEC & dphi = (gradFields.find(ELECTRIC_POTENTIAL))->second;
  flux[0] = -permittivity_* dphi[0];
  flux[1] = -permittivity_* dphi[1];
  flux[2] = -permittivity_* dphi[2];
}
//---------------------------------------------------------------------
} // end namespace
