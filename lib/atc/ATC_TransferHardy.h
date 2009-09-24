/** ATC_TransferHardy : Hardy smoothing  */
#ifndef ATC_TRANSFER_HARDY_H
#define ATC_TRANSFER_HARDY_H

// ATC_Transfer headers
#include "ATC_Transfer.h"
#include "LammpsInterface.h"
#include "ATC_HardyKernel.h"
#include "TimeFilter.h"

// Other headers
#include <map>
#include <list>
using std::list;

namespace ATC {

  enum hardyNormalization { 
    NO_NORMALIZATION=0, 
    VOLUME_NORMALIZATION, NUMBER_NORMALIZATION, MASS_NORMALIZATION
  };

  enum hardyFieldName { 
    HARDY_DENSITY=0, 
    HARDY_DISPLACEMENT, 
    HARDY_MOMENTUM, 
    HARDY_VELOCITY, 
    HARDY_PROJECTED_VELOCITY, 
    HARDY_TEMPERATURE, 
    HARDY_KINETIC_TEMPERATURE, 
    HARDY_STRESS, 
    HARDY_HEAT_FLUX, 
    HARDY_ENERGY,
    HARDY_NUMBER_DENSITY,
    HARDY_ESHELBY_STRESS, 
    HARDY_CAUCHY_BORN_STRESS, 
    HARDY_TRANSFORMED_STRESS,
    HARDY_VACANCY_CONCENTRATION,
    HARDY_TYPE_CONCENTRATION,
    NUM_HARDY_FIELDS 
  };

  /** string to field enum */
  static bool string_to_hardy_field(const string & name, hardyFieldName & index) {
    if      (name=="density")
      index = HARDY_DENSITY;
    else if (name=="displacement")
      index = HARDY_DISPLACEMENT;
    else if (name=="momentum")
      index = HARDY_MOMENTUM;
    else if (name=="velocity")
      index = HARDY_VELOCITY;
    else if (name=="projected_velocity")
      index = HARDY_PROJECTED_VELOCITY;
    else if (name=="temperature")
      index = HARDY_TEMPERATURE;
    else if (name=="kinetic_temperature")
      index = HARDY_KINETIC_TEMPERATURE; // temperature from full KE
    else if (name=="stress")
      index = HARDY_STRESS;
    else if (name=="eshelby_stress")
      index = HARDY_ESHELBY_STRESS;
    else if (name=="cauchy_born_stress")
      index = HARDY_CAUCHY_BORN_STRESS;
    else if (name=="heat_flux")
      index = HARDY_HEAT_FLUX;
    else if (name=="energy")
      index = HARDY_ENERGY;
    else if (name=="number_density")
      index = HARDY_NUMBER_DENSITY;
    else if (name=="transformed_stress")
      index = HARDY_TRANSFORMED_STRESS;
    else if (name=="vacancy_concentration")
      index = HARDY_VACANCY_CONCENTRATION;
    else if (name=="type_concentration")
      index = HARDY_TYPE_CONCENTRATION;
    else
      return false;

    return true;
  };

  /** string to field enum */
  static bool hardy_field_to_string(const int & index,string & name)
  {
    if      (index==HARDY_DENSITY)
      name = "density";
    else if (index==HARDY_DISPLACEMENT)
      name = "displacement";
    else if (index==HARDY_MOMENTUM)
      name = "momentum";
    else if (index == HARDY_VELOCITY)
      name="velocity";
    else if (index == HARDY_PROJECTED_VELOCITY)
      name="projected_velocity";
    else if (index == HARDY_TEMPERATURE)
      name="temperature";
    else if (index == HARDY_KINETIC_TEMPERATURE) 
      name="kinetic_temperature";
    else if (index == HARDY_STRESS)
      name="stress";
    else if (index == HARDY_ESHELBY_STRESS)
      name="eshelby_stress";
    else if (index == HARDY_CAUCHY_BORN_STRESS)
      name="cauchy_born_stress";
    else if (index == HARDY_HEAT_FLUX)
      name="heat_flux";
    else if (index == HARDY_ENERGY)
      name="energy";
    else if (index == HARDY_NUMBER_DENSITY)
      name="number_density";
    else if (index == HARDY_TRANSFORMED_STRESS)
      name="transformed_stress";
    else if (index == HARDY_VACANCY_CONCENTRATION)
      name="vacancy_concentration";
    else if (index == HARDY_TYPE_CONCENTRATION)
      name="type_concentration";
    else
      return false;

    return true;
  };

// Forward declarations
class FE_Engine;
class StressCauchyBorn;

class ATC_TransferHardy : public ATC_Transfer {

 public:
  
  // constructor
  ATC_TransferHardy(std::string groupName, std::string matParamFile = "none");

  // destructor
  ~ATC_TransferHardy();

  /** parser/modifier */
  virtual bool modify(int narg, char **arg);
  
  /** pre time integration */
  virtual void initialize();

  /** post time integration */
  virtual void finish();

  /** first time substep routines */
  virtual void pre_init_integrate();
  virtual void init_integrate_velocity(){};
  virtual void mid_init_integrate(){};
  virtual void init_integrate_position(){};
  virtual void post_init_integrate();
  

  /** second time substep routine */
  virtual void pre_final_integrate();
  virtual void final_integrate(){};
  virtual void post_final_integrate(); 

  virtual void set_ghost_atoms() {};

 private:
  /** pointer to position data : either x_reference or x_current */
  double ** xPointer_;

  /** data */
  map <string, DENS_MAT > hardyData_;
  map <string, DENS_MAT > hardyDataOld_;
  map <string, DENS_MAT > filteredHardyData_;
  DENS_MAT atomicTemperature_;
  DENS_MAT atomicKineticTemperature_;
  DENS_MAT atomicDensity_;
  DENS_MAT atomicMomentum_;
  DENS_MAT atomicDisplacement_;
  DENS_MAT atomicVelocity_;
  DENS_MAT atomicStress_;
  DENS_MAT atomicHeat_;
  DENS_MAT atomicEnergy_;

  /** reference data */
  bool setRefPE_;
  DENS_MAT nodalRefPotentialEnergy_;
  double nodalRefPEvalue_;
  bool setRefPEvalue_;

  /** contour/boundary integral data */
  map < pair<string,string>, const set< PAIR > *  > bndyIntegralData_;
  map < pair<string,string>, const set< PAIR > *  > contourIntegralData_;

  /** workset data */
  SPAR_MAT kernelShpFcn_;
  SPAR_MAT atomicBondTable_;
  DENS_MAT vbar_;
  DENS_MAT ubar_;
  DENS_MAT atomicForceTable_;
  DENS_MAT atomicHeatTable_;
  DENS_MAT uVariationVelocity_;
  vector< SPAR_MAT > gradientTable_;

  /** calculation flags */
  Array<bool> fieldFlags_;
  Array<bool> outputFlags_;
  Array<bool> gradFlags_;
  Array<bool> rateFlags_;
  map<string,int> computes_;

  /** calculation flags */
  Array<int> fieldSizes_;

  /** compute nodal quantities */
  void compute_potential_energy(DENS_MAT & nodalPE);
  void compute_number_density(DENS_MAT & density);
  void compute_mass_density(DENS_MAT & density);
  void compute_displacement(DENS_MAT & displacement,
                            const DENS_MAT & density);
  void compute_momentum(DENS_MAT & momentum);
  void compute_projected_velocity(DENS_MAT & velocity);
  void compute_velocity(DENS_MAT & velocity,
    const DENS_MAT & density, const DENS_MAT & momentum);
  void compute_variation_velocity(DENS_MAT & velocity,
                                  const DENS_MAT & vI);
  void compute_temperature(DENS_MAT & temperature);
  void compute_kinetic_temperature(DENS_MAT & temperature);
  void compute_stress(DENS_MAT & stress);
  void compute_heatflux(DENS_MAT & flux);
  void compute_total_energy(DENS_MAT & energy);
  void compute_eshelby_stress(DENS_MAT & eshebly_stress,
    const DENS_MAT &  energy, const DENS_MAT & stress, 
    const DENS_MAT & displacement_gradient);
  void compute_cauchy_born_stress(DENS_MAT & cb_stress,
    const DENS_MAT & displacement_gradient);
  void compute_transformed_stress(DENS_MAT & stress,
    const DENS_MAT & T, const DENS_MAT & displacement_gradient);
  void compute_vacancy_concentration(DENS_MAT & vacancy_concentration,
    const DENS_MAT & displacement_gradient,
    const DENS_MAT & number_density);
  void compute_type_concentration(DENS_MAT & type_concentration);

  /** compute smooth fields */
  void compute_fields(void); 
  void time_filter_pre (double dt, bool output_now = true); 
  void time_filter_post(double dt, bool output_now = true); 

  /** compute boundary integral */
  void compute_boundary_integrals(void);

  /** output function */
  void output(); 

  /** physics specific filter initialization */
  void init_filter();

  /** kernel */
  ATC_HardyKernel* kernelFunction_;

  /** mapping of atomic pairs to pair index value */
  map< pair< int,int >,int > pairMap_;
  int nPairs_;

  /** routine to compute pair map */
  void compute_pair_map();

  /** routine to calculate matrix of Hardy bond functions */
  void compute_bond_matrix();

  /** routine to calculate matrix of kernel shape functions */
  void compute_kernel_matrix();

  /** routine to calculate matrix of gradient of Hardy functions */
  void compute_gradient_matrix();

  /** routine to check pair map */
  void check_pair_map();

  /** routine to calculate matrix of force & position dyads */
  void compute_force_matrix();

  /** routine to calculate matrix of heat flux vector components */
  void compute_heat_matrix();

  /** routine to calculate stress on-the-fly */
  void compute_potential_stress(DENS_MAT& stress);

  /** routine to calculate force part of the heat flux on-the-fly */
  void compute_potential_heatflux(DENS_MAT& flux);

  /** periodicity flags and lengths */
  int periodicity[3];
  double box_bounds[2][3];
  double box_length[3]; 

  /** routine to adjust node-pair distances due to periodicity */
  void periodicity_correction(DENS_VEC& xaI);

  /** routine to set xPointer to xref or xatom */
  void set_xPointer();

  /** number of atom types */
  int nTypes_;

 protected:

  /** project (number density): given w_\alpha,
               w_I = \sum_\alpha N_{I\alpha} w_\alpha */
  void project_count_normalized(const DENS_MAT & atomData,
                     DENS_MAT & nodeData);

  /** hardy_project (volume density): given w_\alpha, 
               w_I = 1/\Omega_I \sum_\alpha N_{I\alpha} w_\alpha 
               where   \Omega_I = \int_{support region of node I} N_{I} dV  */
  void project_volume_normalized(const DENS_MAT & atomData,
                           DENS_MAT & nodeData);

  /** gradient_compute: given w_I, 
               w_J = \sum_I N_I'{xJ} \dyad w_I 
               where N_I'{xJ} is the gradient of the normalized 
                              shape function of node I evaluated at node J */
  void gradient_compute(const DENS_MAT & inNodeData,
                              DENS_MAT & outNodeData);

  int nNodesGlobal_;

  /** compute kernel shape functions on-the-fly w/o storing N_Ia */
  bool kernelOnTheFly_;

  /** compute stress or heat flux on-the-fly w/o storing B_Iab */
  bool bondOnTheFly_;

  /** if false, no coarse velocity is calculated for kernel-based estimates */
  bool useAtomicShapeFunctions_;

  /** a continuum model to compare to and/or estimate Hardy quantities */
  StressCauchyBorn * cauchyBornStress_;

  const LammpsInterface * constLammpsInterface_;

  Array<TimeFilter *> timeFilters_;

  /** check consistency of fieldFlags_ */
  void check_fieldFlags_consistency(); 

};

};

#endif
