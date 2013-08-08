/** ATC_Transfer :  coarse-graining methods  */
#ifndef ATC_TRANSFER_H
#define ATC_TRANSFER_H

// ATC headers
#include "ATC_Method.h"
#include "MoleculeSet.h"
#include "AtomToMoleculeTransfer.h"

// Other headers
#include <map>
#include <list>
using std::list;

namespace ATC {

// Forward declarations
class FE_Engine;
class StressCauchyBorn;
class TimeFilter;

class ATC_Transfer : public ATC_Method {

 public:
  
  // constructor
  ATC_Transfer(std::string groupName, 
               double **& perAtomArray,
               LAMMPS_NS::Fix * thisFix,
               std::string matParamFile = "none");

  // destructor
  virtual ~ATC_Transfer();

  /** parser/modifier */
  virtual bool modify(int narg, char **arg);

  /** pre time integration */
  virtual void initialize();

  /** post time integration */
  virtual void finish();

  /** first time substep routines */
  virtual void pre_init_integrate();

  /** second time substep routine */
  virtual void pre_final_integrate();
  //virtual void final_integrate(){};
  virtual void post_final_integrate();

  /** communication routines */
  virtual void pre_neighbor() {ATC_Method::pre_neighbor(); neighborReset_ = true;};

  /** output function */
  virtual void output(); 

  /** external access to hardy data and other information*/
  const DENS_MAT * hardy_data(string field) { return &hardyData_[field].quantity(); }

 protected:

  /** pointer to position data : either x_reference or x_current */
  double ** xPointer_;

  /** data */
  TAG_FIELDS hardyData_; 
  SmallMoleculeSet * smallMoleculeSet_; // KKM add
  SmallMoleculeCentroid * moleculeCentroid_; // KKM add
  SmallMoleculeDipoleMoment * dipoleMoment_; // KKM add
  SmallMoleculeQuadrupoleMoment * quadrupoleMoment_; // KKM add
  /** container for dependency managed data */
  vector < DENS_MAN * > outputFields_;
  
  map < string, DENS_MAN * > outputFieldsTagged_;

  DENS_MAN * restrictedCharge_; // WIP/TEMP

  /** work space */ 
  DENS_MAT atomicScalar_;
  DENS_MAT atomicVector_;
  DENS_MAT atomicTensor_;

  /** calculation flags */
  Array<bool> fieldFlags_; 
  Array<bool> outputFlags_;
  Array<bool> gradFlags_;
  Array<bool> rateFlags_;
  map<string,int> computes_;
  bool outputStepZero_;

  /** check whether atoms have shifted box or element or neighbors changed */
  bool neighborReset_;

  //---------------------------------------------------------------
  /** initialization routines */
  //---------------------------------------------------------------
  /** gets baseline data from continuum model */
  virtual void set_continuum_data();
  /** sets up all data necessary to define the computational geometry */
  virtual void set_computational_geometry();
  /** constructs all data which is updated with time integration, i.e. fields */
  virtual void construct_time_integration_data();
  /** create methods, e.g. time integrators, filters */
  virtual void construct_methods();
  /** set up data which is dependency managed */
  virtual void construct_transfers();

  /** compute atom to nodal quantities */
// OBSOLETE
  void compute_energy(DENS_MAT & energy);
  void compute_internal_energy(DENS_MAT & energy);
  void compute_stress(DENS_MAT & stress);
  void compute_heatflux(DENS_MAT & flux);
  /** derived quantities: compute nodal to nodal quantities */
  void compute_eshelby_stress(DENS_MAT & eshebly_stress,
    const DENS_MAT &  energy, const DENS_MAT & stress, 
    const DENS_MAT & displacement_gradient);
  void cauchy_born_stress(const DENS_MAT &dudx, DENS_MAT &T, const DENS_MAT *temp=0);
  void cauchy_born_energy(const DENS_MAT &dudx, DENS_MAT &T, const DENS_MAT *temp=0);
  void cauchy_born_entropic_energy(const DENS_MAT &dudx, DENS_MAT &E, const DENS_MAT & T);
  void compute_transformed_stress(DENS_MAT & stress,
    const DENS_MAT & T, const DENS_MAT & displacement_gradient);
  void compute_polar_decomposition(DENS_MAT & rotation,
    DENS_MAT & stretch, const DENS_MAT & displacement_gradient);
  void compute_elastic_deformation_gradient(DENS_MAT & elastic_def_grad,
    const DENS_MAT & stress, const DENS_MAT & displacement_gradient);
  void compute_elastic_deformation_gradient2(DENS_MAT & elastic_def_grad,
    const DENS_MAT & stress, const DENS_MAT & displacement_gradient);
  /** hybrid computes? */
  void compute_electric_potential(DENS_MAT & electric_potential);
  void compute_vacancy_concentration(DENS_MAT & vacancy_concentration,
    const DENS_MAT & displacement_gradient,
    const DENS_MAT & number_density);
  /** calculate kinetic part of stress */
  virtual void compute_kinetic_stress(DENS_MAT& stress);
  /** calculate stress on-the-fly */
  virtual void compute_potential_stress(DENS_MAT& stress) = 0;
  /** calculate kinetic part of heat flux */
  virtual void compute_kinetic_heatflux(DENS_MAT& flux);
  /** calculate force part of the heat flux on-the-fly */
  virtual void compute_potential_heatflux(DENS_MAT& flux) = 0;
  /** compute molecule to nodal quantities */
  void compute_dipole_moment(DENS_MAT & dipole_moment);
  void compute_quadrupole_moment(DENS_MAT & quadrupole_moment);

  /** calculate dislocation density tensor from DXA output */
  virtual void compute_dislocation_density(DENS_MAT & dislocation_density) = 0;

  /** compute smooth fields */
  void compute_fields(void); 
  void time_filter_pre (double dt); 
  void time_filter_post(double dt); 

  /** mapping of atomic pairs to pair index value */
  class PairMap * pairMap_; 
  class BondMatrix * bondMatrix_; 
  class PairVirial * pairVirial_; 
  class PairPotentialHeatFlux * pairHeatFlux_; 

  /** routine to calculate matrix of force & position dyads */
  void compute_force_matrix();

  /** routine to calculate matrix of heat flux vector components */
  void compute_heat_matrix();

  /** routine to calculate matrix of kernel functions */
  virtual void compute_kernel_matrix_molecule() = 0; //KKM add

  /** calculate projection on the fly*/
  // REFACTOR use AtfKernelFunctionRestriction and derivatives
  virtual void compute_projection(const DENS_MAT & atomData,
                                  DENS_MAT & nodeData) = 0;

  /** routine to calculate matrix of bond functions */
  virtual void compute_bond_matrix(); 

  /** routine to set xPointer to xref or xatom */
  void set_xPointer();

  /** number of atom types */
  int nTypes_;

  /** project : given w_\alpha,
               w_I = \sum_\alpha N_{I\alpha} w_\alpha */
  // REFACTOR AtfShapeFunctionRestriction
  void project(const DENS_MAT & atomData,
                     DENS_MAT & nodeData);
  void project_molecule(const DENS_MAT & molData,
                              DENS_MAT & nodeData); //KKM add
  void project_molecule_gradient(const DENS_MAT & molData,
                              DENS_MAT & nodeData); //KKM add

  /** project (number density): given w_\alpha,
               w_I = \sum_\alpha N_{I\alpha} w_\alpha */
  // REFACTOR AtfNodeWeightedShapeFunctionRestriction
  void project_count_normalized(const DENS_MAT & atomData,
                                      DENS_MAT & nodeData);

  /** hardy_project (volume density): given w_\alpha, 
               w_I = 1/\Omega_I \sum_\alpha N_{I\alpha} w_\alpha 
               where   \Omega_I = \int_{support region of node I} N_{I} dV  */
  // REFACTOR AtfNodeWeightedShapeFunctionRestriction
  void project_volume_normalized(const DENS_MAT & atomData,
                                       DENS_MAT & nodeData);
  void project_volume_normalized_molecule(const DENS_MAT & molData,
                                                DENS_MAT & nodeData); // KKM add    
  void project_volume_normalized_molecule_gradient(const DENS_MAT & molData,
                                                DENS_MAT & nodeData); // KKM add    
  
  
  /** gradient_compute: given w_I, 
               w_J = \sum_I N_I'{xJ} \dyad w_I 
               where N_I'{xJ} is the gradient of the normalized 
                              shape function of node I evaluated at node J */
  // REFACTOR MatToGradBySparse
  void gradient_compute(const DENS_MAT & inNodeData,
                              DENS_MAT & outNodeData);

  int nNodesGlobal_;
  int nComputes_;

  /** workset data */
  VectorDependencyManager<SPAR_MAT * > * gradientMatrix_;

  
  SPAR_MAT atomicBondMatrix_;
  DENS_MAT atomicForceMatrix_;
  DENS_MAT atomicHeatMatrix_;

  /** use pair/bond forces */
  bool hasPairs_;
  bool hasBonds_;

  /** need to reset kernel function and bond matrix */
  bool resetKernelFunction_;

  /** use "exact" serial mode when using DXA to compute dislocation densities */
  bool dxaExactMode_;

  /** a continuum model to compare to and/or estimate  quantities */
  StressCauchyBorn * cauchyBornStress_;

  Array<TimeFilter *> timeFilters_;

  /** check consistency of fieldFlags_ */
  void check_field_dependencies(); 

};

};

#endif
