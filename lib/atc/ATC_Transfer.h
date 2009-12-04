// ATC_Transfer : a base class for atom-continuum transfers & control 
// (derived classes are physics dependent)
// note init() must be called before most functions will work
// we expect only one derived class to exist, i.e. don't call initilize/modify/etc more than once

// NOTE:  implement print

#ifndef ATC_TRANSFER_H
#define ATC_TRANSFER_H

// ATC_Transfer headers
#include "PhysicsModel.h"
#include "MatrixLibrary.h"
#include "Array.h"
#include "Array2D.h"
#include "OutputManager.h"
#include "XT_Function.h"
#include "FE_Element.h"
#include "TimeFilter.h"
#include "LammpsInterface.h"
#include "FE_Engine.h"
#include "ExtrinsicModel.h"

// Other headers
#include <vector>
#include <set>

using namespace std;

namespace ATC {

  // Forward declarations
  class PrescribedDataManager;
  class TimeIntegrator;

  /**
   *  @class  ATC_Transfer
   *  @brief  Base class for atom-continuum transfer operators
   */

  class ATC_Transfer {

  public:
    // NOTE should remove friends and use existing ATC hooks
    friend class ExtrinsicModel; // friend is not inherited
    friend class ExtrinsicModelTwoTemperature;
    friend class ExtrinsicModelDriftDiffusion; 
    friend class ExtrinsicModelElectrostatic;
    friend class ExtrinsicModelElectrostaticElastic;

    //---------------------------------------------------------------
    /** \name enumerated types */
    //---------------------------------------------------------------
    /*@{*/
    /** boundary integration */
    enum BoundaryIntegrationType {
      NO_QUADRATURE=0,
      FE_QUADRATURE,
      FE_INTERPOLATION
    };
    /** boundary integration */
    enum IntegrationDomainType {
      FULL_DOMAIN=0,
      ATOM_DOMAIN,
      FE_DOMAIN
    };
    /** atomic weight specification */
    enum AtomicWeightType {
      USER=0,
      LATTICE,
      ELEMENT,
      REGION,
      GROUP,
      MULTISCALE
    };
    /** shape function location with respect to MD domain */
    enum ShapeFunctionType {
      FE_ONLY = 0,
      MD_ONLY,
      BOUNDARY
    };
    /** */
    enum AtomToElementMapType {
      LAGRANGIAN=0, 
      EULERIAN
    };
    /* enumerated type for coupling matrix structure */
    enum MatrixStructure {
      FULL=0,     // contributions from all nodes
      LOCALIZED,  // contributions only from nodes with sources
      LUMPED      // row-sum lumped version of full matrix
    };
    /** */ 
    enum GroupThermostatType {
      RESCALE_TEMPERATURE,
      RESCALE_TEMPERATURE_RATE,
      RESCALE_RATE
    };
    /** */
    enum GroupThermostatMomentumControl {
      NO_MOMENTUM_CONTROL,
      ZERO_MOMENTUM,
      PRESERVE_MOMENTUM
    };
    /** */
    enum ATC_GroupComputes { LAMBDA_POWER };
    /*@}*/

    /** constructor */
    ATC_Transfer(void);

    /** destructor */
    virtual ~ATC_Transfer();

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** pre integration run */
    virtual void initialize();

    /** post integration run : called at end of run or simulation */
    virtual void finish();

    /** Predictor phase, executed before Verlet */
    virtual void pre_init_integrate() = 0;

    /** Predictor phase, Verlet first step for velocity */
    virtual void init_integrate_velocity(); 

    /** Predictor phase, executed between velocity and position Verlet */
    virtual void mid_init_integrate() = 0;

    /** Predictor phase, Verlet first step for position */
    virtual void init_integrate_position(); 

    /** Predictor phase, executed after Verlet */
    virtual void post_init_integrate() = 0;

    /** Corrector phase, executed before Verlet */
    virtual void pre_final_integrate();

    /** Corrector phase, Verlet second step for velocity */
    virtual void final_integrate(); 

    /** Corrector phase, executed after Verlet*/
    virtual void post_final_integrate()=0;

    //---------------------------------------------------------------
    /** \name memory management and processor information exchange */
    //---------------------------------------------------------------
    /*@{*/
    /** pre_exchange  is our indicator that atoms have moved across processors */
    void pre_exchange() {
      atomSwitch_ = true;
    }
    int memory_usage();
    void grow_arrays(int);
    void copy_arrays(int, int);
    int pack_exchange(int, double *);
    int unpack_exchange(int, double *);
    int pack_comm(int , int *, double *, int, int *);
    void unpack_comm(int, int, double *);
    /*@}*/

    /** Get grouping bit for LAMMPS compatability */
    int get_groupbit() { return groupbit_; }

    //---------------------------------------------------------------
    /** \name Wrappers for calls to lammps region.h information */
    //---------------------------------------------------------------
    /*@{*/
    bool get_region_bounds(const char * regionName,
                           double &xmin, double &xmax,
                           double &ymin, double & ymax,
                           double &zmin, double &zmax,
                           double &xscale,
                           double &yscale,
                           double &zscale);

    /** return pointer to PrescribedDataManager */
    PrescribedDataManager * get_prescribed_data_manager() {
      return prescribedDataMgr_;
    }
    
    /** return referece to ExtrinsicModelManager */
    ExtrinsicModelManager & get_extrinsic_model_manager() {
      return extrinsicModelManager_;
    }

    /** compute scalar for output */
    virtual double compute_scalar() {return 0.;}
    
    /** compute vector for output */
    virtual double compute_vector(int n) {return 0.;}
    /*@}*/

    //---------------------------------------------------------------
    /** \name access methods for output computation data */
    //---------------------------------------------------------------
    /*@{*/
    int scalar_flag() const {return scalarFlag_;}
    int vector_flag() const {return vectorFlag_;}
    int size_vector() const {return sizeVector_;}
    int global_freq() const {return globalFreq_;};
    int extscalar() const {return extScalar_;};
    int extvector() const {return extVector_;};
    int * extlist() {return extList_;};
    /*@}*/

    //---------------------------------------------------------------
    /** \name Access methods for data used by various methods */
    //---------------------------------------------------------------
    /*@{*/
    /** access to name FE fields */
    DENS_MAT &get_field(FieldName thisField){return fields_[thisField];};
    /** access to FE field time derivatives */
    DENS_MAT &get_dot_field(FieldName thisField){return dot_fields_[thisField];};
    /** access to name FE source terms */
    DENS_MAT &get_source(FieldName thisField){return sources_[thisField];};
    /** access to name atomic source terms */
    DENS_MAT &get_atomic_source(FieldName thisField){return atomicSources_[thisField];};
    /** access to name extrinsic source terms */
    DENS_MAT &get_extrinsic_source(FieldName thisField){return extrinsicSources_[thisField];};
    /** access to boundary fluxes */
    DENS_MAT &get_boundary_flux(FieldName thisField){return boundaryFlux_[thisField];};
    /** access to nodal fields of atomic variables */
    DENS_MAT &get_atomic_field(FieldName thisField)
      { return fieldNdFiltered_[thisField]; };
    /** access to auxilliary storage */
    DENS_MAT &get_aux_storage(FieldName thisField)
      { return auxStorage_[thisField]; };
    /** access to all fields */
    FIELDS &get_fields() {return fields_;};
    /** access to all fields rates of change (roc) */
    FIELDS &get_fields_roc() {return dot_fields_;};
    /** access to all boundary fluxes */
    FIELDS &get_boundary_fluxes() {return boundaryFlux_;};
    /** add a new field */
    void add_fields(map<FieldName,int> & newFieldSizes);
    /** access to finite element right-hand side data */
    DENS_MAT &get_field_rhs(FieldName thisField)
      { return rhs_[thisField]; };
    /** access to inverse mass matrices */
    DIAG_MAT &get_mass_mat_inv(FieldName thisField)
      { return massMatInv_[thisField];};
    DIAG_MAT &get_mass_mat_md(FieldName thisField)
      { return massMatsMD_[thisField];};
    /** access FE rate of change */
    DENS_MAT &get_field_roc(FieldName thisField)
      { return dot_fields_[thisField]; };
    /** access atomic rate of change contributions to finite element equation */
    DENS_MAT &get_fe_atomic_field_roc(FieldName thisField)
      { return fieldRateNdFiltered_[thisField]; };
    /** access to atomic rate of change */
    DENS_MAT &get_atomic_field_roc(FieldName thisField)
      {return dot_fieldRateNdFiltered_[thisField]; };
    /** access to second time derivative (2roc) */
    DENS_MAT &get_field_2roc(FieldName thisField)
      { return ddot_fields_[thisField]; };
    /** access for field mask */
    Array2D<bool> &get_field_mask() {return fieldMask_;};
    /** access for shape function weigths */
    DIAG_MAT &get_shape_function_weights(){return shpWeight_;};
    /** access to template matrix for control equations */
    SPAR_MAT &get_m_t_template(){return M_T_Template;};
    /** access for mapping from global nodes to nodes overlapping MD */
    Array<int> &get_node_to_overlap_map(){return nodeToOverlapMap_;};
    /** access for mapping from nodes overlapping MD to global nodes */
    Array<int> &get_overlap_to_node_map(){return overlapToNodeMap_;};
    /** number of nodes whose shape function support overlaps with MD */
    int get_nNode_overlap(){return nNodeOverlap_;};
    /** check if atomic quadrature is being used for MD_ONLY nodes */
    bool atom_quadrature_on(){return atomQuadForInternal_;};
    /** check if lambda is localized */
    bool use_localized_lambda(){return useLocalizedLambda_;};
    /** check if matrix should be lumpted for lambda solve */
    bool use_lumped_lambda_solve(){return useLumpedLambda_;};
    /** get scaled shape function matrix */
    SPAR_MAT &get_nhat_overlap() {return NhatOverlap_;};
    /** get scaled shape function matrix weights */
    DENS_VEC &get_nhat_overlap_weights() {return NhatOverlapWeights_;};
    /** get map general atomic shape function matrix to overlap region */
    SPAR_MAT &get_atom_to_overlap_map() {return Trestrict_;};
    /** internal atom to global map */
    Array<int> &get_internal_to_atom_map() {return internalToAtom_;};
    /** ghost atom to global map */
    Array<int> &get_ghost_to_atom_map() {return ghostToAtom_;};
    /** get number of unique FE nodes */
    int get_nNodes() {return nNodes_;};
    /** get number of spatial dimensions */
    int get_nsd() {return nsd_;};
    /** get number of ATC internal atoms on this processor */
    int get_nlocal() {return nLocal_;};
    /** get total number of LAMMPS atoms on this processor */
    int get_nlocal_total() {return nLocalTotal_;};
    /** get number of ATC ghost atoms on this processor */
    int get_nlocal_ghost() {return nLocalGhost_;};
    /** get number of ATC mask atoms on this processor */
    int get_nlocal_mask() {return nLocalMask_;};
    /** get number of ATC atoms used in lambda computation on this processor */
    int get_nlocal_lambda() {return nLocalLambda_;};
    /** get number of ATC internal atoms */
    int get_ninternal() {return nInternal_;}
    /** get number of ATC ghost atoms */
    int get_nghost() {return nGhost_;};
    /** get the current simulation time */
    double get_sim_time() {return simTime_;};
    /** access to lammps atomic positions */
    double ** get_x() {return lammpsInterface_->xatom();};
    /** access to lammps atomic velocities */
    double ** get_v() {return lammpsInterface_->vatom();};
    /** access to lammps atomic forces */
    double ** get_f() {return lammpsInterface_->fatom();};
    /** access to physics model */
    PhysicsModel * get_physics_model() {return physicsModel_; };
    /** access to time filter */
    TimeFilterManager * get_time_filter_manager() {return &timeFilterManager_;};
    /** access to time integrator */
    TimeIntegrator * get_time_integrator() {return timeIntegrator_;};
    /** access to FE engine */
    FE_Engine * get_fe_engine() {return feEngine_;};
    /** access to faceset names */
    const set<PAIR> &get_faceset(const string & name) const {return (feEngine_->get_feMesh())->get_faceset(name);};
    /** access to overlapped ghost shape function */
    SPAR_MAT & get_shape_function_ghost_overlap(){return shpFcnGhostOverlap_;};
    /** return fixed node flag for a field */
    Array2D<bool> & get_fixed_node_flags(FieldName thisField) {return isFixedNode_[thisField];};
    Array<int> & get_node_type() {return nodeType_;};
    void get_field(/*const*/ char ** args, int &argIndex,
                   FieldName &thisField, int &thisIndex);
    DIAG_MAT & get_atomic_weights() {return atomicWeights_;};
    bool track_charge() {return trackCharge_;};
    /** access to set of DENS_MATs accessed by tagging */
    DENS_MAT & get_tagged_dens_mat(const string & tag) {return taggedDensMats_[tag];};
    /** access to set of SPAR_MATs accessed by tagging */
    SPAR_MAT & get_tagged_spar_mat(const string & tag) {return taggedSparMats_[tag];};
    /*@}*/

    //---------------------------------------------------------------
    /** \name boundary integration */
    //---------------------------------------------------------------
    /*@{*/
    void set_boundary_integration_type(int boundaryIntegrationType) 
      {bndyIntType_ = boundaryIntegrationType;};
    void set_boundary_face_set(const set< pair<int,int> > * boundaryFaceSet) 
      {bndyFaceSet_ = boundaryFaceSet;};
    BoundaryIntegrationType parse_boundary_integration
      (int narg, char **arg, const set< pair<int,int> > * boundaryFaceSet);
    /*@}*/

    //---------------------------------------------------------------
    /** \name FE nodesets/sidesets functions */
    //---------------------------------------------------------------
    /*@{*/
    /** mask for computation of fluxes */
    void set_fixed_nodes();

    /** set initial conditions by changing fields */
    void set_initial_conditions();

    /** calculate and set matrix of sources_ */
    void set_sources();

    /** array indicating fixed nodes for all fields */
    map<FieldName, Array2D<bool> > isFixedNode_;

    /** array indicating if the node is boundary, MD, or FE */
    Array<int> nodeType_;

    /** wrapper for FE_Engine's compute_flux functions */
    void compute_flux(const Array2D<bool> & rhs_mask,
                      const FIELDS &fields, 
                      GRAD_FIELDS &flux,
                      const PhysicsModel * physicsModel=NULL);
    /** wrapper for FE_Engine's compute_boundary_flux functions */
    void compute_boundary_flux(const Array2D<bool> & rhs_mask,
                               const FIELDS &fields, 
                               FIELDS &rhs);

    /** wrapper for FE_Engine's compute_rhs_vector functions */
    void compute_rhs_vector(const Array2D<bool> & rhs_mask,
                            const FIELDS &fields, 
                            FIELDS &rhs,
                            const IntegrationDomainType domain, // = FULL_DOMAIN
                            const PhysicsModel * physicsModel=NULL);

    /** evaluate rhs on a specified domain defined by mask and physics model */
    void evaluate_rhs_integral(const Array2D<bool> & rhs_mask,
                           const FIELDS &fields, 
                           FIELDS &rhs,
                           const IntegrationDomainType domain,
                           const PhysicsModel * physicsModel=NULL);
    
    /** assemble various contributions to the heat flux in the atomic region */
    void compute_atomic_sources(const Array2D<bool> & rhs_mask,
                                const FIELDS &fields, 
                                FIELDS &atomicSources);

    /** multiply inverse mass matrix times given data in place */
    // NOTE change massMatInv to map of pointers and allow for consistent mass matrices
    //      inverted using CG
    void apply_inverse_mass_matrix(MATRIX & data, FieldName thisField)
    { 
      data = massMatInv_[thisField]*data;
    };
    /** multiply inverse mass matrix times given data and return result */
    void apply_inverse_mass_matrix(const MATRIX & data_in, MATRIX & data_out,
                                   FieldName thisField)
    {
      data_out = massMatInv_[thisField]*data_in;
    };

    void apply_inverse_md_mass_matrix(MATRIX & data, FieldName thisField)
    { data = massMatMDInv_[thisField]*data; };
    void apply_inverse_md_mass_matrix(const MATRIX & data_in, MATRIX & data_out,
                                   FieldName thisField)
    { data_out = massMatMDInv_[thisField]*data_in; };
    /*@}*/

    //----------------------------------------------------------------
    /** \name FE mapping operations */
    //----------------------------------------------------------------
    /*@{*/
    /** Mapping between unique nodes and nodes overlapping MD region */
    void map_unique_to_overlap(const MATRIX & uniqueData,
                               MATRIX & overlapData);

    /** Mapping between nodes overlapping MD region to unique nodes */
    void map_overlap_to_unique(const MATRIX & overlapData,
                               MATRIX & uniqueData);
    /*@}*/
                
    //----------------------------------------------------------------
    /** \name Interscale operators */
    //----------------------------------------------------------------
    /*@{*/
    /** Restrict (number density) : given w_\alpha, w_I = 1/V_I \sum_\alpha N_{I\alpha} w_\alpha */
    void restrict(const MATRIX &atomData,
                  MATRIX &nodeData);

    /** Restrict based on atomic volume integration : given w_\alpha, w_I = \sum_\alpha N_{I\alpha} w_\alpha V_\alpha */
    void restrict_unscaled(const MATRIX &atomData,
                           MATRIX &nodeData);

    /** Restrict based on atomic volume integration for volumetric quantities : given w_\alpha, w_I = \sum_\alpha N_{I\alpha} w_\alpha */
    void restrict_volumetric_quantity(const MATRIX &atomData,
                                      MATRIX &nodeData);
    void restrict_volumetric_quantity(const MATRIX &atomData,
                                      MATRIX &nodeData,
                                      const SPAR_MAT &shpFcn);

    /** Project based on an atomic volume integration :  given w_\alpha, \sum_\alpha M_\alpha w_I = \sum_\alpha N_{I\alpha} w_\alpha V_\alpha (note mass matrix has V_\alpha in it) */
    void project(const MATRIX &atomData,
                 MATRIX &nodeData,
                 FieldName thisField);
    void project_md(const MATRIX &atomData,
                    MATRIX &nodeData,
                    FieldName thisField);

    /** Project based on an atomic volume integration for volumetric quantities :  given w_\alpha, \sum_\alpha M_\alpha w_I = \sum_\alpha N_{I\alpha} w_\alpha */
    void project_volumetric_quantity(const MATRIX &atomData,
                                     MATRIX &nodeData,
                                     FieldName thisField);
    void project_volumetric_quantity(const MATRIX &atomData,
                                     MATRIX &nodeData,
                                     const SPAR_MAT &shpFcn,
                                     FieldName thisField);
    void project_md_volumetric_quantity(const MATRIX &atomData,
                                        MATRIX &nodeData,
                                        FieldName thisField);
    void project_md_volumetric_quantity(const MATRIX &atomData,
                                        MATRIX &nodeData,
                                        const SPAR_MAT &shpFcn,
                                        FieldName thisField);

    /** Prolong : given w_I,  w_\alpha = \sum_I N_{I\alpha} w_I */
    void prolong(const MATRIX &nodeData,
                 MATRIX &atomData);

    /** Prolong based on scaled shape functions : given w_I,  w_\alpha = \sum_I 1/V_I N_{I\alpha} w_I */
    void prolong_scaled(const MATRIX &nodeData,
                        MATRIX &atomData);

    /** Prolong onto ghost atoms*/
    void prolong_ghost(const MATRIX &nodeData,
                       MATRIX &atomData);
    /*@}*/
        
    //----------------------------------------------------------------
    /** \name Interscale physics constructors */
    //----------------------------------------------------------------
    /*@{*/
    /** Compute atomic masses */
    void compute_atomic_mass(MATRIX &atomicMasses);

    /** Compute atomic charges */
    void compute_atomic_charge(MATRIX &atomicCharges);

    /** Compute atomic temperature 
        possibly using dot product of two velocities */
    void compute_atomic_temperature(MATRIX &T,
                                    const double * const* v,
                                    double ** v2 = NULL);
    /** Compute atomic kinetic energy 
        possibly using dot product of two velocities */
    void compute_atomic_kinetic_energy(MATRIX &T,
                                       const double * const* v,
                                       double ** v2 = NULL);
    /** Compute atomic power : in units dT/dt */
    void compute_atomic_power(MATRIX &dot_T,
                              const double * const* v,
                              const double * const* f);
    void compute_atomic_temperature_roc(MATRIX &dot_T,
                                        const double * const* v,
                                        const double * const* f);
    void compute_atomic_power(MATRIX &dot_T,
                              const double * const* v,
                              const MATRIX & f);
    void compute_atomic_temperature_roc(MATRIX &dot_T,
                                        const double * const* v,
                                        const MATRIX & f);
    
    // NOTE change these when physical mass matrix is used
    /** Compute magnitude of atomic force */
    void compute_atomic_force_strength(MATRIX &forceStrength,
                                       const double * const* f);
    void compute_atomic_force_strength(MATRIX &forceStrength,
                                       const MATRIX & f);
    void compute_atomic_force_dot(MATRIX &forceDot,
                                  const double * const* f1,
                                  const MATRIX & f2);

    /** Compute lambda power at atoms : in units dT/dt  */
    void compute_atomic_temperature_roc(MATRIX &atomicVdotflam,
                                        const MATRIX &lambda,
                                        const double * const* v);
    void compute_atomic_power(MATRIX &atomicVdotflam,
                              const MATRIX &lambda,
                              const double * const* v);

    /** Compute lambda power at atoms : in units dT/dt  */
    void compute_atomic_lambda_power(MATRIX &atomicPower,
                                     const MATRIX &force,
                                     const double * const* v);
    void compute_atomic_temperature_lambda_roc(MATRIX &atomicPower,
                                               const MATRIX &force,
                                               const double * const* v);
                                                                 
    /** Compute lambda power at atoms with explicit application */
    void compute_lambda_power_explicit(MATRIX &lambdaPower,
                                       const MATRIX &lambda,
                                       const double * const* v,
                                       const double dt);

    /** Compute atomic position */
    void compute_atomic_position(DENS_MAT &atomicDisplacement,
                                 const double * const* x);

    /** Compute atomic center of mass */
    void compute_atomic_centerOfMass(DENS_MAT &atomicCom,
                                     const double * const* x);

    /** Compute atomic displacement */
    void compute_atomic_displacement(DENS_MAT &atomicDisplacement,
                                     const double * const* x);

    /** Compute atomic mass displacement */
    void compute_atomic_centerOfMass_displacement(DENS_MAT &atomicMd,
                                                  const double * const* x);

    /** Compute atomic velocity */
    void compute_atomic_velocity(DENS_MAT &atomicVelocity,
                                 const double * const* v);

    /** Compute atomic momentum */
    void compute_atomic_momentum(DENS_MAT &atomicMomentum,
                                 const double * const* v);

    /** Compute atomic acceleration */
    void compute_atomic_acceleration(DENS_MAT &atomicAcceleration,
                                     const double * const* f);

    /** Compute atomic force */
    void compute_atomic_force(DENS_MAT &atomicForce,
			      const double * const* f);
    /*@}*/

    /** allow FE_Engine to construct ATC structures after mesh is constructed */
    void initialize_mesh_data(void); 

  protected:
  
    /** pointer to lammps interface class */
    LammpsInterface * lammpsInterface_;

    /** pointer to physics model */
    PhysicsModel * physicsModel_;

    /** manager for extrinsic models */
    ExtrinsicModelManager extrinsicModelManager_;

    /** method to create physics model */
    void create_physics_model(const PhysicsType & physicsType,
                              string matFileName);

    /** global flag to indicate atoms have changed processor 
        keyed off of atomSwitch (see pre_exchange() which sets atomSwitch =true)
    */
    int globalSwitch_ ;

    /** a global flag to trigger reset of shape function at will */
    bool resetSwitch_;

    /** flag on if initialization has been performed */
    bool initialized_;

    /** flag to determine if charge is tracked */
    bool trackCharge_;

    TimeFilterManager timeFilterManager_;
    TimeIntegrator * timeIntegrator_;

    /** finite element handler */
    FE_Engine * feEngine_;

    /** prescribed data handler */
    PrescribedDataManager * prescribedDataMgr_;

    /** reference atomic coordinates */
    DENS_MAT atomicCoords_;
    DENS_MAT atomicCoordsMask_;

    /** number of unique FE nodes */
    int nNodes_;

    /** Number of Spatial Dimensions */
    int nsd_;

    /** data for handling atoms crossing processors */
    bool atomSwitch_;

    /** reference position of the atoms */
    double ** xref_;
    double Xprd_,Yprd_,Zprd_;
    double XY_,YZ_,XZ_;
    void set_xref();

    /** current time in simulation */
    double simTime_;

    /** re-read reference positions */
    bool readXref_;
    string xRefFile_;

    //---------------------------------------------------------------
    /** \name FE nodesets/sidesets data */
    //---------------------------------------------------------------
    /*@{*/
    /** mask for computation of fluxes */
    Array2D<bool> fieldMask_;

    /** sources */
    FIELDS sources_;
    FIELDS atomicSources_;
    FIELDS extrinsicSources_;

    /** boundary flux quadrature */
    int bndyIntType_;
    const set< pair<int,int> > * bndyFaceSet_;
    set<string> boundaryFaceNames_;
    /*@}*/

    //---------------------------------------------------------------
    /** \name output data */
    //---------------------------------------------------------------
    /*@{*/
    /** base name for output files */
    string outputPrefix_;
  
    /** output frequency */
    int outputFrequency_;
  
    /** sample frequency */
    int sampleFrequency_;
  
    /** step counter */
    int stepCounter_;

    /** sample counter */
    int sampleCounter_;

    /** atomic output */
    /** base name for output files */
    string outputPrefixAtom_;
  
    /** output frequency */
    int outputFrequencyAtom_;

    /** output object */
    OutputManager mdOutputManager_;
    set<string> atomicOutputMask_;
    /*@}*/
    //---------------------------------------------------------------
    /** \name output functions */
    //---------------------------------------------------------------
    /*@{*/
    void output();
    void atomic_output();
    /*@}*/


    //---------------------------------------------------------------
    /** \name member data related to compute_scalar() and compute_vector() */
    //---------------------------------------------------------------
    /*@{*/
    int scalarFlag_;              // 0/1 if compute_scalar() function exists
    int vectorFlag_;              // 0/1 if compute_vector() function exists
    int sizeVector_;              // N = size of global vector
    int globalFreq_;              // frequency global data is available at
    int extScalar_;               // 0/1 if scalar is intensive/extensive
    int extVector_;               // 0/1/-1 if vector is all int/ext/extlist
    int *extList_;                // list of 0/1 int/ext for each vec component
    /*@}*/

    //---------------------------------------------------------------
    /** \name fields and necessary data for FEM */
    //---------------------------------------------------------------
    /*@{*/
    map<FieldName,int> fieldSizes_;
    FIELDS fields_;
    map<FieldName,DIAG_MAT > massMats_;
    map<FieldName,DIAG_MAT > massMatInv_;
    map<FieldName,DIAG_MAT > massMatsMD_;
    map<FieldName,DIAG_MAT > massMatMDInv_;
    virtual void compute_md_mass_matrix(FieldName thisField,
                                        map<FieldName,DIAG_MAT> & massMats) {};
    DENS_MAT consistentMassInverse_;
    FIELDS rhs_;  // for pde
    FIELDS rhsAtomDomain_; // for thermostat
    FIELDS boundaryFlux_; // for thermostat & rhs pde
    /*@}*/


    //---------------------------------------------------------------
    /** \name time integration and filtering fields */
    //---------------------------------------------------------------
    /*@{*/
    
    FIELDS dot_fields_;
    FIELDS ddot_fields_;
    FIELDS dddot_fields_;
    FIELDS dot_fieldsMD_;
    FIELDS ddot_fieldsMD_;
    FIELDS dot_dot_fieldsMD_;
  
    /** Restricted Fields */
    FIELDS fieldNdOld_;
    FIELDS fieldNdFiltered_;
    FIELDS fieldRateNdOld_;
    FIELDS fieldRateNdFiltered_;
    FIELDS dot_fieldRateNdOld_;
    FIELDS dot_fieldRateNdFiltered_;

    /** auxilliary storage */
    FIELDS auxStorage_;
    /*@}*/

  
    //---------------------------------------------------------------
    /** \name quadrature weights */
    //---------------------------------------------------------------
    /*@{*/
    DIAG_MAT NodeVolumes_;
    DIAG_MAT invNodeVolumes_;
    /** atomic quadrature integration weights (V_\alpha) */
    DIAG_MAT atomicWeights_;
    DIAG_MAT atomicWeightsMask_;
    double atomicVolume_;  // global atomic volume for homogeneous set of atoms
    map<int,double> Valpha_;
    AtomicWeightType atomWeightType_;
    /** weighting factor per shape function: 
        shpWeight_(I,I) =  1/N_I = 1/(\sum_\alpha N_{I\alpha}) */
    DIAG_MAT shpWeight_;
    DIAG_MAT fluxMask_;
    DIAG_MAT fluxMaskComplement_;
    /*@}*/
    //---------------------------------------------------------------
    /** \name quadrature weight function */
    //---------------------------------------------------------------
    /*@{*/
    /** determine weighting method for atomic integration */
    void reset_atomicWeightsLattice();
    void reset_atomicWeightsElement();
    void reset_atomicWeightsRegion();
    void reset_atomicWeightsGroup();
    void reset_atomicWeightsMultiscale(const SPAR_MAT & shapeFunctionMatrix,
                                       DIAG_MAT & atomicVolumeMatrix);

    void compute_consistent_md_mass_matrix(const SPAR_MAT & shapeFunctionMatrix,
                                           SPAR_MAT & mdMassMatrix);

    /** resets shape function matrices based on atoms on this processor */
    virtual void reset_nlocal();
    void reset_coordinates();
    void set_atomic_weights();
    virtual void reset_shape_functions();
    void reset_NhatOverlap();
    /*@}*/

    //---------------------------------------------------------------
    /** \name atom data  */
    //---------------------------------------------------------------
    /*@{*/
    /** bitwise comparisons for boundary (ghost) atoms */
    int groupbit_;
    int groupbitGhost_;
    set<int> igroups_;
    set<int> igroupsGhost_;
  
    /** number of atoms of correct type,
        ghosts are atoms outside our domain of interest
        boundary are atoms contributing to boundary flux terms */
    /** Number of "internal" atoms on this processor */
    int nLocal_;
    /** Number of atoms on this processor */
    int nLocalTotal_;
    int nLocalGhost_;
    int nLocalMask_;
    int nLocalLambda_;
    int nInternal_;
    int nGhost_;
    Array<int> internalToAtom_;
    std::map<int,int> atomToInternal_;
    Array<int> ghostToAtom_;
    DENS_MAT ghostAtomCoords_;
    /*@}*/
    //----------------------------------------------------------------
    /** \name  maps and masks */
    //----------------------------------------------------------------
    /*@{*/
    AtomToElementMapType atomToElementMapType_;
    int atomToElementMapFrequency_;
    Array<int> atomToElementMap_;
    Array<int> ghostAtomToElementMap_;
    /** overlap map, from shapeWeights */
    // -1 is no overlap, otherwise entry is overlap index
    Array<int> nodeToOverlapMap_;
    // mapping from overlap nodes to unique nodes
    Array<int> overlapToNodeMap_;
    int nNodeOverlap_;
    Array<bool> elementMask_;
    Array<int> elementToMaterialMap_;
    Array< set <int> > atomMaterialGroups_;
    int regionID_;
    bool atomQuadForInternal_;
    bool useLocalizedLambda_;
    bool useLumpedLambda_;
    /*@}*/
    //----------------------------------------------------------------
    /** \name  Map related functions */
    //----------------------------------------------------------------
    /*@{*/
    bool check_internal(int eltIdx);
    int check_shape_function_type(int nodeIdx);
    bool intersect_ghost(int eltIdx);
    virtual void set_ghost_atoms() = 0;
    /*@}*/

    //----------------------------------------------------------------
    /** \name shape function matrices */
    //----------------------------------------------------------------
    /*@{*/
    // sparse matrix where columns correspond to global node numbering
    // dimensions are numAtoms X numNodes (the transpose of N_{I\alpha} )
    /** shpFcn_ is N_{I\alpha} the un-normalized shape function evaluated at the atoms */
    SPAR_MAT shpFcn_;
    vector<SPAR_MAT > shpFcnDerivs_;
    SPAR_MAT shpFcnGhost_;
    SPAR_MAT shpFcnGhostOverlap_;
    vector<SPAR_MAT > shpFcnDerivsGhost_;
    SPAR_MAT shpFcnMasked_;
    vector<SPAR_MAT > shpFcnDerivsMask_;
    Array<bool> atomMask_;

    /** map from species string tag to the species density */
    map<string,DENS_MAT> taggedDensMats_;
    /** map from species string tag to shape function and weight matrices */
    map<string,SPAR_MAT> taggedSparMats_;

    /** weighted shape function matrices at overlap nodes
        for use with thermostats */
    // dimensions are numAtoms X numNodesOverlap
    SPAR_MAT NhatOverlap_;
    SPAR_MAT Trestrict_;
    DENS_VEC NhatOverlapWeights_;
    /*@}*/


    //---------------------------------------------------------------
    /** \name thermostat data */
    //---------------------------------------------------------------
    /*@{*/
    /** sparse matrix to store elements needed for CG solve */
    SPAR_MAT M_T_Template;
    DIAG_MAT maskMat_;

    bool equilibriumStart_;
    //---------------------------------------------------------------
    /** \name time filtering */
    //---------------------------------------------------------------
    /*@{*/
    /** allocate memory for time filter */
    void init_filter();
    void update_filter(MATRIX &filteredQuantity,
                       const MATRIX &unfilteredQuantity,
                       MATRIX &unfilteredQuantityOld,
                       const double dt);

    double get_unfiltered_coef(const double dt);
                             
    void update_filter_implicit(MATRIX &filteredQuantity,
                                const MATRIX &unfilteredQuantity,
                                const double dt);
    /*@}*/

    /** group computes : type, group_id -> value */
    map< pair < int, int > , double> groupCompute_;

    /** group computes : type, group_id -> value */
    map< pair < string, FieldName > , double> nsetCompute_;

    /** allow FE_Engine to construct data manager after mesh is constructed */
    void construct_prescribed_data_manager (void); 

    //---------------------------------------------------------------
    /** \name restart procedures */
    //---------------------------------------------------------------
    bool useRestart_;
    string restartFileName_;
    virtual void read_restart_data(string fileName_, OUTPUT_LIST & data);
    virtual void write_restart_data(string fileName_, OUTPUT_LIST & data);
    void pack_fields(OUTPUT_LIST & data);

    //---------------------------------------------------------------
    /** \name neighbor reset frequency */
    //---------------------------------------------------------------
    int neighborResetFrequency_;
  };

};

#endif
