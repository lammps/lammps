#ifndef ATC_COUPLING_H
#define ATC_COUPLING_H

#include <set>
#include <map>
#include <string>
#include <utility>

// ATC headers
#include "ATC_Method.h"
#include "ExtrinsicModel.h"

namespace ATC {

  // Forward declarations
  class PrescribedDataManager;
  class AtomicRegulator;
  class TimeIntegrator;
  class ReferencePositions;

  /**
   *  @class  ATC_Coupling
   *  @brief  Base class for atom-continuum coupling 
   */

  class ATC_Coupling : public ATC_Method {

  public: /** methods */
    
    friend class ExtrinsicModel; // friend is not inherited
    friend class ExtrinsicModelTwoTemperature;
    friend class ExtrinsicModelDriftDiffusion; 
    friend class ExtrinsicModelDriftDiffusionConvection;
    friend class ExtrinsicModelElectrostatic;
    friend class ExtrinsicModelElectrostaticMomentum;
    friend class SchrodingerPoissonSolver;
    friend class SliceSchrodingerPoissonSolver;

    /** constructor */
    ATC_Coupling(std::string groupName, double **& perAtomArray, LAMMPS_NS::Fix * thisFix);

    /** destructor */
    virtual ~ATC_Coupling();

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** pre neighbor */
    virtual void pre_neighbor();

    /** pre exchange */
    virtual void pre_exchange();
    virtual void reset_atoms(){};

    /** pre force */
    virtual void pre_force() {};

    /** post force */
    virtual void post_force();

    /** pre integration run */
    virtual void initialize();

    /** post integration run : called at end of run or simulation */
    virtual void finish();

    /** first time, before atomic integration */
    virtual void pre_init_integrate();

    /** Predictor phase, executed between velocity and position Verlet */
    virtual void mid_init_integrate();

    /** Predictor phase, executed after Verlet */
    virtual void post_init_integrate();

    /** Corrector phase, executed before Verlet */
    virtual void pre_final_integrate(){};

    /** Corrector phase, executed after Verlet*/
    
    virtual void post_final_integrate() {lammpsInterface_->computes_addstep(lammpsInterface_->ntimestep()+1);};

    /** pre/post atomic force calculation in minimize */
    virtual void min_pre_force(){};
    virtual void min_post_force(){};

    // data access
    /** get map general atomic shape function matrix to overlap region */
    SPAR_MAT &get_atom_to_overlap_mat() {return atomToOverlapMat_.set_quantity();};
    /** get map general atomic shape function matrix to overlap region */
    SPAR_MAN &atom_to_overlap_mat() {return atomToOverlapMat_;};
    /** check if atomic quadrature is being used for MD_ONLY nodes */
    bool atom_quadrature_on(){return atomQuadForInternal_;};
    const std::set<std::string> & boundary_face_names() {return boundaryFaceNames_;};
    /** access to boundary integration method */
    int boundary_integration_type() {return bndyIntType_;};
    void set_boundary_integration_type(int boundaryIntegrationType) 
    {bndyIntType_ = boundaryIntegrationType;};
    void set_boundary_face_set(const std::set< std::pair<int,int> > * boundaryFaceSet) 
    {bndyFaceSet_ = boundaryFaceSet;};
    BoundaryIntegrationType parse_boundary_integration
      (int narg, char **arg, const std::set< std::pair<int,int> > * boundaryFaceSet);
    TemperatureDefType temperature_def() const {return temperatureDef_;};
    void set_temperature_def(TemperatureDefType tdef) {temperatureDef_ = tdef;};

//--------------------------------------------------------

    /** access to all boundary fluxes */
    FIELDS &boundary_fluxes() {return boundaryFlux_;}; 
    /** wrapper for FE_Engine's compute_boundary_flux functions */
    void compute_boundary_flux(const Array2D<bool> & rhs_mask,
                               const FIELDS &fields, 
                               FIELDS &rhs,
                               const Array< std::set <int> > atomMaterialGroups,
                               const VectorDependencyManager<SPAR_MAT * > * shpFcnDerivs,
                               const SPAR_MAN * shpFcn = NULL,
                               const DIAG_MAN * atomicWeights = NULL,
                               const MatrixDependencyManager<DenseMatrix, bool> * elementMask = NULL,
                               const RegulatedNodes * nodeSet = NULL);
    /** access to full right hand side / forcing vector */
    FIELDS &rhs() {return rhs_;};

    DENS_MAN &field_rhs(FieldName thisField) { return rhs_[thisField]; };
    /** allow FE_Engine to construct ATC structures after mesh is constructed */
    virtual void initialize_mesh_data(void);  

// public for FieldIntegrator
    bool source_atomic_quadrature(FieldName field)  
      { return (sourceIntegration_ == FULL_DOMAIN_ATOMIC_QUADRATURE_SOURCE); }
    ATC::IntegrationDomainType source_integration() 
      { return sourceIntegration_; }

    /** wrapper for FE_Engine's compute_sources */
    void compute_sources_at_atoms(const RHS_MASK & rhsMask,
                                  const FIELDS & fields,
                                  const PhysicsModel * physicsModel,
                                  FIELD_MATS & atomicSources);
    /** computes tangent matrix using atomic quadrature near FE region */
   void masked_atom_domain_rhs_tangent(const std::pair<FieldName,FieldName> row_col,
                                       const RHS_MASK & rhsMask,      
                                       const FIELDS & fields,                  
                                       SPAR_MAT & stiffness,
                                       const PhysicsModel * physicsModel);
    /** wrapper for FE_Engine's compute_rhs_vector functions */
    void compute_rhs_vector(const RHS_MASK & rhs_mask,
                            const FIELDS &fields, 
                            FIELDS &rhs,
                            const IntegrationDomainType domain, // = FULL_DOMAIN
                            const PhysicsModel * physicsModel=NULL);
   /** wrapper for FE_Engine's compute_tangent_matrix */
   void compute_rhs_tangent(const std::pair<FieldName,FieldName> row_col,
                            const RHS_MASK & rhsMask,      
                            const FIELDS & fields,                  
                            SPAR_MAT & stiffness,
                            const IntegrationDomainType integrationType,
                            const PhysicsModel * physicsModel=NULL);

   /** PDE type */
   WeakEquation::PDE_Type pde_type(const FieldName fieldName) const;
   /** is dynamic PDE */
   bool is_dynamic(const FieldName fieldName) const;

// public for ImplicitSolveOperator
    /** return pointer to PrescribedDataManager */
    PrescribedDataManager * prescribed_data_manager() 
      { return prescribedDataMgr_; }
// public for Kinetostat
    DIAG_MAT &get_mass_mat(FieldName thisField)
      { return massMats_[thisField].set_quantity();};
    /** */
    DENS_MAN &atomic_source(FieldName thisField){return atomicSources_[thisField];};


    //---------------------------------------------------------------
    /** \name materials  */
    //---------------------------------------------------------------
    /*@{*/
    /** access to element to material map */
    Array<int> &element_to_material_map(void){return elementToMaterialMap_;}
    /*@}*/

    /** check if method is tracking charge */
    bool track_charge() {return trackCharge_;};

    void set_mass_mat_time_filter(FieldName thisField,TimeFilterManager::FilterIntegrationType filterIntegrationType);

    /** return referece to ExtrinsicModelManager */
    ExtrinsicModelManager & extrinsic_model_manager() 
      { return extrinsicModelManager_; }
    /** access to time integrator */
    const TimeIntegrator * time_integrator(const FieldName & field) const {
      _ctiIt_ = timeIntegrators_.find(field);
      if (_ctiIt_ == timeIntegrators_.end()) return NULL;
      return _ctiIt_->second;
    };

    //---------------------------------------------------------------
    /** \name managers */
    //---------------------------------------------------------------
    /*@{*/
    /** allow FE_Engine to construct data manager after mesh is constructed */
    void construct_prescribed_data_manager (void); 
    /** method to create physics model */
    void create_physics_model(const PhysicsType & physicsType,
                              std::string matFileName);
    /** access to physics model */
    PhysicsModel * physics_model() {return physicsModel_; };
    /*@}*/

    //---------------------------------------------------------------
    /** \name creation */
    //---------------------------------------------------------------
    /*@{*/
    /** set up atom to material identification */
    virtual void reset_atom_materials();
    /** */
    void reset_node_mask();
    /** */
    void reset_overlap_map();
    /*@}*/

    //---------------------------------------------------------------
    /** \name output/restart */
    //---------------------------------------------------------------
    /*@{*/
    void pack_fields(RESTART_LIST & data);
    void output() { ATC_Method::output(); }
    /*@}*/

    //---------------------------------------------------------------
    /** \name initial & boundary conditions  */
    //---------------------------------------------------------------
    /*@{*/
    /** mask for computation of fluxes */
    void set_fixed_nodes();
    /** set initial conditions by changing fields */
    void set_initial_conditions();
    /*@}*/

    //---------------------------------------------------------------
    /** \name sources  */
    //---------------------------------------------------------------
    /** calculate and set matrix of sources_ */
    void set_sources();
    /** assemble various contributions to the heat flux in the atomic region */
    void compute_atomic_sources(const Array2D<bool> & rhs_mask,
                                const FIELDS &fields, 
                                FIELDS &atomicSources);

    DENS_MAT &get_source(FieldName thisField){return sources_[thisField].set_quantity();};
    DENS_MAN &source(FieldName thisField){return sources_[thisField];};
    FIELDS & sources(){return sources_;};
    /** access to name atomic source terms */
    DENS_MAT &get_atomic_source(FieldName thisField){return atomicSources_[thisField].set_quantity();};
    /** access to name extrinsic source terms */
    DENS_MAT &get_extrinsic_source(FieldName thisField){return extrinsicSources_[thisField].set_quantity();};
    DENS_MAN &extrinsic_source(FieldName thisField){return extrinsicSources_[thisField];};

    /** nodal projection of a field through the physics model */

    void nodal_projection(const FieldName & fieldName,
                          const PhysicsModel * physicsModel,
                          FIELD & field);
    /*@}*/

    //---------------------------------------------------------------
    /** \name fluxes  */
    //---------------------------------------------------------------
    /*@{*/
    /** access for field mask */
    Array2D<bool> &field_mask() {return fieldMask_;}; 
    /** create field mask */
    void reset_flux_mask();
    /** wrapper for FE_Engine's compute_flux functions */
    void compute_flux(const Array2D<bool> & rhs_mask,
                      const FIELDS &fields, 
                      GRAD_FIELD_MATS &flux,
                      const PhysicsModel * physicsModel=NULL);
    /** evaluate rhs on the atomic domain which is near the FE region */
    void masked_atom_domain_rhs_integral(const Array2D<bool> & rhs_mask,
                                         const FIELDS &fields, 
                                         FIELDS &rhs,
                                         const PhysicsModel * physicsModel);
    /** evaluate rhs on a specified domain defined by mask and physics model */
    void evaluate_rhs_integral(const Array2D<bool> & rhs_mask,
                           const FIELDS &fields, 
                           FIELDS &rhs,
                           const IntegrationDomainType domain,
                           const PhysicsModel * physicsModel=NULL);
    /** access to boundary fluxes */
    DENS_MAT &get_boundary_flux(FieldName thisField){return boundaryFlux_[thisField].set_quantity();};
    DENS_MAN &boundary_flux(FieldName thisField){return boundaryFlux_[thisField];};
    /** access to finite element right-hand side data */
    DENS_MAT &get_field_rhs(FieldName thisField)
      { return rhs_[thisField].set_quantity(); };
    /*@}*/

    //---------------------------------------------------------------
    /** \name mass matrices  */
    //---------------------------------------------------------------
    /*@{*/
    // atomic field time derivative filtering
    virtual void init_filter(void);
    // mass matrix filtering
    void delete_mass_mat_time_filter(FieldName thisField);
    /** compute mass matrix for requested field */
    void compute_mass_matrix(FieldName thisField, PhysicsModel * physicsModel = NULL);
    /** updates filtering of MD contributions */
    void update_mass_matrix(FieldName thisField);

  private: /** methods */
    ATC_Coupling(); // do not define

  protected: /** data */

    //---------------------------------------------------------------
    /** initialization routines */
    //---------------------------------------------------------------
    /** sets up all data necessary to define the computational geometry */
    virtual void set_computational_geometry();
    /** constructs all data which is updated with time integration, i.e. fields */
    virtual void construct_time_integration_data();
    /** create methods, e.g. time integrators, filters */
    virtual void construct_methods();
    /** set up data which is dependency managed */
    virtual void construct_transfers();
    /** sets up mol transfers */
    virtual void construct_molecule_transfers();
    /** sets up accumulant & interpolant */
    virtual void construct_interpolant();

    //---------------------------------------------------------------
    /** status */
    //---------------------------------------------------------------
    /*@{*/
    /** flag on if FE nodes in MD region should be initialized to projected MD values */
    bool consistentInitialization_;
    bool equilibriumStart_;
    bool useFeMdMassMatrix_;
    /** flag to determine if charge is tracked */
    bool trackCharge_;
    /** temperature definition model */
    TemperatureDefType temperatureDef_;
    /*@}*/

    //---------------------------------------------------------------
    /** \name managers */
    //---------------------------------------------------------------
    /*@{*/
    /** prescribed data handler */
    PrescribedDataManager * prescribedDataMgr_;
    /** pointer to physics model */
    PhysicsModel * physicsModel_;
    /** manager for extrinsic models */
    ExtrinsicModelManager extrinsicModelManager_;
    /** manager for regulator */
    AtomicRegulator * atomicRegulator_;
    /** managers for time integrators per field */
    std::map<FieldName,TimeIntegrator * > timeIntegrators_;
    /** time integrator iterator */
    mutable std::map<FieldName,TimeIntegrator * >::iterator _tiIt_;
    /** time integrator const iterator */
    mutable std::map<FieldName,TimeIntegrator * >::const_iterator _ctiIt_;
    /*@}*/

    //---------------------------------------------------------------
    /** materials */
    //---------------------------------------------------------------
    /*@{*/
    Array<int> elementToMaterialMap_;  // ATOMIC_TAG * elementToMaterialMap_;
    /** atomic ATC material tag */
    
    
    Array< std::set <int> > atomMaterialGroups_;  // ATOMIC_TAG*atomMaterialGroups_;
    Array< std::set <int> > atomMaterialGroupsMask_;  // ATOMIC_TAG*atomMaterialGroupsMask_;
    /*@}*/

    //---------------------------------------------------------------
    /** computational geometry */
    //---------------------------------------------------------------
    /*@{*/
    bool atomQuadForInternal_;  
    MatrixDependencyManager<DenseMatrix, bool> * elementMask_;
    MatrixDependencyManager<DenseMatrix, bool> * elementMaskMass_;
    MatrixDependencyManager<DenseMatrix, bool> * elementMaskMassMd_;
    MatrixDependencyManager<DenseMatrix, bool> * create_full_element_mask();
    MatrixDependencyManager<DenseMatrix, int> * create_element_set_mask(const std::string & elementSetName);
    LargeToSmallAtomMap * internalToMask_;
    MatrixDependencyManager<DenseMatrix, int> * internalElement_;
    MatrixDependencyManager<DenseMatrix, int> * ghostElement_;
    DenseMatrixTransfer<int> * nodalGeometryType_;
    /*@}*/

    /** \name boundary integration */
    /*@{*/
    /** boundary flux quadrature */
    int bndyIntType_;
    const std::set< std::pair<int,int> > * bndyFaceSet_;
    std::set<std::string> boundaryFaceNames_;
    /*@}*/

    //----------------------------------------------------------------
    /** \name shape function matrices */
    //----------------------------------------------------------------
    /*@{*/
    DIAG_MAN * atomicWeightsMask_;
    SPAR_MAN * shpFcnMask_;
    VectorDependencyManager<SPAR_MAT * > * shpFcnDerivsMask_;
    Array<bool> atomMask_;  
    SPAR_MAN atomToOverlapMat_;
    DIAG_MAN nodalMaskMat_;
    /*@}*/

    //---------------------------------------------------------------
    /** \name PDE data */
    //---------------------------------------------------------------
    /*@{*/
    /** mask for computation of fluxes */
    Array2D<bool> fieldMask_; 
    DIAG_MAT fluxMask_;
    DIAG_MAT fluxMaskComplement_;
    /** sources */
    FIELDS sources_;
    FIELDS atomicSources_;
    FIELDS extrinsicSources_;
    ATC::IntegrationDomainType sourceIntegration_;
    SPAR_MAT stiffnessAtomDomain_;
    /** rhs/forcing terms */
    FIELDS rhs_;  // for pde
    FIELDS rhsAtomDomain_; // for thermostat
    FIELDS boundaryFlux_; // for thermostat & rhs pde
    /*@}*/

    // workspace variables
    mutable DENS_MAT _deltaQuantity_;
  };
};

#endif
