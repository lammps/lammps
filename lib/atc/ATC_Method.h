#ifndef ATC_METHOD_H
#define ATC_METHOD_H

// ATC_Method headers
#include "ATC_TypeDefs.h"
#include "PhysicsModel.h"
#include "MatrixLibrary.h"
#include "Array.h"
#include "Array2D.h"
#include "OutputManager.h"
#include "Function.h"
#include "FE_Element.h"
#include "TimeFilter.h"
#include "LammpsInterface.h"
#include "FE_Engine.h"
#include "ExtrinsicModel.h"
#include "InterscaleOperators.h"
#include "TransferLibrary.h"
#include "GhostManager.h"

// Other headers
#include <vector>
#include <set>
#include <utility>
#include <string>
#include <map>



namespace ATC {

  // forward declarations
  class AtomTimeIntegrator;

  /**
   *  @class  ATC_Method
   *  @brief  Base class for atom-continuum coupling or transfer operators
   */

  class ATC_Method {

  public: /** methods */

    /** constructor */
    ATC_Method(std::string groupName, double **& perAtomArray, LAMMPS_NS::Fix * thisFix);

    /** destructor */
    virtual ~ATC_Method();

    std::string version() {return "2.0";}

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);
    void parse_field(/*const*/ char ** args, int &argIndex,
                   FieldName &thisField, int &thisIndex);

    /** initialize any computes that will be needed prior to the first timestep */
    virtual void init_computes() {
      lammpsInterface_->computes_addstep(lammpsInterface_->ntimestep());
    };

    /** pre integration run */
    virtual void initialize();

    /** Predictor phase, executed before Verlet */
    virtual void pre_init_integrate() {
      feEngine_->partition_mesh(); 
      update_step();
    };

    /** Predictor phase, Verlet first step for velocity */
    virtual void init_integrate_velocity(); 

    /** Predictor phase, executed between velocity and position Verlet */
    virtual void mid_init_integrate(){};

    /** Predictor phase, Verlet first step for position */
    virtual void init_integrate_position(); 

    /** Predictor phase, executed after Verlet */
    virtual void post_init_integrate();

    /** Corrector phase, executed before Verlet */
    virtual void pre_final_integrate(){};

    /** Corrector phase, Verlet second step for velocity */
    virtual void final_integrate();

    /** Corrector phase, executed after Verlet*/
    virtual void post_final_integrate();

    /** post integration run : called at end of run or simulation */
    virtual void finish();

    /** pre/post atomic force calculation */
    virtual void pre_force(){};

    /** pre/post atomic force calculation in minimize */
    virtual void min_pre_force(){};
    virtual void min_post_force();

    /** called at end of step for run or minimize */
    virtual void end_of_step();

    //---------------------------------------------------------------
    /** \name memory management and processor information exchange */
    //---------------------------------------------------------------
    /*@{*/
    /** pre_exchange  is our indicator that atoms have moved across processors */
    virtual void pre_exchange();
    void setup_pre_exchange();
    virtual void pre_neighbor();
    virtual void post_force();
    int doubles_per_atom() const;
    virtual int memory_usage();
    virtual void grow_arrays(int);
    void copy_arrays(int, int);
    int pack_exchange(int, double *);
    int unpack_exchange(int, double *);
    int comm_forward(void) {return sizeComm_;}
    int pack_comm(int , int *, double *, int, int *);
    void unpack_comm(int, int, double *);
    /*@}*/

    //---------------------------------------------------------------
    /** \name managers */
    //---------------------------------------------------------------
    /*@{*/
    /** access to FE engine */
    const FE_Engine * fe_engine() const {return feEngine_;};
    /** access to interscale manager */
    InterscaleManager & interscale_manager() {return interscaleManager_;};
    /** access to lammps interface */
    
    LammpsInterface const * lammps_interface() const {return lammpsInterface_;};
    /** access to time filter */
    TimeFilterManager * time_filter_manager() {return &timeFilterManager_;};
    /*@}*/

    //---------------------------------------------------------------
    /** \name access methods for output computation data */
    //---------------------------------------------------------------
    /*@{*/
    /** compute scalar for output */
    virtual double compute_scalar() {return 0.;}
    /** compute vector for output */
    virtual double compute_vector(int n) {return 0.;}
    /** compute vector for output */
    virtual double compute_array(int irow, int icol) {return 0.;}; 
    int scalar_flag() const {return scalarFlag_;}
    int vector_flag() const {return vectorFlag_;}
    int size_vector() const {return sizeVector_;}
    int peratom_flag() const {return sizePerAtomCols_ > 0;}
    int size_peratom_cols() const {return sizePerAtomCols_;}
    int peratom_freq() const {return 1;}
    void set_peratom_pointer(double ** & ptr) { ptr = perAtomOutput_; }
    int global_freq() const {return scalarVectorFreq_;};
    int extscalar() const {return extScalar_;};
    int extvector() const {return extVector_;};
    int * extlist() {return extList_;};
    int thermo_energy_flag() const {return thermoEnergyFlag_;};
    bool parallel_consistency() const {return parallelConsistency_;};

    /** access to step number */
    int step() const {return stepCounter_;};
    double time() const {return simTime_;};
    double dt() const {return lammpsInterface_->dt();}

    /** time/step functions */
    bool sample_now(void) const
    { 
      int s = step();
      bool now =  ( (sampleFrequency_ > 0) && (s % sampleFrequency_ == 0));
      return now;
    }
    bool output_now(void) const
    { 
      int s = step();

      bool now = ( (outputFrequency_ > 0) && (s == 1 || s % outputFrequency_ == 0) ); 
      now = now || outputNow_;
      return now;
    }
    double output_index(void) const
    {
      if (outputTime_) return time();
      else             return step();
    }

    /** print tracked types and groups */
    int print_tracked() const
    {
      std::string msg = "species:\n";
      for(unsigned int i = 0; i < typeList_.size(); i++) {
        msg+=" type:"+ATC_Utility::to_string(typeList_[i])+" name: "+ typeNames_[i]+"\n";      }
      for(unsigned int i = 0; i < groupList_.size(); i++) {
        msg+=" group (bit):"+ATC_Utility::to_string(groupList_[i])+" name: "+ groupNames_[i]+"\n";
      }
      ATC::LammpsInterface::instance()->print_msg_once(msg);
      return typeList_.size()+groupList_.size();
    }
    std::vector<std::string>  tracked_names() const
    {
      std::vector<std::string> names(typeList_.size()+groupList_.size());
      int j = 0;
      for(unsigned int i = 0; i < typeList_.size(); i++) {
        names[j++] = typeNames_[i];
      }
      for(unsigned int i = 0; i < groupList_.size(); i++) {
        names[j++] = groupNames_[i];
      }
      return names;
    }
    int tag_to_type(std::string tag) const {
      for(unsigned int i = 0; i < typeList_.size(); i++) {
        if (tag == typeNames_[i]) return typeList_[i];
      }
      return -1;
    }
    int type_index(int t) const {
      for(unsigned int i = 0; i < typeList_.size(); i++) {
        if (t == typeList_[i]) return i;
      }
      return -1;
    }
    /*@}*/

    //---------------------------------------------------------------
    /** \name Access methods for sizes */
    //---------------------------------------------------------------
    /*@{*/
    /** get number of unique FE nodes */
    int num_nodes() const {return nNodes_;};
    /** get number of spatial dimensions */
    int nsd() const {return nsd_;};
    /** get number of ATC internal atoms on this processor */
    int nlocal() const {return nLocal_;};
    /** get total number of LAMMPS atoms on this processor */
    int nlocal_total() const {return nLocalTotal_;};
    /** get number of ATC ghost atoms on this processor */
    int nlocal_ghost() const {return nLocalGhost_;};
    /** get the number of all LAMMPS real and parallel ghost atoms on this processor */
    int nproc_ghost() const {return nLocalTotal_ + lammpsInterface_->nghost();};
    /** match group bits */
    bool is_ghost_group(int grpbit) { return (grpbit == groupbitGhost_); }
    bool is_internal_group(int grpbit) { return (grpbit == groupbit_); }
    unsigned int ntracked() { return typeList_.size()+groupList_.size(); }
    bool has_tracked_species() { return typeList_.size()+groupList_.size() > 0; }
    /*@}*/

    virtual void initialize_mesh_data(void){meshDataInitialized_=true;}

    //---------------------------------------------------------------
    /** \name Access methods for data used by various methods */
    //---------------------------------------------------------------
    /*@{*/



    /** access to name FE fields */
    DENS_MAN &field(FieldName thisField){return fields_[thisField];};
    /** access to FE field time derivatives */
    DENS_MAT &get_dot_field(FieldName thisField){return dot_fields_[thisField].set_quantity();};
    DENS_MAN &dot_field(FieldName thisField){return dot_fields_[thisField];};
    /** access to nodal fields of atomic variables */
    DENS_MAT &get_atomic_field(FieldName thisField)
      { return nodalAtomicFields_[thisField].set_quantity(); };
    DENS_MAN &nodal_atomic_field(FieldName thisField)
      { return nodalAtomicFields_[thisField]; };
    /** access to all fields */
    FIELDS &fields() {return fields_;};
    /** access to all fields rates of change (roc) */
    FIELDS &fields_roc() {return dot_fields_;};
    /** add a new field */
    void add_fields(std::map<FieldName,int> & newFieldSizes);
    /** access FE rate of change */
    DENS_MAT &get_field_roc(FieldName thisField)
      { return dot_fields_[thisField].set_quantity(); };
    DENS_MAN &field_roc(FieldName thisField)
      { return dot_fields_[thisField]; };
    /** access atomic rate of change contributions to finite element equation */

    DENS_MAT &get_nodal_atomic_field_roc(FieldName thisField)
      { return nodalAtomicFieldsRoc_[thisField].set_quantity(); };
    DENS_MAN &nodal_atomic_field_roc(FieldName thisField)
      { return nodalAtomicFieldsRoc_[thisField]; };
    /** access to second time derivative (2roc) */

    DENS_MAT &get_field_2roc(FieldName thisField)
      { return ddot_fields_[thisField].set_quantity(); };
    DENS_MAN &field_2roc(FieldName thisField)
      { return ddot_fields_[thisField]; };
    /** access to third time derivative (3roc) */
    DENS_MAT &get_field_3roc(FieldName thisField)
      { return dddot_fields_[thisField].set_quantity(); };
    DENS_MAN &field_3roc(FieldName thisField)
      { return dddot_fields_[thisField]; };
    /** group bit */
    int groupbit() {return groupbit_;};
    /** group bit for ghosts */
    int groupbit_ghost() {return groupbitGhost_;};
    /** internal atom to global map */
    const Array<int> &internal_to_atom_map() {return internalToAtom_;};
    /** ghost atom to global map */
    const Array<int> &ghost_to_atom_map() {return ghostToAtom_;};
    const std::map<int,int> & atom_to_internal_map() {return atomToInternal_;};
    /** access to xref */
    double ** xref() {return xref_;};
    /** access to faceset names */
    const std::set<PAIR> &faceset(const std::string & name) const {return (feEngine_->fe_mesh())->faceset(name);};
    DENS_VEC copy_nodal_coordinates(int i) const { return feEngine_->fe_mesh()->nodal_coordinates(i); }
    /** access to set of DENS_MANs accessed by tagging */

    DENS_MAN & tagged_dens_man(const std::string & tag) {return taggedDensMan_[tag];};
    /** access to atom to element type map */
    AtomToElementMapType atom_to_element_map_type() {return atomToElementMapType_;};
    /** access to atom to element update frequency */
    int atom_to_element_map_frequency() {return atomToElementMapFrequency_;};
    /** flag on whether atc is initialized */
    bool is_initialized() const {return initialized_;};
    /** step number within a run or minimize */
    int local_step() const {return localStep_;};
    /** flags whether a methods reset is required */
    virtual bool reset_methods() const {return (!initialized_) || timeFilterManager_.need_reset() || timeFilterManager_.end_equilibrate() || ghostManager_.need_reset();};
    /** sizes of each field being considered */
    const std::map<FieldName,int> & field_sizes() {return fieldSizes_;};
    /*@}*/
   /** compute the consistent MD mass matrix */
   void compute_consistent_md_mass_matrix(const SPAR_MAT & shapeFunctionMatrix,
                                          SPAR_MAT & mdMassMatrix) const;
    /** access to molecule ids */
    const std::map<std::string,std::pair<MolSize,int> > & molecule_ids() const {return moleculeIds_;};
    /** access to internal element set */
    const std::string & internal_element_set() {return internalElementSet_;};
    

    //----------------------------------------------------------------
    /** \name mass matrix operations */
    //----------------------------------------------------------------
    //      inverted using GMRES
    void apply_inverse_mass_matrix(MATRIX & data, FieldName thisField)
    {
      
      if (useConsistentMassMatrix_(thisField)) {
        //data = consistentMassInverse_*data;
        data = (consistentMassMatsInv_[thisField].quantity())*data;
        return;
      }
      data = (massMatsInv_[thisField].quantity())*data;
    };
    /** multiply inverse mass matrix times given data and return result */
    void apply_inverse_mass_matrix(const MATRIX & data_in, MATRIX & data_out,
                                   FieldName thisField)
    {
      if (useConsistentMassMatrix_(thisField)) { 
        //data_out = consistentMassInverse_*data_in;
        data_out = (consistentMassMatsInv_[thisField].quantity())*data_in;
        return;
      }
      data_out = (massMatsInv_[thisField].quantity())*data_in;
    };
    void apply_inverse_md_mass_matrix(const MATRIX & data_in, MATRIX & data_out,
                                   FieldName thisField)
    { data_out = (massMatsMdInv_[thisField].quantity())*data_in; };


    DIAG_MAN &mass_mat(FieldName thisField)
      { return massMats_[thisField];};

    //---------------------------------------------------------------
    /** \name mass matrices  */
    //---------------------------------------------------------------
    /*@{*/
    /** access to  mass matrices */
    /** access to inverse mass matrices */
    DIAG_MAT &get_mass_mat_inv(FieldName thisField)
      { return massMatsInv_[thisField].set_quantity();};
    DIAG_MAN &mass_mat_inv(FieldName thisField)
      { return massMatsInv_[thisField];};
    /** nodal volumes associated with the atoms, used for the atomic mass matrix */
    AdmtfShapeFunctionRestriction * nodalAtomicVolume_;

    void register_mass_matrix_dependency(DependencyManager * dependent,
                                         FieldName thisField)
    {
      if (useConsistentMassMatrix_(thisField)) { 
        consistentMassMatsInv_[thisField].register_dependence(dependent);
        return;
      }
      massMatsInv_[thisField].register_dependence(dependent);
    };
    
    void apply_inverse_md_mass_matrix(MATRIX & data, FieldName thisField)
    { data = (massMatsMdInv_[thisField].quantity())*data; };
    void register_md_mass_matrix_dependency(DependencyManager * dependent,
                                            FieldName thisField)
    {massMatsMdInv_[thisField].register_dependence(dependent);}
//    /** determine weighting method for atomic integration */
//    void compute_consistent_md_mass_matrix(const SPAR_MAT & shapeFunctionMatrix,
 //                                          SPAR_MAT & mdMassMatrix);
    virtual void compute_md_mass_matrix(FieldName thisField,
                                        DIAG_MAT & massMat) {};
/** access to md mass matrices */

    DIAG_MAN &mass_mat_md_inv(FieldName thisField)
      { return massMatsMdInv_[thisField];};
    DIAG_MAN &set_mass_mat_md(FieldName thisField)
     { return massMatsMd_[thisField]; };
    const DIAG_MAN &mass_mat_md(FieldName thisField) const
     {
       MASS_MATS::const_iterator man =  massMatsMd_.find(thisField);
       if (man == massMatsMd_.end() ) {
         std::string msg = " MD mass for " + field_to_string(thisField) +  " does not exist";
         throw ATC_Error(msg);
       }
       return man->second;
     };
    /*@}*/
                
    //----------------------------------------------------------------
    /** \name Interscale operators */
    //----------------------------------------------------------------
    bool use_md_mass_normalization() const { return mdMassNormalization_;}  
    bool kernel_based() { return kernelBased_; }
    bool kernel_on_the_fly() const { return kernelOnTheFly_;}
    bool has_kernel_function() { return kernelFunction_ != NULL; }
    KernelFunction * kernel_function() { return kernelFunction_; }
    std::vector<int> & type_list() { return typeList_; }
    std::vector<int> & group_list() { return groupList_; }
    SPAR_MAN* interpolant() { return shpFcn_; }
    SPAR_MAN* accumulant() { return accumulant_; }
    DIAG_MAN* accumulant_weights() { return  accumulantWeights_;}
    DIAG_MAN* accumulant_inverse_volumes() { return accumulantInverseVolumes_; }
    int accumulant_bandwidth() const { if (accumulantBandwidth_) return accumulantBandwidth_; else return feEngine_->num_nodes(); };

    PerAtomQuantity<double> * atom_coarsegraining_positions() { return  atomCoarseGrainingPositions_; }
    PerAtomQuantity<double> * atom_reference_positions() { return  atomReferencePositions_; }
    PerAtomQuantity<int> * atom_to_element_map() { return atomElement_;}

    double ke_scale() { return keScale_; }
    double pe_scale() { return peScale_; }
    
    /** from a atom group, find the nodes that have non-zero shape function contributions */

    bool nodal_influence(const int groupbit, std::set<int>& nset, std::set<int>& aset, double tol =1.e-8);
    int nodal_influence(const int groupbit, std::set<int>& nset, std::set<int>& aset, 
      bool ghost,
      double tol =1.e-8);
    /*@{*/

    

 /** Restrict based on atomic volume integration for volumetric quantities : given w_\alpha, w_I = \sum_\alpha N_{I\alpha} w_\alpha */
    void restrict_volumetric_quantity(const MATRIX &atomData,
                                      MATRIX &nodeData); 
    void restrict_volumetric_quantity(const MATRIX &atomData,
                                      MATRIX &nodeData,
                                      const SPAR_MAT & shpFcn);

    /** Prolong : given w_I,  w_\alpha = \sum_I N_{I\alpha} w_I */
    void prolong(const MATRIX &nodeData, MATRIX &atomData);

    //---------------------------------------------------------------
    /** \name quadrature weights */
    //---------------------------------------------------------------
    PerAtomDiagonalMatrix<double> * create_atom_volume();

    //---------------------------------------------------------------
    /** \name access to potential energy reference */
    //---------------------------------------------------------------
    /*@{*/
    DENS_MAN * nodal_ref_potential_energy(void) { return nodalRefPotentialEnergy_; }

  protected: /** methods */
    /** time functions */
    void set_time(double t=0) {simTime_=t;};
    void update_time(double alpha = 1.0)  
    {
      double dt = lammpsInterface_->dt();
      simTime_ += alpha*dt;
      if (dt == 0.0) simTime_ = stepCounter_;
    }
    // note step counter different than lammps step e.g. min
    void update_step(void)  { ++stepCounter_; }

    //---------------------------------------------------------------
    /** initialization routines */
    //---------------------------------------------------------------
    /** gets baseline data from continuum model */
    virtual void set_continuum_data();
    /** sets up all data necessary to define the computational geometry */
    virtual void set_computational_geometry();
    /** constructs all data which is updated with time integration, i.e. fields */
    virtual void construct_time_integration_data() = 0;
    /** create methods, e.g. time integrators, filters */
    virtual void construct_methods();
    /** set up data which is dependency managed */
    virtual void construct_transfers();
    /** sets up accumulant & interpolant */
    virtual void construct_interpolant()=0;
    /** sets up mol transfers */
    virtual void construct_molecule_transfers()=0;

    /** update the peratom output pointers */
    void update_peratom_output(void);

    virtual void  read_restart_data(std::string fileName_, RESTART_LIST & data);
    virtual void write_restart_data(std::string fileName_, RESTART_LIST & data);
    void pack_fields(RESTART_LIST & data); 

   /** mass matrices */
    MASS_MATS massMats_;
    MASS_MATS massMatsInv_;
    MASS_MATS massMatsMd_;
    MASS_MATS massMatsMdInstantaneous_;
    MASS_MATS massMatsMdInv_;
    MASS_MATS massMatsFE_;
    MASS_MATS massMatsAq_;
    MASS_MATS massMatsAqInstantaneous_;
    Array<bool> useConsistentMassMatrix_;
    std::map<FieldName,SPAR_MAN> consistentMassMats_;
    std::map<FieldName,DENS_MAN> consistentMassMatsInv_;
    std::map<FieldName,TimeFilter * > massMatTimeFilters_;

    //---------------------------------------------------------------
    /** \name quadrature weight function */
    //---------------------------------------------------------------
    /*@{*/
    void write_atomic_weights(const std::string filename,const DIAG_MAT & atomicVolumeMatrix);

    /** resets shape function matrices based on atoms on this processor */
    virtual void reset_nlocal();
    virtual void reset_coordinates();
    /*@}*/

    /** re-read reference positions */
    bool read_atomic_ref_positions(const char * filename);
    void remap_ghost_ref_positions(void);
    void adjust_xref_pbc();

    //---------------------------------------------------------------
    /** \name output functions */
    //---------------------------------------------------------------
    /*@{*/
    virtual void output();
    void compute_nodeset_output(void);
    void compute_faceset_output(void);
    void compute_elementset_output(void);
    /*@}*/

    //---------------------------------------------------------------
    /** \name types, groups, and molecules */
    //---------------------------------------------------------------
    /*@{*/
    /** map from species string tag to LAMMPS type id or group bit */
    std::map<std::string,std::pair<MolSize,int> > moleculeIds_;
    /** a list of lammps types & groups ATC tracks */
    std::vector<std::string> typeNames_;
    std::vector<std::string> groupNames_;
    std::vector<int> typeList_;
    std::vector<int> groupList_;
    /*@}*/

    void reset_fields();

  private: /** methods */ 
    ATC_Method(); // do not define

  protected: /** data */
    
    /* parsed input requires changes */
    bool needReset_;

    // managers
    /** pointer to lammps interface class */
    LammpsInterface * lammpsInterface_;

    /** manager for atomic quantities and interscale operations */
    InterscaleManager interscaleManager_;

    TimeFilterManager timeFilterManager_;

    /** check to see if we are integrating the atoms */
    bool integrateInternalAtoms_;
    /** object which integrates atoms */
    AtomTimeIntegrator * atomTimeIntegrator_;
    /** objects which handles integration and modification of ghost atoms */
    GhostManager ghostManager_;

    /** finite element handler */
    FE_Engine * feEngine_;

    // status flags
    /** flag on if initialization has been performed */
    bool initialized_;
    bool meshDataInitialized_;

    /** counter for steps of a run or minimize */
    int localStep_;

    // sizes
    /** size of per atom communication */
    int sizeComm_;

    /** atomic coordinates for coarse graining */
    PerAtomQuantity<double> * atomCoarseGrainingPositions_;
    PerAtomQuantity<double> * atomGhostCoarseGrainingPositions_;
    PerAtomQuantity<double> * atomProcGhostCoarseGrainingPositions_;
    PerAtomQuantity<double> * atomReferencePositions_;

    
    /** number of unique FE nodes */
    int nNodes_;

    /** Number of Spatial Dimensions */
    int nsd_;

#ifdef EXTENDED_ERROR_CHECKING
    /** data for handling atoms crossing processors */
    bool atomSwitch_;
#endif

    /** reference position of the atoms */
    double ** xref_;
    bool readXref_;
    bool needXrefProcessorGhosts_;
    std::string xRefFile_;

    /** flag for tracking displacements or not, depending on physics */
    bool trackDisplacement_;

    /** map from reference positions to element id, pointer is to internal only */
    
    bool needsAtomToElementMap_;
    PerAtomQuantity<int> * atomElement_;
    PerAtomQuantity<int> * atomGhostElement_;

    /* use element sets to define internal and/or ghost regions */
    std::string internalElementSet_;
    
    /** atomic ATC material tag */
    
    
    double Xprd_,Yprd_,Zprd_; // lengths of periodic box in reference frame
    double XY_,YZ_,XZ_;
    double boxXlo_,boxXhi_; // lo/hi bounds of periodic box in reference frame
    double boxYlo_,boxYhi_; // lo/hi bounds of periodic box in reference frame
    double boxZlo_,boxZhi_; // lo/hi bounds of periodic box in reference frame
    // next data members are for consistency with existing ATC_Transfer, but are redundant and do not
    // conform to naming standards, and should be accessible through the mesh
    /** periodicity flags and lengths */
    int periodicity[3];
    double box_bounds[2][3];
    double box_length[3];

    /** pointers to needed atom quantities and transfers */
    FundamentalAtomQuantity * atomMasses_;
    FundamentalAtomQuantity * atomPositions_;
    FundamentalAtomQuantity * atomVelocities_;
    FundamentalAtomQuantity * atomForces_;

    //---------------------------------------------------------------
    /** \name output data */
    //---------------------------------------------------------------
    /*@{*/

  //private:
    bool parallelConsistency_;

    /** base name for output files */
    std::string outputPrefix_;

    /** output flag */ 
    
    bool outputNow_;
    /** output time or step (for lammps compatibility) */
    bool outputTime_;
  
    /** output frequency */
    int outputFrequency_;
  
    /** sample frequency */
    int sampleFrequency_;
  
    /** sample counter */
    int sampleCounter_;

    TAG_FIELDS filteredData_; 

    double peScale_,keScale_;

  //protected:
    /*@}*/

    //---------------------------------------------------------------
    /** \name member data related to compute_scalar() and compute_vector() */
    //---------------------------------------------------------------
    /*@{*/
    int scalarFlag_;              // 0/1 if compute_scalar() function exists
    int vectorFlag_;              // 0/1 if compute_vector() function exists
    int sizeVector_;              // N = size of global vector
    int scalarVectorFreq_;        // frequency compute s/v data is available at
    int sizePerAtomCols_;         // N = size of per atom vector to dump
    
    double **perAtomOutput_;      // per atom data
    double **&perAtomArray_;      // per atom data
    int extScalar_;               // 0/1 if scalar is intensive/extensive
    int extVector_;               // 0/1/-1 if vector is all int/ext/extlist
    int *extList_;                // list of 0/1 int/ext for each vec component
    int thermoEnergyFlag_;        // 0/1 if fix adds to overall energy
    /*@}*/

    //---------------------------------------------------------------
    /** \name fields and necessary data for FEM */
    //---------------------------------------------------------------
    /*@{*/
    std::map<FieldName,int> fieldSizes_;
    FIELDS fields_;
    /*@}*/


    //---------------------------------------------------------------
    /** \name time integration and filtering fields */
    //---------------------------------------------------------------
    /*@{*/
    
    FIELDS dot_fields_; 
    FIELDS ddot_fields_;
    FIELDS dddot_fields_;

    /** Restricted Fields */

    FIELDS nodalAtomicFields_;  // replaces fieldNdFiltered_
    FIELDS nodalAtomicFieldsRoc_; 

    /*@}*/

    //---------------------------------------------------------------
    /** \name quadrature weights */
    //---------------------------------------------------------------
    /*@{*/
    
    DIAG_MAT NodeVolumes_; 
    DIAG_MAN invNodeVolumes_; 
    /** atomic quadrature integration weights (V_\alpha) */
    ProtectedAtomDiagonalMatrix<double> * atomVolume_;
    std::string atomicWeightsFile_;
    bool atomicWeightsWriteFlag_;
    int atomicWeightsWriteFrequency_;
    double atomicVolume_;  // global atomic volume for homogeneous set of atoms
    std::map<int,double> Valpha_;
    AtomicWeightType atomWeightType_;
    /*@}*/

    //---------------------------------------------------------------
    /** \name domain decomposition */
    //---------------------------------------------------------------
    /*@{*/
    DomainDecompositionType domainDecomposition_;
    /*@}*/

    //---------------------------------------------------------------
    /** \name atom data  */
    //---------------------------------------------------------------
    /*@{*/
    /** bitwise comparisons for boundary (ghost) atoms */
    int groupbit_;
    int groupbitGhost_;
    bool needProcGhost_;
    std::string groupTag_;
    std::string groupTagGhost_;
  
    /** number of atoms of correct type,
        ghosts are atoms outside our domain of interest
        boundary are atoms contributing to boundary flux terms */
    /** Number of "internal" atoms on this processor */
    int nLocal_;
    /** Number of atoms on this processor */
    int nLocalTotal_;
    int nLocalGhost_;
    Array<int> internalToAtom_;
    std::map<int,int> atomToInternal_;
    Array<int> ghostToAtom_;
    /*@}*/
    //----------------------------------------------------------------
    /** \name  maps and masks */
    //----------------------------------------------------------------
    /*@{*/

    AtomToElementMapType atomToElementMapType_;
    int atomToElementMapFrequency_;
    int regionID_;
    /*@}*/

    //----------------------------------------------------------------
    /** \name shape function matrices */
    //----------------------------------------------------------------
    /*@{*/
    // sparse matrix where columns correspond to global node numbering
    SPAR_MAN * shpFcn_;
    VectorDependencyManager<SPAR_MAT * > * shpFcnDerivs_;
    /** map from species std::string tag to the species density */
    std::map<std::string,DENS_MAN> taggedDensMan_;

    /** weighted shape function matrices at overlap nodes
        for use with thermostats */
    SPAR_MAN NhatOverlap_;
    /*@}*/

    //----------------------------------------------------------------
    /** \name accumulant matrices */
    //----------------------------------------------------------------
    /*@{*/
    /** compute kernel shape functions on-the-fly w/o storing N_Ia */
    bool mdMassNormalization_;
    bool kernelBased_;
    bool kernelOnTheFly_;
    class KernelFunction * kernelFunction_;
    bool bondOnTheFly_;
    SPAR_MAN* accumulant_;
    SPAR_MAN* accumulantMol_; // KKM add
    SPAR_MAN* accumulantMolGrad_; // KKM add
    SPAR_MAN  kernelAccumulantMol_; // KKM add
    SPAR_MAN  kernelAccumulantMolGrad_; // KKM add
    DIAG_MAN* accumulantWeights_;
    DIAG_MAN* accumulantInverseVolumes_;  
    int accumulantBandwidth_;
    /*@}*/


    //---------------------------------------------------------------
    /** \name restart procedures */
    //---------------------------------------------------------------
    bool useRestart_;
    std::string restartFileName_;

    //---------------------------------------------------------------
    /** \name data specific to node/faceset for global output */
    //---------------------------------------------------------------
    /** group computes : type, group_id -> value */
    std::map< std::pair<std::string, FieldName > , NodesetOperationType> nsetData_;
    std::map< std::pair<std::string,std::string>, FacesetIntegralType  >       fsetData_;
    std::map< std::pair<std::string, FieldName>,ElementsetOperationType > esetData_;


    //---------------------------------------------------------------
    /** \name reference data */ 
    //---------------------------------------------------------------
    bool hasRefPE_;
    bool setRefPE_;
    bool setRefPEvalue_;
    double refPEvalue_;
    bool readRefPE_;
    std::string nodalRefPEfile_;
    DENS_MAN* nodalRefPotentialEnergy_;
    void set_reference_potential_energy(void);

  private: /** data */
    /** current time in simulation */
    double simTime_;

    /** step counter */
    int stepCounter_;

  };

};

#endif
