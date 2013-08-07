/** Atomic Regulator : a base class class for atom-continuum control */

#ifndef ATOMICREGULATOR_H
#define ATOMICREGULATOR_H

// ATC headers
#include "ATC_TypeDefs.h"

// other headers
#include <map>
#include <set>
#include <vector>

namespace ATC {

  static const int myMaxIterations = 0;
  static const double myTolerance = 1.e-10;

  // forward declarations
  class TimeFilter;
  class RegulatorMethod;
  class LambdaMatrixSolver;
  class ATC_Coupling;
  class NodeToSubset;
  class SubsetToNode;
  class RegulatedNodes;
  class ElementMaskNodeSet;
  class LargeToSmallAtomMap;
  template <typename T>
    class PerAtomQuantity;
  template <typename T>
    class ProtectedAtomQuantity;
  template <typename T>
    class PerAtomSparseMatrix;

  /**
   *  @class  AtomicRegulator
   *  @brief  Base class for atom-continuum control
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomicRegulator
  //--------------------------------------------------------
  //--------------------------------------------------------

  class AtomicRegulator {
  
  public:

    /** linear solver types */
    enum LinearSolverType {
      NO_SOLVE=0,
      CG_SOLVE, // conjugate gradient
      RSL_SOLVE   // row-sum lumping solution
    };

    /** regulator target variable */
    enum RegulatorTargetType {
      NONE=0,
      FIELD,
      DERIVATIVE,
      DYNAMICS
    };

    enum RegulatorCouplingType {
      UNCOUPLED=0,
      FLUX,
      GHOST_FLUX,
      FIXED
    };
  
    // constructor
    AtomicRegulator(ATC_Coupling * atc,
                    const string & regulatorPrefix = "");
        
    // destructor
    virtual ~AtomicRegulator();
        
    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** instantiate up the desired method(s) */
    virtual void construct_methods() = 0;

    /** method(s) create all necessary transfer operators */
    virtual void construct_transfers();

    /** initialization of method data */
    virtual void initialize();
        
    /** add output information */
    virtual void output(OUTPUT_LIST & outputData) const;
    virtual double compute_vector(int n) const {return 0;}
    
    /** final work at the end of a run */
    virtual void finish();

    /** reset number of local atoms, as well as atomic data */
    virtual void reset_nlocal();

    /** set up atom to material identification */
    virtual void reset_atom_materials(const Array<int> & elementToMaterialMap,
                                      const MatrixDependencyManager<DenseMatrix, int> * atomElement);
        
    // application steps
    /** apply the regulator in the pre-predictor phase */
    virtual void apply_pre_predictor(double dt, int timeStep);
    /** apply the regulator in the mid-predictor phase */
    virtual void apply_mid_predictor(double dt, int timeStep);
    /** apply the regulator in the post-predictor phase */
    virtual void apply_post_predictor(double dt, int timeStep);
    /** apply the regulator in the pre-correction phase */
    virtual void apply_pre_corrector(double dt, int timeStep);
    /** apply the regulator in the post-correction phase */
    virtual void apply_post_corrector(double dt, int timeStep);

    /** prior to exchanges */
    virtual void pre_force();
    /** prior to exchanges */
    virtual void pre_exchange();

    /** pack fields for restart */
    virtual void pack_fields(RESTART_LIST & data);

    /** thermo output */
    virtual int size_vector(int s) const {return 0;};

    // coupling to FE state
    /** FE state variable regulator is applied to */
    virtual RegulatorTargetType regulator_target() const {return regulatorTarget_;};
    /** type of boundary coupling */
    //TEMP_JAT field variable should be removed
    virtual RegulatorCouplingType coupling_mode(const FieldName field=NUM_TOTAL_FIELDS) const {return couplingMode_;};
    /** compute the thermal boundary flux, must be consistent with regulator */
    virtual void compute_boundary_flux(FIELDS & fields);
    /** add contributions (if any) to the finite element right-hand side */
    virtual void add_to_rhs(FIELDS & rhs);
        
    // data access, intended for method objects
    /** returns a pointer to the DENS_MAN associated with the tag, creates a new data member if necessary */
    DENS_MAN * regulator_data(const string tag, int nCols);
    /** can externally set regulator dynamic contributions */
    virtual void reset_lambda_contribution(const DENS_MAT & target, const FieldName field=NUM_TOTAL_FIELDS) {};
    /** returns a const pointer to the DENS_MAN associated with the tag, or NULL */
    const DENS_MAN * regulator_data(const string tag) const;
    /** return the maximum number of iterations */
    int max_iterations() {return maxIterations_;};
    /** return the solver tolerance */
    double tolerance() {return tolerance_;};
    /** access for ATC transfer */
    ATC_Coupling * atc_transfer() {return atc_;};
    /** access for time filter */
    TimeFilter * time_filter() {return timeFilter_;};
    /** access for number of nodes */
    int num_nodes() {return nNodes_;};
    /** access for number of spatial dimensions */
    int nsd() {return nsd_;};
    /** access for number of local atoms */
    int nlocal() {return nLocal_;}; 
    /** access for boundary integration methods */
    BoundaryIntegrationType boundary_integration_type()
      {return boundaryIntegrationType_;};
    /** access for boundary face sets */
    const set< pair<int,int> > * face_sets()
      { return boundaryFaceSet_;};
    /** access for needing a reset */
    bool need_reset() const {return needReset_;};
    /** force a reset to occur */
    void force_reset() {needReset_ = true;};
    /** check if lambda is localized */
    bool use_localized_lambda() const {return useLocalizedLambda_;};
    /** check if matrix should be lumpted for lambda solve */
    bool use_lumped_lambda_solve() const {return useLumpedLambda_;};
    /** check to see if this direction is being used */
    bool apply_in_direction(int i) const {return applyInDirection_[i];};

    
    /** checks if there are any fixed nodes in the MD region */
    bool md_fixed_nodes(FieldName fieldName = NUM_TOTAL_FIELDS) const;

    /** checks if there are any flux nodes in the MD region */
    bool md_flux_nodes(FieldName fieldName = NUM_TOTAL_FIELDS) const;

    /** returns prefix tag for regulator */
    const string & regulator_prefix() const {return regulatorPrefix_;};
        
  protected:

    // methods
    /** deletes the current regulator method */
    void delete_method();

    /** deletes all unused data */
    void delete_unused_data();

    /** sets all data to be unused */
    void set_all_data_to_unused();

    /** sets all data to be used */
    void set_all_data_to_used();
        
    // data
    /** point to atc_transfer object */
    ATC_Coupling * atc_;
        
    /** how often in number of time steps regulator is applied */
    int howOften_;

    // reset/reinitialize flags
    /** flag to reset data */
    bool needReset_;
    /** reinitialize method */
    void reset_method();

    // regulator data
    /** container for all data, string is tag, bool is false if currently in use */
    map<string, pair<bool,DENS_MAN * > > regulatorData_;
    /** maximum number of iterations used in solving for lambda */
    int maxIterations_;
    /** tolerance used in solving for lambda */
    double tolerance_;

    /** regulator target flag */
    RegulatorTargetType regulatorTarget_;
    /** regulator fe coupling type flag */
    RegulatorCouplingType couplingMode_;
                
    /** number of nodes */
    int nNodes_;
    /** number of spatial dimensions */
    int nsd_;
    /** number of local atoms */
    int nLocal_;

    /** use of localization techniques */
    bool useLocalizedLambda_;
    bool useLumpedLambda_;

    /** restrict application in certain directions */
    vector<bool> applyInDirection_;
        
    // method pointers
    /** time filtering object */
    TimeFilter * timeFilter_;
    /** sets up and solves the regulator equations */
    RegulatorMethod * regulatorMethod_;

    // boundary flux information
    BoundaryIntegrationType boundaryIntegrationType_;
    const set< pair<int,int> > * boundaryFaceSet_;

    /** prefix string for registering data */
    const string regulatorPrefix_;

  private:
    
    // DO NOT define this
    AtomicRegulator();
        
  };

  /**
   *  @class  RegulatorMethod
   *  @brief  Base class for implementation of control algorithms
   */
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class RegulatorMethod
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class RegulatorMethod {
  
  public:
  
    RegulatorMethod(AtomicRegulator * atomicRegulator,
                    const string & regulatorPrefix = "");
        
    virtual ~RegulatorMethod(){};

    /** instantiate all needed data */
    virtual void construct_transfers(){};

    /** pre-"run" initialization */
    virtual void initialize(){};

    /** reset number of local atoms, as well as atomic data */
    virtual void reset_nlocal(){};

    /** set up atom to material identification */
    virtual void reset_atom_materials(const Array<int> & elementToMaterialMap,
                                      const MatrixDependencyManager<DenseMatrix, int> * atomElement){};
        
    /** applies regulator to atoms in the pre-predictor phase */
    virtual void apply_pre_predictor(double dt){};

    /** applies regulator to atoms in the mid-predictor phase */
    virtual void apply_mid_predictor(double dt){};

    /** applies regulator to atoms in the post-predictor phase */
    virtual void apply_post_predictor(double dt){};

    /** applies regulator to atoms in the pre-corrector phase */
    virtual void apply_pre_corrector(double dt){};

    /** applies regulator to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt){};

    /** applies regulator to atoms in the pre-corrector phase */
    virtual void apply_pre_force(double dt){};

    /** applies regulator to atoms in the post-corrector phase */
    virtual void apply_post_force(double dt){};

    /** applies regulator in pre-force phase */
    virtual void pre_force(){};

    /** applies regulator in pre-exchange phase */
    virtual void pre_exchange(){};

    /** applies regulator in post-exchange phase */
    virtual void post_exchange(){};

    /** compute boundary flux, requires regulator input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields);

    /** add contributions (if any) to the finite element right-hand side */
    virtual void add_to_rhs(FIELDS & rhs){};

    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData){};
    virtual double compute_vector(int n) const {return 0;}

    /** final work at the end of a run */
    virtual void finish(){};

    /** pack fields for restart */
    virtual void pack_fields(RESTART_LIST & data){};
        
  protected:

    //data
    /** pointer to atomic regulator object for data */
    AtomicRegulator * atomicRegulator_;

    /** pointer to ATC_transfer object */
    ATC_Coupling * atc_;

    /** boundary flux */
    FIELDS & boundaryFlux_;

    /** field mask for specifying boundary flux */
    Array2D<bool> fieldMask_;

    /** number of nodes */
    int nNodes_;

    /** prefix string for registering data */
    const string regulatorPrefix_;

    /** mapping for atom materials for atomic FE quadrature */
    Array<set<int> > atomMaterialGroups_;

    /** shape function derivative matrices for boundary atoms */
    VectorDependencyManager<SPAR_MAT * > * shpFcnDerivs_;

  private:

    // DO NOT define this
    RegulatorMethod();
  
  };

  /**
   *  @class  RegulatorShapeFunction
   *  @brief  Base class for implementation of regulation algorithms using the shape function matrices
   */
  // DESIGN each regulator handles only one lambda, but solvers and data are added later
  //        add a new function to set the linear solver based on enum CG_SOLVE or RSL_SOLVE and shape function matrix
  //        followed by call to compute sparsity pattern
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class RegulatorShapeFunction
  //    base class for all regulators of general form
  //    of N^T w N lambda = rhs
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class RegulatorShapeFunction : public RegulatorMethod {
  
  public:
  
    RegulatorShapeFunction(AtomicRegulator * atomicRegulator,
                           const string & regulatorPrefix = "");
        
    virtual ~RegulatorShapeFunction();

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-"run" initialization */
    virtual void initialize();

    /** reset number of local atoms, as well as atomic data */
    virtual void reset_nlocal();

    /** set up atom to material identification */
    virtual void reset_atom_materials(const Array<int> & elementToMaterialMap,
                                      const MatrixDependencyManager<DenseMatrix, int> * atomElement);

    /** compute boundary flux, requires regulator input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields);

    /** determine if local shape function matrices are needed */
    virtual bool use_local_shape_functions() const {return false;};

  protected:

    // methods
    /** compute sparsity for matrix */
    void compute_sparsity(void);

    /** solve matrix equation */
    void solve_for_lambda(const DENS_MAT & rhs,
                          DENS_MAT & lambda);

    /** set weighting factor for in matrix Nhat^T * weights * Nhat */
    virtual void set_weights() = 0;

    /** Mapping between unique nodes and nodes overlapping MD region */
    void map_unique_to_overlap(const MATRIX & uniqueData,
                               MATRIX & overlapData);

    /** Mapping between nodes overlapping MD region to unique nodes */
    void map_overlap_to_unique(const MATRIX & overlapData,
                               MATRIX & uniqueData);

    /** sets up the transfer which is the set of nodes being regulated */
    virtual void construct_regulated_nodes();

    /** creates data structure needed for all to regulated node maps */
    virtual void create_node_maps();

    // member data
    /** lambda coupling parameter */
    DENS_MAN * lambda_;

    /** lambda at atomic locations */
    ProtectedAtomQuantity<double> * atomLambdas_;

    /** shape function matrix for use in GLC solve */
    PerAtomSparseMatrix<double> * shapeFunctionMatrix_;

    /** algorithm being used for the linear solver */
    AtomicRegulator::LinearSolverType linearSolverType_;

    /** pre-templated sparsity pattern for N^T * T * N */
    SPAR_MAN matrixTemplate_;

    /** maximum number of iterations used in solving for lambda */
    int maxIterations_;

    /** tolerance used in solving for lambda */
    double tolerance_;

    /** matrix solver object */
    LambdaMatrixSolver * matrixSolver_;

    /** set of nodes used to construct matrix */
    RegulatedNodes * regulatedNodes_;

    /** set of nodes on which lambda is non-zero */
    RegulatedNodes * applicationNodes_;

    /** set of nodes needed for localized boundary quadrature */
    RegulatedNodes * boundaryNodes_;

    /** mapping from all nodes to overlap nodes: -1 is no overlap, otherwise entry is overlap index */
    NodeToSubset * nodeToOverlapMap_;

    /** mapping from overlap nodes to unique nodes */
    SubsetToNode * overlapToNodeMap_;

    /** shape function matrix for boundary atoms */
    SPAR_MAN * shpFcn_;

    /** atomic weights for boundary atoms */
    DIAG_MAN * atomicWeights_;

    /** element mask for boundary elements corresponding to nodeToOverlapMap_ */
    ElementMaskNodeSet * elementMask_;

    /** maps atoms from atc indexing to regulator indexing */
    LargeToSmallAtomMap * lambdaAtomMap_;

    /** weight per-atom transfer */
    PerAtomQuantity<double> * weights_;

    /** number of spatial dimensions */
    int nsd_;

    /** number of ATC internal atoms on this processor */
    int nLocal_;

  private:

    // DO NOT define this
    RegulatorShapeFunction();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class LambdaMatrixSolver
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class LambdaMatrixSolver {
  
  public:
        
    LambdaMatrixSolver(SPAR_MAN & matrixTemplate, SPAR_MAN * shapeFunctionMatrix, int maxIterations, double tolerance);
        
    virtual ~LambdaMatrixSolver(){};

    /** assemble the matrix */
    virtual void assemble_matrix(DIAG_MAT & weights);
        
    /** execute the solver */
    virtual void execute(VECTOR & rhs, VECTOR & lambda)=0;

  protected:

    /** sparse template for the matrix */
    SPAR_MAN & matrixTemplate_;

    /** non-symmetric part of the matrix */
    SPAR_MAN * shapeFunctionMatrix_;

    /** matrix used to solve for lambda */
    SPAR_MAT lambdaMatrix_;
        
    /** maximum number of iterations */
    int maxIterations_;

    /** relative tolerance to solve to */
    double tolerance_;

  private:

    // DO NOT define this
    LambdaMatrixSolver();
  
  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class LambdaMatrixSolverLumped
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class LambdaMatrixSolverLumped : public LambdaMatrixSolver {
  
  public:
        
    LambdaMatrixSolverLumped(SPAR_MAN & matrixTemplate, SPAR_MAN * shapeFunctionMatrix, int maxIterations, double tolerance, const RegulatedNodes * applicationNodes, const NodeToSubset * nodeToOverlapMap);
        
    virtual ~LambdaMatrixSolverLumped(){};

    /** assemble the matrix */
    virtual void assemble_matrix(DIAG_MAT & weights);
        
    /** execute the solver */
    virtual void execute(VECTOR & rhs, VECTOR & lambda); 
        
  protected:

    /** lumped version of the matrix governing lamda */
    DIAG_MAT lumpedMatrix_;

    /** set of regulated nodes */
    const RegulatedNodes * applicationNodes_;

    /** mapping from all nodes to overlap nodes: -1 is no overlap, otherwise entry is overlap index */
    const NodeToSubset * nodeToOverlapMap_;

  private:

    // DO NOT define this
    LambdaMatrixSolverLumped();
  
  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class LambdaMatrixSolverCg
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class LambdaMatrixSolverCg : public LambdaMatrixSolver {
  
  public:
        
    LambdaMatrixSolverCg(SPAR_MAN & matrixTemplate, SPAR_MAN * shapeFunctionMatrix, int maxIterations, double tolerance);
        
    virtual ~LambdaMatrixSolverCg(){};
        
    /** execute the solver */
    virtual void execute(VECTOR & rhs, VECTOR & lambda); 
        
  protected:


  private:

    // DO NOT define this
    LambdaMatrixSolverCg();
  
  };

};

#endif
