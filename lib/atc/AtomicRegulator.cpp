// ATC Headers
#include "AtomicRegulator.h"
#include "ATC_Error.h"
#include "ATC_Coupling.h"
#include "PrescribedDataManager.h"
#include "TimeIntegrator.h"
#include "LinearSolver.h"

using std::map;
using std::string;
using std::set;
using std::pair;

namespace ATC {

  
  // only one regulator method at time, i.e. fixed & flux, thermo & elastic
  // regulator manages lambda variables, creates new ones when requested with dimensions and zero ics (map of tag to lambda)
  // regulator keeps track of which lambda are being used, unused lambdas deleted (map of tag to bool), all tags set to unused on start of initialization
  // method requests needed lambda from regulator
  // method sets up all needed linear solvers, null linear solver does nothing
  // regulator adds nodes to fixed or fluxed lists it owns, based on localization and type
  // method gets lists of fixed nodes and fluxed nodes
  // method lumps fluxed lambdas and truncates fixed lambdas based on single localized bool in regulator
  // inherited methods should be fixed, fluxed, combined

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomicRegulator
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomicRegulator::AtomicRegulator(ATC_Coupling * atc,
                                   const string & regulatorPrefix) :
    atc_(atc),
    howOften_(1),
    needReset_(true),
    maxIterations_(myMaxIterations),
    tolerance_(myTolerance),
    regulatorTarget_(NONE),
    couplingMode_(UNCOUPLED),
    nNodes_(0),
    nsd_(atc_->nsd()),
    nLocal_(0),
    useLocalizedLambda_(false),
    useLumpedLambda_(false),
    timeFilter_(NULL),
    regulatorMethod_(NULL),
    boundaryIntegrationType_(NO_QUADRATURE),
    regulatorPrefix_(regulatorPrefix)
  {
    applyInDirection_.resize(atc_->nsd(),true);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  AtomicRegulator::~AtomicRegulator()
  {
    delete_method();
    set_all_data_to_unused();
    delete_unused_data();
  }

  //--------------------------------------------------------
  //  delete_method:
  //    deletes the method
  //--------------------------------------------------------
  void AtomicRegulator::delete_method()
  {
    if (regulatorMethod_)
      delete regulatorMethod_;
  }

  //--------------------------------------------------------
  //  delete_unused_data:
  //    deletes all data that is currently not in use
  //--------------------------------------------------------
  void AtomicRegulator::delete_unused_data()
  {
    map<string, pair<bool,DENS_MAN * > >::iterator it;
    for (it = regulatorData_.begin(); it != regulatorData_.end(); it++) {
      if (((it->second).first)) {
        delete (it->second).second;
        regulatorData_.erase(it);
      }
    }
  }

  //--------------------------------------------------------
  //  get_regulator_data:
  //    gets a pointer to the requested data, is crated if
  //    if doesn't exist
  //--------------------------------------------------------
  DENS_MAN * AtomicRegulator::regulator_data(const string tag, int nCols)
  {
    DENS_MAN * data(NULL);
    map<string, pair<bool,DENS_MAN * > >::iterator it = regulatorData_.find(tag);
    if (it == regulatorData_.end()) {
      data = new DENS_MAN(nNodes_,nCols);
      regulatorData_.insert(pair<string, pair<bool,DENS_MAN * > >(tag,pair<bool,DENS_MAN * >(false,data)));
    }
    else {
      data = (it->second).second;
      if ((data->nRows() != nNodes_) || (data->nCols() != nCols)) {
        data->reset(nNodes_,nCols);
      }
      (it->second).first = false;
    }
    return data;
  }

  //--------------------------------------------------------
  //  get_regulator_data:
  //    gets a pointer to the requested data, or NULL if
  //    if doesn't exist
  //--------------------------------------------------------
  const DENS_MAN * AtomicRegulator::regulator_data(const string tag) const
  {
    map<string, pair<bool,DENS_MAN * > >::const_iterator it = regulatorData_.find(tag);
    if (it == regulatorData_.end()) {
      return NULL;
    }
    else {
      return const_cast<DENS_MAN * >((it->second).second);
    }
  }

  //--------------------------------------------------------
  //  set_all_data_to_unused:
  //    sets bool such that all data is unused
  //--------------------------------------------------------
  void AtomicRegulator::set_all_data_to_unused()
  {
    map<string, pair<bool,DENS_MAN * > >::iterator it;
    for (it = regulatorData_.begin(); it != regulatorData_.end(); it++) {
      (it->second).first = true;
    }
  }

  //--------------------------------------------------------
  //  set_all_data_to_used:
  //    sets bool such that all data is used
  //--------------------------------------------------------
  void AtomicRegulator::set_all_data_to_used()
  {
    map<string, pair<bool,DENS_MAN * > >::iterator it;
    for (it = regulatorData_.begin(); it != regulatorData_.end(); it++) {
      (it->second).first = false;
    }
  }

  //--------------------------------------------------------
  //  modify:
  //    parses and adjusts controller state based on
  //    user input, in the style of LAMMPS user input
  //--------------------------------------------------------
  bool AtomicRegulator::modify(int /* narg */, char **arg)
  {
    bool foundMatch = false;

        // set parameters for numerical matrix solutions
    /*! \page man_control fix_modify AtC control
      \section syntax
      fix_modify AtC control <physics_type> <solution_parameter> <value>\n
        - physics_type (string) = thermal | momentum\n
        - solution_parameter (string) = max_iterations | tolerance\n
      
      fix_modify AtC transfer <physics_type> control max_iterations <max_iterations>\n
        - max_iterations (int) = maximum number of iterations that will be used by iterative matrix solvers\n

      fix_modify AtC transfer <physics_type> control tolerance <tolerance> \n
        - tolerance (float) = relative tolerance to which matrix equations will be solved\n

      \section examples
      <TT> fix_modify AtC control thermal max_iterations 10 </TT> \n
      <TT> fix_modify AtC control momentum tolerance 1.e-5 </TT> \n
      \section description
      Sets the numerical parameters for the matrix solvers used in the specified control algorithm.  Many solution approaches require iterative solvers, and these methods enable users to provide the maximum number of iterations and the relative tolerance.
      \section restrictions
      only for be used with specific controllers :
      thermal, momentum \n
      They are ignored if a lumped solution is requested
      \section related
      \section default
      max_iterations is the number of rows in the matrix\n
      tolerance is 1.e-10
    */
    int argIndex = 0;
    if (strcmp(arg[argIndex],"max_iterations")==0) {
      argIndex++;
      maxIterations_ = atoi(arg[argIndex]);
      if (maxIterations_ < 1) {
        throw ATC_Error("Bad maximum iteration count");
      }
      needReset_ = true;
      foundMatch = true;
    }
    else if (strcmp(arg[argIndex],"tolerance")==0) {
      argIndex++;
      tolerance_ = atof(arg[argIndex]);
      if (tolerance_ < 0.) {
        throw ATC_Error("Bad tolerance value");
      }
      needReset_ = true;
      foundMatch = true;
    }

    /*! \page man_localized_lambda fix_modify AtC control localized_lambda 
      \section syntax
      fix_modify AtC control localized_lambda <on|off> 
      \section examples
       <TT> fix_modify atc control localized_lambda on </TT> \n
      \section description
      Turns on localization algorithms for control algorithms to restrict the influence of FE coupling or boundary conditions to a region near the boundary of the MD region.  Control algorithms will not affect atoms in elements not possessing faces on the boundary of the region.  Flux-based control is localized via row-sum lumping while quantity control is done by solving a truncated matrix equation.
      \section restrictions 
      \section related
      \section default
      Default is off.
    */
    else if (strcmp(arg[argIndex],"localized_lambda")==0) {
      argIndex++;
      if (strcmp(arg[argIndex],"on")==0) {
        useLocalizedLambda_ = true;
        foundMatch = true;
      }
      else if (strcmp(arg[argIndex],"off")==0) {
        useLocalizedLambda_ = false;
        foundMatch = true;
      }
    }

    
    
    /*! \page man_lumped_lambda_solve fix_modify AtC control lumped_lambda_solve 
      \section syntax
      fix_modify AtC control lumped_lambda_solve <on|off> 
      \section examples
       <TT> fix_modify atc control lumped_lambda_solve on </TT> \n
      \section description
      Command to use or not use lumped matrix for lambda solve
      \section restrictions 
      \section related
      \section default
    */
    else if (strcmp(arg[argIndex],"lumped_lambda_solve")==0) {
      argIndex++;
      if (strcmp(arg[argIndex],"on")==0) {
        useLumpedLambda_ = true;
        foundMatch = true;
      }
      else if (strcmp(arg[argIndex],"off")==0) {
        useLumpedLambda_ = false;
        foundMatch = true;
      }
    }

    /*! \page man_mask_direction fix_modify AtC control mask_direction
      \section syntax
      fix_modify AtC control mask_direction <direction> <on|off> 
      \section examples
       <TT> fix_modify atc control mask_direction 0 on </TT> \n
      \section description
      Command to mask out certain dimensions from the atomic regulator
      \section restrictions 
      \section related
      \section default
    */
    else if (strcmp(arg[argIndex],"mask_direction")==0) {
      argIndex++;
      int dir = atoi(arg[argIndex]);
      argIndex++;
      if (strcmp(arg[argIndex],"on")==0) {
        applyInDirection_[dir] = false;
        foundMatch = true;
      }
      else if (strcmp(arg[argIndex],"off")==0) {
        applyInDirection_[dir] = true;
        foundMatch = true;
      }
    }

    return foundMatch;
  }

  //--------------------------------------------------------
  //  reset_nlocal:
  //    resizes lambda force if necessary
  //--------------------------------------------------------
  void AtomicRegulator::reset_nlocal()
  {
    nLocal_ = atc_->nlocal();
    if (regulatorMethod_)
      regulatorMethod_->reset_nlocal();
  }
  //--------------------------------------------------------
  //  reset_atom_materials:
  //    resets the localized atom to material map
  //--------------------------------------------------------
  void AtomicRegulator::reset_atom_materials(const Array<int> & elementToMaterialMap,
                                             const MatrixDependencyManager<DenseMatrix, int> * atomElement)
  {
    if (regulatorMethod_)
      regulatorMethod_->reset_atom_materials(elementToMaterialMap,
                                             atomElement);
  }
  //--------------------------------------------------------
  //  reset_method:
  //    sets up methods, if necessary
  //--------------------------------------------------------
  void AtomicRegulator::reset_method()
  {
    // set up defaults for anything that didn't get set
    if (!regulatorMethod_)
      regulatorMethod_ = new RegulatorMethod(this);
    if (!timeFilter_)
      timeFilter_ = (atc_->time_filter_manager())->construct();
  }
  //--------------------------------------------------------
  //  md_fixed_nodes:
  //    determines if any fixed nodes overlap the MD region
  //--------------------------------------------------------
  bool AtomicRegulator::md_fixed_nodes(FieldName fieldName) const
  {
    FixedNodes fixedNodes(atc_,fieldName);
    const set<int> & myNodes(fixedNodes.quantity());
    if (myNodes.size() == 0) {
      return false;
    }
    else {
      return true;
    }
  }
  //--------------------------------------------------------
  //  md_flux_nodes:
  //    determines if any nodes with fluxes overlap the MD region
  //--------------------------------------------------------
  bool AtomicRegulator::md_flux_nodes(FieldName fieldName) const
  {
    FluxNodes fluxNodes(atc_,fieldName);
    const set<int> & myNodes(fluxNodes.quantity());
    if (myNodes.size() == 0) {
      return false;
    }
    else {
      return true;
    }
  }
  //--------------------------------------------------------
  //  construct_methods:
  //    sets up methods before a run
  //--------------------------------------------------------
  void AtomicRegulator::construct_methods()
  {
    // get base-line data that was set in stages 1 & 2 of ATC_Method::initialize
    // computational geometry
    nNodes_ = atc_->num_nodes();

    // make sure consistent boundary integration is being used
    atc_->set_boundary_integration_type(boundaryIntegrationType_);
  }

  //--------------------------------------------------------
  //  construct_transfers:
  //    pass through to appropriate transfer constuctors
  //--------------------------------------------------------
  void AtomicRegulator::construct_transfers()
  {
    regulatorMethod_->construct_transfers();
  }

  //--------------------------------------------------------
  //  initialize:
  //    sets up methods before a run
  //--------------------------------------------------------
  void AtomicRegulator::initialize()
  {
    regulatorMethod_->initialize();
    needReset_ = false;
  }

  //--------------------------------------------------------
  //  output:
  //    pass through to appropriate output methods
  //--------------------------------------------------------
  void AtomicRegulator::output(OUTPUT_LIST & outputData) const
  {
    regulatorMethod_->output(outputData);
  }

  //--------------------------------------------------------
  //  finish:
  //    pass through to appropriate end-of-run methods
  //--------------------------------------------------------
  void AtomicRegulator::finish()
  {
    regulatorMethod_->finish();
    set_all_data_to_unused();
  }

  //--------------------------------------------------------
  //  apply_pre_predictor:
  //    applies the controller in the pre-predictor
  //    phase of the time integrator  
  //--------------------------------------------------------
  void AtomicRegulator::apply_pre_predictor(double dt, int timeStep)
  {
    if (timeStep % howOften_==0) // apply full integration scheme, including filter
      regulatorMethod_->apply_pre_predictor(dt);
  }

  //--------------------------------------------------------
  //  apply_mid_predictor:
  //    applies the controller in the mid-predictor
  //    phase of the time integrator  
  //--------------------------------------------------------
  void AtomicRegulator::apply_mid_predictor(double dt, int timeStep)
  {
    if (timeStep % howOften_==0) // apply full integration scheme, including filter
      regulatorMethod_->apply_mid_predictor(dt);
  }

  //--------------------------------------------------------
  //  apply_post_predictor:
  //    applies the controller in the post-predictor
  //    phase of the time integrator  
  //--------------------------------------------------------
  void AtomicRegulator::apply_post_predictor(double dt, int timeStep)
  {
    if (timeStep % howOften_==0) // apply full integration scheme, including filter
      regulatorMethod_->apply_post_predictor(dt);
  }
 
  //--------------------------------------------------------
  //  apply_pre_corrector:
  //    applies the controller in the pre-corrector phase
  //    of the time integrator
  //--------------------------------------------------------
  void AtomicRegulator::apply_pre_corrector(double dt, int timeStep)
  {
    if (timeStep % howOften_==0) // apply full integration scheme, including filter   
      regulatorMethod_->apply_pre_corrector(dt);
  }

  //--------------------------------------------------------
  //  apply_post_corrector:
  //    applies the controller in the post-corrector phase
  //    of the time integrator
  //--------------------------------------------------------
  void AtomicRegulator::apply_post_corrector(double dt, int timeStep)
  {
    if (timeStep % howOften_==0) // apply full integration scheme, including filter
      regulatorMethod_->apply_post_corrector(dt);
  }

  //--------------------------------------------------------
  //  pre_exchange
  //--------------------------------------------------------
  void AtomicRegulator::pre_exchange()
  {
    regulatorMethod_->pre_exchange();
  }

  //--------------------------------------------------------
  //  pre_force
  //--------------------------------------------------------
  void AtomicRegulator::pre_force()
  {
    regulatorMethod_->post_exchange();
  }

  //--------------------------------------------------
  // pack_fields
  //   bundle all allocated field matrices into a list
  //   for output needs
  //--------------------------------------------------
  void AtomicRegulator::pack_fields(RESTART_LIST & data)
  {
    map<string, pair<bool,DENS_MAN * > >::iterator it;
    for (it = regulatorData_.begin(); it != regulatorData_.end(); it++) {
      data[(it->first)] = &(((it->second).second)->set_quantity());
    }
  }

  //--------------------------------------------------------
  //  compute_boundary_flux:
  //    computes the boundary flux to be consistent with
  //    the controller
  //--------------------------------------------------------
  void AtomicRegulator::compute_boundary_flux(FIELDS & fields)
  {
    regulatorMethod_->compute_boundary_flux(fields);
  }

  //--------------------------------------------------------
  //  add_to_rhs:
  //    adds any controller contributions to the FE rhs 
  //--------------------------------------------------------
  void AtomicRegulator::add_to_rhs(FIELDS & rhs)
  {
    regulatorMethod_->add_to_rhs(rhs);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class RegulatorMethod
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  RegulatorMethod::RegulatorMethod(AtomicRegulator * atomicRegulator,
                                   const string & regulatorPrefix) :
    atomicRegulator_(atomicRegulator),
    atc_(atomicRegulator_->atc_transfer()),
    boundaryFlux_(atc_->boundary_fluxes()),
    fieldMask_(NUM_FIELDS,NUM_FLUX),
    nNodes_(atomicRegulator_->num_nodes()),
    regulatorPrefix_(atomicRegulator->regulator_prefix()+regulatorPrefix),
    shpFcnDerivs_(NULL)
  {
    fieldMask_ = false;
  }

  //--------------------------------------------------------
  //  compute_boundary_flux
  //    default computation of boundary flux based on
  //    finite
  //--------------------------------------------------------
  void RegulatorMethod::compute_boundary_flux(FIELDS & fields)
  {
    atc_->compute_boundary_flux(fieldMask_,
                                fields,
                                boundaryFlux_,
                                atomMaterialGroups_,
                                shpFcnDerivs_);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class RegulatorShapeFunction
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  RegulatorShapeFunction::RegulatorShapeFunction(AtomicRegulator * atomicRegulator,
                                                 const string & regulatorPrefix) :
    RegulatorMethod(atomicRegulator,regulatorPrefix),
    lambda_(NULL),
    atomLambdas_(NULL),
    shapeFunctionMatrix_(NULL),
    linearSolverType_(AtomicRegulator::NO_SOLVE),
    maxIterations_(atomicRegulator->max_iterations()),
    tolerance_(atomicRegulator->tolerance()),
    matrixSolver_(NULL),
    regulatedNodes_(NULL),
    applicationNodes_(NULL),
    boundaryNodes_(NULL),
    shpFcn_(NULL),
    atomicWeights_(NULL),
    elementMask_(NULL),
    lambdaAtomMap_(NULL),
    weights_(NULL),
    nsd_(atomicRegulator_->nsd()),
    nLocal_(atomicRegulator_->nlocal())
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  RegulatorShapeFunction::~RegulatorShapeFunction()
  {
    if (matrixSolver_)
      delete matrixSolver_;
  }

  //--------------------------------------------------------
  //  create_node_maps
  //  - creates the node mappings between all nodes and the
  //    subset which are regulated
  //--------------------------------------------------------
  void RegulatorShapeFunction::create_node_maps()
  {
    this->construct_regulated_nodes();

    InterscaleManager & interscaleManager(atc_->interscale_manager());
    nodeToOverlapMap_ = static_cast<NodeToSubset * >(interscaleManager.dense_matrix_int(regulatorPrefix_+"NodeToOverlapMap"));
    if (!nodeToOverlapMap_) {
      nodeToOverlapMap_ = new NodeToSubset(atc_,regulatedNodes_);
      interscaleManager.add_dense_matrix_int(nodeToOverlapMap_,
                                             regulatorPrefix_+"NodeToOverlapMap");
    }
    overlapToNodeMap_ = static_cast<SubsetToNode * >(interscaleManager.dense_matrix_int(regulatorPrefix_+"OverlapToNodeMap"));
    if (!overlapToNodeMap_) {
      overlapToNodeMap_ = new SubsetToNode(nodeToOverlapMap_);
      interscaleManager.add_dense_matrix_int(overlapToNodeMap_,
                                             regulatorPrefix_+"OverlapToNodeMap");
    }
    
  }

  //--------------------------------------------------------
  //  construct_transfers
  //  - create all the needed transfer operators, in this
  //    case weights for the lambda matrix
  //--------------------------------------------------------
  void RegulatorShapeFunction::construct_transfers()
  {
    this->set_weights(); // construct specific weighting matrix transfer

    // specialized quantities for boundary flux integration if the lambda atom map exists
    if (lambdaAtomMap_ && (atomicRegulator_->boundary_integration_type() == FE_INTERPOLATION)) {
      InterscaleManager & interscaleManager(atc_->interscale_manager());

      // atomic weights
      PerAtomDiagonalMatrix<double> * atomWeights(interscaleManager.per_atom_diagonal_matrix("AtomVolume"));
      atomicWeights_ = new MappedDiagonalMatrix(atc_,
                                                atomWeights,
                                                lambdaAtomMap_);
      interscaleManager.add_diagonal_matrix(atomicWeights_,
                                            regulatorPrefix_+"RegulatorAtomWeights");

      // shape function
      shpFcn_ = new RowMappedSparseMatrix(atc_,
                                          interscaleManager.per_atom_sparse_matrix("Interpolant"),
                                          lambdaAtomMap_);
      interscaleManager.add_sparse_matrix(shpFcn_,
                                          regulatorPrefix_+"RegulatorShapeFunction");

      // shape function derivatives
      VectorDependencyManager<SPAR_MAT * > * interpolantGradient = interscaleManager.vector_sparse_matrix("InterpolantGradient");
      if (!interpolantGradient) {
        interpolantGradient = new PerAtomShapeFunctionGradient(atc_);
        interscaleManager.add_vector_sparse_matrix(interpolantGradient,
                                                   "InterpolantGradient");
      }
      shpFcnDerivs_ = new RowMappedSparseMatrixVector(atc_,
                                                      interpolantGradient,
                                                      lambdaAtomMap_);
      interscaleManager.add_vector_sparse_matrix(shpFcnDerivs_,
                                                 regulatorPrefix_+"RegulatorShapeFunctionGradient");
    }
  }

  //--------------------------------------------------------
  //  initialize
  //  - pre-run work, in this cases constructs the linear
  //    solver
  //--------------------------------------------------------
  void RegulatorShapeFunction::initialize()
  {
    if (!shapeFunctionMatrix_) {
      throw ATC_Error("RegulatorShapeFunction::initialize - shapeFunctionMatrix_ must be created before the initialize phase");
    }
    if (matrixSolver_)
      delete matrixSolver_;

    if (linearSolverType_ == AtomicRegulator::RSL_SOLVE) {
      matrixSolver_ = new LambdaMatrixSolverLumped(matrixTemplate_,
                                                   shapeFunctionMatrix_,
                                                   maxIterations_,
                                                   tolerance_,
                                                   applicationNodes_,
                                                   nodeToOverlapMap_);
    }
    else if (linearSolverType_ == AtomicRegulator::CG_SOLVE) {
      matrixSolver_ = new LambdaMatrixSolverCg(matrixTemplate_,
                                               shapeFunctionMatrix_,
                                               maxIterations_,
                                               tolerance_);
    }
    else {
      throw ATC_Error("RegulatorShapeFunction::initialize - unsupported solver type");
    }

    compute_sparsity();
  }

  //--------------------------------------------------------
  //  compute_sparsity
  //  - creates sparsity template
  //--------------------------------------------------------
  void RegulatorShapeFunction::compute_sparsity(void)
  {
    
    // first get local pattern from N N^T
    int nNodeOverlap = nodeToOverlapMap_->size();
    DENS_MAT tmpLocal(nNodeOverlap,nNodeOverlap);
    DENS_MAT tmp(nNodeOverlap,nNodeOverlap);
    const SPAR_MAT & myShapeFunctionMatrix(shapeFunctionMatrix_->quantity());
    if (myShapeFunctionMatrix.nRows() > 0) {
      tmpLocal = myShapeFunctionMatrix.transMat(myShapeFunctionMatrix);
    }
    
    // second accumulate total pattern across processors
    LammpsInterface::instance()->allsum(tmpLocal.ptr(), tmp.ptr(), tmp.size());
    // third extract non-zero entries & construct sparse template
    SPAR_MAT & myMatrixTemplate(matrixTemplate_.set_quantity());
    myMatrixTemplate.reset(nNodeOverlap,nNodeOverlap);
    for (int i = 0; i < nNodeOverlap; i++) {
      for (int j = 0; j < nNodeOverlap; j++) {
        if (abs(tmp(i,j))>0) {
          myMatrixTemplate.add(i,j,0.);
        }
      }
    }
    myMatrixTemplate.compress();
  }

  //--------------------------------------------------------
  //  solve_for_lambda
  //    solves matrix equation for lambda using given rhs
  //--------------------------------------------------------
  void RegulatorShapeFunction::solve_for_lambda(const DENS_MAT & rhs,
                                                DENS_MAT & lambda)
  {

    // assemble N^T W N with appropriate weighting matrix
    
    DIAG_MAT weights;
    if (shapeFunctionMatrix_->nRows() > 0) {
      weights.reset(weights_->quantity());
    }
    matrixSolver_->assemble_matrix(weights);
    
    // solve on overlap nodes
    int nNodeOverlap = nodeToOverlapMap_->size();
    DENS_MAT rhsOverlap(nNodeOverlap,rhs.nCols());
    map_unique_to_overlap(rhs, rhsOverlap);
    DENS_MAT lambdaOverlap(nNodeOverlap,lambda.nCols());

    for (int i = 0; i < rhs.nCols(); i++) {
      CLON_VEC tempLambda(lambdaOverlap,CLONE_COL,i);
      if (atomicRegulator_->apply_in_direction(i)) {
        CLON_VEC tempRHS(rhsOverlap,CLONE_COL,i);
        matrixSolver_->execute(tempRHS,tempLambda);
      }
      else {
        tempLambda = 0.;
      }
    }
    
    // map solution back to all nodes
    map_overlap_to_unique(lambdaOverlap,lambda);
  }

  //--------------------------------------------------------
  //  reset_nlocal:
  //    resets data dependent on local atom count
  //--------------------------------------------------------
  void RegulatorShapeFunction::reset_nlocal()
  {
    RegulatorMethod::reset_nlocal();
    nLocal_ = atomicRegulator_->nlocal();

    
    
    //compute_sparsity();
  }

  //--------------------------------------------------------
  //  reset_atom_materials:
  //    resets the localized atom to material map
  //--------------------------------------------------------
  void RegulatorShapeFunction::reset_atom_materials(const Array<int> & elementToMaterialMap,
                                                    const MatrixDependencyManager<DenseMatrix, int> * atomElement)
  {
    // specialized quantities for boundary flux integration if the lambda atom map exists
    if (lambdaAtomMap_ && (atomicRegulator_->boundary_integration_type() == FE_INTERPOLATION)) {
      int nMaterials = (atc_->physics_model())->nMaterials();
      atomMaterialGroups_.reset(nMaterials);
      const INT_ARRAY & atomToElementMap(atomElement->quantity());
      const INT_ARRAY & map(lambdaAtomMap_->quantity());
      int idx;
      for (int i = 0; i < nLocal_; i++) {
        idx = map(i,0);
        if (idx > -1) {
          atomMaterialGroups_(elementToMaterialMap(atomToElementMap(i,0))).insert(idx);
        }
      }
    }
  }

  //--------------------------------------------------------
  //  map_unique_to_overlap:
  //    maps unique node data to overlap node data
  //--------------------------------------------------------
  void RegulatorShapeFunction::map_unique_to_overlap(const MATRIX & uniqueData,
                                                     MATRIX & overlapData)
  {
    const INT_ARRAY & nodeToOverlapMap(nodeToOverlapMap_->quantity());
    for (int i = 0; i < nNodes_; i++) {
      if (nodeToOverlapMap(i,0) > -1) {
        for (int j = 0; j < uniqueData.nCols(); j++) {
          overlapData(nodeToOverlapMap(i,0),j) = uniqueData(i,j);
        }
      }
    }
  }

  //--------------------------------------------------------
  //  map_overlap_to_unique:
  //    maps overlap node data to unique node data
  //--------------------------------------------------------
  void RegulatorShapeFunction::map_overlap_to_unique(const MATRIX & overlapData,
                                                     MATRIX & uniqueData)
  {
    const INT_ARRAY & overlapToNodeMap(overlapToNodeMap_->quantity());
    uniqueData.resize(nNodes_,overlapData.nCols());
    for (int i = 0; i < overlapToNodeMap.size(); i++) {
      for (int j = 0; j < overlapData.nCols(); j++) {
        uniqueData(overlapToNodeMap(i,0),j) = overlapData(i,j);
      }
    }
  }
  //--------------------------------------------------------
  //  construct_regulated_nodes:
  //    constructs the set of nodes being regulated
  //--------------------------------------------------------
  void RegulatorShapeFunction::construct_regulated_nodes()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());
    regulatedNodes_ = interscaleManager.set_int("RegulatedNodes");

    if (!regulatedNodes_) {
      if (!(atomicRegulator_->use_localized_lambda())) {
        regulatedNodes_ = new RegulatedNodes(atc_);
      }
      else {
        regulatedNodes_ = new AllRegulatedNodes(atc_);
      }
      interscaleManager.add_set_int(regulatedNodes_,
                                    regulatorPrefix_+"RegulatedNodes");
    }

    // application and regulated are same, unless specified
    applicationNodes_ = regulatedNodes_;
    // boundary and regulated nodes are same, unless specified
    boundaryNodes_ = regulatedNodes_;

    // special set of boundary elements
    if (atomicRegulator_->use_localized_lambda()) {
      elementMask_ = interscaleManager.dense_matrix_bool(regulatorPrefix_+"BoundaryElementMask");
      if (!elementMask_) {
        elementMask_ = new ElementMaskNodeSet(atc_,boundaryNodes_);
        interscaleManager.add_dense_matrix_bool(elementMask_,
                                                regulatorPrefix_+"BoundaryElementMask");
      }
    }
  }

  //--------------------------------------------------------
  //  compute_boundary_flux
  //    default computation of boundary flux based on
  //    finite
  //--------------------------------------------------------
  void RegulatorShapeFunction::compute_boundary_flux(FIELDS & fields)
  {
    atc_->compute_boundary_flux(fieldMask_,
                                fields,
                                boundaryFlux_,
                                atomMaterialGroups_,
                                shpFcnDerivs_,
                                shpFcn_,
                                atomicWeights_,
                                elementMask_,
                                boundaryNodes_);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class LambdaMatrixSolver
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //         Grab references to necessary data
  //--------------------------------------------------------
  LambdaMatrixSolver::LambdaMatrixSolver(SPAR_MAN & matrixTemplate, SPAR_MAN * shapeFunctionMatrix, int maxIterations, double tolerance) :
    matrixTemplate_(matrixTemplate),
    shapeFunctionMatrix_(shapeFunctionMatrix),
    maxIterations_(maxIterations),
    tolerance_(tolerance)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  assemble_matrix
  //        Assemble the matrix using the shape function
  //        matrices and weights.  This improves efficiency
  //        when multiple solves or iterations are required.
  //--------------------------------------------------------
  void LambdaMatrixSolver::assemble_matrix(DIAG_MAT & weights)
  {
    // form matrix : sum_a N_Ia * W_a * N_Ja
    
    SPAR_MAT lambdaMatrixLocal(matrixTemplate_.quantity());
    if (weights.nRows()>0)
      lambdaMatrixLocal.weighted_least_squares(shapeFunctionMatrix_->quantity(),weights);

    // swap contributions
    lambdaMatrix_ = matrixTemplate_.quantity();
    LammpsInterface::instance()->allsum(lambdaMatrixLocal.ptr(),
                                        lambdaMatrix_.ptr(), lambdaMatrix_.size());
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class LambdaMatrixSolverLumped
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //         Grab references to necessary data
  //--------------------------------------------------------
  LambdaMatrixSolverLumped::LambdaMatrixSolverLumped(SPAR_MAN & matrixTemplate, SPAR_MAN * shapeFunctionMatrix, int maxIterations, double tolerance, const SetDependencyManager<int> * applicationNodes, const NodeToSubset * nodeToOverlapMap) :
    LambdaMatrixSolver(matrixTemplate,shapeFunctionMatrix,maxIterations,tolerance),
    applicationNodes_(applicationNodes),
    nodeToOverlapMap_(nodeToOverlapMap)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  assemble_matrix
  //        Assemble the matrix using the shape function
  //        matrices and weights.  This improves efficiency
  //        when multiple solves or iterations are required.
  //--------------------------------------------------------
  void LambdaMatrixSolverLumped::assemble_matrix(DIAG_MAT & weights)
  {
    LambdaMatrixSolver::assemble_matrix(weights);
    
    lumpedMatrix_ = lambdaMatrix_.row_sum_lump();
  }

  void LambdaMatrixSolverLumped::execute(VECTOR & rhs, VECTOR & lambda) 
  {
    
    // solve lumped equation
    const set<int> & applicationNodes(applicationNodes_->quantity());
    const INT_ARRAY & nodeToOverlapMap(nodeToOverlapMap_->quantity());
    lambda = 0.;
    set<int>::const_iterator iset;
    for (iset = applicationNodes.begin(); iset != applicationNodes.end(); iset++) {
      int node = nodeToOverlapMap(*iset,0);
      lambda(node) = rhs(node)/lumpedMatrix_(node,node);
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class LambdaMatrixSolverCg
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //         Grab references to necessary data
  //--------------------------------------------------------
  LambdaMatrixSolverCg::LambdaMatrixSolverCg(SPAR_MAN & matrixTemplate, SPAR_MAN * shapeFunctionMatrix, int maxIterations, double tolerance) :
    LambdaMatrixSolver(matrixTemplate,shapeFunctionMatrix,maxIterations,tolerance)
  {
    // do nothing
  }

  void LambdaMatrixSolverCg::execute(VECTOR & rhs, VECTOR & lambda)
  {
    if (lambdaMatrix_.size()<1)
      throw ATC_Error("solver given zero size matrix in LambdaMatrixSolverCg::execute()");


    LinearSolver solver(lambdaMatrix_, ATC::LinearSolver::ITERATIVE_SOLVE_SYMMETRIC, true);
    int myMaxIt = maxIterations_ > 0 ? maxIterations_ : 2*lambdaMatrix_.nRows();
    solver.set_max_iterations(myMaxIt);
    solver.set_tolerance(tolerance_);
    solver.solve(lambda,rhs);
  }
};
