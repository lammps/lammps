// ATC_Transfer Headers
#include "AtomicRegulator.h"
#include "CG.h"
#include "ATC_Error.h"
#include "PrescribedDataManager.h"
#include "TimeIntegrator.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomicRegulator
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomicRegulator::AtomicRegulator(ATC_Transfer * atcTransfer) :
    atcTransfer_(atcTransfer),
    howOften_(1),
    timeFilter_(NULL),
    regulatorMethod_(NULL),
    boundaryIntegrationType_(ATC_Transfer::NO_QUADRATURE),
    nNodes_(0),
    nsd_(0),
    nLocal_(0),
    needReset_(true),
    resetData_(true)
  {
    // nothing to do
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  AtomicRegulator::~AtomicRegulator()
  {
    if (timeFilter_)
      delete timeFilter_;
    destroy();
  }

  //--------------------------------------------------------
  //  destroy:
  //    deallocates all memory
  //--------------------------------------------------------
  void AtomicRegulator::destroy()
  {
    if (regulatorMethod_)
      delete regulatorMethod_;
  }

  //--------------------------------------------------------
  //  modify:
  //    parses and adjusts controller state based on
  //    user input, in the style of LAMMPS user input
  //--------------------------------------------------------
  bool AtomicRegulator::modify(int narg, char **arg)
  {
    bool foundMatch = false;

    return foundMatch;
  }

  //--------------------------------------------------------
  //  reset_nlocal:
  //    resizes lambda force if necessary
  //--------------------------------------------------------
  void AtomicRegulator::reset_nlocal()
  {
    nLocal_ = atcTransfer_->get_nlocal();
    if (nLocal_ > 0)
      lambdaForce_.reset(nLocal_,nsd_);
    if (regulatorMethod_)
      regulatorMethod_->reset_nlocal();
  }

  //--------------------------------------------------------
  //  reset_data:
  //    sets up storage for all data structures
  //--------------------------------------------------------
  void AtomicRegulator::reset_data()
  {
    nNodes_ = atcTransfer_->get_nNodes();
    nsd_    = atcTransfer_->get_nsd();

    if (timeFilter_)
        delete timeFilter_;
    timeFilter_ = NULL;

    resetData_ = false;
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
      timeFilter_ = (atcTransfer_->get_time_filter_manager())->construct();
     
    needReset_ = false;
  }
  
  //--------------------------------------------------------
  //  initialize:
  //    sets up methods before a run
  //--------------------------------------------------------
  void AtomicRegulator::initialize()
  {
    // make sure consistent boundary integration is being used
    atcTransfer_->set_boundary_integration_type(boundaryIntegrationType_);

    // reset data related to local atom count
    reset_nlocal();
  }

  //--------------------------------------------------------
  //  output:
  //    pass through to appropriate output methods
  //--------------------------------------------------------
  void AtomicRegulator::output(double dt, OUTPUT_LIST & outputData) const
  {
    regulatorMethod_->output(dt,outputData);
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
  RegulatorMethod::RegulatorMethod(AtomicRegulator * atomicRegulator) :
    atomicRegulator_(atomicRegulator),
    atcTransfer_(atomicRegulator->get_atc_transfer()),
    fieldMask_(NUM_FIELDS,NUM_FLUX),
    boundaryFlux_(atcTransfer_->get_boundary_fluxes()),
    nNodes_(atomicRegulator_->get_nNodes())
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
    atcTransfer_->compute_boundary_flux(fieldMask_,
                                        fields,
                                        boundaryFlux_);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class RegulatorShapeFunction
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  RegulatorShapeFunction::RegulatorShapeFunction(AtomicRegulator * atomicRegulator) :
    RegulatorMethod(atomicRegulator),
    maxIterations_(50),
    tolerance_(1.e-10),
    nNodeOverlap_(atcTransfer_->get_nNode_overlap()),
    nsd_(atomicRegulator_->get_nsd()),
    lambda_(atomicRegulator_->get_lambda()),
    shapeFunctionMatrix_(atcTransfer_->get_nhat_overlap()),
    glcMatrixTemplate_(atcTransfer_->get_m_t_template()),
    shapeFunctionGhost_(atcTransfer_->get_shape_function_ghost_overlap()),
    internalToAtom_(atcTransfer_->get_internal_to_atom_map()),
    internalToOverlapMap_(atcTransfer_->get_atom_to_overlap_map()),
    ghostToAtom_(atcTransfer_->get_ghost_to_atom_map()),
    nLocal_(0),
    nLocalLambda_(0),
    nLocalGhost_(0)
  {
    if (atcTransfer_->use_lumped_lambda_solve())
      matrixSolver_ = new LambdaMatrixSolverLumped(glcMatrixTemplate_,
                                                   shapeFunctionMatrix_,
                                                   maxIterations_,
                                                   tolerance_);
    else
      matrixSolver_ = new LambdaMatrixSolverCg(glcMatrixTemplate_,
                                               shapeFunctionMatrix_,
                                               maxIterations_,
                                               tolerance_);
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
  //  solve_for_lambda
  //    solves matrix equation for lambda using given rhs
  //--------------------------------------------------------
  void RegulatorShapeFunction::solve_for_lambda(const DENS_MAT & rhs)
  {
    // set up weighting matrix
    DIAG_MAT weights;
    if (nLocalLambda_>0)
      set_weights(weights);
    
    // solve on overlap nodes
    DENS_MAT rhsOverlap(nNodeOverlap_,rhs.nCols());
    atcTransfer_->map_unique_to_overlap(rhs, rhsOverlap);
    DENS_MAT lambdaOverlap(nNodeOverlap_,lambda_.nCols());

    for (int i = 0; i < rhs.nCols(); i++) {
      CLON_VEC tempRHS(rhsOverlap,CLONE_COL,i);
      CLON_VEC tempLambda(lambdaOverlap,CLONE_COL,i);
      matrixSolver_->execute(tempRHS,tempLambda,weights,atcTransfer_);
    }
    
    // map solution back to all nodes
    atcTransfer_->map_overlap_to_unique(lambdaOverlap,lambda_);
  }

  //--------------------------------------------------------
  //  reset_nlocal:
  //    resets data dependent on local atom count
  //--------------------------------------------------------
  void RegulatorShapeFunction::reset_nlocal()
  {
    RegulatorMethod::reset_nlocal();
    nLocal_ = atomicRegulator_->get_nLocal();
    nLocalLambda_ = atcTransfer_->get_nlocal_lambda();
    nLocalGhost_ = atcTransfer_->get_nlocal_ghost();
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
  LambdaMatrixSolver::LambdaMatrixSolver(SPAR_MAT & matrixTemplate, SPAR_MAT & shapeFunctionMatrix, int maxIterations, double tolerance) :
    matrixTemplate_(matrixTemplate),
    shapeFunctionMatrix_(shapeFunctionMatrix),
    maxIterations_(maxIterations),
    tolerance_(tolerance)
  {
    // do nothing
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
  LambdaMatrixSolverLumped::LambdaMatrixSolverLumped(SPAR_MAT & matrixTemplate, SPAR_MAT & shapeFunctionMatrix, int maxIterations, double tolerance) :
    LambdaMatrixSolver(matrixTemplate,shapeFunctionMatrix,maxIterations,tolerance)
  {
    // do nothing
  }

  void LambdaMatrixSolverLumped::execute(VECTOR & rhs, VECTOR & lambda, DIAG_MAT & weights, ATC_Transfer * atcTransfer)
  {
    // form matrix : sum_a N_Ia * W_a * N_Ja
    SPAR_MAT myMatrixLocal(matrixTemplate_);
    if (weights.nRows()>0)
      myMatrixLocal.WeightedLeastSquares(shapeFunctionMatrix_,weights);

    // swap contributions
    SPAR_MAT myMatrix(matrixTemplate_);
    LammpsInterface::instance()->allsum(myMatrixLocal.get_ptr(),
                                        myMatrix.get_ptr(), myMatrix.size());

    DIAG_MAT lumpedMatrix(myMatrix.nRows(),myMatrix.nCols());
    for (int i = 0; i < myMatrix.nRows(); i++)
      for (int j = 0; j < myMatrix.nCols(); j++)
        lumpedMatrix(i,i) += myMatrix(i,j);

    // solve lumped equation
    for (int i = 0; i < rhs.size(); i++)
      lambda(i) = rhs(i)/lumpedMatrix(i,i);
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
  LambdaMatrixSolverCg::LambdaMatrixSolverCg(SPAR_MAT & matrixTemplate, SPAR_MAT & shapeFunctionMatrix, int maxIterations, double tolerance) :
    LambdaMatrixSolver(matrixTemplate,shapeFunctionMatrix,maxIterations,tolerance)
  {
    // do nothing
  }

  void LambdaMatrixSolverCg::execute(VECTOR & rhs, VECTOR & lambda, DIAG_MAT & weights, ATC_Transfer * atcTransfer)
  {
    // form matrix : sum_a N_Ia * W_a * N_Ja
    SPAR_MAT myMatrixLocal(matrixTemplate_);
    if (weights.nRows()>0)
      myMatrixLocal.WeightedLeastSquares(shapeFunctionMatrix_,weights);

    // swap contributions
    SPAR_MAT myMatrix(matrixTemplate_);
    LammpsInterface::instance()->allsum(myMatrixLocal.get_ptr(),
                                        myMatrix.get_ptr(), myMatrix.size());
           
    
    DIAG_MAT preConditioner = myMatrix.get_diag();
    int myMaxIt = 2*myMatrix.nRows(); // note could also use the fixed parameter
    double myTol = tolerance_;

    int convergence = CG(myMatrix, lambda, rhs, preConditioner, myMaxIt, myTol);

    // error if didn't converge
    if (convergence>0)
      throw ATC_Error(0,"CG solver did not converge in LambdaMatrixSolverCg::execute()");   
  }

};
