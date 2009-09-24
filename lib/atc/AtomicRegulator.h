/** Atomic Regulator : a base class class for atom-continuum control */

#ifndef ATOMICREGULATOR_H
#define ATOMICREGULATOR_H

// ATC_Transfer headers
#include "ATC_Transfer.h"

// other headers
#include <map>
#include <set>

namespace ATC {

  // forward declarations
  class TimeFilter;
  class RegulatorMethod;
  class LambdaMatrixSolver;

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
  
    // constructor
    AtomicRegulator(ATC_Transfer * atcTransfer);
        
    // destructor
    ~AtomicRegulator();
        
    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** pre time integration */
    virtual void initialize();
        
    /** add output information */
    virtual void output(double dt, OUTPUT_LIST & outputData) const;

    /** reset number of local atoms, as well as atomic data */
    virtual void reset_nlocal();
        
    // application steps
    /** apply the thermostat in the pre-predictor phase */
    virtual void apply_pre_predictor(double dt, int timeStep);
    /** apply the thermostat in the mid-predictor phase */
    virtual void apply_mid_predictor(double dt, int timeStep);
    /** apply the thermostat in the post-predictor phase */
    virtual void apply_post_predictor(double dt, int timeStep);
    /** apply the thermostat in the pre-correction phase */
    virtual void apply_pre_corrector(double dt, int timeStep);
    /** apply the thermostat in the post-correction phase */
    virtual void apply_post_corrector(double dt, int timeStep);

    // coupling to FE state
    /** compute the thermal boundary flux, must be consistent with thermostat */
    virtual void compute_boundary_flux(FIELDS & fields);
    /** add contributions (if any) to the finite element right-hand side */
    virtual void add_to_rhs(FIELDS & rhs);
        
    // data access, intended for method objects
    /** return value of lambda */
    DENS_MAT & get_lambda() {return lambda_;};
    /** return the atomic force defined by lambda */
    DENS_MAT & get_lambda_force() {return lambdaForce_;};
    /** access for ATC transfer */
    ATC_Transfer * get_atc_transfer() {return atcTransfer_;};
    /** access for time filter */
    TimeFilter * get_time_filter() {return timeFilter_;};
    /** access for number of nodes */
    int get_nNodes() {return nNodes_;};
    /** access for number of spatial dimensions */
    int get_nsd() {return nsd_;};
    /** access for number of local atoms */
    int get_nLocal() {return nLocal_;};
    /** access for boundary integration methods */
    ATC_Transfer::BoundaryIntegrationType get_boundary_integration_type()
      {return boundaryIntegrationType_;};
    /** access for boundary face sets */
    const set< pair<int,int> > * get_face_sets()
      { return boundaryFaceSet_;};
        
  protected:

    void destroy();
        
    /** point to atc_transfer object */
    ATC_Transfer * atcTransfer_;
        
    /** how often in number of time steps thermostat is applied */
    int howOften_;

    // reset/reinitialize flags
    /** flag to see if data requires a reset */
    bool resetData_;
    /** flag to reset data */
    bool needReset_;
    /** resets data structures */
    void reset_data();
    /** reinitialize method */
    void reset_method();

    // regulator data
    /** control parameter */
    DENS_MAT lambda_;
    /** lambda force computed by controller */
    DENS_MAT lambdaForce_;
                
    /** number of nodes */
    int nNodes_;
    /** number of spatial dimensions */
    int nsd_;
    /** number of local atoms */
    int nLocal_;
        
    // method pointers
    /** time filtering object */
    TimeFilter * timeFilter_;
    /** sets up and solves the thermostat equations */
    RegulatorMethod * regulatorMethod_;

    // boundary flux information
    ATC_Transfer::BoundaryIntegrationType boundaryIntegrationType_;
    const set< pair<int,int> > * boundaryFaceSet_;

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
  
    RegulatorMethod(AtomicRegulator * atomicRegulator);
        
    ~RegulatorMethod(){};

    /** reset number of local atoms, as well as atomic data */
    virtual void reset_nlocal(){};
        
    /** applies thermostat to atoms in the pre-predictor phase */
    virtual void apply_pre_predictor(double dt){};

    /** applies thermostat to atoms in the mid-predictor phase */
    virtual void apply_mid_predictor(double dt){};

    /** applies thermostat to atoms in the post-predictor phase */
    virtual void apply_post_predictor(double dt){};

    /** applies thermostat to atoms in the pre-corrector phase */
    virtual void apply_pre_corrector(double dt){};

    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt){};

    /** compute boundary flux, requires thermostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields);

    /** add contributions (if any) to the finite element right-hand side */
    virtual void add_to_rhs(FIELDS & rhs){};

    /** get data for output */
    virtual void output(double dt, OUTPUT_LIST & outputData){};
        
  protected:

    /** pointer to ATC_transfer object */
    ATC_Transfer * atcTransfer_;

    /** pointer to atomic regulator object for data */
    AtomicRegulator * atomicRegulator_;

    /** boundary flux */
    FIELDS & boundaryFlux_;

    /** field mask for specifying boundary flux */
    Array2D<bool> fieldMask_;

    /** number of nodes */
    int nNodes_;

  private:

    // DO NOT define this
    RegulatorMethod();
  
  };

  /**
   *  @class  RegulatorShapeFunction
   *  @brief  Base class for implementation of regulation algorithms using the shape function matrices
   */
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class RegulatorShapeFunction
  //    base class for all regulators of general form
  //    of N^T w N lambda = rhs
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class RegulatorShapeFunction : public RegulatorMethod {
  
  public:
  
    RegulatorShapeFunction(AtomicRegulator * atomicRegulator);
        
    ~RegulatorShapeFunction();

    /** reset number of local atoms, as well as atomic data */
    virtual void reset_nlocal();

  protected:

    // methods
    /** solve matrix equation */
    void solve_for_lambda(const DENS_MAT & rhs);

    /** set weighting factor for in matrix Nhat^T * weights * Nhat */
    virtual void set_weights(DIAG_MAT & weights){};

    // member data
    /** lambda coupling parameter */
    DENS_MAT & lambda_;

    /** shape function matrix for use in GLC solve */
    SPAR_MAT & shapeFunctionMatrix_;

    /** pre-templated sparsity pattern for N^T * T * N */
    SPAR_MAT & glcMatrixTemplate_;

    /**  reference to ATC unity shape function on ghost atoms */
    SPAR_MAT & shapeFunctionGhost_;

    /** maximum number of iterations used in solving for lambda */
    int maxIterations_;

    /** tolerance used in solving for lambda */
    double tolerance_;

    /** matrix solver object */
    LambdaMatrixSolver * matrixSolver_;

    /** maps internal atom ids to LAMMPS atom ids */
    Array<int> & internalToAtom_;

    /** maps internal atoms to overlap atoms */
    SPAR_MAT & internalToOverlapMap_;

    /** maps ghost atom and LAMMPS atom ids */
    Array<int> & ghostToAtom_;

    /** number of overlapping nodes */
    int nNodeOverlap_;

    /** number of spatial dimensions */
    int nsd_;

    /** number of ATC internal atoms on this processor */
    int nLocal_;

    /** number of thermostatted ATC internal atoms on this processor */
    int nLocalLambda_;

    /** number of ATC ghost atoms on this processor */
    int nLocalGhost_;

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
        
    LambdaMatrixSolver(SPAR_MAT & matrixTemplate, SPAR_MAT & shapeFunctionMatrix, int maxIterations, double tolerance);
        
    ~LambdaMatrixSolver(){};
        
    /** execute the solver */
    virtual void execute(VECTOR & rhs, VECTOR & lambda, DIAG_MAT & weights,ATC_Transfer * atcTransfer_=NULL) = 0;
        
  protected:

    /** sparse template for the matrix */
    SPAR_MAT & matrixTemplate_;

    /** non-symmetric part of the matrix */
    SPAR_MAT & shapeFunctionMatrix_;
        
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
        
    LambdaMatrixSolverLumped(SPAR_MAT & matrixTemplate, SPAR_MAT & shapeFunctionMatrix, int maxIterations, double tolerance);
        
    ~LambdaMatrixSolverLumped(){};
        
    /** execute the solver */
    virtual void execute(VECTOR & rhs, VECTOR & lambda, DIAG_MAT & weights,ATC_Transfer * atcTransfer_=NULL);
        
  protected:
        

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
        
    LambdaMatrixSolverCg(SPAR_MAT & matrixTemplate, SPAR_MAT & shapeFunctionMatrix, int maxIterations, double tolerance);
        
    ~LambdaMatrixSolverCg(){};
        
    /** execute the solver */
    virtual void execute(VECTOR & rhs, VECTOR & lambda, DIAG_MAT & weights,ATC_Transfer * atcTransfer_=NULL);
        
  protected:
        


  private:

    // DO NOT define this
    LambdaMatrixSolverCg();
  
  };

};

#endif
