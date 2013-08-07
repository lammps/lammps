// ATC Headers
#include "ExtrinsicModelDriftDiffusion.h"
#include "ATC_Error.h"
#include "FieldEulerIntegrator.h"
#include "ATC_Coupling.h"
#include "LammpsInterface.h"
#include "PrescribedDataManager.h"
#include "PhysicsModel.h"
#include "LinearSolver.h"
#include "PoissonSolver.h"
#include "SchrodingerSolver.h"

// timer
#include "Utility.h"

const double tol = 1.e-8; 
const double zero_tol = 1.e-12; 
const double f_tol = 1.e-8; 

namespace ATC {

enum oneDconservationEnum {ONED_DENSITY=0, ONED_FLUX};

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModelDriftDiffusion
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ExtrinsicModelDriftDiffusion::ExtrinsicModelDriftDiffusion
                                (ExtrinsicModelManager * modelManager,
                                 ExtrinsicModelType modelType,
                                 string matFileName) :
    ExtrinsicModelTwoTemperature(modelManager,modelType,matFileName),
    continuityIntegrator_(NULL),
    
    poissonSolverType_(DIRECT), // ITERATIVE | DIRECT
    poissonSolver_(NULL),
    baseSize_(0),
    electronDensityEqn_(ELECTRON_CONTINUITY),
    fluxUpdateFreq_(1),
    schrodingerSolverType_(DIRECT), // ITERATIVE | DIRECT
    schrodingerSolver_(NULL),
    schrodingerPoissonMgr_(),
    schrodingerPoissonSolver_(NULL),
    maxConsistencyIter_(0), maxConstraintIter_(1), 
    safe_dEf_(0.1), Ef_shift_(0.0),
    oneD_(false), oneDcoor_(0), oneDconserve_(ONED_DENSITY)
  {
     // delete base class's version of the physics model
     if (physicsModel_) delete physicsModel_; 
     if (modelType == DRIFT_DIFFUSION_EQUILIBRIUM) {
       physicsModel_ = new PhysicsModelDriftDiffusionEquilibrium(matFileName); 
       electronDensityEqn_ = ELECTRON_EQUILIBRIUM;
     }
     else if (modelType == DRIFT_DIFFUSION_SCHRODINGER) {
       physicsModel_ = new PhysicsModelDriftDiffusionSchrodinger(matFileName); 
       electronDensityEqn_ = ELECTRON_SCHRODINGER;
       maxConsistencyIter_ = 1;
     }
     else if (modelType == DRIFT_DIFFUSION_SCHRODINGER_SLICE) {
       physicsModel_ = new PhysicsModelDriftDiffusionSchrodingerSlice(matFileName); 
       electronDensityEqn_ = ELECTRON_SCHRODINGER;
       maxConsistencyIter_ = 1;
     }
     else {
       physicsModel_ = new PhysicsModelDriftDiffusion(matFileName); 
     }
     atc_->useConsistentMassMatrix_(ELECTRON_DENSITY) = true;
  }
  
  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ExtrinsicModelDriftDiffusion::~ExtrinsicModelDriftDiffusion()
  {
    if(continuityIntegrator_) delete continuityIntegrator_;
    if(poissonSolver_) delete poissonSolver_;
    if(schrodingerSolver_) delete schrodingerSolver_;
    if(schrodingerPoissonSolver_) delete schrodingerPoissonSolver_;
  }

  //--------------------------------------------------------
  // modify 
  //--------------------------------------------------------
  bool ExtrinsicModelDriftDiffusion::modify(int narg, char **arg)
  {
    bool match = false;
    int argIndx = 0;
    if (!match) {
      match = ExtrinsicModelTwoTemperature::modify(narg, arg);
    }
    return match;
  }

  //--------------------------------------------------------
  //  initialize
  //--------------------------------------------------------
  void ExtrinsicModelDriftDiffusion::initialize()
  {
    // xTTM sets rhsMaskIntrinsic_
    ExtrinsicModelTwoTemperature::initialize();
    nNodes_ = atc_->num_nodes();
    rhs_[ELECTRON_DENSITY].reset(nNodes_,1);
    rhs_[ELECTRIC_POTENTIAL].reset(nNodes_,1);
      
    // set up electron continuity integrator
    Array2D <bool> rhsMask(NUM_TOTAL_FIELDS,NUM_FLUX); 
    rhsMask = false;
    for (int i = 0; i < NUM_FLUX; i++) {
      rhsMask(ELECTRON_DENSITY,i) = atc_->fieldMask_(ELECTRON_DENSITY,i);
    }
    // need to create the bcs for the solver to configure properly
    atc_->set_fixed_nodes();

    if (continuityIntegrator_) delete continuityIntegrator_;
    if (electronTimeIntegration_ == TimeIntegrator::IMPLICIT) { 
      continuityIntegrator_ = new FieldImplicitEulerIntegrator(ELECTRON_DENSITY,
        physicsModel_, atc_->feEngine_, atc_, rhsMask);
    }
    else {
      continuityIntegrator_ = new FieldExplicitEulerIntegrator(ELECTRON_DENSITY,
        physicsModel_, atc_->feEngine_, atc_, rhsMask);
    }


    atc_->compute_mass_matrix(ELECTRON_DENSITY,physicsModel_);
    //(atc_->consistentMassMats_[ELECTRON_DENSITY].quantity()).print("PHYS MASS MAT");
    //DENS_MAT temp = atc_->consistentMassInverse_ - atc_->consistentMassMatInv_[ELECTRON_DENSITY];
    //temp.print("DIFF In MATS");

    // set up poisson solver
    rhsMask = false;
    for (int i = 0; i < NUM_FLUX; i++) {
      rhsMask(ELECTRIC_POTENTIAL,i) = atc_->fieldMask_(ELECTRIC_POTENTIAL,i);
    }
    int type = ATC::LinearSolver::ITERATIVE_SOLVE_SYMMETRIC;
    if (poissonSolverType_ == DIRECT) {  
      type = ATC::LinearSolver::DIRECT_SOLVE;
    }
    if (poissonSolver_) delete poissonSolver_;
    poissonSolver_ = new PoissonSolver(ELECTRIC_POTENTIAL,
      physicsModel_, atc_->feEngine_, atc_->prescribedDataMgr_, atc_,
      rhsMask,type, true);
    poissonSolver_->initialize();

    // set up schrodinger solver
    if ( electronDensityEqn_ == ELECTRON_SCHRODINGER ) {
      int type = ATC::LinearSolver::ITERATIVE_SOLVE_SYMMETRIC;
      if (schrodingerSolverType_ == DIRECT) {  
        type = ATC::LinearSolver::DIRECT_SOLVE;
      }
      if ( schrodingerSolver_ )  delete schrodingerSolver_;
      if ( oneD_ ) { 
        EfHistory_.reset(oneDslices_.size(),2);
        schrodingerSolver_ = new SliceSchrodingerSolver(ELECTRON_DENSITY,
          physicsModel_, atc_->feEngine_, atc_->prescribedDataMgr_, atc_,
          oneDslices_, type, true);
      } 
      else {
        schrodingerSolver_ = new SchrodingerSolver(ELECTRON_DENSITY,
          physicsModel_, atc_->feEngine_, atc_->prescribedDataMgr_, atc_,
          type, true);
      }
      schrodingerSolver_->initialize();

      if ( schrodingerPoissonSolver_ )  delete schrodingerPoissonSolver_;
      schrodingerPoissonSolver_ = schrodingerPoissonMgr_.initialize(
      atc_, schrodingerSolver_, poissonSolver_, physicsModel_);

    }

    if (electronDensityEqn_ == ELECTRON_SCHRODINGER && !(atc_->is_initialized())) {
      ((atc_->fields())[ELECTRON_WAVEFUNCTION].set_quantity()).reset(nNodes_,1);
      ((atc_->fields())[ELECTRON_WAVEFUNCTIONS].set_quantity()).reset(nNodes_,nNodes_);
      ((atc_->fields())[ELECTRON_WAVEFUNCTION_ENERGIES].set_quantity()).reset(nNodes_,1);
    }

  }

  //--------------------------------------------------------
  //  pre initial integration
  //--------------------------------------------------------
  void ExtrinsicModelDriftDiffusion::pre_init_integrate()
  {
    double dt = atc_->lammpsInterface_->dt();
    double time = atc_->time();
    int step = atc_->step();
    if (step % fluxUpdateFreq_ != 0) return; 

    // set Dirchlet data
    atc_->set_fixed_nodes();

    // set Neumann data (atc does not set these until post_final)
    atc_->set_sources();

    // subcyle integration of fast electron variable/s
    
    double idt = dt/nsubcycle_;
    for (int i = 0; i < nsubcycle_ ; ++i) {
      if (electronDensityEqn_ == ELECTRON_CONTINUITY) {
        // update continuity eqn
        if (! atc_->prescribedDataMgr_->all_fixed(ELECTRON_DENSITY) ) 
          continuityIntegrator_->update(idt,time,atc_->fields_,rhs_);
        atc_->set_fixed_nodes(); 
        // solve poisson eqn for electric potential
        if (! atc_->prescribedDataMgr_->all_fixed(ELECTRIC_POTENTIAL) )
          poissonSolver_->solve(atc_->fields(),rhs_);
      } 
      else if (electronDensityEqn_ == ELECTRON_SCHRODINGER) {
        schrodingerPoissonSolver_->solve(rhs_,fluxes_);
      }
      // update electron temperature
      if (! atc_->prescribedDataMgr_->all_fixed(ELECTRON_TEMPERATURE) ) 
        temperatureIntegrator_->update(idt,time,atc_->fields_,rhs_);
      atc_->set_fixed_nodes(); 
    }

    
  }

  //--------------------------------------------------------
  //  output
  //--------------------------------------------------------
  void ExtrinsicModelDriftDiffusion::output(OUTPUT_LIST & outputData)
  {
    ExtrinsicModelTwoTemperature::output(outputData);
    // fields

    outputData["dot_electron_density"] = & (atc_->dot_field(ELECTRON_DENSITY)).set_quantity();
    outputData["joule_heating"]        = & rhs_[ELECTRON_TEMPERATURE].set_quantity();
    Array2D <bool> rhsMask(NUM_FIELDS,NUM_FLUX); rhsMask = false;
    rhsMask(ELECTRON_DENSITY,FLUX) = true;
    rhsMask(ELECTRIC_POTENTIAL,FLUX) = true;
    atc_->compute_flux(rhsMask,atc_->fields_,fluxes_,physicsModel_);
//(fluxes_[ELECTRON_DENSITY][0]).print("J_x");
    outputData["electron_flux_x"]      = & fluxes_[ELECTRON_DENSITY][0];
    outputData["electron_flux_y"]      = & fluxes_[ELECTRON_DENSITY][1];
    outputData["electron_flux_z"]      = & fluxes_[ELECTRON_DENSITY][2];
    outputData["electric_field_x"]     = & fluxes_[ELECTRIC_POTENTIAL][0];
    outputData["electric_field_y"]     = & fluxes_[ELECTRIC_POTENTIAL][1];
    outputData["electric_field_z"]     = & fluxes_[ELECTRIC_POTENTIAL][2];
    if (electronDensityEqn_ == ELECTRON_SCHRODINGER ) {
      SPAR_MAT K;
      Array2D <bool> rhsMask(NUM_FIELDS,NUM_FLUX);
      rhsMask = false;
      rhsMask(ELECTRON_WAVEFUNCTION,FLUX) = true;
      pair<FieldName,FieldName> row_col(ELECTRON_WAVEFUNCTION,
                                        ELECTRON_WAVEFUNCTION);
      atc_->feEngine_->compute_tangent_matrix(
        rhsMask, row_col, atc_->fields(), physicsModel_,
        atc_->element_to_material_map(), K);
      phiTotal_.reset(K.nRows(),1);
      const DIAG_MAT & inv_dV = (atc_->invNodeVolumes_).quantity();
      for (int i = 0; i < K.nRows() ; i++) {
        phiTotal_(i,0) = 0.0;
        for (int j = 0; j < K.nCols() ; j++) {
          phiTotal_(i,0) += K(i,j);
        }
        phiTotal_(i,0) *= inv_dV(i,i);
      }
      outputData["V_total"]      = & phiTotal_;
    }
    // globals
    double nSum = ((atc_->field(ELECTRON_DENSITY)).quantity()).col_sum();
    atc_->feEngine_->add_global("total_electron_density",nSum);
  }

  //--------------------------------------------------------
  //  size_vector
  //--------------------------------------------------------
  int ExtrinsicModelDriftDiffusion::size_vector(int intrinsicSize)
  {
    int xSize = ExtrinsicModelTwoTemperature::size_vector(intrinsicSize);
    baseSize_ = intrinsicSize  + xSize;
    xSize += 1;
    return xSize;
  }

  //--------------------------------------------------------
  //  compute_vector
  //--------------------------------------------------------
  bool ExtrinsicModelDriftDiffusion::compute_vector(int n, double & value)
  {
    // output[1] = total electron density

    bool match = ExtrinsicModelTwoTemperature::compute_vector(n,value);
    if (match) return match;

    if (n == baseSize_) {
      double nSum = ((atc_->field(ELECTRON_DENSITY)).quantity()).col_sum();
      value = nSum;
      return true;
    }
    return false;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModelDriftDiffusionConvection
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ExtrinsicModelDriftDiffusionConvection::ExtrinsicModelDriftDiffusionConvection
                                (ExtrinsicModelManager * modelManager,
                                 ExtrinsicModelType modelType,
                                 string matFileName) :
    ExtrinsicModelDriftDiffusion(modelManager,modelType,matFileName),
    cddmPoissonSolver_(NULL),
    baseSize_(0)
  {
     // delete base class's version of the physics model
     if (physicsModel_) delete physicsModel_; 
     if (modelType == CONVECTIVE_DRIFT_DIFFUSION_SCHRODINGER) {
       physicsModel_ = new PhysicsModelDriftDiffusionConvectionSchrodinger(matFileName); 
       electronDensityEqn_ = ELECTRON_SCHRODINGER;
     }
     else {
       physicsModel_ = new PhysicsModelDriftDiffusionConvection(matFileName); 
     }
     atc_->useConsistentMassMatrix_(ELECTRON_VELOCITY) = true;
     atc_->useConsistentMassMatrix_(ELECTRON_TEMPERATURE) = true;
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ExtrinsicModelDriftDiffusionConvection::~ExtrinsicModelDriftDiffusionConvection()
  {
    if (cddmPoissonSolver_) delete cddmPoissonSolver_;
    for (vector<LinearSolver * >::const_iterator iter=velocitySolvers_.begin();
                                   iter != velocitySolvers_.end(); iter++)
      if (*iter) delete *iter;
  }

  //--------------------------------------------------------
  //  initialize
  //--------------------------------------------------------
  void ExtrinsicModelDriftDiffusionConvection::initialize()
  {
    ExtrinsicModelDriftDiffusion::initialize();

    // change temperature integrator to be Crank-Nicolson
    if (electronTimeIntegration_ == TimeIntegrator::IMPLICIT) {
      if (temperatureIntegrator_) delete temperatureIntegrator_;
      Array2D <bool> rhsMask(NUM_FIELDS,NUM_FLUX);
      rhsMask = false;
      for (int i = 0; i < NUM_FLUX; i++) {
        rhsMask(ELECTRON_TEMPERATURE,i) = atc_->fieldMask_(ELECTRON_TEMPERATURE,i);
      }
      temperatureIntegrator_ = new FieldImplicitEulerIntegrator(ELECTRON_TEMPERATURE,
                                                                physicsModel_,
                                                                atc_->feEngine_, atc_, 
                                                                rhsMask);
    }

    nNodes_ = atc_->num_nodes();
    nsd_ = atc_->nsd();
    rhs_[ELECTRON_VELOCITY].reset(nNodes_,nsd_);

    
    atc_->set_fixed_nodes(); // needed to correctly set BC data
    // initialize Poisson solver
    if (cddmPoissonSolver_) delete cddmPoissonSolver_;
    Array2D <bool> rhsMask(NUM_FIELDS,NUM_FLUX);
    rhsMask = false;
    rhsMask(ELECTRIC_POTENTIAL,FLUX) = true;
    pair<FieldName,FieldName> row_col(ELECTRIC_POTENTIAL,ELECTRIC_POTENTIAL);
    SPAR_MAT stiffness;
    (atc_->feEngine_)->compute_tangent_matrix(rhsMask,row_col, atc_->fields(), physicsModel_, 
                                                atc_->element_to_material_map(), stiffness);
    
    const BC_SET & bcs = (atc_->prescribedDataMgr_->bcs(ELECTRIC_POTENTIAL))[0];

    cddmPoissonSolver_ = new LinearSolver(stiffness, bcs, poissonSolverType_,
                                          -1, true);

    // initialize velocity solver
    const BCS & velocityBcs = atc_->prescribedDataMgr_->bcs(ELECTRON_VELOCITY);
    DENS_MAT velocityRhs(nNodes_,nsd_);
    atc_->compute_mass_matrix(ELECTRON_VELOCITY,physicsModel_);
    SPAR_MAT & velocityMassMat = (atc_->consistentMassMats_[ELECTRON_VELOCITY]).set_quantity();

    for (int i = 0; i < nsd_; i++ ) {
      LinearSolver * myVelocitySolver =
        new LinearSolver(velocityMassMat, velocityBcs[i],
                         LinearSolver::AUTO_SOLVE, -1, true);
      velocitySolvers_.push_back(myVelocitySolver);
    }
  }

  //--------------------------------------------------------
  //  pre initial integration
  //--------------------------------------------------------
  void ExtrinsicModelDriftDiffusionConvection::pre_init_integrate()
  {
    double dt = atc_->lammpsInterface_->dt();
    double time = atc_->time();
    int step = atc_->step();
    if (step % fluxUpdateFreq_ != 0) return; 

    // set Dirchlet data
    atc_->set_fixed_nodes();

    // set Neumann data (atc does not set these until post_final)
    atc_->set_sources();

    // subcyle integration of fast electron variable/s
    
    double idt = dt/nsubcycle_;
    for (int i = 0; i < nsubcycle_ ; ++i) {
      // update electron temperature mass matrix
      atc_->compute_mass_matrix(ELECTRON_VELOCITY,physicsModel_);
      // update electron velocity
      if (!(atc_->prescribedDataMgr_)->all_fixed(ELECTRON_VELOCITY))  {
        //const BCS & bcs 
        //  = atc_->prescribedDataMgr_->bcs(ELECTRON_VELOCITY);
        Array2D <bool> rhsMask(NUM_FIELDS,NUM_FLUX); 
        rhsMask = false;
        rhsMask(ELECTRON_VELOCITY,SOURCE) = atc_->fieldMask_(ELECTRON_VELOCITY,SOURCE);
        rhsMask(ELECTRON_VELOCITY,FLUX) = atc_->fieldMask_(ELECTRON_VELOCITY,FLUX);
        FIELDS rhs;
        rhs[ELECTRON_VELOCITY].reset(nNodes_,nsd_);
        atc_->compute_rhs_vector(rhsMask, atc_->fields_, rhs, atc_->source_integration(), physicsModel_);
        const DENS_MAT & velocityRhs =  rhs[ELECTRON_VELOCITY].quantity();
        // add a solver for electron momentum  
        DENS_MAT & velocity = (atc_->field(ELECTRON_VELOCITY)).set_quantity();
        for (int j = 0; j < nsd_; ++j) {
          if (! atc_->prescribedDataMgr_->all_fixed(ELECTRON_VELOCITY,j) ) {
            CLON_VEC v = column(velocity,j);
            const CLON_VEC r = column(velocityRhs,j);
            (velocitySolvers_[j])->solve(v,r);
          }
        }
      }
      
      //atc_->set_fixed_nodes();
      
      if (electronDensityEqn_ == ELECTRON_CONTINUITY) {
        // update continuity eqn
        if (! atc_->prescribedDataMgr_->all_fixed(ELECTRON_DENSITY) ) 
          continuityIntegrator_->update(idt,time,atc_->fields_,rhs_);
        atc_->set_fixed_nodes(); 
        // solve poisson eqn for electric potential
        
        if (! atc_->prescribedDataMgr_->all_fixed(ELECTRIC_POTENTIAL) ) {
        //poissonSolver_->solve(atc_->fields_,rhs_);
          Array2D <bool> rhsMask(NUM_FIELDS,NUM_FLUX); 
          rhsMask = false;
          rhsMask(ELECTRIC_POTENTIAL,SOURCE) = atc_->fieldMask_(ELECTRIC_POTENTIAL,SOURCE);
          rhsMask(ELECTRIC_POTENTIAL,PRESCRIBED_SOURCE) = atc_->fieldMask_(ELECTRIC_POTENTIAL,PRESCRIBED_SOURCE);
          FIELDS rhs;
          rhs[ELECTRIC_POTENTIAL].reset(nNodes_,1);
          atc_->compute_rhs_vector(rhsMask, atc_->fields_, rhs, atc_->source_integration(), physicsModel_);
          CLON_VEC x =column((atc_->field(ELECTRIC_POTENTIAL)).set_quantity(),0);
          const CLON_VEC r =column(rhs[ELECTRIC_POTENTIAL].quantity(),0);
          cddmPoissonSolver_->solve(x,r);
        }
      } 
      else if (electronDensityEqn_ == ELECTRON_SCHRODINGER) {
        schrodingerPoissonSolver_->solve(rhs_,fluxes_);
      }
      
      atc_->set_fixed_nodes(); 
      // update electron temperature mass matrix
      atc_->compute_mass_matrix(ELECTRON_TEMPERATURE,physicsModel_);
      // update electron temperature
      if (! atc_->prescribedDataMgr_->all_fixed(ELECTRON_TEMPERATURE) ) 
        temperatureIntegrator_->update(idt,time,atc_->fields_,rhs_);
      atc_->set_fixed_nodes(); 
      
    }

  }

  //--------------------------------------------------------
  //  output
  //--------------------------------------------------------
  void ExtrinsicModelDriftDiffusionConvection::output(OUTPUT_LIST & outputData)
  {
    ExtrinsicModelDriftDiffusion::output(outputData);
    
    
    //FIELD jouleHeating(atc_->num_nodes(),1);
    //set_kinetic_energy_source(atc_->fields(),jouleHeating);
    outputData["joule_heating"] = & (atc_->extrinsic_source(TEMPERATURE)).set_quantity();

    // globals
    DENS_MAT nodalKineticEnergy;
    compute_nodal_kinetic_energy(nodalKineticEnergy);
    double kineticEnergy = nodalKineticEnergy.sum();
    atc_->feEngine_->add_global("total_electron_kinetic_energy",kineticEnergy);
  }

  //--------------------------------------------------------
  //  size_vector
  //--------------------------------------------------------
  int ExtrinsicModelDriftDiffusionConvection::size_vector(int intrinsicSize)
  {
    int xSize = ExtrinsicModelDriftDiffusion::size_vector(intrinsicSize);
    baseSize_ = intrinsicSize  + xSize;
    xSize += 1;
    return xSize;
  }

  //--------------------------------------------------------
  //  compute_vector
  //--------------------------------------------------------
  bool ExtrinsicModelDriftDiffusionConvection::compute_vector(int n, double & value)
  {
    // output[1] = total electron kinetic energy

    bool match = ExtrinsicModelDriftDiffusion::compute_vector(n,value);
    if (match) return match;

    if (n == baseSize_) {
      DENS_MAT nodalKineticEnergy;
      compute_nodal_kinetic_energy(nodalKineticEnergy);
      value = nodalKineticEnergy.sum();
      return true;
    }
    return false;
  }

  //--------------------------------------------------------
  //  compute_kinetic_energy
  //--------------------------------------------------------
  void ExtrinsicModelDriftDiffusionConvection::compute_nodal_kinetic_energy(DENS_MAT & kineticEnergy)
  {
    DENS_MAT & velocity((atc_->field(ELECTRON_VELOCITY)).set_quantity());
    SPAR_MAT & velocityMassMat = (atc_->consistentMassMats_[ELECTRON_VELOCITY]).set_quantity();
    kineticEnergy.reset(nNodes_,1);
    
    for (int j = 0; j < nsd_; j++) {
      CLON_VEC myVelocity(velocity,CLONE_COL,j);
      DENS_MAT velocityMat(nNodes_,1);
      for (int i = 0; i < nNodes_; i++)
        velocityMat(i,0) = myVelocity(i);
      kineticEnergy += velocityMat.mult_by_element(myVelocity);
    }

    kineticEnergy = 0.5*velocityMassMat*kineticEnergy;
  }
};
