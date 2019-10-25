// ATC transfer headers
#include "ElasticTimeIntegrator.h"
#include "ATC_Coupling.h"
#include "TimeFilter.h"
#include "ATC_Error.h"
#include "PerAtomQuantityLibrary.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class MomentumTimeIntegrator
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //-------------------------------------------------------- 
  MomentumTimeIntegrator::MomentumTimeIntegrator(ATC_Coupling * atc,
                                                 TimeIntegrationType timeIntegrationType) :
    TimeIntegrator(atc, timeIntegrationType)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  modify
  //    parses inputs and modifies state of the integrator
  //--------------------------------------------------------
  bool MomentumTimeIntegrator::modify(int /* narg */, char **arg)
  {
    bool foundMatch = false;
    int argIndex = 0;
    // time integration scheme
    /*! \page man_momentum_time_integration fix_modify AtC time_integration (momentum)
      \section syntax
      fix_modify AtC time_integration <descriptor> \n
      - descriptor (string) = time integration type  \n
      
      various time integration methods for the finite elements\n
      \section description
      verlet - atomic velocity update with 2nd order Verlet, nodal temperature update with 2nd order Verlet, kinetostats based on controlling force \n
      fractional_step - atomic velocity update with 2nd order Verlet, mixed nodal momentum update, 2nd order Verlet for continuum and exact 2nd order Verlet for atomic contributions, kinetostats based on controlling discrete momentum changes\n
      gear - atomic velocity update with 2nd order Verlet, nodal temperature update with 3rd or 4th order Gear, kinetostats based on controlling power \n
      \section examples
      <TT> fix_modify atc time_integration verlet </TT> \n
      <TT> fix_modify atc time_integration fractional_step </TT> \n
      \section description
      \section related
      see \ref man_fix_atc
      \section default
      none 
    */
    if (strcmp(arg[argIndex],"verlet")==0) {
      timeIntegrationType_ = VERLET;
      needReset_ = true;
      foundMatch = true;
    }
    else if (strcmp(arg[argIndex],"fractional_step")==0) {
      timeIntegrationType_ = FRACTIONAL_STEP;
      needReset_ = true;
      foundMatch = true;
    }
    else if (strcmp(arg[argIndex],"gear")==0) {
      timeIntegrationType_ = GEAR;
      needReset_ = true;
      foundMatch = true;
    }
    return foundMatch;
  }

  //--------------------------------------------------------
  //  construct_methods
  //    creates algorithm objects
  //--------------------------------------------------------
  void MomentumTimeIntegrator::construct_methods()
  {
    if (atc_->reset_methods()) {
      if (timeIntegrationMethod_)
        delete timeIntegrationMethod_;
          
      if (timeFilterManager_->need_reset()) {
        switch (timeIntegrationType_) {
          case VERLET:
            timeFilter_ = timeFilterManager_->construct(TimeFilterManager::IMPLICIT);
            atc_->set_mass_mat_time_filter(MOMENTUM,TimeFilterManager::IMPLICIT);
            break;
          case FRACTIONAL_STEP:
          case GEAR:
            timeFilter_ = timeFilterManager_->construct(TimeFilterManager::EXPLICIT_IMPLICIT);
            atc_->set_mass_mat_time_filter(MOMENTUM,TimeFilterManager::EXPLICIT_IMPLICIT);
            break;
          default:
            throw ATC_Error("Uknown time integration type in ThermalTimeIntegrator::Initialize()");
        }
      }

      if (timeFilterManager_->filter_dynamics()) {
        switch (timeIntegrationType_) {
          case VERLET: {
            timeIntegrationMethod_ = new ElasticTimeIntegratorVerletFiltered(this);
            break;
          }
        default:
          throw ATC_Error("Uknown time integration type in MomentumTimeIntegrator::Initialize()");
        }
      }
      else {
        switch (timeIntegrationType_) {
          case VERLET: {
            timeIntegrationMethod_ = new ElasticTimeIntegratorVerlet(this);
            break;
          }
          case FRACTIONAL_STEP: {
            timeIntegrationMethod_ = new ElasticTimeIntegratorFractionalStep(this);
            break;
          }
          case GEAR: {
            timeIntegrationMethod_ = new FluidsTimeIntegratorGear(this);
            break;
          }
        default:
          throw ATC_Error("Uknown time integration type in MomentumTimeIntegrator::Initialize()");
        }
      }   
    }
  }

  //--------------------------------------------------------
  //  pack_fields
  //    add persistent variables to data list
  //--------------------------------------------------------
  void MomentumTimeIntegrator::pack_fields(RESTART_LIST & data)
  {
    data["NodalAtomicForceFiltered"] = & nodalAtomicForceFiltered_.set_quantity();
    data["NodalAtomicMomentumFiltered"] = & nodalAtomicMomentumFiltered_.set_quantity();
    TimeIntegrator::pack_fields(data);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class MomentumIntegrationMethod
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //        Grab data from ATC
  //-------------------------------------------------------- 
  MomentumIntegrationMethod::MomentumIntegrationMethod(MomentumTimeIntegrator * momentumTimeIntegrator) :
    TimeIntegrationMethod(momentumTimeIntegrator),
    timeFilter_(timeIntegrator_->time_filter()),
    velocity_(atc_->field(VELOCITY)),
    acceleration_(atc_->field_roc(VELOCITY)),
    nodalAtomicVelocityOut_(atc_->nodal_atomic_field(VELOCITY)),
    velocityRhs_(atc_->field_rhs(VELOCITY)),
    nodalAtomicForceOut_(atc_->nodal_atomic_field_roc(VELOCITY))
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  construct_transfers
  //        Grab existing managed quantities,
  //        create the rest
  //--------------------------------------------------------
  void MomentumIntegrationMethod::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());
    nodalAtomicVelocity_ = interscaleManager.dense_matrix("NodalAtomicVelocity");
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ElasticTimeIntegratorVerlet
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //-------------------------------------------------------- 
  ElasticTimeIntegratorVerlet::ElasticTimeIntegratorVerlet(MomentumTimeIntegrator * momentumTimeIntegrator) :
    MomentumIntegrationMethod(momentumTimeIntegrator),
    displacement_(atc_->field(DISPLACEMENT)),
    nodalAtomicDisplacementOut_(atc_->nodal_atomic_field(DISPLACEMENT)),
    nodalAtomicForceFiltered_(momentumTimeIntegrator->nodal_atomic_force_filtered()),
    nodalAtomicDisplacement_(NULL),
    nodalAtomicForce_(NULL)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  construct_transfers
  //        Grab existing managed quantities,
  //        create the rest
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::construct_transfers()
  {
    MomentumIntegrationMethod::construct_transfers();
    InterscaleManager & interscaleManager = atc_->interscale_manager();
    nodalAtomicDisplacement_ = interscaleManager.dense_matrix("NodalAtomicDisplacement");
    nodalAtomicForce_ = interscaleManager.dense_matrix("NodalAtomicForce");
  }

  //--------------------------------------------------------
  //  initialize
  //        initialize all data
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::initialize()
  {
    MomentumIntegrationMethod::initialize();

    // sets up time filter for cases where variables temporally filtered
    // this time integrator should use an implicit filter
    TimeFilterManager * timeFilterManager = (timeIntegrator_->atc())->time_filter_manager();
    if (timeFilterManager->need_reset()) {
      timeFilter_->initialize(nodalAtomicForce_->quantity());
    }
    
    if (!(timeFilterManager->end_equilibrate())) {
      nodalAtomicForceFiltered_.reset(atc_->num_nodes(),atc_->nsd());
    }

    if (!(timeFilterManager->filter_dynamics())){
      //post_process();
      //compute_nodal_forces(velocityRhs_.set_quantity());
      velocityRhs_ = nodalAtomicForce_->quantity();
    }
  }

  //--------------------------------------------------------
  //  pre_initial_integrate1
  //    time integration before Verlet step 1
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::pre_initial_integrate1(double dt)
  {
    explicit_1(velocity_.set_quantity(),acceleration_.quantity(),.5*dt);
  }

  
  //--------------------------------------------------------
  //  post_initial_integrate1
  //    time integration after Verlet step 1
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::post_initial_integrate1(double dt)
  {
    
    //      for improved accuracy, but this would be inconsistent with
    //      the atomic integration scheme
    explicit_1(displacement_.set_quantity(),velocity_.quantity(),dt);
  }

  //--------------------------------------------------------
  //  pre_final_integrate1
  //    first time integration computations 
  //    before Verlet step 2
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::pre_final_integrate1(double dt)
  {
    // integrate filtered atomic force
    timeFilter_->apply_post_step2(nodalAtomicForceFiltered_.set_quantity(),
                                  nodalAtomicForce_->quantity(),dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate2
  //    second time integration computations
  //    after Verlet step 2
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::post_final_integrate2(double dt)
  {
    atc_->apply_inverse_mass_matrix(velocityRhs_.quantity(),
                                    acceleration_.set_quantity(),
                                    VELOCITY);
    explicit_1(velocity_.set_quantity(),acceleration_.quantity(),.5*dt);
  }

  //--------------------------------------------------------
  //  add_to_rhs
  //    add integrated atomic force contributions
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::add_to_rhs()
  {
    // Compute MD contribution to FEM equation
    velocityRhs_ += nodalAtomicForce_->quantity();
  }

  //--------------------------------------------------------
  //  post_process
  //    post processing of variables before output
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::post_process()
  {
    nodalAtomicDisplacementOut_ = nodalAtomicDisplacement_->quantity();
    nodalAtomicVelocityOut_ = nodalAtomicVelocity_->quantity();
  }

  //--------------------------------------------------------
  //  output
  //    add variables to output list
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::output(OUTPUT_LIST & outputData)
  {
    DENS_MAT & nodalAtomicForce(nodalAtomicForce_->set_quantity());
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData["NodalAtomicForce"] = &nodalAtomicForce;
    }
  }

  //--------------------------------------------------------
  //  finish
  //    finalize state of nodal atomic quantities
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::finish()
  {
    post_process();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ElasticTimeIntegratorVerletFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //-------------------------------------------------------- 
  ElasticTimeIntegratorVerletFiltered::ElasticTimeIntegratorVerletFiltered(MomentumTimeIntegrator * momentumTimeIntegrator) :
    ElasticTimeIntegratorVerlet(momentumTimeIntegrator),
    nodalAtomicAcceleration_(atc_->nodal_atomic_field_roc(VELOCITY))
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  pre_initial_integrate1
  //    time integration before Verlet step 1
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerletFiltered::pre_initial_integrate1(double dt)
  {
    explicit_1(velocity_.set_quantity(),acceleration_.quantity(),.5*dt);
    explicit_1(nodalAtomicVelocityOut_.set_quantity(),nodalAtomicAcceleration_.quantity(),.5*dt);
  }

  //--------------------------------------------------------
  //  post_initial_integrate1
  //    time integration after Verlet step 1
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerletFiltered::post_initial_integrate1(double dt)
  {
    
    //      for improved accuracy, but this would be inconsistent with
    //      the atomic integration scheme
    explicit_1(displacement_.set_quantity(),velocity_.quantity(),dt);
    explicit_1(nodalAtomicDisplacementOut_.set_quantity(),nodalAtomicVelocityOut_.quantity(),dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate2
  //    second time integration after Verlet step 2
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerletFiltered::post_final_integrate2(double dt)
  {
    DENS_MAT velocityRoc(velocityRhs_.nRows(),velocityRhs_.nCols());
    atc_->apply_inverse_mass_matrix(velocityRhs_.quantity(),
                                    acceleration_.set_quantity(),
                                    VELOCITY);
    explicit_1(velocity_.set_quantity(),acceleration_.quantity(),.5*dt);
    
    atc_->apply_inverse_md_mass_matrix(nodalAtomicForceFiltered_.quantity(),
                                       nodalAtomicAcceleration_.set_quantity(),
                                       VELOCITY);
    explicit_1(nodalAtomicVelocityOut_.set_quantity(),nodalAtomicAcceleration_.quantity(),.5*dt);
  }

  //--------------------------------------------------------
  //  add_to_rhs
  //    add integrated atomic force contributions
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerletFiltered::add_to_rhs()
  {
    // MD contributions to FE equations
    velocityRhs_ += nodalAtomicForceFiltered_.set_quantity();
  }

  //--------------------------------------------------------
  //  output
  //    add variables to output list
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerletFiltered::output(OUTPUT_LIST & outputData)
  {
    DENS_MAT & nodalAtomicForce(nodalAtomicForceFiltered_.set_quantity());
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData["NodalAtomicForce"] = &nodalAtomicForce;
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ElasticTimeIntegratorFractionalStep
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //-------------------------------------------------------- 
  ElasticTimeIntegratorFractionalStep::ElasticTimeIntegratorFractionalStep(MomentumTimeIntegrator * momentumTimeIntegrator) :
    MomentumIntegrationMethod(momentumTimeIntegrator),
    displacement_(atc_->field(DISPLACEMENT)),
    nodalAtomicDisplacementOut_(atc_->nodal_atomic_field(DISPLACEMENT)),
    nodalAtomicForceFiltered_(momentumTimeIntegrator->nodal_atomic_force_filtered()),
    nodalAtomicMomentum_(NULL),
    nodalAtomicMomentumFiltered_(momentumTimeIntegrator->nodal_atomic_momentum_filtered()),
    nodalAtomicDisplacement_(NULL),
    nodalAtomicMomentumOld_(atc_->num_nodes(),atc_->nsd()), 
    nodalAtomicVelocityOld_(atc_->num_nodes(),atc_->nsd())
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  construct_transfers
  //        Grab existing managed quantities,
  //        create the rest
  //--------------------------------------------------------
  void ElasticTimeIntegratorFractionalStep::construct_transfers()
  {
    MomentumIntegrationMethod::construct_transfers();
    InterscaleManager & interscaleManager = atc_->interscale_manager();
    nodalAtomicMomentum_ = interscaleManager.dense_matrix("NodalAtomicMomentum");
    nodalAtomicDisplacement_ = interscaleManager.dense_matrix("NodalAtomicDisplacement");
  }

  //--------------------------------------------------------
  //  initialize
  //        initialize all data
  //--------------------------------------------------------
  void ElasticTimeIntegratorFractionalStep::initialize()
  {
    MomentumIntegrationMethod::initialize();

    // initial force to zero
    nodalAtomicForce_.reset(atc_->num_nodes(),atc_->nsd());

    // sets up time filter for cases where variables temporally filtered
    // this time integrator should use Crank-Nicholson filter for 2nd order accuracy
    TimeFilterManager * timeFilterManager = (timeIntegrator_->atc())->time_filter_manager();
    if (timeFilterManager->need_reset()) {
      // the form of this integrator implies no time filters that require history data can be used
      timeFilter_->initialize();
    }
    
    // sets up time filter for post-processing the filtered power
    // this time integrator should use an explicit-implicit filter
    // to mirror the 2nd order Verlet integration scheme
    // It requires no history information so initial value just sizes arrays
    if (!(timeFilterManager->end_equilibrate())) {
      // implies an initial condition of the instantaneous atomic energy
      nodalAtomicMomentumFiltered_ = nodalAtomicMomentum_->quantity();
      nodalAtomicForceFiltered_.reset(atc_->num_nodes(),atc_->nsd());
    }
  }

  //--------------------------------------------------------
  //  pre_initial_integrate1
  //    time integration before Verlet step 1, used to
  //    provide the baseline momentum and displacement to
  //    quantify the change over the timestep
  //--------------------------------------------------------
  void ElasticTimeIntegratorFractionalStep::pre_initial_integrate1(double dt)
  {
    // initialize changes in momentum
    const DENS_MAT & myNodalAtomicMomentum(nodalAtomicMomentum_->quantity());
    // updated filtered energy using explicit-implicit scheme
    timeFilter_->apply_pre_step1(nodalAtomicMomentumFiltered_.set_quantity(),
                                 myNodalAtomicMomentum,dt);
  }

  //--------------------------------------------------------
  //  pre_initial_integrate2
  //    second time integration after kinetostat application
  //    to compute MD contributions to momentum change
  //--------------------------------------------------------
  void ElasticTimeIntegratorFractionalStep::pre_initial_integrate2(double dt)
  {
    // used for updating change in velocity from mass matrix change
    this->compute_old_time_data();

    // update filtered nodal atomic force
    timeFilter_->apply_pre_step1(nodalAtomicForceFiltered_.set_quantity(),
                                 nodalAtomicForce_,dt);

    // store current force for use later
    nodalAtomicForce_ = nodalAtomicMomentum_->quantity();
    nodalAtomicForce_ *= -1.;
  }

  //--------------------------------------------------------
  //  post_initial_integrate1
  //    time integration after Verlet step 1
  //--------------------------------------------------------
  void ElasticTimeIntegratorFractionalStep::post_initial_integrate1(double dt)
  {
    // atomic contributions to change in momentum
    // compute change in restricted atomic momentum
    const DENS_MAT & nodalAtomicMomentum(nodalAtomicMomentum_->quantity());
    nodalAtomicForce_ += nodalAtomicMomentum;

    // update FE velocity with change in velocity from MD
    DENS_MAT & atomicVelocityDelta(atomicVelocityDelta_.set_quantity());
    atc_->apply_inverse_mass_matrix(nodalAtomicForce_,
                                    atomicVelocityDelta,
                                    VELOCITY);
    velocity_ += atomicVelocityDelta;
 
    // approximation to force for output
    nodalAtomicForce_ /= 0.5*dt;
    timeFilter_->apply_post_step1(nodalAtomicForceFiltered_.set_quantity(),
                                  nodalAtomicForce_,dt);

    // change to velocity from FE dynamics
    atc_->apply_inverse_mass_matrix(velocityRhs_.quantity(),
                                    acceleration_.set_quantity(),
                                    VELOCITY);
    explicit_1(velocity_.set_quantity(),acceleration_.quantity(),0.5*dt);

    // used for updating change in momentum from mass matrix change
    atc_->apply_inverse_mass_matrix(nodalAtomicMomentum,
                                    nodalAtomicVelocityOld_,
                                    VELOCITY);
    nodalAtomicMomentumOld_ = nodalAtomicMomentum;

    // get nodal momentum for second part of force update
    nodalAtomicForce_ = nodalAtomicMomentum;
    nodalAtomicForce_ *= -1.;

    // update nodal displacements
    explicit_1(displacement_.set_quantity(),velocity_.quantity(),dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate2
  //    second time integration computations
  //    after Verlet step 2
  //--------------------------------------------------------
  void ElasticTimeIntegratorFractionalStep::post_final_integrate2(double dt)
  {
    // atomic contributions to change in momentum
    // compute change in restricted atomic momentum
    nodalAtomicForce_ += nodalAtomicMomentum_->quantity();
    
    // update FE temperature with change in temperature from MD
    compute_velocity_delta(nodalAtomicForce_,dt);
    velocity_ += atomicVelocityDelta_.quantity();
    
    // approximation to power for output
    nodalAtomicForce_ /= 0.5*dt;
    timeFilter_->apply_post_step1(nodalAtomicForceFiltered_.set_quantity(),
                                  nodalAtomicForce_,dt);
    
    // change to velocity from FE dynamics
    atc_->apply_inverse_mass_matrix(velocityRhs_.quantity(),
                                    acceleration_.set_quantity(),
                                    VELOCITY);
    explicit_1(velocity_.set_quantity(),acceleration_.quantity(),0.5*dt);
  }
  //--------------------------------------------------------
  //  post_process
  //    post processing of variables before output
  //--------------------------------------------------------
  void ElasticTimeIntegratorFractionalStep::post_process()
  {
    nodalAtomicDisplacementOut_ = nodalAtomicDisplacement_->quantity();
    nodalAtomicVelocityOut_ = nodalAtomicVelocity_->quantity();
  }

  //--------------------------------------------------------
  //  output
  //    add variables to output list
  //--------------------------------------------------------
  void ElasticTimeIntegratorFractionalStep::output(OUTPUT_LIST & outputData)
  {
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData["NodalAtomicForce"] = & nodalAtomicForce_;
    }
  }

  //--------------------------------------------------------
  //  finish
  //    finalize state of nodal atomic quantities
  //--------------------------------------------------------
  void ElasticTimeIntegratorFractionalStep::finish()
  {
    post_process();
  }

  //--------------------------------------------------------
  //  compute_old_time_data
  //--------------------------------------------------------
  void ElasticTimeIntegratorFractionalStep::compute_old_time_data()
  {
    const DENS_MAT & myNodalAtomicMomentum(nodalAtomicMomentum_->quantity());
    atc_->apply_inverse_mass_matrix(myNodalAtomicMomentum,
                                    nodalAtomicVelocityOld_,
                                    VELOCITY);
    nodalAtomicMomentumOld_ = myNodalAtomicMomentum;
  }

  //--------------------------------------------------------
  //  compute_velocity_delta
  //--------------------------------------------------------
  void ElasticTimeIntegratorFractionalStep::compute_velocity_delta(const DENS_MAT & atomicMomentumDelta,
                                                                   double /* dt */)
  {
    DENS_MAT & myAtomicVelocityDelta(atomicVelocityDelta_.set_quantity());
    myAtomicVelocityDelta = nodalAtomicMomentumOld_ + atomicMomentumDelta;
    atc_->apply_inverse_mass_matrix(myAtomicVelocityDelta,
                                    VELOCITY);
    myAtomicVelocityDelta += -1.*nodalAtomicVelocityOld_;
  }
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FluidsTimeIntegratorGear
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //        Grab data from ATC
  //--------------------------------------------------------

  FluidsTimeIntegratorGear::FluidsTimeIntegratorGear(MomentumTimeIntegrator * momentumTimeIntegrator) :
    MomentumIntegrationMethod(momentumTimeIntegrator),
    nodalAtomicForceFiltered_(momentumTimeIntegrator->nodal_atomic_force_filtered()),
    nodalAtomicMomentum_(NULL),
    nodalAtomicMomentumFiltered_(momentumTimeIntegrator->nodal_atomic_momentum_filtered()),
    atomicVelocityDelta_(atc_->num_nodes(),atc_->nsd()),
    nodalAtomicMomentumOld_(atc_->num_nodes(),atc_->nsd()),
    nodalAtomicVelocityOld_(atc_->num_nodes(),atc_->nsd()),
    velocity2Roc_(atc_->field_2roc(VELOCITY))
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  construct_transfers
  //        Grab existing managed quantities,
  //        create the rest
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::construct_transfers()
  {
    MomentumIntegrationMethod::construct_transfers();
    InterscaleManager & interscaleManager(atc_->interscale_manager());
    nodalAtomicMomentum_ = interscaleManager.dense_matrix("NodalAtomicMomentum");
  }

  //--------------------------------------------------------
  //  initialize
  //        initialize all data
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::initialize()
  {
    MomentumIntegrationMethod::initialize();

    // initial power to zero
    nodalAtomicForce_.reset(atc_->num_nodes(),atc_->nsd());

    // sets up time filter for cases where variables temporally filtered
    // this time integrator should use Crank-Nicholson filter for 2nd order accuracy
    TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
    if (timeFilterManager->need_reset()) {
      // the form of this integrator implies no time filters that require history data can be used
      timeFilter_->initialize();
    }

    // sets up time filter for post-processing the filtered power
    // this time integrator should use an explicit-implicit filter
    // to mirror the 2nd order Verlet integration scheme
    // It requires no history information so initial value just sizes arrays
    if (!timeFilterManager->end_equilibrate()) {
      // implies an initial condition of the instantaneous atomic energy
      // for the corresponding filtered variable, consistent with the temperature
      nodalAtomicMomentumFiltered_ = nodalAtomicMomentum_->quantity();  
      nodalAtomicForceFiltered_.reset(atc_->num_nodes(),atc_->nsd());
    }
  }

  //--------------------------------------------------------
  //  pre_initial_integrate1
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::pre_initial_integrate1(double dt)
  {
    const DENS_MAT & myNodalAtomicMomentum(nodalAtomicMomentum_->quantity());
    // updated filtered momentum using explicit-implicit scheme
    timeFilter_->apply_pre_step1(nodalAtomicMomentumFiltered_.set_quantity(),
                                 myNodalAtomicMomentum,dt);
  }

  //--------------------------------------------------------
  //  pre_initial_integrate2
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::pre_initial_integrate2(double dt)
  {
    // used for updating change in velocity from mass matrix change
    this->compute_old_time_data();

    // update FE contributions
    apply_gear_predictor(dt);

    // update filtered nodal atomic force
    
    //      that way kinetostat and integrator can be consistent
    timeFilter_->apply_pre_step1(nodalAtomicForceFiltered_.set_quantity(),
                                 nodalAtomicForce_,dt);

    // store current momentum for use later
    nodalAtomicForce_ = nodalAtomicMomentum_->quantity();
    nodalAtomicForce_ *= -1.;
  }

  //--------------------------------------------------------
  //  pre_final_integrate1
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::pre_final_integrate1(double dt)
  {
    
    //      before the new rhs is computed but after atomic velocity is updated.
    // compute change in restricted atomic momentum
    nodalAtomicForce_ += nodalAtomicMomentum_->quantity();
    
    // update FE velocity with change in velocity from MD
    compute_velocity_delta(nodalAtomicForce_,dt);
    velocity_ += atomicVelocityDelta_.quantity();
    
    // approximation to force for output
    nodalAtomicForce_ /= dt;
    timeFilter_->apply_post_step1(nodalAtomicForceFiltered_.set_quantity(),
                                  nodalAtomicForce_,dt);

    // make sure nodes are fixed
    atc_->set_fixed_nodes();
  }

  //--------------------------------------------------------
  //  post_final_integrate1
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::post_final_integrate1(double dt)
  {
    // Finish updating temperature with FE contributions
    atc_->apply_inverse_mass_matrix(velocityRhs_.quantity(),
                                    _velocityResidual_,VELOCITY);
    _velocityResidual_ -= acceleration_.quantity();
    _velocityResidual_ *= dt;
    apply_gear_corrector(_velocityResidual_,dt);
  }

  //--------------------------------------------------------
  //  post_process
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::post_final_integrate2(double dt)
  {
    // update filtered atomic energy
    timeFilter_->apply_post_step1(nodalAtomicMomentumFiltered_.set_quantity(),
                                  nodalAtomicMomentum_->quantity(),dt);
  }

  //--------------------------------------------------------
  //  output
  //    add variables to output list
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::post_process()
  {
    nodalAtomicForceOut_ = nodalAtomicForce_;
    nodalAtomicVelocityOut_ = nodalAtomicVelocity_->quantity();
  }

  //--------------------------------------------------------
  //  output
  //    add variables to output list
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::output(OUTPUT_LIST & outputData)
  {
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData["NodalAtomicForce"] = & nodalAtomicForce_;
    }
  }

  //--------------------------------------------------------
  //  finish
  //    finalize state of nodal atomic quantities
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::finish()
  {
    post_process();
  }

  //--------------------------------------------------------
  //  apply_gear_predictor
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::apply_gear_predictor(double dt)
  {
    gear1_3_predict(velocity_.set_quantity(),
                    acceleration_.set_quantity(),
                    velocity2Roc_.quantity(),dt);
  }

  //--------------------------------------------------------
  //  apply_gear_corrector
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::apply_gear_corrector(const DENS_MAT & residual, double dt)
  {
    gear1_3_correct(velocity_.set_quantity(),
                    acceleration_.set_quantity(),
                    velocity2Roc_.set_quantity(),
                    residual,dt);
  }

  //--------------------------------------------------------
  //  compute_old_time_data
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::compute_old_time_data()
  {
    const DENS_MAT & myNodalAtomicMomentum(nodalAtomicMomentum_->quantity());
    atc_->apply_inverse_mass_matrix(myNodalAtomicMomentum,
                                    nodalAtomicVelocityOld_,
                                    VELOCITY);
    nodalAtomicMomentumOld_ = myNodalAtomicMomentum;
  }

  //--------------------------------------------------------
  //  compute_velocity_delta
  //--------------------------------------------------------
  void FluidsTimeIntegratorGear::compute_velocity_delta(const DENS_MAT & atomicMomentumDelta,
                                                        double /* dt */)
  {
    DENS_MAT & myAtomicVelocityDelta(atomicVelocityDelta_.set_quantity());
    myAtomicVelocityDelta = nodalAtomicMomentumOld_ + atomicMomentumDelta;
    atc_->apply_inverse_mass_matrix(myAtomicVelocityDelta,
                                    VELOCITY);
    myAtomicVelocityDelta -= nodalAtomicVelocityOld_;
  }
};
