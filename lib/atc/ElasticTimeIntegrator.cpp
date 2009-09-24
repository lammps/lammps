// ATC transfer headers
#include "ElasticTimeIntegrator.h"
#include "ATC_Transfer.h"
#include "TimeFilter.h"
#include "ATC_Error.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ElasticTimeIntegrator
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //        Grab data from ATC
  //-------------------------------------------------------- 
  ElasticTimeIntegrator::ElasticTimeIntegrator(ATC_Transfer * atcTransfer,
                                               TimeIntegrationType timeIntegrationType) :
    TimeIntegrator(atcTransfer, timeIntegrationType)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  modify
  //    parses inputs and modifies state of the filter
  //--------------------------------------------------------
  bool ElasticTimeIntegrator::modify(int narg, char **arg)
  {
    // currently no parsing for elastic time integration
    return false;
  }

  //--------------------------------------------------------
  //  initialize
  //    sets up all the necessary data
  //--------------------------------------------------------
  void ElasticTimeIntegrator::initialize()
  {
    if (needReset_ || timeFilterManager_->need_reset()) {
      if (timeIntegrationMethod_)
        delete timeIntegrationMethod_;
          
      if (timeFilterManager_->need_reset()) {
        if (timeFilter_)
          delete timeFilter_;
        timeFilter_ = timeFilterManager_->construct(TimeFilterManager::IMPLICIT);
      }

      if (timeFilterManager_->filter_dynamics()) {
        timeIntegrationMethod_ = new ElasticTimeIntegratorVerletFiltered(this);
      }
      else {
        timeIntegrationMethod_ = new ElasticTimeIntegratorVerlet(this);
      }   
    }
    TimeIntegrator::initialize();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ElasticIntegrationMethod
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //        Grab data from ATC
  //-------------------------------------------------------- 
  ElasticIntegrationMethod::ElasticIntegrationMethod(ElasticTimeIntegrator * elasticTimeIntegrator) :
    TimeIntegrationMethod(elasticTimeIntegrator),
    timeFilter_(timeIntegrator_->get_time_filter()),
    displacement_(atcTransfer_->get_field(DISPLACEMENT)),
    velocity_(atcTransfer_->get_field(VELOCITY)),
    acceleration_(atcTransfer_->get_field_roc(VELOCITY)),
    nodalAtomicDisplacement_(atcTransfer_->get_atomic_field(DISPLACEMENT)),
    nodalAtomicVelocity_(atcTransfer_->get_atomic_field(VELOCITY)),
    velocityRhs_(atcTransfer_->get_field_rhs(VELOCITY)),
    nodalAtomicForce_(atcTransfer_->get_fe_atomic_field_roc(VELOCITY)),
    forceFilteringData_(atcTransfer_->get_aux_storage(VELOCITY))
  {
    // do nothing
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ElasticTimeIntegratorVerlet
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //-------------------------------------------------------- 
  ElasticTimeIntegratorVerlet::ElasticTimeIntegratorVerlet(ElasticTimeIntegrator * elasticTimeIntegrator) :
    ElasticIntegrationMethod(elasticTimeIntegrator)
  {
    TimeFilterManager * timeFilterManager = (timeIntegrator_->get_atc_transfer())->get_time_filter_manager();
    if (timeFilterManager->need_reset()) {
      timeFilter_->initialize(nodalAtomicForce_);
    }
    // reset filtering data, if needed
    if (!(timeFilterManager->end_equilibrate())) {
      forceFilteringData_.reset(atcTransfer_->get_nNodes(),atcTransfer_->get_nsd());
    }
  }

  //--------------------------------------------------------
  //  mid_initial_integrate1
  //    time integration at the mid-point of Verlet step 1
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::mid_initial_integrate1(double dt)
  {
    explicit_1(velocity_,acceleration_,.5*dt);
  }

  //--------------------------------------------------------
  //  post_initial_integrate1
  //    time integration after Verlet step 1
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::post_initial_integrate1(double dt)
  {
    // NOTE could use explicit_2 with velocityRhs_ as the 2nd derivative
    //      for improved accuracy, but this would be inconsistent with
    //      the atomic integration scheme
    explicit_1(displacement_,velocity_,dt);
  }

  //--------------------------------------------------------
  //  pre_final_integrate1
  //    first time integration computations 
  //    before Verlet step 2
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::pre_final_integrate1(double dt)
  {
    // Compute MD contribution to FEM equation
    DENS_MAT atomicForces;
    atcTransfer_->compute_atomic_force(atomicForces,atcTransfer_->get_f());
    atcTransfer_->restrict_volumetric_quantity(atomicForces,nodalAtomicForce_);
    timeFilter_->apply_post_step2(forceFilteringData_,nodalAtomicForce_,dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate1
  //    second time integration computations
  //    after Verlet step 2
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::post_final_integrate1(double dt)
  {
    atcTransfer_->apply_inverse_mass_matrix(velocityRhs_,
                                            acceleration_,
                                            VELOCITY);
    explicit_1(velocity_,acceleration_,.5*dt);
  }

  //--------------------------------------------------------
  //  post_process
  //    post processing of variables before output
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerlet::post_process(double dt)
  {
    DENS_MAT atomicQuantity;
    atcTransfer_->compute_atomic_momentum(atomicQuantity,atcTransfer_->get_v());
    atcTransfer_->project_md_volumetric_quantity(atomicQuantity,nodalAtomicVelocity_,VELOCITY);
    
    atcTransfer_->compute_atomic_centerOfMass_displacement(atomicQuantity,atcTransfer_->get_x());
    atcTransfer_->project_md_volumetric_quantity(atomicQuantity,nodalAtomicDisplacement_,VELOCITY);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ElasticTimeIntegratorVerletFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //-------------------------------------------------------- 
  ElasticTimeIntegratorVerletFiltered::ElasticTimeIntegratorVerletFiltered(ElasticTimeIntegrator * elasticTimeIntegrator) :
    ElasticTimeIntegratorVerlet(elasticTimeIntegrator),
    nodalAtomicAcceleration_(atcTransfer_->get_atomic_field_roc(VELOCITY))
  {
    // swap filtered and unfiltered forces
    if ((timeIntegrator_->get_time_filter_manager())->end_equilibrate()) {
      DENS_MAT temp(nodalAtomicForce_);
      nodalAtomicForce_ = forceFilteringData_;
      forceFilteringData_ = temp;
    }
  }

  //--------------------------------------------------------
  //  mid_initial_integrate1
  //    time integration at the mid-point of Verlet step 1
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerletFiltered::mid_initial_integrate1(double dt)
  {
    explicit_1(velocity_,acceleration_,.5*dt);
    explicit_1(nodalAtomicVelocity_,nodalAtomicAcceleration_,.5*dt);
  }

  //--------------------------------------------------------
  //  post_initial_integrate1
  //    time integration after Verlet step 1
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerletFiltered::post_initial_integrate1(double dt)
  {
    // NOTE could use explicit_2 with velocityRhs_ as the 2nd derivative
    //      for improved accuracy, but this would be inconsistent with
    //      the atomic integration scheme
    explicit_1(displacement_,velocity_,dt);
    explicit_1(nodalAtomicDisplacement_,nodalAtomicVelocity_,dt);
  }

  //--------------------------------------------------------
  //  pre_final_integrate1
  //    first time integration computations 
  //    before Verlet step 2
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerletFiltered::pre_final_integrate1(double dt)
  {
    // Compute MD contribution to FEM equation
    DENS_MAT atomicForces;
    atcTransfer_->compute_atomic_force(atomicForces,atcTransfer_->get_f());

    // apply time filtering to instantaneous atomic force for FE equation
    // NOTE would an explicit-implicit time filter be more accurate?
    atcTransfer_->restrict_volumetric_quantity(atomicForces,forceFilteringData_);
    timeFilter_->apply_post_step2(nodalAtomicForce_,forceFilteringData_,dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate1
  //    second time integration computations
  //    after Verlet step 2
  //--------------------------------------------------------
  void ElasticTimeIntegratorVerletFiltered::post_final_integrate1(double dt)
  {
    DENS_MAT velocityRoc(velocityRhs_.nRows(),velocityRhs_.nCols());
    atcTransfer_->apply_inverse_mass_matrix(velocityRhs_,
                                            acceleration_,
                                            VELOCITY);
    explicit_1(velocity_,acceleration_,.5*dt);

    atcTransfer_->apply_inverse_md_mass_matrix(nodalAtomicForce_,
                                               nodalAtomicAcceleration_,
                                               VELOCITY);
    explicit_1(nodalAtomicVelocity_,nodalAtomicAcceleration_,.5*dt);
  }

};
