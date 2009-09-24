// ATC transfer headers
#include "TimeIntegrator.h"
#include "ATC_Transfer.h"
#include "ATC_Error.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeIntegrator
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  TimeIntegrator::TimeIntegrator(ATC_Transfer * atcTransfer,
				 TimeIntegrationType timeIntegrationType) :
    atcTransfer_(atcTransfer),
    timeIntegrationMethod_(NULL),
    timeFilter_(NULL),
    timeFilterManager_(atcTransfer_->get_time_filter_manager()),
    timeIntegrationType_(timeIntegrationType),
    needReset_(true)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  TimeIntegrator::~TimeIntegrator()
  {
    if (timeFilter_)
      delete timeFilter_;
    
    if (timeIntegrationMethod_)
      delete timeIntegrationMethod_;
  }

    //--------------------------------------------------------
  //  pre_initial_integrate1
  //    first time integration computations
  //    before Verlet step 1
  //--------------------------------------------------------
  void TimeIntegrator::pre_initial_integrate1(double dt)
  {
    timeIntegrationMethod_->pre_initial_integrate1(dt);
  }

  //--------------------------------------------------------
  //  pre_initial_integrate2
  //    second time integration computations
  //    before Verlet step 1
  //--------------------------------------------------------
  void TimeIntegrator::pre_initial_integrate2(double dt)
  {
    timeIntegrationMethod_->pre_initial_integrate2(dt);
  }

  //--------------------------------------------------------
  //  mid_initial_integrate1
  //    first time integration computations
  //    at the mid-point of Verlet step 1
  //--------------------------------------------------------
  void TimeIntegrator::mid_initial_integrate1(double dt)
  {
    timeIntegrationMethod_->mid_initial_integrate1(dt);
  }

  //--------------------------------------------------------
  //  mid_initial_integrate2
  //    second time integration computations
  //    at the mid-point of Verlet step 1
  //--------------------------------------------------------
  void TimeIntegrator::mid_initial_integrate2(double dt)
  {
    timeIntegrationMethod_->mid_initial_integrate2(dt);
  }

  //--------------------------------------------------------
  //  post_initial_integrate1
  //    first time integration computations
  //    after Verlet step 1
  //--------------------------------------------------------
  void TimeIntegrator::post_initial_integrate1(double dt)
  {
    timeIntegrationMethod_->post_initial_integrate1(dt);
  }

  //--------------------------------------------------------
  //  post_initial_integrate2
  //    second time integration computations
  //    after Verlet step 1
  //--------------------------------------------------------
  void TimeIntegrator::post_initial_integrate2(double dt)
  {
    timeIntegrationMethod_->post_initial_integrate2(dt);
  }

  //--------------------------------------------------------
  //  pre_final_integrate1
  //    first time integration computations 
  //    before Verlet step 2
  //--------------------------------------------------------
  void TimeIntegrator::pre_final_integrate1(double dt)
  {
    timeIntegrationMethod_->pre_final_integrate1(dt);
  }

  //--------------------------------------------------------
  //  pre_final_integrate2
  //    second time integration computations
  //    before Verlet step 2
  //--------------------------------------------------------
  void TimeIntegrator::pre_final_integrate2(double dt)
  {
    timeIntegrationMethod_->pre_final_integrate2(dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate1
  //    first time integration computations 
  //    after Verlet step 2
  //--------------------------------------------------------
  void TimeIntegrator::post_final_integrate1(double dt)
  {
    timeIntegrationMethod_->post_final_integrate1(dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate2
  //    second time integration computations
  //    after Verlet step 2
  //--------------------------------------------------------
  void TimeIntegrator::post_final_integrate2(double dt)
  {
    timeIntegrationMethod_->post_final_integrate2(dt);
  }

  //--------------------------------------------------------
  //  post_process
  //    perform any post processing calculations
  //--------------------------------------------------------
  void TimeIntegrator::post_process(double dt)
  {
    timeIntegrationMethod_->post_process(dt);
  }

  //--------------------------------------------------------
  //  output
  //    add variables to output list
  //--------------------------------------------------------
  void TimeIntegrator::output(OUTPUT_LIST & outputData)
  {
    timeIntegrationMethod_->output(outputData);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeIntegrationMethod
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //        Grab data from ATC
  //-------------------------------------------------------- 
  TimeIntegrationMethod::TimeIntegrationMethod(TimeIntegrator * timeIntegrator) :
    timeIntegrator_(timeIntegrator),
    atcTransfer_(timeIntegrator_->get_atc_transfer())
  {
    // do nothing
  }

};
