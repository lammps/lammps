#ifndef TIME_INTEGRATOR_H
#define TIME_INTEGRATOR_H

// ATC_Transfer headers
#include "MatrixLibrary.h"
#include "TimeFilter.h"
#include "ATC_TypeDefs.h"

using namespace std;
namespace ATC {

  // forward declarations
  class ATC_Transfer;
  class TimeIntegrationMethod;

  /**
   *  @class  TimeIntegrator
   *  @brief  Base class fo various time integrators for FE quantities
   */
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeIntegrator
  //--------------------------------------------------------
  //--------------------------------------------------------

  class TimeIntegrator {
  
  public:

    /** types of time integration */
    enum TimeIntegrationType {
      STEADY,
      VERLET,
      GEAR,
      FRACTIONAL_STEP,
      EXPLICIT,
      IMPLICIT,
      CRANK_NICOLSON
    };
      
    // constructor
    TimeIntegrator(ATC_Transfer * atcTransfer,
                   TimeIntegrationType timeIntegrationType = STEADY);
        
    // destructor
    virtual ~TimeIntegrator();
        
    /** parser/modifier */
    virtual bool modify(int narg, char **arg){return false;};
        
    /** pre time integration */
    virtual void initialize(){needReset_ = false;};

    /** flag if reset is needed */
    bool need_reset() {return needReset_;};
        
    // time step methods, corresponding to ATC_Transfer
    /** first part of pre_initial_integrate */
    virtual void pre_initial_integrate1(double dt);
    /** second part of pre_initial_integrate */
    virtual void pre_initial_integrate2(double dt);
        
    /** first part of mid_initial_integrate */
    virtual void mid_initial_integrate1(double dt);
    /** second part of mid_initial_integrate */
    virtual void mid_initial_integrate2(double dt);
        
    /** first part of post_initial_integrate */
    virtual void post_initial_integrate1(double dt);
    /** second part of post_initial_integrate */
    virtual void post_initial_integrate2(double dt);
        
    /** first part of pre_final_integrate */
    virtual void pre_final_integrate1(double dt);
    /** second part of pre_final_integrate */
    virtual void pre_final_integrate2(double dt);
        
    /** first part of post_final_integrate */
    virtual void post_final_integrate1(double dt);
    /** second part of post_final_integrate */
    virtual void post_final_integrate2(double dt);

    /** post processing step */
    virtual void post_process(double dt);
    /** add output data */
    virtual void output(OUTPUT_LIST & outputData);

    // Member data access
    /** access to time integration type */
    TimeIntegrationType get_time_integration_type() const
    { return timeIntegrationType_; };

    /** access to ATC Transfer object */
    ATC_Transfer * get_atc_transfer() {return atcTransfer_;};

    /** access to time filter object */
    TimeFilter * get_time_filter() {return timeFilter_;};

    /** access to time filter manager object */
    TimeFilterManager * get_time_filter_manager() {return timeFilterManager_;};

    /** force the integrator to be reset */
    void force_reset() {needReset_ = true;};

    /** force the integrator not to be reset */
    void force_no_reset() {needReset_ = false;};

  protected:

    /** pointer to time integrator method */
    TimeIntegrationMethod * timeIntegrationMethod_;
        
    /** pointer to access ATC methods */
    ATC_Transfer * atcTransfer_;
        
    /** time filter for specific updates */
    TimeFilter * timeFilter_;

    /** time filter manager for getting time filtering info */
    TimeFilterManager * timeFilterManager_;
        
    /** type of integration scheme being used */
    TimeIntegrationType timeIntegrationType_;
  
    /** flat to reset data */
    bool needReset_;

  private:

    // DO NOT define this
    TimeIntegrator();

  };

  /**
   *  @class  TimeIntegrationMethod
   *  @brief  Base class fo various time integration methods
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeIntegrationMethod
  //     Base class for time integration methods which
  //     update the FE quantities
  //--------------------------------------------------------
  //--------------------------------------------------------

  class TimeIntegrationMethod {
  
  public:
  
    // constructor
    TimeIntegrationMethod(TimeIntegrator * timeIntegrator);
        
    // destructor
    virtual ~TimeIntegrationMethod(){};
        
    // time step methods, corresponding to ATC_Transfer and TimeIntegrator
    /** first part of pre_initial_integrate */
    virtual void pre_initial_integrate1(double dt){};
    /** second part of pre_initial_integrate */
    virtual void pre_initial_integrate2(double dt){};
        
    /** first part of mid_initial_integrate */
    virtual void mid_initial_integrate1(double dt){};
    /** second part of mid_initial_integrate */
    virtual void mid_initial_integrate2(double dt){};
        
    /** first part of post_initial_integrate */
    virtual void post_initial_integrate1(double dt){};
    /** second part of post_initial_integrate */
    virtual void post_initial_integrate2(double dt){};
        
    /** first part of pre_final_integrate */
    virtual void pre_final_integrate1(double dt){};
    /** second part of pre_final_integrate */
    virtual void pre_final_integrate2(double dt){};
        
    /** first part of post_final_integrate */
    virtual void post_final_integrate1(double dt){};
    /** second part of post_final_integrate */
    virtual void post_final_integrate2(double dt){};

    /** post processing step */
    virtual void post_process(double dt){};
    /** add output data */
    virtual void output(OUTPUT_LIST & outputData){};
        
  protected:

    /** owning time integrator */
    TimeIntegrator * timeIntegrator_;

    /** associated ATC transfer object */
    ATC_Transfer * atcTransfer_;

  private:

    // DO NOT define this
    TimeIntegrationMethod();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  time integration functions not associated
  //  with any particular class
  //--------------------------------------------------------
  //--------------------------------------------------------

  static void gear1_4_predict(MATRIX & f,
                              MATRIX & dot_f,
                              MATRIX & ddot_f,
                              MATRIX & dddot_f,
                              double dt)
  // 4th order Gear integrator for 1rst order ODE predictor step
  {
    f      = f + dot_f*dt + ddot_f*(1./2.*dt*dt) + dddot_f*(1./6.*dt*dt*dt);
    dot_f  = dot_f + ddot_f*dt+dddot_f*(1./2.*dt*dt);
    ddot_f = ddot_f + dddot_f*dt;
  };

  static void gear1_3_predict(MATRIX & f,
                              MATRIX & dot_f,
                              MATRIX & ddot_f,
                              double dt)
  // 3rd order Gear integrator for 1rst order ODE predictor step
  {
    f      = f + dot_f*dt + ddot_f*(1./2.*dt*dt);
    dot_f  = dot_f + ddot_f*dt;
  };

  static void gear1_4_correct(MATRIX & f,
                              MATRIX & dot_f,
                              MATRIX & ddot_f,
                              MATRIX & dddot_f,
                              const MATRIX & R_f,
                              double dt)
  // 4th order Gear integrator for 1rst order ODE corrector step
  {
    f       = f       + (3./8.)*R_f;
    dot_f   = dot_f   + (1./dt)*R_f;
    ddot_f  = ddot_f  + (3./2./dt/dt)*R_f;
    dddot_f = dddot_f + (1./dt/dt/dt)*R_f;
  };

  static void gear1_3_correct(MATRIX & f,
                              MATRIX & dot_f,
                              MATRIX & ddot_f,
                              const MATRIX & R_f,
                              double dt)
  // 3rd order Gear integrator for 1rst order ODE corrector step
  {
    f      = f      + (5./12.)*R_f;
    dot_f  = dot_f  + (1./dt)*R_f;
    ddot_f = ddot_f + (1./dt/dt)*R_f;
  };
  
  static void explicit_1(MATRIX & f,
                         MATRIX & dot_f,
                         double dt)
  // 1rst order explict ODE update
  {
    f = f + dt*dot_f;
  };

  static void explicit_2(MATRIX & f,
                         MATRIX & dot_f,
                         MATRIX & ddot_f,
                         double dt)
  // 2nd order explict ODE update
  {
    f = f + dt*dot_f + .5*dt*dt*ddot_f;
  };

};

#endif

