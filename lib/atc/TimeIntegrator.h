#ifndef TIME_INTEGRATOR_H
#define TIME_INTEGRATOR_H

#include "MatrixLibrary.h"
#include "TimeFilter.h"
#include "ATC_TypeDefs.h"

namespace ATC {

  // forward declarations
  class ATC_Method;
  class ATC_Coupling;
  class TimeIntegrationMethod;

  /**
   *  @class  AtomTimeIntegrator
   *  @brief  Base class for various time integrators for atomic quantities (replacing other lammps fixes)
   */

  class AtomTimeIntegrator {
  
  public:

    // constructor
    AtomTimeIntegrator(){};

    // destructor
    virtual ~AtomTimeIntegrator(){};

    /** create and get necessary transfer operators */
    virtual void construct_transfers(){};
        
   /** Predictor phase, Verlet first step for velocity */
    virtual void init_integrate_velocity(double dt){};

    /** Predictor phase, Verlet first step for position */
    virtual void init_integrate_position(double dt){};

    /** Corrector phase, Verlet second step for velocity */
    virtual void final_integrate(double dt){};

  };

  /**
   *  @class  AtomTimeIntegratorType
   *  @brief  class for applying velocity-verlet based on atom type
   */

  class AtomTimeIntegratorType : public AtomTimeIntegrator {
  
  public:

    // constructor
    AtomTimeIntegratorType(ATC_Method * atc, AtomType atomType);

    // destructor
    virtual ~AtomTimeIntegratorType(){};

    /** create and get necessary transfer operators */
    virtual void construct_transfers();
        
   /** Predictor phase, Verlet first step for velocity */
    virtual void init_integrate_velocity(double dt);

    /** Predictor phase, Verlet first step for position */
    virtual void init_integrate_position(double dt);

    /** Corrector phase, Verlet second step for velocity */
    virtual void final_integrate(double dt);

  protected:

    /** pointer to atc object */
    ATC_Method * atc_;

    /** atom type this is applied to */
    AtomType atomType_;

    /** atomic masses */
    DENS_MAN * mass_;

    /** atomic positions */
    DENS_MAN * position_;

    /** atomic velocities */
    DENS_MAN * velocity_;

    /** atomic forces */
    DENS_MAN * force_;

    // workspace
    DENS_MAT _deltaQuantity_;

  private:

    // DO NOT define this
    AtomTimeIntegratorType();

  };

  /**
   *  @class  TimeIntegrator
   *  @brief  Base class for various time integrators for FE quantities
   */
  
  class TimeIntegrator {
  
  public:

    /** types of time integration */
    enum TimeIntegrationType {
      NONE=0,
      STEADY,
      VERLET,
      GEAR,
      FRACTIONAL_STEP,
      EXPLICIT,
      IMPLICIT,
      CRANK_NICOLSON,
      DIRECT
    };
      
    // constructor
    TimeIntegrator(ATC_Coupling * atc,
                   TimeIntegrationType timeIntegrationType = STEADY);
        
    // destructor
    virtual ~TimeIntegrator();
        
    /** parser/modifier */
    virtual bool modify(int narg, char **arg){return false;};

    /** create objects to implement requested numerical method */
    virtual void construct_methods() = 0;

    /** create and get necessary transfer operators */
    virtual void construct_transfers();

    /** pre time integration initialization of data */
    virtual void initialize();

    /** flag if reset is needed */
    bool need_reset() const {return needReset_;};
        
    // time step methods, corresponding to ATC_Coupling
    /** first part of pre_initial_integrate */
    virtual void pre_initial_integrate1(double dt);
    /** second part of pre_initial_integrate */
    virtual void pre_initial_integrate2(double dt);

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
    /** third part of post_final_integrate */
    virtual void post_final_integrate3(double dt);

    /** checks to see if first RHS computation is needed */
    virtual bool has_final_predictor();
    /** checks to see if second RHS computation is needed */
    virtual bool has_final_corrector();

    /** adds any contributions from time integrator to RHS */
    virtual void add_to_rhs();
    /** post processing step prior to output */
    virtual void post_process();
    /** add output data */
    virtual void output(OUTPUT_LIST & outputData);
    /** pack persistent fields */
    virtual void pack_fields(RESTART_LIST & data);

    /** finalize any data */
    virtual void finish();

    // Member data access
    /** access to time integration type */
    TimeIntegrationType time_integration_type() const
    { return timeIntegrationType_; };

    /** access to ATC Transfer object */
    ATC_Coupling * atc() {return atc_;};

    /** access to time filter object */
    TimeFilter * time_filter() {return timeFilter_;};

    /** access to time filter manager object */
    TimeFilterManager * time_filter_manager() {return timeFilterManager_;};

    /** force the integrator to be reset */
    void force_reset() {needReset_ = true;};

    /** force the integrator not to be reset */
    void force_no_reset() {needReset_ = false;};

  protected:

    /** pointer to time integrator method */
    TimeIntegrationMethod * timeIntegrationMethod_;
        
    /** pointer to access ATC methods */
    ATC_Coupling * atc_;
        
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
   *  @brief  Base class for time integration methods which update FE quantities
   */

  class TimeIntegrationMethod {
  
  public:
  
    // constructor
    TimeIntegrationMethod(TimeIntegrator * timeIntegrator);
        
    // destructor
    virtual ~TimeIntegrationMethod(){};

    /** create and get necessary transfer operators */
    virtual void construct_transfers(){};
    /** pre time integration */
    virtual void initialize(){};
        
    // time step methods, corresponding to ATC_Coupling and TimeIntegrator
    /** first part of pre_initial_integrate */
    virtual void pre_initial_integrate1(double dt){};
    /** second part of pre_initial_integrate */
    virtual void pre_initial_integrate2(double dt){};

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
    /** third part of post_final_integrate */
    virtual void post_final_integrate3(double dt){};

    /** checks to see if first RHS computation is needed */
    virtual bool has_final_predictor() {return false;};
    /** checks to see if second RHS computation is needed */
    virtual bool has_final_corrector() {return false;};

    /** adds any contributions from time integrator to RHS */
    virtual void add_to_rhs() {};
    /** post processing step */
    virtual void post_process(){};
    /** add output data */
    virtual void output(OUTPUT_LIST & outputData){};
    /** pack persistent fields */
    virtual void pack_fields(RESTART_LIST & data){};

    /** finalize any states */
    virtual void finish(){};
        
  protected:

    /** owning time integrator */
    TimeIntegrator * timeIntegrator_;

    /** associated ATC transfer object */
    ATC_Coupling * atc_;

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

  inline void gear1_4_predict(MATRIX & f,
                              MATRIX & dot_f,
                              MATRIX & ddot_f,
                              const MATRIX & dddot_f,
                              double dt)
  // 4th order Gear integrator for 1rst order ODE predictor step
  {
    f      = f + dot_f*dt + ddot_f*(1./2.*dt*dt) + dddot_f*(1./6.*dt*dt*dt);
    dot_f  = dot_f + ddot_f*dt+dddot_f*(1./2.*dt*dt);
    ddot_f = ddot_f + dddot_f*dt;
  };

  inline void gear1_3_predict(MATRIX & f,
                              MATRIX & dot_f,
                              const MATRIX & ddot_f,
                              double dt)
  // 3rd order Gear integrator for 1rst order ODE predictor step
  {
    f      = f + dot_f*dt + ddot_f*(1./2.*dt*dt);
    dot_f  = dot_f + ddot_f*dt;
  };

  inline void gear1_4_correct(MATRIX & f,
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

  inline void gear1_3_correct(MATRIX & f,
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
  
  inline void explicit_1(MATRIX & f,
                         const MATRIX & dot_f,
                         double dt)
  // 1rst order explict ODE update
  {
    f = f + dt*dot_f;
  };

  inline void explicit_2(MATRIX & f,
                         const MATRIX & dot_f,
                         const MATRIX & ddot_f,
                         double dt)
  // 2nd order explict ODE update
  {
    f = f + dt*dot_f + .5*dt*dt*ddot_f;
  };

};

#endif

