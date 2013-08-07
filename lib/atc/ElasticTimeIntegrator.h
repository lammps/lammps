#ifndef ELASTIC_TIME_INTEGRATOR_H
#define ELASTIC_TIME_INTEGRATOR_H

/** MomentumTimeIntegrator : a class to implement various elasticity integrators for FE quantities */

// ATC headers
#include "TimeIntegrator.h"

using namespace std;
namespace ATC {

    // forward declarations
    class MomentumIntegrationMethod;

    /**
     *  @class  MomentumTimeIntegrator
     *  @brief  Base class for various time integrators for elasticity FE quantities which handles parsing and stores basic data structures
     */

    class MomentumTimeIntegrator : public TimeIntegrator {
  
    public:
  
        // constructor
        MomentumTimeIntegrator(ATC_Coupling * atc,
                               TimeIntegrationType timeIntegrationType);
        
        // destructor
        virtual ~MomentumTimeIntegrator(){};

        /** parser/modifier */
        virtual bool modify(int narg, char **arg);
        
        /** create objects to implement requested numerical method */
        virtual void construct_methods();

        /** pack persistent fields */
        virtual void pack_fields(RESTART_LIST & data);

        // Member data access
        /** access for filtered atomic force */
        DENS_MAN & nodal_atomic_force_filtered(){return nodalAtomicForceFiltered_;};

        /** access for filtered atomic momentum */
        // note:  nodalAtomicMomentum_ should always be reset as it tracks the original momentum + MD evolution
        DENS_MAN & nodal_atomic_momentum_filtered(){return nodalAtomicMomentumFiltered_;};


    protected:

        /** filtered atomic force */
        
        DENS_MAN nodalAtomicForceFiltered_;

        /** filtered atomic momentum due initial conditions and MD updates */
        DENS_MAN nodalAtomicMomentumFiltered_;


    private:

        // DO NOT define this
        MomentumTimeIntegrator();
  
    };

    /**
     *  @class  MomentumIntegrationMethod
     *  @brief  Base class for various time integration methods for the momentum equation
     */

    class MomentumIntegrationMethod : public TimeIntegrationMethod {
  
    public:
  
        // constructor
        MomentumIntegrationMethod(MomentumTimeIntegrator * momentumTimeIntegrator);
        
        // destructor
        virtual ~MomentumIntegrationMethod(){};

        /** create and get necessary transfer operators */
        virtual void construct_transfers();

        /** checks to see if first RHS computation is needed */
        virtual bool has_final_predictor() {return true;};
        
    protected:
 
        /** time filtering object */
        TimeFilter * timeFilter_;
        
        /** finite element velocity field */
        DENS_MAN & velocity_;
        /** finite element acceleration field */
        DENS_MAN & acceleration_;
        
        /** atomic nodal velocity field */
        DENS_MAN & nodalAtomicVelocityOut_;
        /** right-hand side of velocity equation */
        DENS_MAN & velocityRhs_;
        /** force at nodes from atomic quantities */
        DENS_MAN & nodalAtomicForceOut_;

        /** transfer for computing nodal atomic velocity */
        DENS_MAN * nodalAtomicVelocity_;

    private:

        // DO NOT define this
        MomentumIntegrationMethod();
  
    };

    /**
     *  @class  ElasticTimeIntegratorVerlet
     *  @brief  Verlet integration for FE elastic quantities.  Uses the second order Verlet integration to update the finite element velocity and displacement fields, i.e. the same integration used for the atomic velocities and positions.
     */

    class ElasticTimeIntegratorVerlet : public MomentumIntegrationMethod {
  
    public:
  
      // constructor
      ElasticTimeIntegratorVerlet(MomentumTimeIntegrator * momentumTimeIntegrator);
        
      // destructor
      virtual ~ElasticTimeIntegratorVerlet(){};

      /** create and get necessary transfer operators */
      virtual void construct_transfers();
      
      /** pre time integration initialization of data */
      virtual void initialize();
      
      // time step methods, corresponding to ATC_Transfer
      /** first part of mid_initial_integrate */
      virtual void mid_initial_integrate1(double dt);
      /** first part of post_initial_integrate */
      virtual void post_initial_integrate1(double dt);
      /** first part of pre_final_integrate */
      virtual void pre_final_integrate1(double dt);
      /** second part of post_final_integrate */
      virtual void post_final_integrate2(double dt);
      /** adds any contributions from time integrator to RHS */
      virtual void add_to_rhs();
      /** post processing step before output */
      virtual void post_process();
        
      /** add output data */
      virtual void output(OUTPUT_LIST & outputData);
      
      /** operations at end of a run */
      virtual void finish();
      
    protected:
      
      /** finite element displacement field */
      DENS_MAN & displacement_;
      
      /** atomic nodal displacement field */
      DENS_MAN & nodalAtomicDisplacementOut_;
      
      /** filtered atomic force */
      DENS_MAN & nodalAtomicForceFiltered_;
      
      /** transfer for computing atomic displacement */
      DENS_MAN * nodalAtomicDisplacement_;

      /** transfer for computing nodal atomic force */
      DENS_MAN * nodalAtomicForce_;
      
    private:
      
      // DO NOT define this
      ElasticTimeIntegratorVerlet();
      
    };
    
    /**
     *  @class  ElasticTimeIntegratorVerlet
     *  @brief  Verlet integration for FE elastic quantities with time filtering
     */

    class ElasticTimeIntegratorVerletFiltered : public ElasticTimeIntegratorVerlet {
  
    public:
  
        // constructor
        ElasticTimeIntegratorVerletFiltered(MomentumTimeIntegrator * momentumTimeIntegrator);
        
        // destructor
        virtual ~ElasticTimeIntegratorVerletFiltered(){};
        
        // time step methods, corresponding to ATC_Transfer
        /** first part of mid_initial_integrate */
        virtual void mid_initial_integrate1(double dt);
        /** first part of post_initial_integrate */
        virtual void post_initial_integrate1(double dt);
        /** second part of post_final_integrate */
        virtual void post_final_integrate2(double dt);
        /** adds any contributions from time integrator to RHS */
        virtual void add_to_rhs();
        /** post processing step before output */
        virtual void post_process(){};

        /** add output data */
        virtual void output(OUTPUT_LIST & outputData);
        
    protected:

        /** atomic nodal acceleration field */
        DENS_MAN & nodalAtomicAcceleration_;

    private:

        // DO NOT define this
        ElasticTimeIntegratorVerletFiltered();
  
    };

  /**
   *  @class  ElasticTimeIntegratorFractionalStep 
   *  @brief  Class for using 2nd order Verlet integration to update FE contributions to momentum field
   *          (Uses same update for the atomic contributions to the finite 
   *           elements as are used by the LAMMPS integration scheme 
   *           for the atomic velocities and positions, i.e. Verlet.)
   */ 

  class ElasticTimeIntegratorFractionalStep : public MomentumIntegrationMethod {

  public:

    // constructor
    ElasticTimeIntegratorFractionalStep(MomentumTimeIntegrator * momentumTimeIntegrator);
        
    // destructor
    virtual ~ElasticTimeIntegratorFractionalStep() {};
    
    /** create and get necessary transfer operators */
    virtual void construct_transfers();

    /** pre time integration initialization of data */
    virtual void initialize();
        
    // time step methods, corresponding to ATC_Transfer
    /** first part of pre_initial_integrate */
    virtual void pre_initial_integrate1(double dt);
    /** second part of mid_initial_integrate */
    virtual void pre_initial_integrate2(double dt);
    /** first part of mid_initial_integrate */
    virtual void mid_initial_integrate1(double dt);
    /** first part of post_initial_integrate */
    virtual void post_initial_integrate1(double dt);
    /** second part of post_final_integrate */
    virtual void post_final_integrate2(double dt);
    /** post processing step before output */
    virtual void post_process();

    /** finalize state of some unfiltered variables */
    virtual void finish();

    /** add output data */
    virtual void output(OUTPUT_LIST & outputData);

  protected:

    // methods
    /** compute old energy and temperature for use in time integrators */
    virtual void compute_old_time_data();

    /** computes temperature change associated with atomic energy change */
    virtual void compute_velocity_delta(const DENS_MAT & atomicMomentumDelta,
                                        double dt);

    // data
    /** finite element displacement field */
    DENS_MAN & displacement_;
      
    /** atomic nodal displacement field */
    DENS_MAN & nodalAtomicDisplacementOut_;

    /** equivalent nodal force due to atomic momentum change */
    DENS_MAT nodalAtomicForce_;
      
    /** filtered atomic force */
    DENS_MAN & nodalAtomicForceFiltered_;

    /** transfer for computing atomic momentum */
    DENS_MAN * nodalAtomicMomentum_;

    /** filtered atomic momentum */
    DENS_MAN & nodalAtomicMomentumFiltered_;
      
    /** transfer for computing atomic displacement */
    DENS_MAN * nodalAtomicDisplacement_;

    /** change in FE velocity due to atomic motions */
    DENS_MAN atomicVelocityDelta_;

    /** restricted atomic momentum from previous time step */
    DENS_MAT nodalAtomicMomentumOld_;

    /** FE atomic velocity contribution from previous time step */
    DENS_MAT nodalAtomicVelocityOld_;


  private:

    // DO NOT define this
    ElasticTimeIntegratorFractionalStep();

  };

  /**
   *  @class  FluidsTimeIntegratorGear
   *  @brief  Class for using 3rd order Gear integration to update FE contributions to momentum field
   *          and fractional step method for atomic contributions
   */ 

  class FluidsTimeIntegratorGear : public MomentumIntegrationMethod {

  public:

    // constructor
    FluidsTimeIntegratorGear(MomentumTimeIntegrator * MomentumTimeIntegrator);
        
    // destructor
    virtual ~FluidsTimeIntegratorGear() {};
    
    /** create and get necessary transfer operators */
    virtual void construct_transfers();

    /** pre time integration initialization of data */
    virtual void initialize();
        
    // time step methods, corresponding to ATC_Transfer
    /** first part of pre_initial_integrate */
    virtual void pre_initial_integrate1(double dt);
    /** second part of pre_initial_integrate */
    virtual void pre_initial_integrate2(double dt);
    /** first part of pre_final_integrate */
    virtual void pre_final_integrate1(double dt);
    /** first part of post_final_integrate */
    virtual void post_final_integrate1(double dt);
    /** second part of post_final_integrate */
    virtual void post_final_integrate2(double dt);

    /** post processing step before output */
    virtual void post_process();

    /** finalize state of some unfiltered variables */
    virtual void finish();

    /** add output data */
    virtual void output(OUTPUT_LIST & outputData);

  protected:

    // methods
    /** applies Gear predictor */
    virtual void apply_gear_predictor(double dt);

    /** applies Gear corrector */
    virtual void apply_gear_corrector(const DENS_MAT & R_theta,
                                      double dt);

    /** compute old energy and temperature for use in time integrators */
    virtual void compute_old_time_data();

    /** computes temperature change associated with atomic energy change */
    virtual void compute_velocity_delta(const DENS_MAT & atomicMomentumDelta,
                                        double dt);

    // data
    /** equivalent nodal force due to atomic momentum change */
    DENS_MAT nodalAtomicForce_;
      
    /** filtered atomic force */
    DENS_MAN & nodalAtomicForceFiltered_;

    /** transfer for computing atomic momentum */
    DENS_MAN * nodalAtomicMomentum_;

    /** filtered atomic momentum */
    DENS_MAN & nodalAtomicMomentumFiltered_;

    /** change in FE velocity due to atomic motions */
    DENS_MAN atomicVelocityDelta_;

    /** restricted atomic momentum from previous time step */
    DENS_MAT nodalAtomicMomentumOld_;

    /** FE atomic velocity contribution from previous time step */
    DENS_MAT nodalAtomicVelocityOld_;

    /** finite element velocity 2nd time derivative */
    DENS_MAN & velocity2Roc_;

    // workspace for gear integration
    DENS_MAT _velocityResidual_;

  private:

    // DO NOT define this
    FluidsTimeIntegratorGear();

    };

};

#endif

