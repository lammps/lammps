/** ElasticTimeIntegrator : a class to implement various elasticity integrators for FE quantities */

#ifndef ELASTIC_TIME_INTEGRATOR_H
#define ELASTIC_TIME_INTEGRATOR_H

// ATC_Transfer headers
#include "TimeIntegrator.h"

using namespace std;
namespace ATC {

    // forward declarations
    class ElasticIntegrationMethod;

    /**
     *  @class  ElasticTimeIntegrator
     *  @brief  Base class fo various time integrators for elasticity FE quantities
     */

    //--------------------------------------------------------
    //--------------------------------------------------------
    //  Class ElasticTimeIntegrator
    //     Base class for elastic integrators which handles
    //     parsing and stores basic data structures
    //--------------------------------------------------------
    //--------------------------------------------------------

    class ElasticTimeIntegrator : public TimeIntegrator {
  
    public:
  
        // constructor
        ElasticTimeIntegrator(ATC_Transfer * atcTransfer,
                              TimeIntegrationType timeIntegrationType);
        
        // destructor
        virtual ~ElasticTimeIntegrator(){};

        /** parser/modifier */
        virtual bool modify(int narg, char **arg);
        
        /** pre time integration */
        virtual void initialize();

    private:

        // DO NOT define this
        ElasticTimeIntegrator();
  
    };

    /**
     *  @class  ElasticIntegrationMethod
     *  @brief  Base class fo various time integration methods for elasticity FE quantities
     */

    //--------------------------------------------------------
    //--------------------------------------------------------
    //  Class ElasticIntegrationMethod
    //     Base class for elastic integration methods which
    //     update the FE quantities in time
    //--------------------------------------------------------
    //--------------------------------------------------------

    class ElasticIntegrationMethod : public TimeIntegrationMethod {
  
    public:
  
        // constructor
        ElasticIntegrationMethod(ElasticTimeIntegrator * elasticTimeIntegrator);
        
        // destructor
        virtual ~ElasticIntegrationMethod(){};
        
    protected:
 
        /** time filtering application object */
        TimeFilter * timeFilter_;

        /** finite element displacement field */
        DENS_MAT & displacement_;
        /** finite element velocity field */
        DENS_MAT & velocity_;
        /** finite element acceleration field */
        DENS_MAT & acceleration_;

        /** atomic nodal displacement field */
        DENS_MAT & nodalAtomicDisplacement_;
        /** atomic nodal velocity field */
        DENS_MAT & nodalAtomicVelocity_;

        /** right-hand side of velocity equation */
        DENS_MAT & velocityRhs_;

        /** force at nodes from atomic quantities */
        DENS_MAT & nodalAtomicForce_;

        /** filtered power for computation during equilibration */
        DENS_MAT & forceFilteringData_;

    private:

        // DO NOT define this
        ElasticIntegrationMethod();
  
    };

    /**
     *  @class  ElasticTimeIntegratorVerlet
     *  @brief  Verlet integration for FE elastic quantities
     */

    //--------------------------------------------------------
    //--------------------------------------------------------
    //  Class ElasticTimeIntegratorVerlet
    //     Uses the second order Verlet integration to update
    //     the finite element velocity and displacement
    //     fields, i.e. the same integration used for the
    //     atomic velocities and positions.
    //--------------------------------------------------------
    //--------------------------------------------------------

    class ElasticTimeIntegratorVerlet : public ElasticIntegrationMethod {
  
    public:
  
        // constructor
        ElasticTimeIntegratorVerlet(ElasticTimeIntegrator * elasticTimeIntegrator);
        
        // destructor
        virtual ~ElasticTimeIntegratorVerlet(){};
        
        // time step methods, corresponding to ATC_Transfer
        
        /** first part of mid_initial_integrate */
        virtual void mid_initial_integrate1(double dt);
        
        /** first part of post_initial_integrate */
        virtual void post_initial_integrate1(double dt);
        
        /** first part of pre_final_integrate */
        virtual void pre_final_integrate1(double dt);

        /** first part of post_final_integrate */
        virtual void post_final_integrate1(double dt);

        /** post processing step before output */
        virtual void post_process(double dt); 
        
    private:
  
        // DO NOT define this
        ElasticTimeIntegratorVerlet();
  
    };

    /**
     *  @class  ElasticTimeIntegratorVerlet
     *  @brief  Verlet integration for FE elastic quantities with time filtering
     */

    //--------------------------------------------------------
    //--------------------------------------------------------
    //  Class ElasticTimeIntegratorVerletFiltered
    //--------------------------------------------------------
    //--------------------------------------------------------

    class ElasticTimeIntegratorVerletFiltered : public ElasticTimeIntegratorVerlet {
  
    public:
  
        // constructor
        ElasticTimeIntegratorVerletFiltered(ElasticTimeIntegrator * elasticTimeIntegrator);
        
        // destructor
        virtual ~ElasticTimeIntegratorVerletFiltered(){};
        
        // time step methods, corresponding to ATC_Transfer
        
        /** first part of mid_initial_integrate */
        virtual void mid_initial_integrate1(double dt);
        
        /** first part of post_initial_integrate */
        virtual void post_initial_integrate1(double dt);
        
        /** first part of pre_final_integrate */
        virtual void pre_final_integrate1(double dt);
        
        /** first part of post_final_integrate */
        virtual void post_final_integrate1(double dt);

        /** post processing step before output */
        virtual void post_process(double dt){};
        
    protected:

        /** atomic nodal acceleration field */
        DENS_MAT & nodalAtomicAcceleration_;

    private:

        // DO NOT define this
        ElasticTimeIntegratorVerletFiltered();
  
    };

};

#endif
