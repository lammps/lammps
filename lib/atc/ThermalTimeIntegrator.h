#ifndef THERMAL_TIME_INTEGRATOR_H
#define THERMAL_TIME_INTEGRATOR_H

// ATC headers
#include "TimeIntegrator.h"

namespace ATC {

  // forward declarations
  class ThermalIntegrationMethod;
  class AtfShapeFunctionRestriction;

  /**
   *  @class  ThermalTimeIntegrator
   *  @brief  Class for various time integrators for thermal FE quantities
   *          (handles parsing and stores basic data structures)
   */

  class ThermalTimeIntegrator : public TimeIntegrator {
  
  public:
  
    // constructor
    ThermalTimeIntegrator(ATC_Coupling * atc,
                          TimeIntegrationType timeIntegrationType);
        
    // destructor
    virtual ~ThermalTimeIntegrator(){};

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** create objects to implement requested numerical method */
    virtual void construct_methods();

    /** pack persistent fields */
    virtual void pack_fields(RESTART_LIST & data);

    // Member data access
    /** access for filtered atomic power */
    DENS_MAN & nodal_atomic_power_filtered(){return nodalAtomicPowerFiltered_;};
    
    /** access for filtered atomic energy */
    // note:  nodalAtomicEnergy_ should always be reset as it tracks the original energy + MD evolution
    DENS_MAN & nodal_atomic_energy_filtered(){return nodalAtomicEnergyFiltered_;};
        
  protected:

    /** filtered atomic power */
    
    DENS_MAN nodalAtomicPowerFiltered_;

    /** filtered atomic energy due initial conditions and MD updates */
    DENS_MAN nodalAtomicEnergyFiltered_;

  private:

    // DO NOT define this
    ThermalTimeIntegrator();
  
  };

  /**
   *  @class  ThermalIntegrationMethod
   *  @brief  Class for various time integration methods for thermal FE quantities
   */

  class ThermalIntegrationMethod : public TimeIntegrationMethod {
  
  public:
  
    // constructor
    ThermalIntegrationMethod(ThermalTimeIntegrator * thermalTimeIntegrator);
        
    // destructor
    virtual ~ThermalIntegrationMethod() {};

    /** create and get necessary transfer operators */
    virtual void construct_transfers();

    /** checks to see if first RHS computation is needed */
    virtual bool has_final_predictor() {return true;};

  protected:

    /** time filtering object */
    TimeFilter * timeFilter_;

    /** finite element temperature field */
    DENS_MAN & temperature_;
    /** finite element temperature Rate of change (Roc) */
    DENS_MAN & temperatureRoc_;
    /** finite element temperature 2nd time derivative */
    DENS_MAN & temperature2Roc_;

    /** atomic nodal temperature field for output */
    DENS_MAN & nodalAtomicTemperatureOut_;

    /** interscale operator for instantaneous temperature */
    DENS_MAN * nodalAtomicTemperature_;

    /** right-hand side of temperature equation */
    DENS_MAN & temperatureRhs_;

    /** finite element power from atomic quantities for output */
    DENS_MAN & nodalAtomicPowerOut_;

    /** workspace for gear integration */
    DENS_MAT _temperatureResidual_;

  private:

    // DO NOT define this
    ThermalIntegrationMethod();
  
  };

  /**
   *  @class  ThermalTimeIntegratorGear
   *  @brief  Class uses 3rd order Gear integration for time integration of FE temperature field 
   */

  class ThermalTimeIntegratorGear : public ThermalIntegrationMethod {
  
  public:
  
    // constructor
    ThermalTimeIntegratorGear(ThermalTimeIntegrator * ThermalTimeIntegrator);
        
    // destructor
    virtual ~ThermalTimeIntegratorGear() {};

    /** create and get necessary transfer operators */
    virtual void construct_transfers();

    /** pre time integration initialization of data */
    virtual void initialize();
        
    // time step methods, corresponding to ATC_Transfer
    /** second part of pre_initial_integrate */
    virtual void pre_initial_integrate2(double dt);
    /** first part of post_final_integrate */
    virtual void post_final_integrate1(double dt);

    /** parallel post-processing operations pre-output */
    virtual void post_process();

    /** finalize state of some unfiltered variables */
    virtual void finish();

  protected:

    /** filtered atomic power */
    DENS_MAN & nodalAtomicPowerFiltered_;

    /** instantaneous atomic power */
    AtfShapeFunctionRestriction * nodalAtomicPower_;

  private:
  
    // DO NOT define this
    ThermalTimeIntegratorGear();
  
  };

  /**
   *  @class  ThermalTimeIntegratorGearFiltered
   *  @brief  Gear integration for FE thermal quantities with time filtering
   */

  class ThermalTimeIntegratorGearFiltered : public ThermalTimeIntegratorGear {
  
  public:
  
    // constructor
    ThermalTimeIntegratorGearFiltered(ThermalTimeIntegrator * thermalTimeIntegrator);
        
    // destructor
    virtual ~ThermalTimeIntegratorGearFiltered(){};
        
    // time step methods, corresponding to ATC_Transfer
    /** second part of pre_initial_integrate */
    virtual void pre_initial_integrate2(double dt);
    /** first part of post_final_integrate */
    virtual void post_final_integrate1(double dt);
    /** third part of post_final_integrate */
    virtual void post_final_integrate3(double dt);

    /** parallel post-processing operations pre-output */
    virtual void post_process();
        
  protected:

    /** finite element temperature 3rd time derivative */
    DENS_MAN & temperature3Roc_;

  private:

    // DO NOT define this
    ThermalTimeIntegratorGearFiltered();
  
  };

  /**
   *  @class  ThermalTimeIntegratorFractionalStep 
   *  @brief  Class for using 3rd order Gear integration to update FE contributions to temperature field
   *          (Uses same update for the atomic contributions to the finite 
   *           elements as are used by the LAMMPS integration scheme 
   *           for the atomic velocities and positions, i.e. Verlet.)
   */ 

  class ThermalTimeIntegratorFractionalStep : public ThermalIntegrationMethod {

  public:

    // constructor
    ThermalTimeIntegratorFractionalStep(ThermalTimeIntegrator * ThermalTimeIntegrator);
        
    // destructor
    virtual ~ThermalTimeIntegratorFractionalStep() {};
    
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
    /** third part of post_final_integrate */
    virtual void post_final_integrate3(double dt);

    /** checks to see if first RHS computation is needed */
    virtual bool has_final_corrector() {return true;};

    /** post-process data */
    virtual void post_process();

    /** finalize state of some unfiltered variables */
    virtual void finish();

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
    virtual void compute_temperature_delta(const DENS_MAT & atomicEnergyDelta,
                                           double dt);

    // data
    /** filtered restricted atomic energy */
    DENS_MAN & nodalAtomicEnergyFiltered_;

    /** filtered atomic power, for post-processing only */
    DENS_MAN & nodalAtomicPowerFiltered_;

    /** change in FE temperature due to atomic motions */
    DENS_MAN atomicTemperatureDelta_;

    /** fractional step auxiliary storage for restricted atomic energy */
    DENS_MAN * nodalAtomicEnergy_;

    /** power associated with thermostat for post-processing */
    DENS_MAT nodalAtomicPower_;

    /** restricted atomic energy from previous time step */
    DENS_MAN nodalAtomicEnergyOld_;

    /** FE atomic temperature contribution from previous time step */
    DENS_MAN nodalAtomicTemperatureOld_;

  private:

    // DO NOT define this
    ThermalTimeIntegratorFractionalStep();

  };

  /**
   *  @class  ThermalTimeIntegratorFractionalStepFiltered
   *  @brief  Class for using filtered results from a 3rd order Gear integration to update FE contributions to temperature field
   */

  class ThermalTimeIntegratorFractionalStepFiltered : public ThermalTimeIntegratorFractionalStep {

  public:

    // constructor
    ThermalTimeIntegratorFractionalStepFiltered(ThermalTimeIntegrator * ThermalTimeIntegrator);
        
    // destructor
    virtual ~ThermalTimeIntegratorFractionalStepFiltered();
        
    // time step methods, corresponding to ATC_Transfer
    /** first part of pre_initial_integrate */
    virtual void pre_initial_integrate1(double dt);

    /** add output data */
    virtual void output(OUTPUT_LIST & outputData);

    /** finalize state of some unfiltered variables */
    virtual void finish(){};

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
    virtual void compute_temperature_delta(const DENS_MAT & atomicEnergyDelta,
                                           double dt);

    // data
    /** nodal temperature 3rd time derivative */
    DENS_MAN & temperature3Roc_;

  private:

    // DO NOT define this
    ThermalTimeIntegratorFractionalStepFiltered();

  };
};
#endif
