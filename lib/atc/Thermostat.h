#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include "AtomicRegulator.h"
#include "PerAtomQuantityLibrary.h"
#include <map>
#include <set>
#include <string>

namespace ATC {

  static const int myLambdaMaxIterations = 50;

  // forward declarations
  class ThermalTimeIntegrator;
  class AtfShapeFunctionRestriction;
  class FundamentalAtomQuantity;
  class PrescribedDataManager;

  /**
   *  @class  Thermostat
   *  @brief  Manager class for atom-continuum control of thermal energy
   */
  class Thermostat : public AtomicRegulator {
  
  public:

    // constructor
    Thermostat(ATC_Coupling * atc,
               const std::string & regulatorPrefix = "");
        
    // destructor
    virtual ~Thermostat(){};
        
    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** instantiate up the desired method(s) */
    virtual void construct_methods();
        
    // data access, intended for method objects
    /** reset the nodal power to a prescribed value */
    virtual void reset_lambda_contribution(const DENS_MAT & target);

    /** return value for the correction maximum number of iterations */
    int lambda_max_iterations() {return lambdaMaxIterations_;};

  protected:

    // data regarding fixed nodes and applied fluxes
    /** set of all fixed nodes */
    std::set<int> fixedNodes_;
    /** set of all nodes which have a flux applied */
    std::set<int> fluxNodes_;

    /** maximum number of iterations used in iterative solve for lambda */
    int lambdaMaxIterations_;

  private:
    
    // DO NOT define this
    Thermostat();
        
  };

  /**
   *  @class  ThermostatShapeFunction
   *  @brief  Class for thermostat algorithms using the shape function matrices
   *          (thermostats have general for of N^T w N lambda = rhs)
   */
  
  class ThermostatShapeFunction : public RegulatorShapeFunction {

  
  public:
  
    ThermostatShapeFunction(Thermostat * thermostat,
                            const std::string & regulatorPrefix = "");
        
    virtual ~ThermostatShapeFunction() {};

    /** instantiate all needed data */
    virtual void construct_transfers();

  protected:

    // methods
    /** set weighting factor for in matrix Nhat^T * weights * Nhat */
    virtual void set_weights();

    // member data
    /** pointer to thermostat object for data */
    Thermostat * thermostat_;

    
    /** MD mass matrix */
    DIAG_MAN & mdMassMatrix_;  

    /** pointer to atom velocities */
    FundamentalAtomQuantity * atomVelocities_;

    /** workspace variables */
    DENS_VEC _weightVector_, _maskedWeightVector_;

  private:

    // DO NOT define this
    ThermostatShapeFunction();

  };

  /**
   *  @class  ThermostatRescale
   *  @brief  Enforces constraint on atomic kinetic energy based on FE temperature
   */
  
  class ThermostatRescale : public ThermostatShapeFunction {
  
  public:
  
    ThermostatRescale(Thermostat * thermostat);
        
    virtual ~ThermostatRescale() {};

    /** instantiate all needed data */
    virtual void construct_transfers();
        
    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);

    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);
        
  protected:

    /** apply solution to atomic quantities */
    void apply_to_atoms(PerAtomQuantity<double> * atomVelocities);

    /** correct the RHS for complex temperature definitions */
    virtual void correct_rhs(DENS_MAT & rhs) {};  // base class does no correction, assuming kinetic definition
        
    /** FE temperature field */
    DENS_MAN & nodalTemperature_;

    /** construction for prolongation of lambda to atoms */
    AtomicVelocityRescaleFactor * atomVelocityRescalings_;

    /** workspace variables */
    DENS_MAT _rhs_;

  private:

    // DO NOT define this
    ThermostatRescale();
  
  };

  /**
   *  @class  ThermostatRescaleMixedKePe
   *  @brief  Enforces constraint on atomic kinetic energy based on FE temperature
   *          when the temperature is a mix of the KE and PE
   */
  
  class ThermostatRescaleMixedKePe : public ThermostatRescale {
  
  public:
  
    ThermostatRescaleMixedKePe(Thermostat * thermostat);
        
    virtual ~ThermostatRescaleMixedKePe() {};

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-run initialization of method data */
    virtual void initialize();
        
  protected:

    /** correct the RHS for inclusion of the PE */
    virtual void correct_rhs(DENS_MAT & rhs);

    /** nodal fluctuating potential energy */
    DENS_MAN * nodalAtomicFluctuatingPotentialEnergy_;

    /** fraction of temperature from KE */
    double keMultiplier_;

    /** fraction of temperature from PE */
    double peMultiplier_;

  private:

    // DO NOT define this
    ThermostatRescaleMixedKePe();
  
  };

  /**
   *  @class  ThermostatGlcFs
   *  @brief  Class for thermostat algorithms based on Gaussian least constraints (GLC) for fractional step (FS) algorithsm
   */
  
  class ThermostatGlcFs : public ThermostatShapeFunction {
  
  public:
  
    ThermostatGlcFs(Thermostat * thermostat,
                    const std::string & regulatorPrefix = "");
        
    virtual ~ThermostatGlcFs() {};

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-run initialization of method data */
    virtual void initialize();

    /** applies thermostat to atoms in the predictor phase */
    virtual void apply_pre_predictor(double dt);

    /** applies thermostat to atoms in the pre-corrector phase */
    virtual void apply_pre_corrector(double dt);

    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);
    
    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);

    /* flag for performing the full lambda prediction calculation */
    bool full_prediction();

  protected:

    // methods
    /** determine mapping from all nodes to those to which the thermostat applies */
    void compute_rhs_map();

    /** sets up appropriate rhs for thermostat equations */
    virtual void set_thermostat_rhs(DENS_MAT & rhs,
                                    double dt) = 0;

    /** apply forces to atoms */
    virtual void apply_to_atoms(PerAtomQuantity<double> * atomicVelocity,
                                const DENS_MAN * nodalAtomicEnergy,
                                const DENS_MAT & lambdaForce,
                                DENS_MAT & nodalAtomicLambdaPower,
                                double dt);

    /** add contributions from thermostat to FE energy */
    virtual void add_to_energy(const DENS_MAT & nodalLambdaPower,
                               DENS_MAT & deltaEnergy,
                               double dt) = 0;

    /* sets up and solves the linear system for lambda */
    virtual void compute_lambda(double dt,
                                bool iterateSolution = true);

    /** solves the non-linear equation for lambda iteratively */
    void iterate_lambda(const MATRIX & rhs);

    // member data
    /** reference to AtC FE temperature */
    DENS_MAN & temperature_;

    /** pointer to a time filtering object */
    TimeFilter * timeFilter_;

    /** power induced by lambda */
    DENS_MAN * nodalAtomicLambdaPower_;

    /** filtered lambda power */
    DENS_MAN * lambdaPowerFiltered_;

    /** atomic force induced by lambda */
    AtomicThermostatForce * atomThermostatForces_;

    /** pointer to atom masses */
    FundamentalAtomQuantity * atomMasses_;

    /** pointer to the values of lambda interpolated to atoms */
    DENS_MAN * rhsLambdaSquared_;

    /** hack to determine if first timestep has been passed */
    bool isFirstTimestep_;

    /** maximum number of iterations used in iterative solve for lambda */
    int lambdaMaxIterations_;

    /** nodal atomic energy */
    DENS_MAN * nodalAtomicEnergy_;

    /** local version of velocity used as predicted final veloctiy */
    PerAtomQuantity<double> * atomPredictedVelocities_;

    /** predicted nodal atomic energy */
    AtfShapeFunctionRestriction * nodalAtomicPredictedEnergy_;

    /** pointer for force applied in first time step */
    DENS_MAN * firstHalfAtomForces_;

    /** FE temperature change from thermostat during predictor phase in second half of timestep */
    DENS_MAT deltaEnergy1_;

    /** FE temperature change from thermostat during corrector phase in second half of timestep */
    DENS_MAT deltaEnergy2_;

    /** right-hand side data for thermostat equation */
    DENS_MAT rhs_;

    /** mapping from all to regulated nodes */
    DENS_MAT rhsMap_;

    /** fraction of timestep over which constraint is exactly enforced */
    double dtFactor_;

    // workspace
    DENS_MAT _lambdaPowerOutput_; // power applied by lambda in output format
    DENS_MAT _velocityDelta_; // change in velocity when lambda force is applied
    DENS_VEC _lambdaOverlap_; // lambda in MD overlapping FE nodes
    DENS_MAT _lambdaOld_; // lambda from previous iteration
    DENS_MAT _rhsOverlap_; // normal RHS vector mapped to overlap nodes
    DENS_VEC _rhsTotal_; // normal + 2nd order RHS for the iteration loop

  private:

    // DO NOT define this
    ThermostatGlcFs();

  };

  /**
   *  @class  ThermostatFlux
   *  @brief  Class enforces GLC on atomic forces based on FE power when using fractional step time integration
   */

  class ThermostatFlux : public ThermostatGlcFs {
  
  public:
  
    ThermostatFlux(Thermostat * thermostat,
                   const std::string & regulatorPrefix = "");
        
    virtual ~ThermostatFlux() {};

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-run initialization of method data */
    virtual void initialize();
       
  protected:

    /** sets up appropriate rhs for thermostat equations */
    virtual void set_thermostat_rhs(DENS_MAT & rhs,
                                    double dt);

    /** add contributions from thermostat to FE energy */
    virtual void add_to_energy(const DENS_MAT & nodalLambdaPower,
                               DENS_MAT & deltaEnergy,
                               double dt);

    /** sets up the transfer which is the set of nodes being regulated */
    virtual void construct_regulated_nodes();

    // data
    /** reference to ATC sources coming from prescribed data, AtC coupling, and extrinsic coupling */
    DENS_MAN & heatSource_;

  private:

    // DO NOT define this
    ThermostatFlux();

  };

  /**
   *  @class  ThermostatFixed
   *  @brief  Class enforces GLC on atomic forces based on FE power when using fractional step time integration
   */

  class ThermostatFixed : public ThermostatGlcFs {
  
  public:
  
    ThermostatFixed(Thermostat * thermostat,
                    const std::string & regulatorPrefix = "");
        
    virtual ~ThermostatFixed() {};

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-run initialization of method data */
    virtual void initialize();
        
    /** applies thermostat to atoms in the predictor phase */
    virtual void apply_pre_predictor(double dt);

    /** applies thermostat to atoms in the pre-corrector phase */
    virtual void apply_pre_corrector(double dt);

    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);

    /** compute boundary flux, requires thermostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields)
      {boundaryFlux_[TEMPERATURE] = 0.;};

    /** determine if local shape function matrices are needed */
    virtual bool use_local_shape_functions() const {return atomicRegulator_->use_localized_lambda();};
        
  protected:

    // methods
    /** initialize data for tracking the change in nodal atomic temperature */
    virtual void initialize_delta_nodal_atomic_energy(double dt);

    /** compute the change in nodal atomic temperature */
    virtual void compute_delta_nodal_atomic_energy(double dt);

    /** sets up appropriate rhs for thermostat equations */
    virtual void set_thermostat_rhs(DENS_MAT & rhs,
                                    double dt);

    /** add contributions from thermostat to FE energy */
    virtual void add_to_energy(const DENS_MAT & nodalLambdaPower,
                               DENS_MAT & deltaEnergy,
                               double dt);

    /* sets up and solves the linear system for lambda */
    virtual void compute_lambda(double dt,
                                bool iterateSolution = true);

    /** flag for halving the applied force to mitigate numerical errors */
    bool halve_force();

    /** sets up the transfer which is the set of nodes being regulated */
    virtual void construct_regulated_nodes();

    // data
    /** change in FE energy over a timestep */
    DENS_MAT deltaFeEnergy_;

    /** initial FE energy used to compute change */
    DENS_MAT initialFeEnergy_;

    /** change in restricted atomic FE energy over a timestep */
    DENS_MAT deltaNodalAtomicEnergy_;

    /** intial restricted atomic FE energy used to compute change */
    DENS_MAT initialNodalAtomicEnergy_;

    /** filtered nodal atomic energy */
    DENS_MAN nodalAtomicEnergyFiltered_;

    /** forces depending on predicted velocities for correct updating with fixed nodes */
    AtomicThermostatForce * atomThermostatForcesPredVel_;

    /** coefficient to account for effect of time filtering on rhs terms */
    double filterCoefficient_;

    /** kinetic energy multiplier in total energy (used for temperature expression) */
    double keMultiplier_;

    // workspace
    DENS_MAT _tempNodalAtomicEnergyFiltered_; // stores filtered energy change in atoms for persistence during predictor

  private:

    // DO NOT define this
    ThermostatFixed();

  };

  /**
   *  @class  ThermostatFluxFiltered
   *  @brief  Class enforces GLC on atomic forces based on FE power when using fractional step time integration
   *          in conjunction with time filtering
   */

  class ThermostatFluxFiltered : public ThermostatFlux {
  
  public:
  
    ThermostatFluxFiltered(Thermostat * thermostat,
                           const std::string & regulatorPrefix = "");
        
    virtual ~ThermostatFluxFiltered() {};

    /** pre-run initialization of method data */
    virtual void initialize();

    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);

    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);
       
  protected:

    /** sets up appropriate rhs for thermostat equations */
    virtual void set_thermostat_rhs(DENS_MAT & rhs,
                                    double dt);

    /** add contributions from thermostat to FE energy */
    virtual void add_to_energy(const DENS_MAT & nodalLambdaPower,
                               DENS_MAT & deltaEnergy,
                               double dt);

    // data
    /** heat source time history required to get instantaneous heat sources */
    DENS_MAT heatSourceOld_;
    DENS_MAT instantHeatSource_;
    DENS_MAT timeStepSource_;

  private:

    // DO NOT define this
    ThermostatFluxFiltered();

  };

  /**
   *  @class  ThermostatFixedFiltered
   *  @brief  Class for thermostatting using the temperature matching constraint and is compatible with
 the fractional step time-integration with time filtering
   */
  
  class ThermostatFixedFiltered : public ThermostatFixed {
  
  public:
  
    ThermostatFixedFiltered(Thermostat * thermostat,
                            const std::string & regulatorPrefix = "");
        
    virtual ~ThermostatFixedFiltered() {};

    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);
        
  protected:

    // methods
    /** initialize data for tracking the change in nodal atomic temperature */
    virtual void initialize_delta_nodal_atomic_energy(double dt);

    /** compute the change in nodal atomic temperature */
    virtual void compute_delta_nodal_atomic_energy(double dt);

    /** sets up appropriate rhs for thermostat equations */
    virtual void set_thermostat_rhs(DENS_MAT & rhs,
                                    double dt);

    /** add contributions from thermostat to temperature for uncoupled nodes */
    virtual void add_to_energy(const DENS_MAT & nodalLambdaPower,
                               DENS_MAT & deltaEnergy,
                               double dt);

  private:

    // DO NOT define this
    ThermostatFixedFiltered();

  };

  /**
   *  @class  ThermostatFluxFixed
   *  @brief  Class for thermostatting using the temperature matching constraint one one set of nodes and the flux matching constraint on another
   */

  class ThermostatFluxFixed : public RegulatorMethod {

  public:

    ThermostatFluxFixed(Thermostat * thermostat,
                        bool constructThermostats = true);
        
    virtual ~ThermostatFluxFixed();

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-run initialization of method data */
    virtual void initialize();

    /** applies thermostat to atoms in the predictor phase */
    virtual void apply_pre_predictor(double dt);

    /** applies thermostat to atoms in the pre-corrector phase */
    virtual void apply_pre_corrector(double dt);

    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);
    
    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);

    /** compute boundary flux, requires thermostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields)
      {thermostatBcs_->compute_boundary_flux(fields);};

  protected:

    // data
    /** thermostat for imposing the fluxes */
    ThermostatFlux * thermostatFlux_;

    /** thermostat for imposing fixed nodes */
    ThermostatFixed * thermostatFixed_;

    /** pointer to whichever thermostat should compute the flux, based on coupling method */
    ThermostatGlcFs * thermostatBcs_;

  private:

    // DO NOT define this
    ThermostatFluxFixed();
  };

  /**
   *  @class  ThermostatFluxFixedFiltered
   *  @brief  Class for thermostatting using the temperature matching constraint one one set of nodes and the flux matching constraint on another with time filtering
   */

  class ThermostatFluxFixedFiltered : public ThermostatFluxFixed {

  public:

    ThermostatFluxFixedFiltered(Thermostat * thermostat);
        
    virtual ~ThermostatFluxFixedFiltered(){};

  private:

    // DO NOT define this
    ThermostatFluxFixedFiltered();

  };

  /**
   *  @class  ThermostatGlc
   *  @brief  Class for thermostat algorithms based on Gaussian least constraints (GLC)
   */
  
  class ThermostatGlc : public ThermostatShapeFunction {
  
  public:
  
    ThermostatGlc(Thermostat * thermostat);
        
    virtual ~ThermostatGlc() {};

    /** instantiate all needed data */
    virtual void construct_transfers();

  protected:

    // methods
    /** apply forces to atoms */
    virtual void apply_to_atoms(PerAtomQuantity<double> * atomicVelocity, 
                                const DENS_MAT & lambdaForce,
                                double dt);

    // member data
    /** pointer to a time filtering object */
    TimeFilter * timeFilter_;

    /** filtered lambda power */
    DENS_MAN * lambdaPowerFiltered_;

    /** atomic force induced by lambda */
    PerAtomQuantity<double> * atomThermostatForces_;

    /** pointer to access prescribed data for fixed nodes */
    PrescribedDataManager * prescribedDataMgr_;

    /** pointer to atom masses */
    FundamentalAtomQuantity * atomMasses_;

    /** workspace variables */
    DENS_MAT _velocityDelta_;

  private:

    // DO NOT define this
    ThermostatGlc();

  };

  /**
   *  @class  ThermostatPowerVerlet
   *  @brief  Class for thermostatting using the heat flux matching constraint and is compatible with
 the Gear time-integration
   */
  
  class ThermostatPowerVerlet : public ThermostatGlc {
  
  public:

    ThermostatPowerVerlet(Thermostat * thermostat);
        
    virtual ~ThermostatPowerVerlet() {};

    /** instantiate all needed data */
    virtual void construct_transfers();
    
    /** pre-run initialization of method data */
    virtual void initialize();
        
    /** applies thermostat to atoms in the predictor phase */
    virtual void apply_pre_predictor(double dt);

    /** applies thermostat to atoms in the pre-corrector phase */
    virtual void apply_pre_corrector(double dt);

    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);

    /** final tasks of a run */
    virtual void finish();

    /** determine if local shape function matrices are needed */
    virtual bool use_local_shape_functions() const {return (!(atomicRegulator_->use_lumped_lambda_solve()) && atomicRegulator_->use_localized_lambda());};

  protected:

    /** nodal temperature rate of change */
    DENS_MAN & nodalTemperatureRoc_;

    /** reference to ATC sources coming from prescribed data, AtC coupling, and extrinsic coupling */
    DENS_MAN & heatSource_;

    /** pointer to nodal atomic power */
    DENS_MAN * nodalAtomicPower_;

    /** power applied to each atom by lambda force */
    AtfShapeFunctionRestriction * nodalAtomicLambdaPower_;

    /** workspace variables */
    DENS_MAT _rhs_;

    /** sets up and solves thermostat equations */
    virtual void compute_thermostat(double dt);

    /** sets up appropriate rhs for thermostat equations */
    virtual void set_thermostat_rhs(DENS_MAT & rhs_nodes);

    /** add contributions (if any) to the finite element right-hand side */
    virtual void add_to_rhs(FIELDS & rhs);

    // workspace
    DENS_MAT _nodalAtomicLambdaPowerOut_; // power induced by lambda in output format

  private:

    // DO NOT define this
    ThermostatPowerVerlet();

  };
  
  /**
   *  @class  ThermostatHooverVerlet
   *  @brief  Classfor thermostatting using the temperature matching constraint and is compatible with
 Gear time-integration 
   */
  
  class ThermostatHooverVerlet : public ThermostatPowerVerlet {
  
  public:

    ThermostatHooverVerlet(Thermostat * thermostat);
        
    virtual ~ThermostatHooverVerlet() {};

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** final tasks of a run */
    virtual void finish() {};

    /** compute boundary flux, requires thermostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields)
      {boundaryFlux_[TEMPERATURE] = 0.;};

  protected:

    /** lambda coupling parameter for hoover thermostat */
    DENS_MAN * lambdaHoover_;

    /** workspace variables */
    DENS_MAT _myNodalLambdaPower_;

    /** sets up and solves thermostat equations */
    virtual void compute_thermostat(double dt);

    /** sets up Hoover component of the thermostat */
    void set_hoover_rhs(DENS_MAT & rhs);

    /** add Hoover contributions to lambda power */
    void add_to_lambda_power(const DENS_MAT & myLambdaForce,
                             double dt);

    /** power applied to each atom by hoover lambda force */
    AtfShapeFunctionRestriction * nodalAtomicHooverLambdaPower_;

    /** add contributions (if any) to the finite element right-hand side */
    virtual void add_to_rhs(FIELDS & rhs);

  private:

    // DO NOT implement this
    ThermostatHooverVerlet();

  };

  /**
   *  @class  ThermostatPowerVerletFiltered
   *  @brief  Class for thermostatting using the heat flux matching constraint and is compatible with
 Gear time-integration with time filtering
   */
  
  class ThermostatPowerVerletFiltered : public ThermostatPowerVerlet {
  
  public:

    ThermostatPowerVerletFiltered(Thermostat * thermostat);
        
    virtual ~ThermostatPowerVerletFiltered(){};

    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);

    /** compute boundary flux, requires thermostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields);

  protected:

    /** sets up appropriate rhs for thermostat equations */
    virtual void set_thermostat_rhs(DENS_MAT & rhs_nodes);

    /** add contributions (if any) to the finite element right-hand side */
    virtual void add_to_rhs(FIELDS & rhs);

    /** nodal temperature 2nd rate of change (i.e. second time derivative) */
    DENS_MAN & nodalTemperature2Roc_;

    /** reference to ATC rate of change sources coming from prescribed data, AtC coupling, and extrinsic coupling */
    DENS_MAN heatSourceRoc_;

    /** references to ATC field rates of changing for inverting the filtered heat sources */
    FIELDS & fieldsRoc_;
   
    /** flux rate of changes for inverting filtered fluxes */
    FIELDS fluxRoc_;

    /** time scale for the time filter */
    double filterScale_;

  private:

    // DO NOT define this
    ThermostatPowerVerletFiltered();

  };

  /**
   *  @class  ThermostatHooverVerletFiltered
   *  @brief  Class for thermostatting using the temperature matching constraint and is compatible with
 Gear time-integration with time filtering
   */
  
  class ThermostatHooverVerletFiltered : public ThermostatPowerVerletFiltered {
  
  public:

    ThermostatHooverVerletFiltered(Thermostat * thermostat);
        
    virtual ~ThermostatHooverVerletFiltered() {};

    /** instantiate all needed data */
    virtual void construct_transfers();
    
    /** final tasks of a run */
    virtual void finish() {};
    
    /** compute boundary flux, requires thermostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields)
      {boundaryFlux_[TEMPERATURE] = 0.;};

  protected:

    /** lambda coupling parameter for hoover thermostat */
    DENS_MAN * lambdaHoover_;

    /** workspace variables */
    DENS_MAT _myNodalLambdaPower_;

    /** sets up and solves thermostat equations */
    virtual void compute_thermostat(double dt);

    /** sets up Hoover component of the thermostat */
    void set_hoover_rhs(DENS_MAT & rhs);

    /** add Hoover contributions to lambda power */
    void add_to_lambda_power(const DENS_MAT & myLambdaForce,
                             double dt);

    /** power applied to each atom by hoover lambda force */
    DENS_MAN * nodalAtomicHooverLambdaPower_;

    /** add contributions (if any) to the finite element right-hand side */
    virtual void add_to_rhs(FIELDS & rhs);

  private:

    // DO NOT implement this
    ThermostatHooverVerletFiltered();

  };

};

#endif
