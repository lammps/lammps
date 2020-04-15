#ifndef KINETOTHERMOSTAT_H
#define KINETOTHERMOSTAT_H

#include "AtomicRegulator.h"
#include "PerAtomQuantityLibrary.h"
#include "Kinetostat.h"
#include "Thermostat.h"
#include <map>
#include <set>
#include <string>

namespace ATC {

  static const int myCouplingMaxIterations = 50;

  // forward declarations
  class MomentumTimeIntegrator;
  class ThermalTimeIntegrator;
  class AtfShapeFunctionRestriction;
  class FundamentalAtomQuantity;
  class PrescribedDataManager;

  /**
   *  @class  KinetoThermostat
   *  @brief  Manager class for atom-continuum simulataneous control of momentum and thermal energy
   */
  class KinetoThermostat : public AtomicRegulator {
  
  public:

    // constructor
    KinetoThermostat(ATC_Coupling * atc,
                     const std::string & regulatorPrefix = "");
        
    // destructor
    virtual ~KinetoThermostat(){};
        
    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** instantiate up the desired method(s) */
    virtual void construct_methods();

        
    // data access, intended for method objects
    /** reset the nodal power to a prescribed value */
    virtual void reset_lambda_contribution(const DENS_MAT & target,
                                           const FieldName field);

    /** return value for the max number of mechanical/thermal coupling iterations */
    int coupling_max_iterations() const {return couplingMaxIterations_;};

  protected:

    /** maximum number of iterations used to solved coupled thermo/mechanical problem */
    int couplingMaxIterations_;

  private:
    
    // DO NOT define this
    KinetoThermostat();
        
  };

  /**
   *  @class  KinetoThermostatShapeFunction
   *  @brief  Class for kinetostat/thermostat algorithms using the shape function matrices
   *          (thermostats have general for of N^T w N lambda = rhs)
   */
  
  class KinetoThermostatShapeFunction : public RegulatorMethod {

  
  public:
  
    KinetoThermostatShapeFunction(AtomicRegulator * kinetoThermostat,
                                  int couplingMaxIterations,
                                  const std::string & /* regulatorPrefix */) : RegulatorMethod(kinetoThermostat),
      couplingMaxIterations_(couplingMaxIterations) {};
    KinetoThermostatShapeFunction(AtomicRegulator * kinetoThermostat,
                                  int couplingMaxIterations)
      : RegulatorMethod(kinetoThermostat), couplingMaxIterations_(couplingMaxIterations) {};
        
    virtual ~KinetoThermostatShapeFunction() {};

    /** instantiate all needed data */
    virtual void construct_transfers() = 0;

    /** initialize all data */
    virtual void initialize() {tolerance_ = atomicRegulator_->tolerance();};

  protected:

    /** maximum number of iterations between energy and momentum regulators */
    int couplingMaxIterations_;

    /** tolerance */
    double tolerance_;

  private:

    // DO NOT define this
    KinetoThermostatShapeFunction();

  };

  /**
   *  @class  VelocityRescaleCombined
   *  @brief  Enforces constraints on atomic velocity based on FE temperature and velocity
   */
  
  class VelocityRescaleCombined : public VelocityGlc {
  
  public:

    friend class KinetoThermostatRescale; // since this is basically a set of member functions for friend
  
    VelocityRescaleCombined(AtomicRegulator * kinetostat);
        
    virtual ~VelocityRescaleCombined(){};

    /** pre-run initialization of method data */
    virtual void initialize();
    
    /** applies kinetostat to atoms */
    virtual void apply_mid_predictor(double /* dt */){};
    /** applies kinetostat to atoms */
    virtual void apply_post_corrector(double /* dt */){};
    
    /** local shape function matrices are incompatible with this mode */
    virtual bool use_local_shape_functions() const {return false;};

  protected:

    // data
    /** reference to AtC FE velocity */
    DENS_MAN & velocity_;
    
    /** RHS correct based on thermostat */
    DENS_MAN * thermostatCorrection_;
  
    // methods
    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs, double dt);

    // disable un-needed functionality
    /** does initial filtering operations before main computation */
    virtual void apply_pre_filtering(double /* dt */){};
    /** applies kinetostat correction to atoms */
    virtual void apply_kinetostat(double /* dt */) {};
    /** computes the nodal FE force applied by the kinetostat */
    virtual void compute_nodal_lambda_force(double /* dt */){};
    /** apply any required corrections for localized kinetostats */
    virtual void apply_localization_correction(const DENS_MAT & /* source */,
                                               DENS_MAT & /* nodalField */,
                                               double /* weight */){};
    virtual void apply_localization_correction(const DENS_MAT & /* source */,
                                               DENS_MAT & /* nodalField */){};

  private:

    // DO NOT define this
    VelocityRescaleCombined();
  
  };

  /**
   *  @class  ThermostatRescaleCombined
   *  @brief  Enforces constraint on atomic kinetic energy based on FE temperature and velocity
   */
  
  class ThermostatRescaleCombined : public ThermostatRescale {
  
  public:
  
    ThermostatRescaleCombined(AtomicRegulator * thermostat);
        
    virtual ~ThermostatRescaleCombined() {};

    /** pre-run initialization of method data */
    virtual void initialize();
        
    // deactivate un-needed methods
    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double /* dt */){};
        
  protected:

    // data
    /** RHS correct based on kinetostat */
    DENS_MAN * kinetostatCorrection_;

    // deactivate un-needed methods
    /** apply solution to atomic quantities */
    virtual void apply_to_atoms(PerAtomQuantity<double> * /* atomVelocities */){};

    /** construct the RHS vector */
    virtual void set_rhs(DENS_MAT & rhs);

  private:

    // DO NOT define this
    ThermostatRescaleCombined();
  
  };

  /**
   *  @class  KinetoThermostatRescale
   *  @brief  Enforces constraints on atomic kinetic energy and velocity based on FE temperature and velocity
   */
  
  class KinetoThermostatRescale : public KinetoThermostatShapeFunction {
  
  public:
  
    KinetoThermostatRescale(AtomicRegulator * kinetoThermostat,
                            int couplingMaxIterations);
        
    virtual ~KinetoThermostatRescale();

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-run initialization of method data */
    virtual void initialize();
        
    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);

    /** compute boundary flux, requires thermostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & /* fields */)
    {boundaryFlux_[TEMPERATURE] = 0.; boundaryFlux_[VELOCITY] = 0.;};

    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);
        
  protected:

    // methods
    /** apply solution to atomic quantities */
    virtual void apply_to_atoms(PerAtomQuantity<double> * atomVelocities);

    /** creates the appropriate rescaling thermostat */
    virtual ThermostatRescale * construct_rescale_thermostat();

    // data
    /** pointer to atom velocities */
    FundamentalAtomQuantity * atomVelocities_;

    /** clone of FE velocity field */
    DENS_MAN & nodalVelocities_;

    /** lambda coupling parameter for momentum */
    DENS_MAN * lambdaMomentum_;

    /** lambda coupling parameter for energy */
    DENS_MAN * lambdaEnergy_;

    /** pointer to rescaled velocity fluctuations */
    PerAtomQuantity<double> * atomicFluctuatingVelocityRescaled_;

    /** pointer to streaming velocity */
    PerAtomQuantity<double> * atomicStreamingVelocity_;

    /** rescaling thermostat */
    ThermostatRescale * thermostat_;

    /** velocity regulating kinetostat */
    VelocityRescaleCombined * kinetostat_;

    // workspace
    DENS_MAT _lambdaEnergyOld_, _lambdaMomentumOld_, _diff_;

  private:

    // DO NOT define this
    KinetoThermostatRescale();
  
  };

  /**
   *  @class  ThermostatRescaleMixedKePeCombined
   *  @brief  Enforces constraint on atomic kinetic energy based on FE temperature and velocity when the temperature is comprised of both KE and PE contributions
   */
  
  class ThermostatRescaleMixedKePeCombined : public ThermostatRescaleMixedKePe {
  
  public:
  
    ThermostatRescaleMixedKePeCombined(AtomicRegulator * thermostat);
        
    virtual ~ThermostatRescaleMixedKePeCombined() {};

    /** pre-run initialization of method data */
    virtual void initialize();
        
    // deactivate un-needed methods
    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double /* dt */){};
        
  protected:

    // data
    /** RHS correct based on kinetostat */
    DENS_MAN * kinetostatCorrection_;

    // deactivate un-needed methods
    /** apply solution to atomic quantities */
    virtual void apply_to_atoms(PerAtomQuantity<double> * /* atomVelocities */){};

    /** construct the RHS vector */
    virtual void set_rhs(DENS_MAT & rhs);

  private:

    // DO NOT define this
    ThermostatRescaleMixedKePeCombined();
  
  };

  /**
   *  @class  KinetoThermostatRescaleMixedKePe
   *  @brief  Enforces constraint on atomic kinetic energy based on FE temperature
   *          when the temperature is a mix of the KE and PE
   */
  
  class KinetoThermostatRescaleMixedKePe : public KinetoThermostatRescale {
  
  public:
  
    KinetoThermostatRescaleMixedKePe(AtomicRegulator * kinetoThermostat,
                                     int couplingMaxIterations);

    virtual ~KinetoThermostatRescaleMixedKePe() {};
        
  protected:

    /** creates the appropriate rescaling thermostat */
    virtual ThermostatRescale * construct_rescale_thermostat();

  private:

    // DO NOT define this
    KinetoThermostatRescaleMixedKePe();
  
  };

  /**
   *  @class  KinetoThermostatGlcFs
   *  @brief  Class for regulation algorithms based on Gaussian least constraints (GLC) for fractional step (FS) algorithsm
   */
  
  class KinetoThermostatGlcFs : public KinetoThermostatShapeFunction {
  
  public:
  
    KinetoThermostatGlcFs(AtomicRegulator * kinetoThermostat,
                          int couplingMaxIterations,
                          const std::string & regulatorPrefix = "");
        
    virtual ~KinetoThermostatGlcFs() {};

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

    /** apply forces to atoms */
    virtual void apply_to_atoms(PerAtomQuantity<double> * atomicVelocity,
                                const DENS_MAN * nodalAtomicMomentum,
                                const DENS_MAN * nodalAtomicEnergy,
                                const DENS_MAT & lambdaForce,
                                DENS_MAT & nodalAtomicLambdaForce,
                                DENS_MAT & nodalAtomicLambdaPower,
                                double dt);

    // USE BASE CLASSES FOR THESE
    /** add contributions from regulator to FE momentum */
    virtual void add_to_momentum(const DENS_MAT & nodalLambdaForce,
                               DENS_MAT & deltaForce,
                               double dt) = 0;

    /** add contributions from regulator to FE energy */
    virtual void add_to_energy(const DENS_MAT & nodalLambdaPower,
                               DENS_MAT & deltaEnergy,
                               double dt) = 0;

    /* sets up and solves the linear system for lambda */
    virtual void compute_lambda(double dt,
                                bool iterateSolution = true);

    // member data
    /** reference to AtC FE velocity */
    DENS_MAN & velocity_;

    /** reference to AtC FE temperature */
    DENS_MAN & temperature_;

    /** pointer to a time filtering object */
    TimeFilter * timeFilter_;

    /** force induced by lambda */
    DENS_MAN * nodalAtomicLambdaForce_;

    /** filtered lambda force */
    DENS_MAN * lambdaForceFiltered_;

    /** power induced by lambda */
    DENS_MAN * nodalAtomicLambdaPower_;

    /** filtered lambda power */
    DENS_MAN * lambdaPowerFiltered_;

    /** atomic force induced by lambda */
    AtomicThermostatForce * atomRegulatorForces_;

    /** atomic force induced by thermostat lambda */
    AtomicThermostatForce * atomThermostatForces_;

    /** pointer to atom masses */
    FundamentalAtomQuantity * atomMasses_;

    /** pointer to atom velocities */
    FundamentalAtomQuantity * atomVelocities_;

    /** hack to determine if first timestep has been passed */
    bool isFirstTimestep_;

    /** nodal atomic momentum */
    DENS_MAN * nodalAtomicMomentum_;

    /** nodal atomic energy */
    DENS_MAN * nodalAtomicEnergy_;

    /** local version of velocity used as predicted final veloctiy */
    PerAtomQuantity<double> * atomPredictedVelocities_;

    /** predicted nodal atomic momentum */
    AtfShapeFunctionRestriction * nodalAtomicPredictedMomentum_;

    /** predicted nodal atomic energy */
    AtfShapeFunctionRestriction * nodalAtomicPredictedEnergy_;

    /** pointer for force applied in first time step */
    DENS_MAN * firstHalfAtomForces_;

    /** FE momentum change from regulator during predictor phase in second half of timestep */
    DENS_MAT deltaMomentum_;

    /** FE temperature change from regulator during predictor phase in second half of timestep */
    DENS_MAT deltaEnergy1_;

    /** FE temperature change from regulator during corrector phase in second half of timestep */
    DENS_MAT deltaEnergy2_;

    /** fraction of timestep over which constraint is exactly enforced */
    double dtFactor_;

    // workspace
    DENS_MAT _lambdaForceOutput_; // force applied by lambda in output format
    DENS_MAT _lambdaPowerOutput_; // power applied by lambda in output format
    DENS_MAT _velocityDelta_; // change in velocity when lambda force is applied

  private:

    // DO NOT define this
    KinetoThermostatGlcFs();

  };
/* #ifdef WIP_JAT */
/*   /\** */
/*    *  @class  ThermostatFlux */
/*    *  @brief  Class enforces GLC on atomic forces based on FE power when using fractional step time integration */
/*    *\/ */

/*   class ThermostatFlux : public ThermostatGlcFs { */
  
/*   public: */
  
/*     ThermostatFlux(Thermostat * thermostat, */
/*                    const std::string & regulatorPrefix = ""); */
        
/*     virtual ~ThermostatFlux() {}; */

/*     /\** instantiate all needed data *\/ */
/*     virtual void construct_transfers(); */

/*     /\** pre-run initialization of method data *\/ */
/*     virtual void initialize(); */
       
/*   protected: */

/*     /\** sets up appropriate rhs for thermostat equations *\/ */
/*     virtual void set_thermostat_rhs(DENS_MAT & rhs, */
/*                                     double dt); */

/*     /\** add contributions from thermostat to FE energy *\/ */
/*     virtual void add_to_energy(const DENS_MAT & nodalLambdaPower, */
/*                                DENS_MAT & deltaEnergy, */
/*                                double dt); */

/*     /\** sets up the transfer which is the set of nodes being regulated *\/ */
/*     virtual void construct_regulated_nodes(); */

/*     // data */
/*     /\** reference to ATC sources coming from prescribed data, AtC coupling, and extrinsic coupling *\/ */
/*     DENS_MAN & heatSource_; */

/*   private: */

/*     // DO NOT define this */
/*     ThermostatFlux(); */

/*   }; */

/*   /\** */
/*    *  @class  ThermostatFixed */
/*    *  @brief  Class enforces GLC on atomic forces based on FE power when using fractional step time integration */
/*    *\/ */

/*   class ThermostatFixed : public ThermostatGlcFs { */
  
/*   public: */
  
/*     ThermostatFixed(Thermostat * thermostat, */
/*                     const std::string & regulatorPrefix = ""); */
        
/*     virtual ~ThermostatFixed() {}; */

/*     /\** instantiate all needed data *\/ */
/*     virtual void construct_transfers(); */

/*     /\** pre-run initialization of method data *\/ */
/*     virtual void initialize(); */
        
/*     /\** applies thermostat to atoms in the predictor phase *\/ */
/*     virtual void apply_pre_predictor(double dt); */

/*     /\** applies thermostat to atoms in the pre-corrector phase *\/ */
/*     virtual void apply_pre_corrector(double dt); */

/*     /\** applies thermostat to atoms in the post-corrector phase *\/ */
/*     virtual void apply_post_corrector(double dt); */

/*     /\** compute boundary flux, requires thermostat input since it is part of the coupling scheme *\/ */
/*     virtual void compute_boundary_flux(FIELDS & fields) */
/*       {boundaryFlux_[TEMPERATURE] = 0.;}; */

/*     /\** determine if local shape function matrices are needed *\/ */
/*     virtual bool use_local_shape_functions() const {return atomicRegulator_->use_localized_lambda();}; */
        
/*   protected: */

/*     // methods */
/*     /\** initialize data for tracking the change in nodal atomic temperature *\/ */
/*     virtual void initialize_delta_nodal_atomic_energy(double dt); */

/*     /\** compute the change in nodal atomic temperature *\/ */
/*     virtual void compute_delta_nodal_atomic_energy(double dt); */

/*     /\** sets up appropriate rhs for thermostat equations *\/ */
/*     virtual void set_thermostat_rhs(DENS_MAT & rhs, */
/*                                     double dt); */

/*     /\** add contributions from thermostat to FE energy *\/ */
/*     virtual void add_to_energy(const DENS_MAT & nodalLambdaPower, */
/*                                DENS_MAT & deltaEnergy, */
/*                                double dt); */

/*     /\* sets up and solves the linear system for lambda *\/ */
/*     virtual void compute_lambda(double dt, */
/*                                 bool iterateSolution = true); */

/*     /\** flag for halving the applied force to mitigate numerical errors *\/ */
/*     bool halve_force(); */

/*     /\** sets up the transfer which is the set of nodes being regulated *\/ */
/*     virtual void construct_regulated_nodes(); */

/*     // data */
/*     /\** change in FE energy over a timestep *\/ */
/*     DENS_MAT deltaFeEnergy_; */

/*     /\** initial FE energy used to compute change *\/ */
/*     DENS_MAT initialFeEnergy_; */

/*     /\** change in restricted atomic FE energy over a timestep *\/ */
/*     DENS_MAT deltaNodalAtomicEnergy_; */

/*     /\** initial restricted atomic FE energy used to compute change *\/ */
/*     DENS_MAT initialNodalAtomicEnergy_; */

/*     /\** filtered nodal atomic energy *\/ */
/*     DENS_MAN nodalAtomicEnergyFiltered_; */

/*     /\** forces depending on predicted velocities for correct updating with fixed nodes *\/ */
/*     AtomicThermostatForce * atomThermostatForcesPredVel_; */

/*     /\** coefficient to account for effect of time filtering on rhs terms *\/ */
/*     double filterCoefficient_; */

/*     /\** kinetic energy multiplier in total energy (used for temperature expression) *\/ */
/*     double keMultiplier_; */

/*     // workspace */
/*     DENS_MAT _tempNodalAtomicEnergyFiltered_; // stores filtered energy change in atoms for persistence during predictor */

/*   private: */

/*     // DO NOT define this */
/*     ThermostatFixed(); */

/*   }; */

/*   /\** */
/*    *  @class  ThermostatFluxFiltered */
/*    *  @brief  Class enforces GLC on atomic forces based on FE power when using fractional step time integration */
/*    *          in conjunction with time filtering */
/*    *\/ */

/*   class ThermostatFluxFiltered : public ThermostatFlux { */
  
/*   public: */
  
/*     ThermostatFluxFiltered(Thermostat * thermostat, */
/*                            const std::string & regulatorPrefix = ""); */
        
/*     virtual ~ThermostatFluxFiltered() {}; */

/*     /\** pre-run initialization of method data *\/ */
/*     virtual void initialize(); */

/*     /\** applies thermostat to atoms in the post-corrector phase *\/ */
/*     virtual void apply_post_corrector(double dt); */

/*     /\** get data for output *\/ */
/*     virtual void output(OUTPUT_LIST & outputData); */
       
/*   protected: */

/*     /\** sets up appropriate rhs for thermostat equations *\/ */
/*     virtual void set_thermostat_rhs(DENS_MAT & rhs, */
/*                                     double dt); */

/*     /\** add contributions from thermostat to FE energy *\/ */
/*     virtual void add_to_energy(const DENS_MAT & nodalLambdaPower, */
/*                                DENS_MAT & deltaEnergy, */
/*                                double dt); */

/*     // data */
/*     /\** heat source time history required to get instantaneous heat sources *\/ */
/*     DENS_MAT heatSourceOld_; */
/*     DENS_MAT instantHeatSource_; */
/*     DENS_MAT timeStepSource_; */

/*   private: */

/*     // DO NOT define this */
/*     ThermostatFluxFiltered(); */

/*   }; */

/*   /\** */
/*    *  @class  ThermostatFixedFiltered */
/*    *  @brief  Class for thermostatting using the temperature matching constraint and is compatible with */
/*  the fractional step time-integration with time filtering */
/*    *\/ */
  
/*   class ThermostatFixedFiltered : public ThermostatFixed { */
  
/*   public: */
  
/*     ThermostatFixedFiltered(Thermostat * thermostat, */
/*                             const std::string & regulatorPrefix = ""); */
        
/*     virtual ~ThermostatFixedFiltered() {}; */

/*     /\** get data for output *\/ */
/*     virtual void output(OUTPUT_LIST & outputData); */
        
/*   protected: */

/*     // methods */
/*     /\** initialize data for tracking the change in nodal atomic temperature *\/ */
/*     virtual void initialize_delta_nodal_atomic_energy(double dt); */

/*     /\** compute the change in nodal atomic temperature *\/ */
/*     virtual void compute_delta_nodal_atomic_energy(double dt); */

/*     /\** sets up appropriate rhs for thermostat equations *\/ */
/*     virtual void set_thermostat_rhs(DENS_MAT & rhs, */
/*                                     double dt); */

/*     /\** add contributions from thermostat to temperature for uncoupled nodes *\/ */
/*     virtual void add_to_energy(const DENS_MAT & nodalLambdaPower, */
/*                                DENS_MAT & deltaEnergy, */
/*                                double dt); */

/*   private: */

/*     // DO NOT define this */
/*     ThermostatFixedFiltered(); */

/*   }; */

/*   /\** */
/*    *  @class  ThermostatFluxFixed */
/*    *  @brief  Class for thermostatting using the temperature matching constraint one one set of nodes and the flux matching constraint on another */
/*    *\/ */

/*   class ThermostatFluxFixed : public RegulatorMethod { */

/*   public: */

/*     ThermostatFluxFixed(Thermostat * thermostat, */
/*                         bool constructThermostats = true); */
        
/*     virtual ~ThermostatFluxFixed(); */

/*     /\** instantiate all needed data *\/ */
/*     virtual void construct_transfers(); */

/*     /\** pre-run initialization of method data *\/ */
/*     virtual void initialize(); */

/*     /\** applies thermostat to atoms in the predictor phase *\/ */
/*     virtual void apply_pre_predictor(double dt); */

/*     /\** applies thermostat to atoms in the pre-corrector phase *\/ */
/*     virtual void apply_pre_corrector(double dt); */

/*     /\** applies thermostat to atoms in the post-corrector phase *\/ */
/*     virtual void apply_post_corrector(double dt); */
    
/*     /\** get data for output *\/ */
/*     virtual void output(OUTPUT_LIST & outputData); */

/*     /\** compute boundary flux, requires thermostat input since it is part of the coupling scheme *\/ */
/*     virtual void compute_boundary_flux(FIELDS & fields) */
/*       {thermostatBcs_->compute_boundary_flux(fields);}; */

/*   protected: */

/*     // data */
/*     /\** thermostat for imposing the fluxes *\/ */
/*     ThermostatFlux * thermostatFlux_; */

/*     /\** thermostat for imposing fixed nodes *\/ */
/*     ThermostatFixed * thermostatFixed_; */

/*     /\** pointer to whichever thermostat should compute the flux, based on coupling method *\/ */
/*     ThermostatGlcFs * thermostatBcs_; */

/*   private: */

/*     // DO NOT define this */
/*     ThermostatFluxFixed(); */
/*   }; */

/*   /\** */
/*    *  @class  ThermostatFluxFixedFiltered */
/*    *  @brief  Class for thermostatting using the temperature matching constraint one one set of nodes and the flux matching constraint on another with time filtering */
/*    *\/ */

/*   class ThermostatFluxFixedFiltered : public ThermostatFluxFixed { */

/*   public: */

/*     ThermostatFluxFixedFiltered(Thermostat * thermostat); */
        
/*     virtual ~ThermostatFluxFixedFiltered(){}; */

/*   private: */

/*     // DO NOT define this */
/*     ThermostatFluxFixedFiltered(); */

/*   }; */
/* #endif */

};

#endif
