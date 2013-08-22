#ifndef KINETOTHERMOSTAT_H
#define KINETOTHERMOSTAT_H

#include "AtomicRegulator.h"
#include "PerAtomQuantityLibrary.h"
//TEMP_JAT - transitional headers until we have a new method
#include "Kinetostat.h"
#include "Thermostat.h"
#include <map>
#include <set>
#include <string>

namespace ATC {
/* #ifdef WIP_JAT */
/*   static const int myLambdaMaxIterations = 50; */
/* #endif */
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

    //TEMP_JAT - required temporarily while we have two regulators
    
    /** initialization of method data */
    virtual void initialize();
    /** method(s) create all necessary transfer operators */
    virtual void construct_transfers();
    /** reset number of local atoms, as well as atomic data */
    virtual void reset_nlocal();
    /** set up atom to material identification */
    virtual void reset_atom_materials(const Array<int> & elementToMaterialMap,
                                      const MatrixDependencyManager<DenseMatrix, int> * atomElement);
    /** apply the regulator in the pre-predictor phase */
    virtual void apply_pre_predictor(double dt, int timeStep);
    /** apply the regulator in the mid-predictor phase */
    virtual void apply_mid_predictor(double dt, int timeStep);
    /** apply the regulator in the post-predictor phase */
    virtual void apply_post_predictor(double dt, int timeStep);
    /** apply the regulator in the pre-correction phase */
    virtual void apply_pre_corrector(double dt, int timeStep);
    /** apply the regulator in the post-correction phase */
    virtual void apply_post_corrector(double dt, int timeStep);
    /** prior to exchanges */
    virtual void pre_exchange();
    /** prior to force calculation */
    virtual void pre_force();
    /** force a reset to occur */
    void force_reset() {kinetostat_.force_reset();thermostat_.force_reset();};
    /** compute the thermal boundary flux, must be consistent with regulator */
    virtual void compute_boundary_flux(FIELDS & fields);
    /** type of boundary coupling */
    virtual RegulatorCouplingType coupling_mode(const FieldName) const;
    virtual void output(OUTPUT_LIST & outputData) const;
    /** final work at the end of a run */
    virtual void finish();
    /** add contributions (if any) to the finite element right-hand side */
    virtual void add_to_rhs(FIELDS & rhs);

    /** pack fields for restart */
    virtual void pack_fields(RESTART_LIST & data);
        
    // data access, intended for method objects
    /** reset the nodal power to a prescribed value */
    virtual void reset_lambda_contribution(const DENS_MAT & target,
                                           const FieldName field);
    //TEMP_JAT - will use real accessor
    /** return value for the correction maximum number of iterations */
    int lambda_max_iterators() {return thermostat_.lambda_max_iterations();};

  protected:

    //TEMP_JAT - individual controllers for now
    Kinetostat kinetostat_;
    Thermostat thermostat_;
  private:
    
    // DO NOT define this
    KinetoThermostat();
        
  };
/* #ifdef WIP_JAT */
/*   /\** */
/*    *  @class  ThermostatShapeFunction */
/*    *  @brief  Class for thermostat algorithms using the shape function matrices */
/*    *          (thermostats have general for of N^T w N lambda = rhs) */
/*    *\/ */
  
/*   class ThermostatShapeFunction : public RegulatorShapeFunction { */

  
/*   public: */
  
/*     ThermostatShapeFunction(Thermostat * thermostat, */
/*                             const std::string & regulatorPrefix = ""); */
        
/*     virtual ~ThermostatShapeFunction() {}; */

/*     /\** instantiate all needed data *\/ */
/*     virtual void construct_transfers(); */

/*   protected: */

/*     // methods */
/*     /\** set weighting factor for in matrix Nhat^T * weights * Nhat *\/ */
/*     virtual void set_weights(); */

/*     // member data */
/*     /\** pointer to thermostat object for data *\/ */
/*     Thermostat * thermostat_; */

/*     
/*     /\** MD mass matrix *\/ */
/*     DIAG_MAN & mdMassMatrix_;  

/*     /\** pointer to atom velocities *\/ */
/*     FundamentalAtomQuantity * atomVelocities_; */

/*     /\** workspace variables *\/ */
/*     DENS_VEC _weightVector_, _maskedWeightVector_; */

/*   private: */

/*     // DO NOT define this */
/*     ThermostatShapeFunction(); */

/*   }; */

/*   /\** */
/*    *  @class  ThermostatRescale */
/*    *  @brief  Enforces constraint on atomic kinetic energy based on FE temperature */
/*    *\/ */
  
/*   class ThermostatRescale : public ThermostatShapeFunction { */
  
/*   public: */
  
/*     ThermostatRescale(Thermostat * thermostat); */
        
/*     virtual ~ThermostatRescale() {}; */

/*     /\** instantiate all needed data *\/ */
/*     virtual void construct_transfers(); */
        
/*     /\** applies thermostat to atoms in the post-corrector phase *\/ */
/*     virtual void apply_post_corrector(double dt); */
/* #ifdef WIP_JAT */
/*     /\** compute boundary flux, requires thermostat input since it is part of the coupling scheme *\/ */
/*     virtual void compute_boundary_flux(FIELDS & fields) */
/*       {boundaryFlux_[TEMPERATURE] = 0.;}; */
/* #endif */

/*     /\** get data for output *\/ */
/*     virtual void output(OUTPUT_LIST & outputData); */
        
/*   protected: */

/*     /\** apply solution to atomic quantities *\/ */
/*     void apply_to_atoms(PerAtomQuantity<double> * atomVelocities); */

/*     /\** correct the RHS for complex temperature definitions *\/ */
/*     virtual void correct_rhs(DENS_MAT & rhs) {};  // base class does no correction, assuming kinetic definition */
        
/*     /\** FE temperature field *\/ */
/*     DENS_MAN & nodalTemperature_; */

/*     /\** construction for prolongation of lambda to atoms *\/ */
/*     AtomicVelocityRescaleFactor * atomVelocityRescalings_; */

/*     /\** workspace variables *\/ */
/*     DENS_MAT _rhs_; */

/*   private: */

/*     // DO NOT define this */
/*     ThermostatRescale(); */
  
/*   }; */

/*   /\** */
/*    *  @class  ThermostatRescaleMixedKePe */
/*    *  @brief  Enforces constraint on atomic kinetic energy based on FE temperature */
/*    *          when the temperature is a mix of the KE and PE */
/*    *\/ */
  
/*   class ThermostatRescaleMixedKePe : public ThermostatRescale { */
  
/*   public: */
  
/*     ThermostatRescaleMixedKePe(Thermostat * thermostat); */
        
/*     virtual ~ThermostatRescaleMixedKePe() {}; */

/*     /\** instantiate all needed data *\/ */
/*     virtual void construct_transfers(); */

/*     /\** pre-run initialization of method data *\/ */
/*     virtual void initialize(); */
        
/*   protected: */

/*     /\** correct the RHS for inclusion of the PE *\/ */
/*     virtual void correct_rhs(DENS_MAT & rhs); */

/*     /\** nodal fluctuating potential energy *\/ */
/*     DENS_MAN * nodalAtomicFluctuatingPotentialEnergy_; */

/*     /\** fraction of temperature from KE *\/ */
/*     double keMultiplier_; */

/*     /\** fraction of temperature from PE *\/ */
/*     double peMultiplier_; */

/*   private: */

/*     // DO NOT define this */
/*     ThermostatRescaleMixedKePe(); */
  
/*   }; */

/*   /\** */
/*    *  @class  ThermostatGlcFs */
/*    *  @brief  Class for thermostat algorithms based on Gaussian least constraints (GLC) for fractional step (FS) algorithsm */
/*    *\/ */
  
/*   class ThermostatGlcFs : public ThermostatShapeFunction { */
  
/*   public: */
  
/*     ThermostatGlcFs(Thermostat * thermostat, */
/*                     const std::string & regulatorPrefix = ""); */
        
/*     virtual ~ThermostatGlcFs() {}; */

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

/*     /\* flag for performing the full lambda prediction calculation *\/ */
/*     bool full_prediction(); */

/*   protected: */

/*     // methods */
/*     /\** determine mapping from all nodes to those to which the thermostat applies *\/ */
/*     void compute_rhs_map(); */

/*     /\** sets up appropriate rhs for thermostat equations *\/ */
/*     virtual void set_thermostat_rhs(DENS_MAT & rhs, */
/*                                     double dt) = 0; */

/*     /\** apply forces to atoms *\/ */
/*     virtual void apply_to_atoms(PerAtomQuantity<double> * atomicVelocity, */
/*                                 const DENS_MAN * nodalAtomicEnergy, */
/*                                 const DENS_MAT & lambdaForce, */
/*                                 DENS_MAT & nodalAtomicLambdaPower, */
/*                                 double dt); */

/*     /\** add contributions from thermostat to FE energy *\/ */
/*     virtual void add_to_energy(const DENS_MAT & nodalLambdaPower, */
/*                                DENS_MAT & deltaEnergy, */
/*                                double dt) = 0; */

/*     /\* sets up and solves the linear system for lambda *\/ */
/*     virtual void compute_lambda(double dt, */
/*                                 bool iterateSolution = true); */

/*     /\** solves the non-linear equation for lambda iteratively *\/ */
/*     void iterate_lambda(const MATRIX & rhs); */

/*     // member data */
/*     /\** reference to AtC FE temperature *\/ */
/*     DENS_MAN & temperature_; */

/*     /\** pointer to a time filtering object *\/ */
/*     TimeFilter * timeFilter_; */

/*     /\** power induced by lambda *\/ */
/*     DENS_MAN * nodalAtomicLambdaPower_; */

/*     /\** filtered lambda power *\/ */
/*     DENS_MAN * lambdaPowerFiltered_; */

/*     /\** atomic force induced by lambda *\/ */
/*     AtomicThermostatForce * atomThermostatForces_; */

/*     /\** pointer to atom masses *\/ */
/*     FundamentalAtomQuantity * atomMasses_; */

/*     /\** pointer to the values of lambda interpolated to atoms *\/ */
/*     DENS_MAN * rhsLambdaSquared_; */

/*     /\** hack to determine if first timestep has been passed *\/ */
/*     bool isFirstTimestep_; */

/*     /\** maximum number of iterations used in iterative solve for lambda *\/ */
/*     int lambdaMaxIterations_; */

/*     /\** nodal atomic energy *\/ */
/*     DENS_MAN * nodalAtomicEnergy_; */

/*     /\** local version of velocity used as predicted final veloctiy *\/ */
/*     PerAtomQuantity<double> * atomPredictedVelocities_; */

/*     /\** predicted nodal atomic energy *\/ */
/*     AtfShapeFunctionRestriction * nodalAtomicPredictedEnergy_; */

/*     /\** pointer for force applied in first time step *\/ */
/*     DENS_MAN * firstHalfAtomForces_; */

/*     /\** FE temperature change from thermostat during predictor phase in second half of timestep *\/ */
/*     DENS_MAT deltaEnergy1_; */

/*     /\** FE temperature change from thermostat during corrector phase in second half of timestep *\/ */
/*     DENS_MAT deltaEnergy2_; */

/*     /\** right-hand side data for thermostat equation *\/ */
/*     DENS_MAT rhs_; */

/*     /\** mapping from all to regulated nodes *\/ */
/*     DENS_MAT rhsMap_; */

/*     /\** fraction of timestep over which constraint is exactly enforced *\/ */
/*     double dtFactor_; */

/*     // workspace */
/*     DENS_MAT _lambdaPowerOutput_; // power applied by lambda in output format */
/*     DENS_MAT _velocityDelta_; // change in velocity when lambda force is applied */
/*     DENS_VEC _lambdaOverlap_; // lambda in MD overlapping FE nodes */
/*     DENS_MAT _lambdaOld_; // lambda from previous iteration */
/*     DENS_MAT _rhsOverlap_; // normal RHS vector mapped to overlap nodes */
/*     DENS_VEC _rhsTotal_; // normal + 2nd order RHS for the iteration loop */

/*   private: */

/*     // DO NOT define this */
/*     ThermostatGlcFs(); */

/*   }; */

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

/*     /\** intial restricted atomic FE energy used to compute change *\/ */
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
