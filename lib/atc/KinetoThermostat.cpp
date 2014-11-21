#include "KinetoThermostat.h"
#include "ATC_Coupling.h"
#include "ATC_Error.h"
#include "PrescribedDataManager.h"
#include "ElasticTimeIntegrator.h"
#include "ThermalTimeIntegrator.h"
#include "TransferOperator.h"

using namespace std;

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class KinetoThermostat
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  KinetoThermostat::KinetoThermostat(ATC_Coupling * atc,
                                     const string & regulatorPrefix) :
    AtomicRegulator(atc,regulatorPrefix),
    couplingMaxIterations_(myCouplingMaxIterations)
  {
    // nothing to do
  }

  //--------------------------------------------------------
  //  modify:
  //    parses and adjusts thermostat state based on
  //    user input, in the style of LAMMPS user input
  //--------------------------------------------------------
  bool KinetoThermostat::modify(int narg, char **arg)
  {
    bool foundMatch = false;
    return foundMatch;
  }

  //--------------------------------------------------------
  //  reset_lambda_contribution:
  //    resets the thermostat generated power to a
  //    prescribed value
  //--------------------------------------------------------
  void KinetoThermostat::reset_lambda_contribution(const DENS_MAT & target,
                                                   const FieldName field)
  {
    if (field==VELOCITY) {
      DENS_MAN * lambdaForceFiltered = regulator_data("LambdaForceFiltered",atc_->nsd());
      *lambdaForceFiltered = target;
    }
    else if (field == TEMPERATURE) {
      DENS_MAN * lambdaPowerFiltered = regulator_data("LambdaPowerFiltered",1);
      *lambdaPowerFiltered = target;
    }
    else {
      throw ATC_Error("KinetoThermostat::reset_lambda_contribution - invalid field given");
    }
  }

  //--------------------------------------------------------
  //  construct_methods:
  //    instantiations desired regulator method(s)
  
  //    dependence, but in general there is also a
  //    time integrator dependence.  In general the 
  //    precedence order is:
  //    time filter -> time integrator -> thermostat
  //    In the future this may need to be added if
  //    different types of time integrators can be
  //    specified.
  //--------------------------------------------------------
  void KinetoThermostat::construct_methods()
  {
    // get data associated with stages 1 & 2 of ATC_Method::initialize
    AtomicRegulator::construct_methods();

    if (atc_->reset_methods()) {
      // eliminate existing methods
      delete_method();

      // error check time integration methods
      TimeIntegrator::TimeIntegrationType myEnergyIntegrationType = (atc_->time_integrator(TEMPERATURE))->time_integration_type();
      TimeIntegrator::TimeIntegrationType myMomentumIntegrationType = (atc_->time_integrator(VELOCITY))->time_integration_type();
      if (myEnergyIntegrationType != TimeIntegrator::FRACTIONAL_STEP || myMomentumIntegrationType != TimeIntegrator::FRACTIONAL_STEP) {
        throw ATC_Error("KinetoThermostat::construct_methods - this scheme only valid with fractional step integration");
      }

      // update time filter
      TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
      if (timeFilterManager->need_reset() ) {
        timeFilter_ = timeFilterManager->construct(TimeFilterManager::EXPLICIT_IMPLICIT);
      }
      
      if (timeFilterManager->filter_dynamics()) {
        switch (regulatorTarget_) {
        case NONE: {
          regulatorMethod_ = new RegulatorMethod(this);
          break;
        }
        case FIELD: { // error check, rescale and filtering not supported together
          throw ATC_Error("KinetoThermostat::construct_methods - Cannot use rescaling thermostat with time filtering");
          break;
        }
        case DYNAMICS: {
        }
        default:
          throw ATC_Error("Unknown thermostat type in Thermostat::initialize");
        }
      }
      else {
        switch (regulatorTarget_) {
        case NONE: {
          regulatorMethod_ = new RegulatorMethod(this);
          break;
        }
        case FIELD: {
          if (atc_->temperature_def()==KINETIC) {
            regulatorMethod_ = new KinetoThermostatRescale(this,couplingMaxIterations_);
          }
          else if (atc_->temperature_def()==TOTAL) {
            regulatorMethod_ = new KinetoThermostatRescaleMixedKePe(this,couplingMaxIterations_);
          }
          else
            throw ATC_Error("Unknown temperature definition");
          break;
        }
        default:
          throw ATC_Error("Unknown thermostat target in Thermostat::initialize");
        }
      }
      
      AtomicRegulator::reset_method();
    }
    else {
      set_all_data_to_used();
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class VelocityRescaleCombined
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and kinetostat data
  //--------------------------------------------------------
  VelocityRescaleCombined::VelocityRescaleCombined(AtomicRegulator * kinetostat) :
    VelocityGlc(kinetostat),
    velocity_(atc_->field(VELOCITY)),
    thermostatCorrection_(NULL)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all data
  //--------------------------------------------------------
  void VelocityRescaleCombined::initialize()
  {
    VelocityGlc::initialize();
    thermostatCorrection_ = (atc_->interscale_manager()).dense_matrix("NodalAtomicFluctuatingMomentumRescaled");
  }

  //--------------------------------------------------------
  //  set_kinetostat_rhs
  //            sets up the right-hand side of the
  //            kinetostat equations
  //--------------------------------------------------------
  void VelocityRescaleCombined::set_kinetostat_rhs(DENS_MAT & rhs, double dt)
  {
    rhs = ((atc_->mass_mat_md(VELOCITY)).quantity())*(velocity_.quantity());
    rhs -= thermostatCorrection_->quantity();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatRescaleCombined
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ThermostatRescaleCombined::ThermostatRescaleCombined(AtomicRegulator * thermostat) :
    ThermostatRescale(thermostat),
    kinetostatCorrection_(NULL)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all data
  //--------------------------------------------------------
  void ThermostatRescaleCombined::initialize()
  {
    ThermostatRescale::initialize();
    kinetostatCorrection_ = (atc_->interscale_manager()).dense_matrix("NodalAtomicCombinedRescaleThermostatError");
  }

  //--------------------------------------------------------
  //  set_rhs:
  //    constructs the RHS vector with the target
  //    temperature
  //--------------------------------------------------------
  void ThermostatRescaleCombined::set_rhs(DENS_MAT & rhs)
  {
    ThermostatRescale::set_rhs(rhs);
    rhs -= kinetostatCorrection_->quantity();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class KinetoThermostatRescale
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  KinetoThermostatRescale::KinetoThermostatRescale(AtomicRegulator * kinetoThermostat,
                                                   int couplingMaxIterations) :
    KinetoThermostatShapeFunction(kinetoThermostat,couplingMaxIterations),
    atomVelocities_(NULL),
    nodalVelocities_(atc_->field(VELOCITY)),
    lambdaMomentum_(NULL),
    lambdaEnergy_(NULL),
    atomicFluctuatingVelocityRescaled_(NULL),
    atomicStreamingVelocity_(NULL),
    thermostat_(NULL),
    kinetostat_(NULL)
  {
    thermostat_ = this->construct_rescale_thermostat();
    kinetostat_ = new VelocityRescaleCombined(kinetoThermostat);
    // data associated with stage 3 in ATC_Method::initialize
    lambdaMomentum_ = atomicRegulator_->regulator_data(regulatorPrefix_+"LambdaMomentum",atc_->nsd());
    lambdaEnergy_ = atomicRegulator_->regulator_data(regulatorPrefix_+"LambdaEnergy",1); 
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  KinetoThermostatRescale::~KinetoThermostatRescale()
  {
    if (thermostat_) delete thermostat_;
    if (kinetostat_) delete kinetostat_;
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void KinetoThermostatRescale::construct_transfers()
  {
    // construct independent transfers first
    thermostat_->construct_transfers();
    kinetostat_->construct_transfers();

    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // get atom velocity data from manager
    atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY);

    // transfers requiring terms from both regulators
    // rescaled velocity fluctuations
    atomicFluctuatingVelocityRescaled_ = new AtomicFluctuatingVelocityRescaled(atc_);
    interscaleManager.add_per_atom_quantity(atomicFluctuatingVelocityRescaled_,
                                            "AtomFluctuatingVelocityRescaled");

    // streaming velocity component
    atomicStreamingVelocity_ = interscaleManager.per_atom_quantity("AtomLambdaMomentum");

    // rescaled momentum fluctuations, error term for kinetostat rhs
    PerAtomQuantity<double> * tempAtom = new AtomicMomentum(atc_,
                                                            atomicFluctuatingVelocityRescaled_);
    interscaleManager.add_per_atom_quantity(tempAtom,"AtomFluctuatingMomentumRescaled");
    DENS_MAN * tempNodes = new AtfShapeFunctionRestriction(atc_,
                                                           tempAtom,
                                                           interscaleManager.per_atom_sparse_matrix("Interpolant"));
    interscaleManager.add_dense_matrix(tempNodes,
                                       "NodalAtomicFluctuatingMomentumRescaled");

    // error term for thermostat rhs
    tempAtom = new AtomicCombinedRescaleThermostatError(atc_);
    interscaleManager.add_per_atom_quantity(tempAtom,"AtomCombinedRescaleThermostatError");
    tempNodes = new AtfShapeFunctionRestriction(atc_,
                                                tempAtom,
                                                interscaleManager.per_atom_sparse_matrix("Interpolant"));
    interscaleManager.add_dense_matrix(tempNodes,
                                       "NodalAtomicCombinedRescaleThermostatError");
  }

  //--------------------------------------------------------
  //  construct_rescale_thermostat
  //    constructs the appropriate rescaling thermostat
  //    varied through inheritance
  //--------------------------------------------------------
  ThermostatRescale * KinetoThermostatRescale::construct_rescale_thermostat()
  {
    return new ThermostatRescaleCombined(atomicRegulator_);
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all data
  //--------------------------------------------------------
  void KinetoThermostatRescale::initialize()
  {
    KinetoThermostatShapeFunction::initialize();
    thermostat_->initialize();
    kinetostat_->initialize();
  }

  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the thermostat in the post corrector phase
  //--------------------------------------------------------
  void KinetoThermostatRescale::apply_post_corrector(double dt)
  {
    // initial guesses
    lambdaMomentum_->set_quantity() = nodalVelocities_.quantity();
    lambdaEnergy_->set_quantity() = 1.;

    int iteration = 0;
    double eErr, pErr;
    while (iteration < couplingMaxIterations_) {
      _lambdaMomentumOld_ = lambdaMomentum_->quantity();
      _lambdaEnergyOld_ = lambdaEnergy_->quantity();

      // update thermostat
      thermostat_->compute_thermostat(dt);
      
      // update kinetostat
      kinetostat_->compute_kinetostat(dt);

      // check convergence
      _diff_ = lambdaEnergy_->quantity() - _lambdaEnergyOld_;
      eErr = _diff_.col_norm()/_lambdaEnergyOld_.col_norm();
      _diff_ = lambdaMomentum_->quantity() - _lambdaMomentumOld_;
      pErr = _diff_.col_norm()/_lambdaMomentumOld_.col_norm();
      if (eErr < tolerance_ && pErr < tolerance_) {
        break;
      }
      iteration++;
    }
    if (iteration == couplingMaxIterations_) {
      stringstream message;
      message << "WARNING: Iterative solve for lambda failed to converge after " << couplingMaxIterations_ << " iterations, final tolerance was " << std::max(eErr,pErr) << "\n";
      ATC::LammpsInterface::instance()->print_msg(message.str());
    }

    // application of rescaling lambda due
    apply_to_atoms(atomVelocities_);
  }

  //--------------------------------------------------------
  //  apply_lambda_to_atoms:
  //    applies the velocity rescale with an existing lambda
  //    note oldAtomicQuantity and dt are not used
  //--------------------------------------------------------
  void KinetoThermostatRescale::apply_to_atoms(PerAtomQuantity<double> * atomVelocities)
  {
    *atomVelocities = atomicFluctuatingVelocityRescaled_->quantity() + atomicStreamingVelocity_->quantity();
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void KinetoThermostatRescale::output(OUTPUT_LIST & outputData)
  {
    thermostat_->output(outputData);
    kinetostat_->output(outputData);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatRescaleMixedKePeCombined
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ThermostatRescaleMixedKePeCombined::ThermostatRescaleMixedKePeCombined(AtomicRegulator * thermostat) :
    ThermostatRescaleMixedKePe(thermostat),
    kinetostatCorrection_(NULL)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all data
  //--------------------------------------------------------
  void ThermostatRescaleMixedKePeCombined::initialize()
  {
    ThermostatRescaleMixedKePe::initialize();
    kinetostatCorrection_ = (atc_->interscale_manager()).dense_matrix("NodalAtomicCombinedRescaleThermostatError");
  }

  //--------------------------------------------------------
  //  set_rhs:
  //    constructs the RHS vector with the target
  //    temperature
  //--------------------------------------------------------
  void ThermostatRescaleMixedKePeCombined::set_rhs(DENS_MAT & rhs)
  {
    ThermostatRescaleMixedKePe::set_rhs(rhs);
    rhs -= kinetostatCorrection_->quantity();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class KinetoThermostatRescaleMixedKePe
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  KinetoThermostatRescaleMixedKePe::KinetoThermostatRescaleMixedKePe(AtomicRegulator * kinetoThermostat,
                                                                     int couplingMaxIterations) :
    KinetoThermostatRescale(kinetoThermostat,couplingMaxIterations)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  construct_rescale_thermostat
  //    constructs the appropriate rescaling thermostat
  //    varied through inheritance
  //--------------------------------------------------------
  ThermostatRescale * KinetoThermostatRescaleMixedKePe::construct_rescale_thermostat()
  {
    return new ThermostatRescaleMixedKePeCombined(atomicRegulator_);
  }


  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class KinetoThermostatGlcFs
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  KinetoThermostatGlcFs::KinetoThermostatGlcFs(AtomicRegulator * kinetoThermostat,
                                               int couplingMaxIterations,
                                               const string & regulatorPrefix) :
    KinetoThermostatShapeFunction(kinetoThermostat,couplingMaxIterations,regulatorPrefix),
    velocity_(atc_->field(VELOCITY)),
    temperature_(atc_->field(TEMPERATURE)),
    timeFilter_(atomicRegulator_->time_filter()),
    nodalAtomicLambdaForce_(NULL),
    lambdaForceFiltered_(NULL),
    nodalAtomicLambdaPower_(NULL),
    lambdaPowerFiltered_(NULL),
    atomRegulatorForces_(NULL),
    atomThermostatForces_(NULL),
    atomMasses_(NULL),
    atomVelocities_(NULL),
    isFirstTimestep_(true),
    nodalAtomicMomentum_(NULL),
    nodalAtomicEnergy_(NULL),
    atomPredictedVelocities_(NULL),
    nodalAtomicPredictedMomentum_(NULL),
    nodalAtomicPredictedEnergy_(NULL),
    firstHalfAtomForces_(NULL),
    dtFactor_(0.)
  {
    // construct/obtain data corresponding to stage 3 of ATC_Method::initialize
    //nodalAtomicLambdaPower_ = thermostat->regulator_data(regulatorPrefix_+"NodalAtomicLambdaPower",1);
    //lambdaPowerFiltered_ = thermostat_->regulator_data(regulatorPrefix_+"LambdaPowerFiltered",1);
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void KinetoThermostatGlcFs::construct_transfers()
  {
    // BASES INITIALIZE
    // GRAB ANY COPIES FROM BASES

    // TOTAL REGULATOR FORCE
    // TOTAL FIRST HALF FORCE, IF NECESSARY
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void KinetoThermostatGlcFs::initialize()
  {
    RegulatorMethod::initialize();

    // INITIALIZE BASES

    // MAKE SURE ANY NEEDED POINTERS FROM BASES ARE COPIED BY HERE
  }

  //--------------------------------------------------------
  //  apply_to_atoms:
  //     determines what if any contributions to the
  //     atomic moition is needed for
  //     consistency with the thermostat
  //     and computes the instantaneous induced power
  //--------------------------------------------------------
  void KinetoThermostatGlcFs::apply_to_atoms(PerAtomQuantity<double> * atomicVelocity,
                                             const DENS_MAN * nodalAtomicMomentum,
                                             const DENS_MAN * nodalAtomicEnergy,
                                             const DENS_MAT & lambdaForce,
                                             DENS_MAT & nodalAtomicLambdaForce,
                                             DENS_MAT & nodalAtomicLambdaPower,
                                             double dt)
  {
    // compute initial contributions to lambda force and power
    nodalAtomicLambdaPower = nodalAtomicEnergy->quantity();
    nodalAtomicLambdaPower *= -1.;
    nodalAtomicLambdaForce = nodalAtomicMomentum->quantity();
    nodalAtomicLambdaForce *= -1.;

    // apply lambda force to atoms
    _velocityDelta_ = lambdaForce;
    _velocityDelta_ /= atomMasses_->quantity();
    _velocityDelta_ *= dt;
    (*atomicVelocity) += _velocityDelta_;

    // finalize lambda force and power
    nodalAtomicLambdaForce += nodalAtomicMomentum->quantity();
    nodalAtomicLambdaPower += nodalAtomicEnergy->quantity();
  }

  //--------------------------------------------------------
  //  full_prediction:
  //    flag to perform a full prediction calcalation
  //    for lambda rather than using the old value
  //--------------------------------------------------------
  bool KinetoThermostatGlcFs::full_prediction()
  {
    // CHECK BOTH BASES
    return false;
  }

  //--------------------------------------------------------
  //  apply_predictor:
  //    apply the thermostat to the atoms in the first step
  //    of the Verlet algorithm
  //--------------------------------------------------------
  void KinetoThermostatGlcFs::apply_pre_predictor(double dt)
  {
    DENS_MAT & lambdaForceFiltered(lambdaForceFiltered_->set_quantity());
    DENS_MAT & nodalAtomicLambdaForce(nodalAtomicLambdaForce_->set_quantity());
    DENS_MAT & myLambdaPowerFiltered(lambdaPowerFiltered_->set_quantity());
    DENS_MAT & myNodalAtomicLambdaPower(nodalAtomicLambdaPower_->set_quantity());

    // update filtered forces power, equivalent to measuring changes in momentum and energy
    timeFilter_->apply_pre_step1(lambdaForceFiltered,nodalAtomicLambdaForce,dt);
    timeFilter_->apply_pre_step1(myLambdaPowerFiltered,myNodalAtomicLambdaPower,dt);

    // apply lambda force to atoms and compute instantaneous lambda power for first half of time step
    this->apply_to_atoms(atomVelocities_,nodalAtomicMomentum_,nodalAtomicEnergy_,
                         firstHalfAtomForces_->quantity(),
                         nodalAtomicLambdaForce,myNodalAtomicLambdaPower,0.5*dt);

    // update nodal variables for first half of time step
    // velocity
    this->add_to_momentum(nodalAtomicLambdaForce,deltaMomentum_,0.5*dt);
    atc_->apply_inverse_mass_matrix(deltaMomentum_,VELOCITY);
    velocity_ += deltaMomentum_;
    // temperature
    this->add_to_energy(myNodalAtomicLambdaPower,deltaEnergy1_,0.5*dt);

    // start update of filtered lambda force and power using temporary (i.e., 0 valued) quantities for first part of update
    nodalAtomicLambdaForce = 0.;
    timeFilter_->apply_post_step1(lambdaForceFiltered,nodalAtomicLambdaForce,dt);
    myNodalAtomicLambdaPower = 0.;
    timeFilter_->apply_post_step1(myLambdaPowerFiltered,myNodalAtomicLambdaPower,dt);
  }

  //--------------------------------------------------------
  //  apply_pre_corrector:
  //    apply the thermostat to the atoms in the first part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void KinetoThermostatGlcFs::apply_pre_corrector(double dt)
  {
    // CHECK WHEN CREATING PREDICTED VELOCITIES IN BASE REGULATORS, ONLY NEED ONE
    (*atomPredictedVelocities_) = atomVelocities_->quantity();

    // do full prediction if we just redid the shape functions
    if (full_prediction()) {
      this->compute_lambda(dt);
      
      atomThermostatForces_->unfix_quantity();  // allow update of atomic force applied by lambda
    }

    // apply lambda force to atoms and compute instantaneous lambda power to predict second half of time step
    DENS_MAT & myNodalAtomicLambdaForce(nodalAtomicLambdaForce_->set_quantity());
    DENS_MAT & myNodalAtomicLambdaPower(nodalAtomicLambdaPower_->set_quantity());
    apply_to_atoms(atomPredictedVelocities_,
                   nodalAtomicPredictedMomentum_,nodalAtomicPredictedEnergy_,
                   firstHalfAtomForces_->quantity(),
                   myNodalAtomicLambdaForce,myNodalAtomicLambdaPower,0.5*dt);

    if (full_prediction())
      atomThermostatForces_->fix_quantity();
    
    // SPLIT OUT FUNCTION TO CREATE DELTA VARIABLES IN BASES, ONLY NEED THESE
    // update predicted nodal variables for second half of time step
    // velocity
    this->add_to_momentum(myNodalAtomicLambdaForce,deltaMomentum_,0.5*dt);
    atc_->apply_inverse_mass_matrix(deltaMomentum_,VELOCITY);
    velocity_ += deltaMomentum_;
    // temperature
    this->add_to_energy(myNodalAtomicLambdaPower,deltaEnergy2_,0.5*dt);
    // following manipulations performed this way for efficiency
    deltaEnergy1_ += deltaEnergy2_;
    atc_->apply_inverse_mass_matrix(deltaEnergy1_,TEMPERATURE);
    temperature_ += deltaEnergy1_;
  }

  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the thermostat to the atoms in the second part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void KinetoThermostatGlcFs::apply_post_corrector(double dt)
  {
    // remove predicted power effects
    // velocity
    DENS_MAT & myVelocity(velocity_.set_quantity());
    myVelocity -= deltaMomentum_;
    // temperature
    DENS_MAT & myTemperature(temperature_.set_quantity());
    atc_->apply_inverse_mass_matrix(deltaEnergy2_,TEMPERATURE);
    myTemperature -= deltaEnergy2_;

    // set up equation and update lambda
    this->compute_lambda(dt);

    // apply lambda force to atoms and compute instantaneous lambda power for second half of time step
    DENS_MAT & nodalAtomicLambdaForce(nodalAtomicLambdaForce_->set_quantity());
    DENS_MAT & myNodalAtomicLambdaPower(nodalAtomicLambdaPower_->set_quantity());
    // allow computation of force applied by lambda using current velocities
    atomThermostatForces_->unfix_quantity();
    atomThermostatForces_->quantity();
    atomThermostatForces_->fix_quantity();
    apply_to_atoms(atomVelocities_,nodalAtomicMomentum_,nodalAtomicEnergy_,
                   atomRegulatorForces_->quantity(),
                   nodalAtomicLambdaForce,myNodalAtomicLambdaPower,0.5*dt);

    // finalize filtered lambda force and power by adding latest contribution
    timeFilter_->apply_post_step2(lambdaForceFiltered_->set_quantity(),
                                  nodalAtomicLambdaForce,dt);
    timeFilter_->apply_post_step2(lambdaPowerFiltered_->set_quantity(),
                                  myNodalAtomicLambdaPower,dt);

    // update nodal variables for second half of time step
    // velocity
    this->add_to_momentum(nodalAtomicLambdaForce,deltaMomentum_,0.5*dt);
    atc_->apply_inverse_mass_matrix(deltaMomentum_,VELOCITY);
    velocity_ += deltaMomentum_;
    // temperature
    this->add_to_energy(myNodalAtomicLambdaPower,deltaEnergy2_,0.5*dt);
    atc_->apply_inverse_mass_matrix(deltaEnergy2_,TEMPERATURE);
    myTemperature += deltaEnergy2_;
    
    
    isFirstTimestep_ = false;
  }

  //--------------------------------------------------------
  //  compute_lambda:
  //   sets up and solves linear system for lambda, if the
  //   bool is true it iterators to a non-linear solution
  //--------------------------------------------------------
  void KinetoThermostatGlcFs::compute_lambda(double dt,
                                       bool iterateSolution)
  {
    // ITERATIVE SOLUTION
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void KinetoThermostatGlcFs::output(OUTPUT_LIST & outputData)
  {
    // DO NOT CALL INDIVIDUAL REGULATORS
    // OUTPUT TOTAL FORCE AND TOTAL POWER 
    // OUTPUT EACH LAMBDA
  }


};
