// ATC transfer headers
#include "ThermalTimeIntegrator.h"
#include "TransferOperator.h"
#include "ATC_Coupling.h"
#include "TimeFilter.h"
#include "ATC_Error.h"
#include "PerAtomQuantityLibrary.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermalTimeIntegrator
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //-------------------------------------------------------- 
  ThermalTimeIntegrator::ThermalTimeIntegrator(ATC_Coupling * atc,
                                               TimeIntegrationType timeIntegrationType) :
    TimeIntegrator(atc, timeIntegrationType)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  modify
  //    parses inputs and modifies state of the integrator
  //--------------------------------------------------------
  bool ThermalTimeIntegrator::modify(int /* narg */, char **arg)
  {
    bool foundMatch = false;
    int argIndex = 0;

    // time integration scheme
    /*! \page man_thermal_time_integration fix_modify AtC time_integration (thermal)
      \section syntax
      fix_modify AtC time_integration <descriptor> \n
      - descriptor (string) = time integration type  \n
      
      various time integration methods for the finite elements\n
      \section description
      gear - atomic velocity update with 2nd order Verlet, nodal temperature update with 3rd or 4th order Gear, thermostats based on controlling power \n
      fractional_step - atomic velocity update with 2nd order Verlet, mixed nodal temperature update, 3/4 Gear for continuum and 2 Verlet for atomic contributions, thermostats based on controlling discrete energy changes\n
      \section examples
      <TT> fix_modify atc time_integration gear </TT> \n
      <TT> fix_modify atc time_integration fractional_step </TT> \n
      \section description
      \section related
      see \ref man_fix_atc
      \section default
      none 
    */
    if (strcmp(arg[argIndex],"gear")==0) {
      timeIntegrationType_ = GEAR;
      needReset_ = true;
      foundMatch = true;
    }
    else if (strcmp(arg[argIndex],"fractional_step")==0) {
      timeIntegrationType_ = FRACTIONAL_STEP;
      needReset_ = true;
      foundMatch = true;
    }
    return foundMatch;
  }

  //--------------------------------------------------------
  //  construct_methods
  //    creates algorithm objects
  //--------------------------------------------------------
  void ThermalTimeIntegrator::construct_methods()
  {
    if (atc_->reset_methods()) {
      if (timeIntegrationMethod_) delete timeIntegrationMethod_;

      if (timeFilterManager_->need_reset()) {
        switch (timeIntegrationType_) {
          case GEAR: {
            timeFilter_ = timeFilterManager_->construct(TimeFilterManager::IMPLICIT);
            atc_->set_mass_mat_time_filter(TEMPERATURE,TimeFilterManager::EXPLICIT);
            break;
          }
          case FRACTIONAL_STEP: {
            timeFilter_ = timeFilterManager_->construct(TimeFilterManager::EXPLICIT_IMPLICIT);
            atc_->set_mass_mat_time_filter(TEMPERATURE,TimeFilterManager::EXPLICIT_IMPLICIT);
            break;
          }
          default:
            throw ATC_Error("Unknown time integration type in ThermalTimeIntegrator::Initialize()");
        }
      }
      
      if (timeFilterManager_->filter_dynamics()) {
        switch (timeIntegrationType_) {
          case GEAR: {
            timeIntegrationMethod_ = new ThermalTimeIntegratorGearFiltered(this);
            break; 
          }
          case FRACTIONAL_STEP: {
            timeIntegrationMethod_ = new ThermalTimeIntegratorFractionalStepFiltered(this);
            break;
          }
        default:
          throw ATC_Error("Unknown time integration type in ThermalTimeIntegrator::Initialize()");
        }
      }
      else {
        switch (timeIntegrationType_) {
          case GEAR: {
            timeIntegrationMethod_ = new ThermalTimeIntegratorGear(this);
            break;
          }
          case FRACTIONAL_STEP: {
            timeIntegrationMethod_ = new ThermalTimeIntegratorFractionalStep(this);
            break;
          }
        default:
          throw ATC_Error("Unknown time integration type in ThermalTimeIntegrator::Initialize()");
        }
      }
    }
  }

  //--------------------------------------------------------
  //  pack_fields
  //    add persistent variables to data list
  //--------------------------------------------------------
  void ThermalTimeIntegrator::pack_fields(RESTART_LIST & data)
  {
    data["NodalAtomicPowerFiltered"] = & nodalAtomicPowerFiltered_.set_quantity();
    data["NodalAtomicEnergyFiltered"] = & nodalAtomicEnergyFiltered_.set_quantity();
    TimeIntegrator::pack_fields(data);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermalIntegrationMethod
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //        Grab data from ATC
  //--------------------------------------------------------
  ThermalIntegrationMethod::ThermalIntegrationMethod(ThermalTimeIntegrator * thermalTimeIntegrator) :
    TimeIntegrationMethod(thermalTimeIntegrator),
    timeFilter_(thermalTimeIntegrator->time_filter()),
    temperature_(atc_->field(TEMPERATURE)),
    temperatureRoc_(atc_->field_roc(TEMPERATURE)),
    temperature2Roc_(atc_->field_2roc(TEMPERATURE)),
    nodalAtomicTemperatureOut_(atc_->nodal_atomic_field(TEMPERATURE)),
    nodalAtomicTemperature_(NULL),
    temperatureRhs_(atc_->field_rhs(TEMPERATURE)),
    nodalAtomicPowerOut_(atc_->nodal_atomic_field_roc(TEMPERATURE))
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  construct_transfers
  //        Grab existing managed quantities,
  //        create the rest
  //--------------------------------------------------------
  void ThermalIntegrationMethod::construct_transfers()
  {
    nodalAtomicTemperature_ =
      (atc_->interscale_manager()).dense_matrix("NodalAtomicTemperature");
  }
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermalIntegratorGear
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ThermalTimeIntegratorGear::ThermalTimeIntegratorGear(ThermalTimeIntegrator * thermalTimeIntegrator) :
    ThermalIntegrationMethod(thermalTimeIntegrator),
    nodalAtomicPowerFiltered_(thermalTimeIntegrator->nodal_atomic_power_filtered())
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  construct_transfers
  //        Grab existing managed quantities,
  //        create the rest
  //--------------------------------------------------------
  void ThermalTimeIntegratorGear::construct_transfers()
  {
    ThermalIntegrationMethod::construct_transfers();
    InterscaleManager & interscaleManager = atc_->interscale_manager();
        
    // add in power computation
    DotTwiceKineticEnergy * dotTwiceKineticEnergy = 
      new DotTwiceKineticEnergy(atc_);
    interscaleManager.add_per_atom_quantity(dotTwiceKineticEnergy,"DotTwiceKineticEnergy");
    nodalAtomicPower_ = new AtfShapeFunctionRestriction(atc_,
                                                        dotTwiceKineticEnergy,
                                                        interscaleManager.per_atom_sparse_matrix("Interpolant"));
    interscaleManager.add_dense_matrix(nodalAtomicPower_,"NodalAtomicPower");
  }

  //--------------------------------------------------------
  //  initialize
  //        initialize all data
  //--------------------------------------------------------
  void ThermalTimeIntegratorGear::initialize()
  {
    ThermalIntegrationMethod::initialize();

    // sets up time filter for cases where variables temporally filtered
    // this time integrator should use an implicit filter
    TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
    if (timeFilterManager->need_reset()) {
      // Some time filters need the old value for the power
      timeFilter_->initialize(nodalAtomicPower_->quantity());
    }

    if (!timeFilterManager->end_equilibrate()) {
      nodalAtomicPowerFiltered_.reset(atc_->num_nodes(),1);
    }

    if (!timeFilterManager->filter_dynamics()) {
      temperatureRhs_ = nodalAtomicPower_->quantity();
    }
  }

  //--------------------------------------------------------
  //  pre_initial_integrate2
  ///   time integration before Verlet step 1
  //--------------------------------------------------------
  void ThermalTimeIntegratorGear::pre_initial_integrate2(double dt)
  {
    // Predict nodal temperatures and time derivatives based on FE data
    // use 3rd order Gear
    gear1_3_predict(temperature_.set_quantity(),
                    temperatureRoc_.set_quantity(),
                    temperature2Roc_.quantity(),dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate1
  //    time integration after Verlet step 2
  //--------------------------------------------------------
  void ThermalTimeIntegratorGear::post_final_integrate1(double dt)
  {
    const DENS_MAT & myNodalAtomicPower(nodalAtomicPower_->quantity());
    timeFilter_->apply_post_step2(nodalAtomicPowerFiltered_.set_quantity(),
                                  myNodalAtomicPower,dt);
    temperatureRhs_ += myNodalAtomicPower;
    
    // Finish updating temperature
    _temperatureResidual_.resize(atc_->num_nodes(),1);
    atc_->apply_inverse_mass_matrix(temperatureRhs_.quantity(),
                                    _temperatureResidual_,
                                    TEMPERATURE);
    _temperatureResidual_ -= temperatureRoc_.quantity();
    _temperatureResidual_ *= dt;
    
    gear1_3_correct(temperature_.set_quantity(),
                    temperatureRoc_.set_quantity(),
                    temperature2Roc_.set_quantity(),
                    _temperatureResidual_,dt);
  }

  //--------------------------------------------------------
  //  post_process
  //    do any post-processing calculations required for
  //    output phase
  //--------------------------------------------------------
  void ThermalTimeIntegratorGear::post_process()
  {
    nodalAtomicPowerOut_ = nodalAtomicPower_->quantity();
    nodalAtomicTemperatureOut_ = nodalAtomicTemperature_->quantity();
  }

  //--------------------------------------------------------
  //  finish
  //    finalize state of nodal atomic quantities
  //--------------------------------------------------------
  void ThermalTimeIntegratorGear::finish()
  {
    post_process();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermalTimeIntegratorGearFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //        Grab data from ATC
  //--------------------------------------------------------

  ThermalTimeIntegratorGearFiltered::ThermalTimeIntegratorGearFiltered(ThermalTimeIntegrator * thermalTimeIntegrator) :
    ThermalTimeIntegratorGear(thermalTimeIntegrator),
    temperature3Roc_(atc_->field_3roc(TEMPERATURE))
  {
    // do nothing
    
    //      specifically if history data is required and we need another time filter object for the fields
  }

  //--------------------------------------------------------
  //  pre_initial_integrate2
  //    time integration before Verlet step 1
  //--------------------------------------------------------
  void ThermalTimeIntegratorGearFiltered::pre_initial_integrate2(double dt)
  {
    // Predict nodal temperatures and time derivatives based on FE data
    // use 3rd order Gear
    gear1_4_predict(temperature_.set_quantity(),
                    temperatureRoc_.set_quantity(),
                    temperature2Roc_.set_quantity(),
                    temperature3Roc_.quantity(),dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate1
  //    first time integration computations 
  //    after Verlet step 2
  //--------------------------------------------------------
  void ThermalTimeIntegratorGearFiltered::post_final_integrate1(double dt)
  {
    DENS_MAT & myNodalAtomicPowerFiltered(nodalAtomicPowerFiltered_.set_quantity());
    timeFilter_->apply_post_step2(myNodalAtomicPowerFiltered,nodalAtomicPower_->quantity(),dt);
    temperatureRhs_ += myNodalAtomicPowerFiltered;
    
    // Finish updating temperature
    _temperatureResidual_.resize(atc_->num_nodes(),1);
    atc_->apply_inverse_mass_matrix(temperatureRhs_.quantity(),
                                            _temperatureResidual_,
                                            TEMPERATURE);
    _temperatureResidual_ -= temperatureRoc_.quantity();
    _temperatureResidual_ *= dt;
    gear1_4_correct(temperature_.set_quantity(),
                    temperatureRoc_.set_quantity(),
                    temperature2Roc_.set_quantity(),
                    temperature3Roc_.set_quantity(),
                    _temperatureResidual_,dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate3
  //    third time integration computations
  //    after Verlet step 2
  //--------------------------------------------------------
  void ThermalTimeIntegratorGearFiltered::post_final_integrate3(double dt)
  {
    // update filtered atomic temperature
    timeFilter_->apply_post_step2(nodalAtomicTemperatureOut_.set_quantity(),
                                  nodalAtomicTemperature_->quantity(),dt);
  }

  //--------------------------------------------------------
  //  post_process
  //    do any post-processing calculations required for
  //    output phase
  //--------------------------------------------------------
  void ThermalTimeIntegratorGearFiltered::post_process()
  {
    nodalAtomicPowerOut_ = nodalAtomicPowerFiltered_.quantity();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermalIntegratorFractionalStep
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //        Grab data from ATC
  //--------------------------------------------------------

  ThermalTimeIntegratorFractionalStep::ThermalTimeIntegratorFractionalStep(ThermalTimeIntegrator * thermalTimeIntegrator) :
    ThermalIntegrationMethod(thermalTimeIntegrator),
    nodalAtomicEnergyFiltered_(thermalTimeIntegrator->nodal_atomic_energy_filtered()),
    nodalAtomicPowerFiltered_(thermalTimeIntegrator->nodal_atomic_power_filtered()),
    atomicTemperatureDelta_(atc_->num_nodes(),1),
    nodalAtomicEnergy_(NULL),
    nodalAtomicEnergyOld_(atc_->num_nodes(),1),
    nodalAtomicTemperatureOld_(atc_->num_nodes(),1)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  construct_transfers
  //        Grab existing managed quantities,
  //        create the rest
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStep::construct_transfers()
  {
    ThermalIntegrationMethod::construct_transfers();
    InterscaleManager & interscaleManager(atc_->interscale_manager());
    nodalAtomicEnergy_ = interscaleManager.dense_matrix("NodalAtomicEnergy");
  }

  //--------------------------------------------------------
  //  initialize
  //        initialize all data
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStep::initialize()
  {
    ThermalIntegrationMethod::initialize();

    // initial power to zero
    nodalAtomicPower_.reset(atc_->num_nodes(),1);

    // sets up time filter for cases where variables temporally filtered
    // this time integrator should use Crank-Nicholson filter for 2nd order accuracy
    TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
    if (timeFilterManager->need_reset()) {
      // the form of this integrator implies no time filters that require history data can be used
      timeFilter_->initialize();
    }

    // sets up time filter for post-processing the filtered power
    // this time integrator should use an explicit-implicit filter
    // to mirror the 2nd order Verlet integration scheme
    // It requires no history information so initial value just sizes arrays
    if (!timeFilterManager->end_equilibrate()) {
      // implies an initial condition of the instantaneous atomic energy
      // for the corresponding filtered variable, consistent with the temperature
      nodalAtomicEnergyFiltered_ = nodalAtomicEnergy_->quantity();  
      nodalAtomicPowerFiltered_.reset(atc_->num_nodes(),1);
    }
  }

  //--------------------------------------------------------
  //  pre_initial_integrate1
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStep::pre_initial_integrate1(double dt)
  {
    const DENS_MAT & myNodalAtomicEnergy(nodalAtomicEnergy_->quantity());
    // updated filtered energy using explicit-implicit scheme
    timeFilter_->apply_pre_step1(nodalAtomicEnergyFiltered_.set_quantity(),
                                 myNodalAtomicEnergy,dt);
  }

  //--------------------------------------------------------
  //  pre_initial_integrate2
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStep::pre_initial_integrate2(double dt)
  {
    // used for updating change in temperature from mass matrix change
    this->compute_old_time_data();

    // update FE contributions
    apply_gear_predictor(dt);

    // update filtered nodal atomic power
    
    //      that way thermostat and integrator can be consistent
    timeFilter_->apply_pre_step1(nodalAtomicPowerFiltered_.set_quantity(),
                                 nodalAtomicPower_,dt);

    // store current energy for use later
    nodalAtomicPower_ = nodalAtomicEnergy_->quantity();
    nodalAtomicPower_ *= -1.;
  }

  //--------------------------------------------------------
  //  pre_final_integrate1
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStep::pre_final_integrate1(double dt)
  {
    
    //      before the new rhs is computed but after atomic velocity is updated
    //      to allow for general notions of temperature beyond kinetic.
    // compute change in restricted atomic energy
    nodalAtomicPower_ += nodalAtomicEnergy_->quantity();
    
    // update FE temperature with change in temperature from MD
    compute_temperature_delta(nodalAtomicPower_,dt);
    temperature_ += atomicTemperatureDelta_.quantity();
    
    // approximation to power for output
    nodalAtomicPower_ /= dt;
    timeFilter_->apply_post_step1(nodalAtomicPowerFiltered_.set_quantity(),
                                  nodalAtomicPower_,dt);

    // make sure nodes are fixed
    atc_->set_fixed_nodes();
  }

  //--------------------------------------------------------
  //  post_final_integrate1
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStep::post_final_integrate1(double dt)
  {
    // Finish updating temperature with FE contributions
    atc_->apply_inverse_mass_matrix(temperatureRhs_.quantity(),
                                    _temperatureResidual_,TEMPERATURE);
    _temperatureResidual_ -= temperatureRoc_.quantity();
    _temperatureResidual_ *= dt;
    apply_gear_corrector(_temperatureResidual_,dt);
  }

  //--------------------------------------------------------
  //  post_final_integrate3
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStep::post_final_integrate3(double dt)
  {
    // update filtered atomic energy
    timeFilter_->apply_post_step1(nodalAtomicEnergyFiltered_.set_quantity(),
                                  nodalAtomicEnergy_->quantity(),dt);
  }

  //--------------------------------------------------------
  //  output
  //    add variables to output list
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStep::post_process()
  {
    nodalAtomicPowerOut_ = nodalAtomicPower_;
    nodalAtomicTemperatureOut_ = nodalAtomicTemperature_->quantity();
  }

  //--------------------------------------------------------
  //  finish
  //    finalize state of nodal atomic quantities
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStep::finish()
  {
    post_process();
  }

  //--------------------------------------------------------
  //  apply_gear_predictor
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStep::apply_gear_predictor(double dt)
  {
    gear1_3_predict(temperature_.set_quantity(),
                    temperatureRoc_.set_quantity(),
                    temperature2Roc_.quantity(),dt);
  }

  //--------------------------------------------------------
  //  apply_gear_corrector
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStep::apply_gear_corrector(const DENS_MAT & R_theta, double dt)
  {
    gear1_3_correct(temperature_.set_quantity(),
                    temperatureRoc_.set_quantity(),
                    temperature2Roc_.set_quantity(),
                    R_theta,dt);
  }

  //--------------------------------------------------------
  //  compute_old_time_data
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStep::compute_old_time_data()
  {
    const DENS_MAT & myNodalAtomicEnergy(nodalAtomicEnergy_->quantity());
    atc_->apply_inverse_mass_matrix(myNodalAtomicEnergy,
                                    nodalAtomicTemperatureOld_.set_quantity(),
                                    TEMPERATURE);
    nodalAtomicEnergyOld_ = myNodalAtomicEnergy;
  }

  //--------------------------------------------------------
  //  compute_temperature_delta
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStep::compute_temperature_delta(const DENS_MAT & atomicEnergyDelta,
                                                                      double /* dt */)
  {
    DENS_MAT & myAtomicTemperatureDelta(atomicTemperatureDelta_.set_quantity());
    myAtomicTemperatureDelta = nodalAtomicEnergyOld_.quantity() + atomicEnergyDelta;
    atc_->apply_inverse_mass_matrix(myAtomicTemperatureDelta,
                                    TEMPERATURE);
    myAtomicTemperatureDelta += -1.*(nodalAtomicTemperatureOld_.quantity());
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermalTimeIntegratorFracionalStepFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //        Grab data from ATC
  //--------------------------------------------------------

  ThermalTimeIntegratorFractionalStepFiltered::ThermalTimeIntegratorFractionalStepFiltered(ThermalTimeIntegrator * thermalTimeIntegrator) :
    ThermalTimeIntegratorFractionalStep(thermalTimeIntegrator),
    temperature3Roc_(atc_->field_3roc(TEMPERATURE))
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ThermalTimeIntegratorFractionalStepFiltered::~ThermalTimeIntegratorFractionalStepFiltered()
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  pre_initial_integrate1
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStepFiltered::pre_initial_integrate1(double dt)
  {
    // determine change in temperature if no forces were applied over this timestep

    // relevant coefficients from time filter
    double coefF1 = timeFilter_->filtered_coefficient_pre_s1(dt);
    double coefF2 = timeFilter_->filtered_coefficient_post_s1(dt);
    double coefU1 = timeFilter_->unfiltered_coefficient_pre_s1(dt);
    double coefU2 = timeFilter_->unfiltered_coefficient_post_s1(dt);

    DENS_MAT & myAtomicTemperatureDelta(atomicTemperatureDelta_.set_quantity());
    DENS_MAT & myNodalAtomicEnergyFiltered(nodalAtomicEnergyFiltered_.set_quantity());
    const DENS_MAT & myNodalAtomicEnergy(nodalAtomicEnergy_->quantity());

    // composite from change after two step update of current filtered energy
    myAtomicTemperatureDelta = (coefF1*coefF2-1.)*myNodalAtomicEnergyFiltered;
    // change in filtered temperature from current energy from this and next time levels
    myAtomicTemperatureDelta += (coefU1*coefF2+coefU2)*myNodalAtomicEnergy;

    // updated filtered energy using explicit-implicit scheme
    // nodalAtomicEnergy_ is either set from initialization or from the end of the last timestep
    timeFilter_->apply_pre_step1(myNodalAtomicEnergyFiltered,myNodalAtomicEnergy,dt);
  }

  //--------------------------------------------------------
  //  output
  //    add variables to output list
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStepFiltered::output(OUTPUT_LIST & outputData)
  {
    atc_->apply_inverse_md_mass_matrix(nodalAtomicEnergyFiltered_.quantity(),
                                       nodalAtomicTemperatureOut_.set_quantity(),
                                       TEMPERATURE);
    DENS_MAT & nodalAtomicPower(nodalAtomicPowerFiltered_.set_quantity());
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData["NodalAtomicPower"] = &nodalAtomicPower;
    }
  }

  //--------------------------------------------------------
  //  apply_gear_predictor
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStepFiltered::apply_gear_predictor(double dt)
  {
    gear1_4_predict(temperature_.set_quantity(),
                    temperatureRoc_.set_quantity(),
                    temperature2Roc_.set_quantity(),
                    temperature3Roc_.quantity(),dt);
  }

  //--------------------------------------------------------
  //  apply_gear_corrector
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStepFiltered::apply_gear_corrector(const DENS_MAT & R_theta, double dt)
  {
    gear1_4_correct(temperature_.set_quantity(),
                    temperatureRoc_.set_quantity(),
                    temperature2Roc_.set_quantity(),
                    temperature3Roc_.set_quantity(),
                    R_theta,dt);
  }

  //--------------------------------------------------------
  //  compute_temperature_delta
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStepFiltered::compute_old_time_data()
  {
    const DENS_MAT & myNodalAtomicEnergyFiltered(nodalAtomicEnergyFiltered_.quantity());
    atc_->apply_inverse_mass_matrix(myNodalAtomicEnergyFiltered,
                                    nodalAtomicTemperatureOld_.set_quantity(),
                                    TEMPERATURE);
    nodalAtomicEnergyOld_ = myNodalAtomicEnergyFiltered;
  }

  //--------------------------------------------------------
  //  compute_old_time_data
  //--------------------------------------------------------
  void ThermalTimeIntegratorFractionalStepFiltered::compute_temperature_delta(const DENS_MAT & atomicEnergyDelta,
                                                                              double dt)
  {
    DENS_MAT & myAtomicTemperatureDelta(atomicTemperatureDelta_.set_quantity());
    double coefU2 = timeFilter_->unfiltered_coefficient_post_s1(dt);
    myAtomicTemperatureDelta += nodalAtomicEnergyOld_.quantity() + coefU2*atomicEnergyDelta;
    atc_->apply_inverse_mass_matrix(myAtomicTemperatureDelta,
                                    TEMPERATURE);
    myAtomicTemperatureDelta += -1.*nodalAtomicTemperatureOld_.quantity();
  }
};
