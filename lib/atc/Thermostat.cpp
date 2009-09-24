// ATC_Transfer Headers
#include "Thermostat.h"
#include "CG.h"
#include "ATC_Error.h"
#include "PrescribedDataManager.h"
#include "TimeIntegrator.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class Thermostat
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  Thermostat::Thermostat(ATC_Transfer * atcTransfer) :
    AtomicRegulator(atcTransfer),
    thermostatType_(NONE)
  {
    // nothing to do
  }

  //--------------------------------------------------------
  //  modify:
  //    parses and adjusts thermostat state based on
  //    user input, in the style of LAMMPS user input
  //--------------------------------------------------------
  bool Thermostat::modify(int narg, char **arg)
  {
    bool foundMatch = false;
    
    // thermostat type
    /*! \page man_thermal_control fix_modify AtC transfer thermal control
      \section syntax
      fix_modify AtC transfer thermal control <control_type> <optional args>\n
	- control_type (string) = none | rescale | hoover | flux\n
      
      fix_modify AtC transfer thermal control rescale <frequency>\n
        - frequency (int) = time step frequency for applying velocity rescaling \n

      fix_modify AtC transfer thermal control hoover \n

      fix_modify AtC transfer thermal control flux <boundary_integration_type> <face_set_id(optional)>\n
        - boundary_integration_type (string) = faceset | interpolate\n
        - face_set_id (string), optional = id of boundary face set, if not specified
      (or not possible when the atomic domain does not line up with 
      mesh boundaries) defaults to an atomic-quadrature approximate 
      evaulation, does not work with interpolate\n
      \section examples
      <TT> fix_modify AtC transfer thermal control none </TT> \n
      <TT> fix_modify AtC transfer thermal control rescale 10 </TT> \n
      <TT> fix_modify AtC transfer thermal control hoover </TT> \n
      <TT> fix_modify AtC transfer thermal control flux bndy_faces </TT> \n
      \section description
      Sets the energy exchange mechansim from the finite elements to the atoms, managed through a control algorithm.  Rescale computes a scale factor for each atom to match the finite element temperature.  Hoover is a Gaussian least-constraint isokinetic thermostat enforces that the nodal restricted atomic temperature matches the finite element temperature.  Flux is a similar mode, but rather adds energy to the atoms based on conservation of energy.  Hoover and flux allows the prescription of sources or fixed temperatures on the atoms. 
      \section restrictions
      only for be used with specific transfers :
      thermal (rescale, hoover, flux), two_temperature (flux) \n
      rescale not valid with time filtering activated
      \section related
      \section default
      none\n
      rescale frequency is 1\n
      flux boundary_integration_type is interpolate
    */
    int argIndex = 0;
    if (strcmp(arg[argIndex],"none")==0) { // restore defaults
      thermostatType_ = NONE;
      howOften_ = 1;
      boundaryIntegrationType_ = ATC_Transfer::NO_QUADRATURE;
      foundMatch = true;
    }
    else if (strcmp(arg[argIndex],"rescale")==0) {
      argIndex++;
      howOften_ = atoi(arg[argIndex]);
      if (howOften_ < 1) {
	throw ATC_Error(0,"Bad rescaling thermostat frequency");
      }
      else {
	thermostatType_ = RESCALE;
	boundaryIntegrationType_ = ATC_Transfer::NO_QUADRATURE;
	foundMatch = true;
      }
    }
    else if (strcmp(arg[argIndex],"hoover")==0) {
      thermostatType_ = HOOVER;
      howOften_ = 1;
      boundaryIntegrationType_ = ATC_Transfer::NO_QUADRATURE;
      foundMatch = true;
    }
    else if (strcmp(arg[argIndex],"flux")==0) {
      thermostatType_ = FLUX;
      howOften_ = 1;
      argIndex++;
      boundaryIntegrationType_ = atcTransfer_->parse_boundary_integration(narg-argIndex,&arg[argIndex],boundaryFaceSet_);
      foundMatch = true;
    }

    if (!foundMatch)
      foundMatch = AtomicRegulator::modify(narg,arg);
    if (foundMatch)
      needReset_ = true;
    return foundMatch;
  }

  //--------------------------------------------------------
  //  reset_nodal_lambda_power:
  //    resets the thermostat generated power to a
  //    prescribed value
  //--------------------------------------------------------
  void Thermostat::reset_lambda_power(DENS_MAT & target)
  {
    for (int i = 0; i < nNodes_; ++i)
      lambdaPowerFiltered_(i) = target(i,0);
  }

  //--------------------------------------------------------
  //  initialize:
  //    sets up methods before a run
  //    NOTE we currently only include time filter
  //    dependence, but in general there is also a
  //    time integrator dependence.  In general the 
  //    precedence order is:
  //    time filter -> time integrator -> kinetostat
  //    In the future this may need to be added if
  //    different types of time integrators can be
  //    specified.
  //--------------------------------------------------------
  void Thermostat::initialize()
  {   
    //TimeIntegrationType myIntegrationType = (atcTransfer_->get_time_integrator())->get_time_integration_type();
    // HACK until thermal time integrator is implemented
    TimeIntegrator::TimeIntegrationType myIntegrationType = TimeIntegrator::VERLET;
    TimeFilterManager * timeFilterManager = atcTransfer_->get_time_filter_manager();

    // reset data if needed and perform any error/conflict checking
    if (resetData_) {
      AtomicRegulator::reset_data();

      // set up storage
      lambda_.reset(nNodes_,1);
      nodalAtomicLambdaPower_.reset(nNodes_);
      lambdaPowerFiltered_.reset(nNodes_);
    }

    if (needReset_ || timeFilterManager->need_reset() || timeFilterManager->end_equilibrate()) {
      // eliminate existing methods
      destroy();
      if (timeFilterManager->need_reset()) {
        if (timeFilter_)
          delete timeFilter_;
        if (myIntegrationType == TimeIntegrator::VERLET)
         timeFilter_ = timeFilterManager->construct(TimeFilterManager::EXPLICIT);
        else if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP)
         timeFilter_ = timeFilterManager->construct(TimeFilterManager::EXPLICIT_IMPLICIT);
      }
      
      if (timeFilterManager->filter_dynamics()) {
	switch (thermostatType_) {
	case NONE: {
	  break;
	}
	case RESCALE: { // error check, rescale and filtering not supported together
	  throw ATC_Error(0,"Cannot use rescaling thermostat with time filtering");
	  break;
	}
	case HOOVER: {
          if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP)
            //regulatorMethod_ = new ThermostatHooverFiltered(this);
            cout << "temporary line\n";
          else 
            regulatorMethod_ = new ThermostatHooverVerletFiltered(this);
	  break;
        }
        case FLUX: {
          if (atcTransfer_->use_localized_lambda())
            throw ATC_Error(0,"Cannot use flux thermostating with localized lambda");
          if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP)
            regulatorMethod_ = new ThermostatPowerFiltered(this);
          else
            regulatorMethod_ = new ThermostatPowerVerletFiltered(this); 
	  break;
	}
	default:
	  throw ATC_Error(0,"Unknown thermostat type in Thermostat::initialize");
	}
      }
      else {
	switch (thermostatType_) {
	case NONE: {
	  break;
	}
	case RESCALE: {
          regulatorMethod_ = new ThermostatRescale(this);
	  break;
	}
	case HOOVER: {
          if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP)
            regulatorMethod_ = new ThermostatHoover(this);
          else
            regulatorMethod_ = new ThermostatHooverVerlet(this); 
	  break;
	}
        case FLUX: {
          if (atcTransfer_->use_localized_lambda())
            throw ATC_Error(0,"Cannot use flux thermostating with localized lambda");
          if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP)
            regulatorMethod_ = new ThermostatPower(this);
          else
            regulatorMethod_ = new ThermostatPowerVerlet(this); 
	  break;
	}
	default:
	  throw ATC_Error(0,"Unknown thermostat type in Thermostat::initialize");
	}
      }
      
      AtomicRegulator::reset_method();
    }

    AtomicRegulator::initialize();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatShapeFunction
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ThermostatShapeFunction::ThermostatShapeFunction(Thermostat * thermostat) :
    RegulatorShapeFunction(thermostat),
    thermostat_(thermostat),
    v_(atcTransfer_->get_v()),
    mdMassMatrix_(atcTransfer_->get_mass_mat_md(TEMPERATURE))
  {
    fieldMask_(TEMPERATURE,FLUX) = true;
  }

  //---------------------------------------------------------
  //  set_weights:
  //    set the diagonal weighting matrix to be the atomic
  //    temperatures
  //---------------------------------------------------------
  void ThermostatShapeFunction::set_weights(DIAG_MAT & weights)
  {
    if (nLocalLambda_>0) {
      DENS_VEC weightVector(nLocal_);
      for (int i = 0; i < nLocal_; i++)
        for (int j = 0; j < nsd_; j++)
          weightVector(i) += v_[internalToAtom_(i)][j]*v_[internalToAtom_(i)][j];
      DENS_VEC maskedWeightVector(nLocalLambda_);
      maskedWeightVector = internalToOverlapMap_*weightVector;
      weights.reset(maskedWeightVector);
    }
  }

  //--------------------------------------------------------
  //  reset_nlocal:
  //    resets data dependent on local atom count
  //--------------------------------------------------------
  void ThermostatShapeFunction::reset_nlocal()
  {
    RegulatorShapeFunction::reset_nlocal();
    if (nLocal_ > 0) {
      atcTransfer_->compute_atomic_mass(atomicMass_);
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatRescale
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ThermostatRescale::ThermostatRescale(Thermostat * thermostat) :
    ThermostatShapeFunction(thermostat),
    nodalTemperature_(atcTransfer_->get_field(TEMPERATURE))
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the thermostat in the post corrector phase
  //--------------------------------------------------------
  void ThermostatRescale::apply_post_corrector(double dt)
  {
    // compute right-hand side
    DENS_MAT rhs(nNodes_,1);
    rhs = mdMassMatrix_*nodalTemperature_;
    
    // solve equations
    solve_for_lambda(rhs);

    // application of rescaling lambda due
    apply_to_atoms(v_,dt);
  }

  //--------------------------------------------------------
  //  apply_lambda_to_atoms:
  //    applies the velocity rescale with an existing lambda
  //    note oldAtomicQuantity and dt are not used
  //--------------------------------------------------------
  void ThermostatRescale::apply_to_atoms(double ** atomicQuantity,
                                         const double dt)
  {
    if (nLocal_>0) {
      DENS_MAT lambdaAtom(nLocal_,1);
      atcTransfer_->prolong(lambda_,lambdaAtom);

      double xi;
      for (int i = 0; i < nLocal_; i++) {
        xi = sqrt(lambdaAtom(i,0)/atomicMass_(i));
        for (int j = 0; j < nsd_; j++) {
          atomicQuantity[internalToAtom_(i)][j] *= xi;
        }
      }
    }
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void ThermostatRescale::output(double dt, OUTPUT_LIST & outputData)
  {
    outputData["lambda"] = &(lambda_);
  }
 
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatGlc
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ThermostatGlc::ThermostatGlc(Thermostat * thermostat) :
    ThermostatShapeFunction(thermostat),
    timeFilter_(atomicRegulator_->get_time_filter()),
    nodalAtomicLambdaPower_(thermostat_->get_nodal_atomic_lambda_power()),
    lambdaPowerFiltered_(thermostat_->get_filtered_lambda_power()),
    lambdaForce_(atomicRegulator_->get_lambda_force()),
    isFixedNode_(atcTransfer_->get_fixed_node_flags(TEMPERATURE))
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  apply_to_atoms:
  //            determines what if any contributions to the
  //            atomic moition is needed for
  //            consistency with the thermostat
  //--------------------------------------------------------
  void ThermostatGlc::apply_to_atoms(double ** atomicQuantity,
                                     const DENS_MAT & lambdaForce,
                                     double dt)
  {
    if (nLocal_>0) {
      // explicit update
      for (int i = 0; i < nLocal_; i++) {
        for (int j = 0; j < nsd_; j++) {
          atomicQuantity[internalToAtom_(i)][j] += dt*lambdaForce(i,j)/atomicMass_(i);
        }
      }
    }
  }

  //---------------------------------------------------------
  //  compute_lambda_force:
  //    computes the atomic force applied by lambda
  //---------------------------------------------------------
  void ThermostatGlc::compute_lambda_force(double * const * atomicQuantity,
                                           DENS_MAT & lambdaForce)
  {
    if (nLocal_>0) {
      // scaled prolongation to (unique) nodes
      DENS_MAT lambdaAtom(nLocal_,1);
      atcTransfer_->prolong(lambda_,lambdaAtom);
      
      // explicit update
      for (int i = 0; i < nLocal_; i++) {
        for (int j = 0; j < nsd_; j++) { 
          lambdaForce(i,j) = -lambdaAtom(i,0)*atomicQuantity[internalToAtom_(i)][j];
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatPower
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatPower::ThermostatPower(Thermostat * thermostat) :
    ThermostatGlc(thermostat),
    nodalTemperature_(atcTransfer_->get_field(TEMPERATURE)),
    nodalAtomicTemperature_(atcTransfer_->get_atomic_field(TEMPERATURE)),
    heatSource_(atcTransfer_->get_atomic_source(TEMPERATURE)),
    nodalWeights_(atcTransfer_->get_shape_function_weights()),
    myAtomicVelocity_(NULL),
    f_(atcTransfer_->get_f()),
    nLocalTotal_(0),
    lambdaMaxIterations_(3),
    lambdaTolerance_(1.e-14),
    filterCoefficient_(1.),
    Nstar_(atcTransfer_->get_shape_function_ghost_overlap()),
    ghostToAtom_(atcTransfer_->get_ghost_to_atom_map())
  {
    deltaTemperature_.reset(nNodes_,1);
    deltaNodalAtomicTemperature_.reset(nNodes_,1);
    nodalTemperatureRoc_.reset(nNodes_,1);

    // sets up time filter for cases where variables temporally filtered
    TimeFilterManager * timeFilterManager = atcTransfer_->get_time_filter_manager();
    if (!timeFilterManager->end_equilibrate()) {
      // NOTE do we need a good initial guess for lambdaForce_?
      nodalAtomicLambdaPower_ = 0.;
      lambdaPowerFiltered_ = 0.;
      timeFilter_->initialize(lambdaPowerFiltered_);
    }
  }

  //--------------------------------------------------------
  //  Destructor
  //         delete thermostat data
  //--------------------------------------------------------
  ThermostatPower::~ThermostatPower()
  {
    destroy();
  }

  //--------------------------------------------------------
  //  destroy
  //         deletes allocated memory
  //--------------------------------------------------------
  void ThermostatPower::destroy()
  {
    if (myAtomicVelocity_) {
      for (int i = 0; i < nLocalTotal_; i++)
        delete [] myAtomicVelocity_[i];
      delete [] myAtomicVelocity_;
    }
  }

  //--------------------------------------------------------
  //  reset_nlocal:
  //    resizes myAtomicVelocity force if necessary
  //--------------------------------------------------------
  void ThermostatPower::reset_nlocal()
  {
    destroy();
    ThermostatGlc::reset_nlocal();
    nLocalTotal_ = atcTransfer_->get_nlocal_total();
    if (nLocal_ > 0) {
      myAtomicVelocity_ = new double*[nLocalTotal_];
      for (int i = 0; i < nLocalTotal_; i++)
        myAtomicVelocity_[i] = new double[nsd_];
    }
  }

  //--------------------------------------------------------
  //  apply_predictor:
  //    apply the thermostat to the atoms in the first step
  //    of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatPower::apply_pre_predictor(double dt)
  {
    // initialize delta temperatures
    deltaTemperature_ = -1.*nodalTemperature_;
    deltaNodalAtomicTemperature_ = -1.*nodalAtomicTemperature_;
    
    // apply lambda force to atoms and compute instantaneous lambda power
    apply_to_atoms(v_,lambdaForce_,dt);
    
    // update nodal variables
    update_nodal_quantities_predictor(dt);
  }

  //--------------------------------------------------------
  //  update_nodal_quantities_predictor:
  //    updates all needed nodal quantities after the
  //    predictor step
  //--------------------------------------------------------
  void ThermostatPower::update_nodal_quantities_predictor(double dt)
  {
    atcTransfer_->apply_inverse_mass_matrix(nodalAtomicLambdaPower_,
                                            nodalTemperatureRoc_,
                                            TEMPERATURE);
    DENS_MAT boundaryFluxRoc(nNodes_,1);
    atcTransfer_->apply_inverse_mass_matrix(boundaryFlux_[TEMPERATURE],
                                            boundaryFluxRoc,
                                            TEMPERATURE);
    nodalTemperatureRoc_ += boundaryFluxRoc;

     // update FE temperature
    nodalTemperature_ += dt*nodalTemperatureRoc_;

    // finalize filtered lambda power
    timeFilter_->apply_pre_step1(lambdaPowerFiltered_,nodalAtomicLambdaPower_,dt);
  }

  //--------------------------------------------------------
  //  apply_pre_corrector:
  //    apply the thermostat to the atoms in the first part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatPower::apply_pre_corrector(double dt)
  {
    // copy velocities over into temporary storage
    if (nLocal_ > 0) {
      for (int i = 0; i < nLocalTotal_; i++)
        for (int j = 0; j < nsd_; j++)
          myAtomicVelocity_[i][j] = v_[i][j];
    }
    
    // compute predicted changes in FE and nodal atomic temperatures
    compute_delta_nodal_atomic_temperature(dt);
    DENS_MAT myDeltaTemperature(deltaTemperature_);
    deltaTemperature_ += nodalTemperature_;

    // set up rhs for lambda equation
    DENS_MAT rhs(nNodes_,1);
    set_thermostat_rhs(rhs,dt);
    
    // solve linear system for lambda guess
    solve_for_lambda(rhs);

    // iterate to solution
    iterate_lambda(rhs,v_,dt);

    // compute force applied by lambda
    compute_lambda_force(v_,lambdaForce_);

    // apply lambda force to atoms and compute instantaneous lambda power
    apply_to_atoms(myAtomicVelocity_,lambdaForce_,dt);

    // update nodal variables
    update_nodal_quantities_pre_corrector(dt);
    
    // reset temporary variables
    deltaTemperature_ = myDeltaTemperature;
  }

  //--------------------------------------------------------
  //  compute_delta_nodal_atomic_temperature:
  //    computes the change in the nodal atomic temperature
  //    that has occured over the past timestep
  //--------------------------------------------------------
  void ThermostatPower::compute_delta_nodal_atomic_temperature(double dt)
  {
    // set delta temperatures based on predicted atomic velocities
    DENS_MAT atomicTemperature(nLocal_,1);
    atcTransfer_->compute_atomic_temperature(atomicTemperature, v_);
    DENS_MAT nodalAtomicTemperature(nNodes_,1);
    // NOTE change this when we use physical mass matrices
    atcTransfer_->restrict(atomicTemperature,nodalAtomicTemperature);
    deltaNodalAtomicTemperature_ += nodalAtomicTemperature;
  }

  //--------------------------------------------------------
  //  update_nodal_quantities_pre_corrector:
  //    updates all needed nodal quantities after the
  //    pre-corrector step
  //--------------------------------------------------------
  void ThermostatPower::update_nodal_quantities_pre_corrector(double dt)
  {
    atcTransfer_->apply_inverse_mass_matrix(nodalAtomicLambdaPower_,
                                            nodalTemperatureRoc_,
                                            TEMPERATURE);
    DENS_MAT boundaryFluxRoc(nNodes_,1);
    atcTransfer_->apply_inverse_mass_matrix(boundaryFlux_[TEMPERATURE],
                                            boundaryFluxRoc,
                                            TEMPERATURE);
    nodalTemperatureRoc_ += boundaryFluxRoc;

    // update FE temperature
    nodalTemperature_ += dt*nodalTemperatureRoc_;
  }

  //--------------------------------------------------------
  //  undo_update_nodal_quantities_pre_corrector:
  //    undoes the nodal quantities update from after the
  //    pre-corrector step
  //--------------------------------------------------------
  void ThermostatPower::undo_update_nodal_quantities_pre_corrector(double dt)
  {
    // remove predicted power effects
    nodalTemperature_ -= dt*nodalTemperatureRoc_;
  }

  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the thermostat to the atoms in the second part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatPower::apply_post_corrector(double dt)
  {
    // remove predicted power effects
    undo_update_nodal_quantities_pre_corrector(dt);

    // set delta FE temperature
    deltaTemperature_ += nodalTemperature_;

    // set up rhs for lambda equation
    DENS_MAT rhs(nNodes_,1);
    set_thermostat_rhs(rhs,dt);
    
    // solve linear system for lambda guess
    solve_for_lambda(rhs);

    // iterate to solution
    iterate_lambda(rhs,myAtomicVelocity_,dt);

    // compute force applied by lambda
    compute_lambda_force(myAtomicVelocity_,lambdaForce_);

    // apply lambda force to atoms and compute instantaneous lambda power
    apply_to_atoms(v_,lambdaForce_,dt);

    // update nodal variables
    update_nodal_quantities_post_corrector(dt);
  }

  //--------------------------------------------------------
  //  update_nodal_quantities_post_corrector:
  //    updates all needed nodal quantities after the
  //    post-corrector step
  //--------------------------------------------------------
  void ThermostatPower::update_nodal_quantities_post_corrector(double dt)
  {
    // finalize filtered lambda power
    timeFilter_->apply_post_step1(lambdaPowerFiltered_,nodalAtomicLambdaPower_,dt);

    atcTransfer_->apply_inverse_mass_matrix(nodalAtomicLambdaPower_,
                                            nodalTemperatureRoc_,
                                            TEMPERATURE);
    DENS_MAT boundaryFluxRoc(nNodes_,1);
    atcTransfer_->apply_inverse_mass_matrix(boundaryFlux_[TEMPERATURE],
                                            boundaryFluxRoc,
                                            TEMPERATURE);
    nodalTemperatureRoc_ += boundaryFluxRoc;

     // update FE temperature
    nodalTemperature_ += dt*nodalTemperatureRoc_;
  }

  //--------------------------------------------------------
  //  apply_to_atoms:
  //    add the thermostat force to the atoms and compute
  //    the instantaneous induced power
  //--------------------------------------------------------
  void ThermostatPower::apply_to_atoms(double ** atomicVelocity,
                                       const DENS_MAT & lambdaForce,
                                       double dt)
  {
    // compute initial contributions to lambda power
    DENS_VEC atomicTemperature(nLocal_);
    DENS_VEC nodalAtomicTemperature(nNodes_);
    atcTransfer_->compute_atomic_temperature(atomicTemperature, atomicVelocity);
    atcTransfer_->restrict_unscaled(atomicTemperature,nodalAtomicTemperature);
    nodalAtomicLambdaPower_ = -1.*nodalAtomicTemperature;

    ThermostatGlc::apply_to_atoms(atomicVelocity,lambdaForce,dt);

    // finalize lambda power
    atcTransfer_->compute_atomic_temperature(atomicTemperature, atomicVelocity);
    atcTransfer_->restrict_unscaled(atomicTemperature,nodalAtomicTemperature);
    nodalAtomicLambdaPower_ += nodalAtomicTemperature;
    nodalAtomicLambdaPower_ = (1./dt)*nodalAtomicLambdaPower_;
  }

  //---------------------------------------------------------
  //  set_weights:
  //    set the diagonal weighting matrix to be the atomic
  //    temperatures
  //---------------------------------------------------------
  void ThermostatPower::set_weights(DIAG_MAT & weights)
  {
    if (nLocalLambda_>0) {
      DENS_VEC weightVector(nLocal_);
      DENS_VEC maskedWeightVector(nLocalLambda_);
      atcTransfer_->compute_atomic_temperature(weightVector,v_,myAtomicVelocity_);
      maskedWeightVector = internalToOverlapMap_*weightVector;
      weights.reset(maskedWeightVector);
    }
  }

  //--------------------------------------------------------
  //  set_thermostat_rhs:
  //    sets up the right-hand side including boundary
  //    fluxes (coupling & prescribed), heat sources, and
  //    fixed (uncoupled) nodes
  //--------------------------------------------------------
  void ThermostatPower::set_thermostat_rhs(DENS_MAT & rhs_nodes,
                                           double dt)
  {
    // NOTE is the scaling right? (a) E/t/kB & (b) T/t/kB
    // only tested with flux != 0 + ess bc = 0

    // (a) for flux based : 
    // form rhs :  2/3kB * W_I^-1 * \int N_I r dV
    // vs  Wagner, CMAME, 2008 eq(24) RHS_I = 2/(3kB) flux_I
    // fluxes are set in ATC transfer
    // NOTE change coef and nodalWeights when physical mass matrices are used
    // and remove rhoCp
    double rhoCp = LammpsInterface::instance()->heat_capacity();
    double k_boltzmann = LammpsInterface::instance()->kBoltzmann();
    double coef = filterCoefficient_*rhoCp*2./(nsd_*k_boltzmann);
    rhs_nodes = coef*heatSource_*nodalWeights_;

    // (b) for essential bcs (fixed nodes) :
    // form rhs : (delThetaV - delTheta)/dt
    DENS_MAT rhs_nodes_prescribed(nNodes_,1);
    rhs_nodes_prescribed = (1./dt)*filterCoefficient_*(deltaNodalAtomicTemperature_-deltaTemperature_);

    // replace rhs for prescribed nodes
    // NOTE can we make is_fixed member data?
    for (int i = 0; i < nNodes_; i++) {
      if (isFixedNode_(i,0)) {
        rhs_nodes(i,0) = rhs_nodes_prescribed(i,0);
      }
    }
  }

  //--------------------------------------------------------
  //  iterate_lambda:
  //    iteratively solves the equations for lambda
  //    for the higher order dt corrections, assuming
  //    an initial guess for lambda
  //--------------------------------------------------------
  void ThermostatPower::iterate_lambda(const MATRIX & rhs,
                                       double * const * atomicVelocity,
                                       double dt)
  {
    DENS_VEC lambdaOverlap(nNodeOverlap_);
    atcTransfer_->map_unique_to_overlap(lambda_, lambdaOverlap);

    DENS_VEC atomicTemperature;
    DENS_MAT lambdaAtom;
    DENS_MAT lambdaSquaredT;
    if (nLocalLambda_>0) {
      lambdaAtom.reset(nLocal_,1);
      lambdaSquaredT.reset(nLocal_,1);
      atomicTemperature.reset(nLocal_);
      atcTransfer_->compute_atomic_temperature(atomicTemperature, atomicVelocity);
    }
  
    DENS_MAT lambdaOld(nNodes_,1);
    DENS_MAT rhsOverlap(nNodeOverlap_,1);
    DENS_MAT rhsTotalLocal(nNodes_,1);
    DENS_MAT rhsTotal(nNodes_,1);
    // solve assuming we get initial guess for lambda
    for (int i = 0; i < lambdaMaxIterations_; ++i) {
      lambdaOld = lambda_;
      
      // assemble quadratic term
      if (nLocalLambda_ > 0) {
        atcTransfer_->prolong_scaled(lambda_, lambdaAtom);
        for (int i = 0; i < nLocal_; i++) {
          lambdaSquaredT(i,0) = lambdaAtom(i,0)*lambdaAtom(i,0)*atomicTemperature(i);
        }
        rhsOverlap = (dt/4.)*shapeFunctionMatrix_.transMat(internalToOverlapMap_*lambdaSquaredT);
        atcTransfer_->map_overlap_to_unique(rhsOverlap,rhsTotalLocal);
      }
      LammpsInterface::instance()->allsum(rhsTotalLocal.get_ptr(),rhsTotal.get_ptr(),rhsTotal.size());

      // the system with the new rhs
      rhsTotal += rhs;
      solve_for_lambda(rhsTotal);

      // check convergence
      DENS_MAT difference = lambda_-lambdaOld;
      if (difference.col_norm()/lambdaOld.col_norm() < lambdaTolerance_) break;
      
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatHoover
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatHoover::ThermostatHoover(Thermostat * thermostat) :
    ThermostatPower(thermostat),
    nodeToOverlapMap_(atcTransfer_->get_node_to_overlap_map())
  {
    lambdaForceHoover_.reset(nLocal_,nsd_);
  }

  //--------------------------------------------------------
  //  reset_nlocal:
  //    resizes Hoover lambda force if necessary
  //--------------------------------------------------------
  void ThermostatHoover::reset_nlocal()
  {
    ThermostatPower::reset_nlocal();
    if (nLocal_ > 0)
      lambdaForceHoover_.reset(nLocal_,nsd_);
  }

  //--------------------------------------------------------
  //  apply_predictor:
  //    apply the thermostat to the atoms in the first step
  //    of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatHoover::apply_pre_predictor(double dt)
  {
    ThermostatPower::apply_pre_predictor(dt);

    // apply lambda force to atoms and compute instantaneous lambda power
    DENS_VEC myNodalAtomicLambdaPower(nodalAtomicLambdaPower_);
    apply_to_atoms(v_,lambdaForceHoover_,dt);

    // compute nodal atomic power from Hoover coupling
    // only add in contribution to uncoupled nodes
    if (atcTransfer_->use_localized_lambda())
      add_to_lambda_power_predictor(dt);
    else
      nodalAtomicLambdaPower_ = 0.;

    // finish updating lambda power
    nodalAtomicLambdaPower_ += myNodalAtomicLambdaPower;
  }

  //--------------------------------------------------------
  //  add_to_lambda_power_predictor:
  //    adds Hoover contributions to the lambda power term
  //--------------------------------------------------------
  void ThermostatHoover::add_to_lambda_power_predictor(double dt)
  {
    atcTransfer_->apply_inverse_mass_matrix(nodalAtomicLambdaPower_,
                                            nodalTemperatureRoc_,
                                            TEMPERATURE);

    // update FE temperature
    for (int i = 0; i < nNodes_; ++i) {
      if (nodeToOverlapMap_(i)==-1)
        nodalTemperature_(i,0) += dt*nodalTemperatureRoc_(i,0);
      else
        nodalAtomicLambdaPower_(i) = 0.;
    }

    // finalize filtered lambda power
    timeFilter_->apply_pre_step2(lambdaPowerFiltered_,nodalAtomicLambdaPower_,dt);
  }

  //--------------------------------------------------------
  //  apply_pre_corrector:
  //    apply the thermostat to the atoms in the first part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatHoover::apply_pre_corrector(double dt)
  {
    ThermostatPower::apply_pre_corrector(dt);
    
    // set up Hoover rhs
    DENS_MAT myDeltaTemperature(deltaTemperature_);
    deltaTemperature_ += nodalTemperature_;
    DENS_MAT rhs(nNodes_,1);
    set_hoover_rhs(rhs,dt);

    // solve linear system for lambda
    solve_for_lambda(rhs);

    // iterate to solution
    iterate_lambda(rhs,v_,dt);

    // compute force applied by lambda
    compute_lambda_force(v_,lambdaForceHoover_);

    // apply the force to the atoms
    DENS_VEC myNodalAtomicLambdaPower(nodalAtomicLambdaPower_);
    apply_to_atoms(myAtomicVelocity_,lambdaForceHoover_,dt);

    // compute nodal atomic power from Hoover coupling
    // only add in contribution to uncoupled nodes
    if (atcTransfer_->use_localized_lambda())
      add_to_lambda_power_pre_corrector(dt);
    else
      nodalAtomicLambdaPower_ = 0.;

    // finish updating lambda power
    nodalAtomicLambdaPower_ += myNodalAtomicLambdaPower;

    // reset temporary variables
    deltaTemperature_ = myDeltaTemperature;
  }

  //--------------------------------------------------------
  //  add_to_lambda_power_pre_corrector:
  //    adds Hoover contributions to the lambda power term
  //--------------------------------------------------------
  void ThermostatHoover::add_to_lambda_power_pre_corrector(double dt)
  {
    DENS_MAT myNodalTemperatureRoc(nodalTemperatureRoc_);
    atcTransfer_->apply_inverse_mass_matrix(nodalAtomicLambdaPower_,
                                            nodalTemperatureRoc_,
                                            TEMPERATURE);

    // update FE temperature
    for (int i = 0; i < nNodes_; ++i) {
      if (nodeToOverlapMap_(i)==-1)
        nodalTemperature_(i,0) += dt*nodalTemperatureRoc_(i,0);
      else {
        nodalTemperatureRoc_(i,0) = 0.;
        nodalAtomicLambdaPower_(i) = 0.;
      }
    }

    // correct nodal temperature rate of change
    nodalTemperatureRoc_ += myNodalTemperatureRoc;
  }

  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the thermostat to the atoms in the second part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatHoover::apply_post_corrector(double dt)
  {
    ThermostatPower::apply_post_corrector(dt);
    
    // set up Hoover rhs
    DENS_MAT rhs(nNodes_,1);
    set_hoover_rhs(rhs,dt);
    
    // solve linear system for lambda
    solve_for_lambda(rhs);

    // iterate to solution
    iterate_lambda(rhs,myAtomicVelocity_,dt);

    // compute force applied by lambda
    compute_lambda_force(myAtomicVelocity_,lambdaForceHoover_);

    // apply the force to the atoms
    DENS_VEC myNodalAtomicLambdaPower(nodalAtomicLambdaPower_);
    apply_to_atoms(v_,lambdaForceHoover_,dt);
    
    // compute nodal atomic power from Hoover coupling
    // only add in contribution to uncoupled nodes
    if (atcTransfer_->use_localized_lambda())
      add_to_lambda_power_post_corrector(dt);
    else
      nodalAtomicLambdaPower_ = 0.;

    // finish updating lambda power
    nodalAtomicLambdaPower_ += myNodalAtomicLambdaPower;
  }

  //--------------------------------------------------------
  //  add_to_lambda_power_post_corrector:
  //    adds Hoover contributions to the lambda power term
  //--------------------------------------------------------
  void ThermostatHoover::add_to_lambda_power_post_corrector(double dt)
  {
    // finalize filtered lambda power
    timeFilter_->apply_post_step2(lambdaPowerFiltered_,nodalAtomicLambdaPower_,dt);

    atcTransfer_->apply_inverse_mass_matrix(nodalAtomicLambdaPower_,
                                            nodalTemperatureRoc_,
                                            TEMPERATURE);

    // update FE temperature
    for (int i = 0; i < nNodes_; ++i) {
      if (nodeToOverlapMap_(i)==-1)
        nodalTemperature_(i,0) += dt*nodalTemperatureRoc_(i,0);
      else
        nodalAtomicLambdaPower_(i) = 0.;
    }
  }

  //--------------------------------------------------------
  //  set_hoover_rhs:
  //    sets up the right-hand side for Hoover coupling with
  //    boundary nodes
  //--------------------------------------------------------
  void ThermostatHoover::set_hoover_rhs(DENS_MAT & rhs_nodes,
                                        double dt)
  {
    // form rhs : (delThetaV - delTheta)/dt
    rhs_nodes = (1./dt)*filterCoefficient_*(deltaNodalAtomicTemperature_-deltaTemperature_);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatPowerFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatPowerFiltered::ThermostatPowerFiltered(Thermostat * thermostat) :
    ThermostatPower(thermostat)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  update_nodal_quantities_predictor:
  //    updates all needed nodal quantities after the
  //    predictor step
  //--------------------------------------------------------
  void ThermostatPowerFiltered::update_nodal_quantities_predictor(double dt)
  {
    // update nodal atomic temperature
    //atcTransfer_->apply_inverse_mass_matrix(lambdaPowerFiltered_,
    //                                            nodalTemperatureRoc_,
    //                                            TEMPERATURE);
    // NOTE temporary until new mass matrices are used
    DIAG_MAT & atomicWeights(atcTransfer_->get_atomic_weights());
    nodalTemperatureRoc_ = (1./atomicWeights(0,0))*nodalWeights_*lambdaPowerFiltered_;
    nodalAtomicTemperature_ += dt*nodalTemperatureRoc_;

    // update FE temperature
    atcTransfer_->apply_inverse_mass_matrix(lambdaPowerFiltered_,
                                            nodalTemperatureRoc_,
                                            TEMPERATURE);
    //DENS_MAT boundaryFluxRoc(nNodes_,1);
    //atcTransfer_->apply_inverse_mass_matrix(boundaryFlux_[TEMPERATURE],
    //                                        boundaryFluxRoc,
    //                                        TEMPERATURE);
    //nodalTemperatureRoc_ += boundaryFluxRoc;
    nodalTemperature_ += 0.5*dt*nodalTemperatureRoc_;

    // finalize filtered lambda power
    timeFilter_->apply_pre_step1(lambdaPowerFiltered_,nodalAtomicLambdaPower_,dt);
  }

  //--------------------------------------------------------
  //  compute_delta_nodal_atomic_temperature:
  //    computes the change in the nodal atomic temperature
  //    that has occured over the past timestep
  //--------------------------------------------------------
  void ThermostatPowerFiltered::compute_delta_nodal_atomic_temperature(double dt)
  {
    // compute next lambda power assuming zero new force
    //DENS_VEC myLambdaPowerFiltered(lambdaPowerFiltered_);
    DENS_VEC zeroNodalLambdaPower(nNodes_);
    timeFilter_->apply_post_step1(lambdaPowerFiltered_,zeroNodalLambdaPower,dt);
            
    // update nodal atomic temperature
    //atcTransfer_->apply_inverse_mass_matrix(lambdaPowerFiltered_,
    //                                            nodalTemperatureRoc_,
    //                                            TEMPERATURE);
    // NOTE temporary until new mass matrices are used
    DIAG_MAT & atomicWeights(atcTransfer_->get_atomic_weights());
    nodalTemperatureRoc_ = (1./atomicWeights(0,0))*nodalWeights_*lambdaPowerFiltered_;
    nodalAtomicTemperature_ += dt*nodalTemperatureRoc_;

    // update FE temperature
    atcTransfer_->apply_inverse_mass_matrix(lambdaPowerFiltered_,
                                            nodalTemperatureRoc_,
                                            TEMPERATURE);
    nodalTemperature_ += .5*dt*nodalTemperatureRoc_;

    deltaNodalAtomicTemperature_ += nodalAtomicTemperature_;
  }

  //--------------------------------------------------------
  //  update_nodal_quantities_pre_corrector:
  //    updates all needed nodal quantities after the
  //    pre-corrector step
  //--------------------------------------------------------
  void ThermostatPowerFiltered::update_nodal_quantities_pre_corrector(double dt)
  {
    atcTransfer_->apply_inverse_mass_matrix(nodalAtomicLambdaPower_,
                                            nodalTemperatureRoc_,
                                            TEMPERATURE);

    double filterCoefficient = timeFilter_->get_unfiltered_coefficient_post_s1(dt);

    //DENS_MAT boundaryFluxRoc(nNodes_,1);
    //atcTransfer_->apply_inverse_mass_matrix(boundaryFlux_[TEMPERATURE],
    //                                        boundaryFluxRoc,
    //                                        TEMPERATURE);
    //nodalTemperatureRoc_ = boundaryFluxRoc + filterCoefficient*nodalTemperatureRoc_;
    nodalTemperatureRoc_ = 0.5*filterCoefficient*nodalTemperatureRoc_;
    // update FE temperature
    nodalTemperature_ += dt*nodalTemperatureRoc_;
  }

//   //--------------------------------------------------------
//   //  undo_update_nodal_quantities_pre_corrector:
//   //    undoes the update of all nodal quantities from after
//   //    the pre-corrector step
//   //--------------------------------------------------------
//   void ThermostatPowerFiltered::undo_update_nodal_quantities_pre_corrector(double dt)
//   {
//     double filterCoefficient = timeFilter_->get_unfiltered_coefficient_pre_s2(dt);

//     // update FE temperature
//     nodalTemperature_ -= dt*filterCoefficient*nodalTemperatureRoc_;
//   }

  //--------------------------------------------------------
  //  update_nodal_quantities_post_corrector:
  //    updates all needed nodal quantities after the
  //    post-corrector step
  //--------------------------------------------------------
  void ThermostatPowerFiltered::update_nodal_quantities_post_corrector(double dt)
  {
    // finalize filtered lambda power
    timeFilter_->apply_post_step2(lambdaPowerFiltered_,nodalAtomicLambdaPower_,dt);
    double filterCoefficient = timeFilter_->get_unfiltered_coefficient_post_s1(dt);

    // update nodal atomic temperature
    //atcTransfer_->apply_inverse_mass_matrix(nodalAtomicLambdaPower_,
    //                                            nodalTemperatureRoc_,
    //                                            TEMPERATURE);
    // NOTE temporary until new mass matrices are used
    DIAG_MAT & atomicWeights(atcTransfer_->get_atomic_weights());
    nodalTemperatureRoc_ = (1./atomicWeights(0,0))*nodalWeights_*nodalAtomicLambdaPower_;
    //nodalTemperatureRoc_.print("NTR");
    nodalAtomicTemperature_ += filterCoefficient*dt*nodalTemperatureRoc_;

    // update FE temperature
    atcTransfer_->apply_inverse_mass_matrix(nodalAtomicLambdaPower_,
                                            nodalTemperatureRoc_,
                                            TEMPERATURE);
    //DENS_MAT boundaryFluxRoc(nNodes_,1);
    //atcTransfer_->apply_inverse_mass_matrix(boundaryFlux_[TEMPERATURE],
    //                                        boundaryFluxRoc,
    //                                        TEMPERATURE);
    //nodalTemperatureRoc_ = boundaryFluxRoc + filterCoefficient*nodalTemperatureRoc_;
    nodalTemperatureRoc_ = .5*filterCoefficient*nodalTemperatureRoc_;
    nodalTemperature_ += dt*nodalTemperatureRoc_;
  }

  //--------------------------------------------------------
  //  set_thermostat_rhs:
  //    sets up the right-hand side including boundary
  //    fluxes (coupling & prescribed), heat sources, and
  //    fixed (uncoupled) nodes
  //--------------------------------------------------------
  void ThermostatPowerFiltered::set_thermostat_rhs(DENS_MAT & rhs_nodes,
                                                   double dt)
  {
    filterCoefficient_ = 1./(timeFilter_->get_unfiltered_coefficient_post_s1(dt));
    // NOTE change the 0.5 when we change how the lambda power is accounted
    // NOTE change rhoCp^-1 when physical mass matrices are used
    double rhoCp = LammpsInterface::instance()->heat_capacity();
    heatSource_ += (0.5)*lambdaPowerFiltered_;
    ThermostatPower::set_thermostat_rhs(rhs_nodes,dt);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatPowerVerlet
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatPowerVerlet::ThermostatPowerVerlet(Thermostat * thermostat) :
    ThermostatGlc(thermostat),
    nodalTemperatureRoc_(atcTransfer_->get_field_roc(TEMPERATURE)),
    heatSource_(atcTransfer_->get_atomic_source(TEMPERATURE)),
    f_(atcTransfer_->get_f())
  {
    // sets up time filter for cases where variables temporally filtered
    TimeFilterManager * timeFilterManager = atcTransfer_->get_time_filter_manager();
    if (!timeFilterManager->end_equilibrate()) {
      nodalAtomicLambdaPower_ = 0.;
      lambdaPowerFiltered_ = 0.;
      timeFilter_->initialize(lambdaPowerFiltered_);
    }
  }

  //--------------------------------------------------------
  //  apply_predictor:
  //    apply the thermostat to the atoms in the first step
  //    of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatPowerVerlet::apply_pre_predictor(double dt)
  {
    compute_thermostat(dt);

    // apply lambda force to atoms
    apply_to_atoms(v_,lambdaForce_,dt);
  }

  //--------------------------------------------------------
  //  apply_pre_corrector:
  //    apply the thermostat to the atoms in the first part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatPowerVerlet::apply_pre_corrector(double dt)
  {
    compute_thermostat(dt);

    // apply lambda force to atoms
    apply_to_atoms(v_,lambdaForce_,dt);
  }

  //--------------------------------------------------------
  //  add_to_rhs:
  //            determines what if any contributions to the
  //            finite element equations are needed for
  //            consistency with the thermostat
  //--------------------------------------------------------
  void ThermostatPowerVerlet::add_to_rhs(FIELDS & rhs)
  {
    rhs[TEMPERATURE] += nodalAtomicLambdaPower_ + boundaryFlux_[TEMPERATURE];
  }

  //--------------------------------------------------------
  //  set_thermostat_rhs:
  //    sets up the right-hand side including boundary
  //    fluxes (coupling & prescribed), heat sources, and
  //    fixed (uncoupled) nodes
  //--------------------------------------------------------
  void ThermostatPowerVerlet::set_thermostat_rhs(DENS_MAT & rhs_nodes)
  {
    // (a) for flux based : 
    // form rhs :  \int N_I r dV
    // vs  Wagner, CMAME, 2008 eq(24) RHS_I = 2/(3kB) flux_I
    // fluxes are set in ATC transfer
    rhs_nodes = heatSource_;
    
    // (b) for ess. bcs
    // form rhs : {sum_a (2 * N_Ia * v_ia * f_ia) - (dtheta/dt)_I}
    DENS_MAT atomicPower;
    atcTransfer_->compute_atomic_power(atomicPower, v_, f_);
    DENS_MAT rhs_nodes_prescribed(nNodes_,1);
    atcTransfer_->restrict_volumetric_quantity(atomicPower, rhs_nodes_prescribed);
    rhs_nodes_prescribed -= 0.5*(mdMassMatrix_*nodalTemperatureRoc_);
    
    // replace rhs for prescribed nodes
    for (int i = 0; i < nNodes_; i++) {
      if (isFixedNode_(i,0)) {
        rhs_nodes(i,0) = rhs_nodes_prescribed(i,0);
      }
    }
  }

  //--------------------------------------------------------
  //  compute_thermostat:
  //    sets up and solves the thermostat equations since
  //    they are the same at different parts of the time
  //    step
  //--------------------------------------------------------
  void ThermostatPowerVerlet::compute_thermostat(double dt)
  {
    // set up rhs
    DENS_MAT rhs(nNodes_,1);
    set_thermostat_rhs(rhs);
    
    // solve linear system for lambda
    solve_for_lambda(rhs);
    
    // compute force applied by lambda
    compute_lambda_force(v_,lambdaForce_);

    // compute nodal atomic power
    compute_nodal_lambda_power(dt);
  }

  //--------------------------------------------------------
  //  compute_nodal_lambda_power:
  //    determines the power exerted by the thermostat at
  //    each FE node
  //--------------------------------------------------------
  void ThermostatPowerVerlet::compute_nodal_lambda_power(double dt)
  {
    DENS_VEC atomicLambdaPower;
    atcTransfer_->compute_atomic_power(atomicLambdaPower,v_,lambdaForce_);
    atcTransfer_->restrict_volumetric_quantity(atomicLambdaPower,nodalAtomicLambdaPower_);
    nodalAtomicLambdaPower_ *= 2.; // accounts for kinetic definition of temperature
    timeFilter_->apply_pre_step1(lambdaPowerFiltered_,nodalAtomicLambdaPower_,dt);
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void ThermostatPowerVerlet::output(double dt, OUTPUT_LIST & outputData)
  {
    outputData["lambda"] = &(lambda_);
    outputData["nodalLambdaPower"] = &(nodalAtomicLambdaPower_);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatHooverVerlet
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatHooverVerlet::ThermostatHooverVerlet(Thermostat * thermostat) :
    ThermostatPowerVerlet(thermostat),
    nodeToOverlapMap_(atcTransfer_->get_node_to_overlap_map())
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  compute_thermostat:
  //    sets up and solves the thermostat equations since
  //    they are the same at different parts of the time
  //    step
  //--------------------------------------------------------
  void ThermostatHooverVerlet::compute_thermostat(double dt)
  {
    // apply prescribed/extrinsic sources and fixed nodes
    ThermostatPowerVerlet::compute_thermostat(dt);
    
    // set up Hoover rhs
    DENS_MAT rhs(nNodes_);
    set_hoover_rhs(rhs);
    
    // solve linear system for lambda
    solve_for_lambda(rhs);
    
    // compute force applied by lambda
    DENS_MAT lambdaForceHoover(nLocal_,nsd_);
    compute_lambda_force(v_,lambdaForceHoover);
    lambdaForce_ += lambdaForceHoover;

    // compute nodal atomic power from Hoover coupling
    // only add in contribution to uncoupled nodes
    if (atcTransfer_->use_localized_lambda())
      add_to_lambda_power(lambdaForceHoover,dt);
  }

  //--------------------------------------------------------
  //  set_hoover_rhs:
  //    sets up the right-hand side for fixed value,
  //    i.e. Hoover coupling
  //--------------------------------------------------------
  void ThermostatHooverVerlet::set_hoover_rhs(DENS_MAT & rhs_nodes)
  {
    // form rhs : sum_a ( N_Ia * v_ia * f_ia) - 0.5*M_MD*(dtheta/dt)_I
    DENS_MAT atomicPower(atcTransfer_->get_nlocal(),1);
    atcTransfer_->compute_atomic_power(atomicPower, v_, f_);
    atcTransfer_->restrict_volumetric_quantity(atomicPower, rhs_nodes);
    rhs_nodes -= 0.5*(mdMassMatrix_*nodalTemperatureRoc_);
  }

  //--------------------------------------------------------
  //  add_to_nodal_lambda_power:
  //    determines the power exerted by the Hoover 
  //    thermostat at each FE node
  //--------------------------------------------------------
  void ThermostatHooverVerlet::add_to_lambda_power(DENS_MAT & myLambdaForce,
                                                   double dt)
  {
    DENS_VEC atomicLambdaPower;
    atcTransfer_->compute_atomic_power(atomicLambdaPower,v_,myLambdaForce);
    DENS_VEC myNodalLambdaPower(nNodes_);
    atcTransfer_->restrict_volumetric_quantity(atomicLambdaPower,myNodalLambdaPower);
    myNodalLambdaPower *= 2.;
    for (int i = 0; i < nNodes_; ++i) {
      if (nodeToOverlapMap_(i)==-1)
        nodalAtomicLambdaPower_(i) += myNodalLambdaPower(i);
      else
        myNodalLambdaPower(i) = 0.;
    }
    timeFilter_->apply_post_step1(lambdaPowerFiltered_,myNodalLambdaPower,dt);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatPowerVerletFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatPowerVerletFiltered::ThermostatPowerVerletFiltered(Thermostat * thermostat) :
    ThermostatPowerVerlet(thermostat),
    nodalTemperature2Roc_(atcTransfer_->get_field_2roc(TEMPERATURE)),
    fieldsRoc_(atcTransfer_->get_fields_roc()),
    filterScale_((atcTransfer_->get_time_filter_manager())->get_filter_scale())
  {
    heatSourceRoc_.reset(nNodes_,1);
    fluxRoc_[TEMPERATURE].reset(nNodes_,1);
  }

  //--------------------------------------------------------
  //  compute_boundary_flux
  //    also sets time derivatives of boundary flux and
  //    heat sources
  //--------------------------------------------------------
  void ThermostatPowerVerletFiltered::compute_boundary_flux(FIELDS & fields)
  {
    ThermostatPowerVerlet::compute_boundary_flux(fields);
    
    // compute boundary flux rate of change
    fluxRoc_[TEMPERATURE] = 0.;
    atcTransfer_->compute_boundary_flux(fieldMask_,
                                        fieldsRoc_,
                                        fluxRoc_);

    // NOTE add time derivatives from sources, currently not implemented

    // compute extrinsic model rate of change
    (atcTransfer_->get_extrinsic_model_manager()).set_sources(fieldsRoc_,fluxRoc_);
    heatSourceRoc_ = fluxRoc_[TEMPERATURE];

    
  }

  //--------------------------------------------------------
  //  add_to_rhs:
  //            determines what if any contributions to the
  //            finite element equations are needed for
  //            consistency with the thermostat
  //--------------------------------------------------------
  void ThermostatPowerVerletFiltered::add_to_rhs(std::map<FieldName,DenseMatrix<double> > & rhs)
  {
    rhs[TEMPERATURE] += lambdaPowerFiltered_ + boundaryFlux_[TEMPERATURE];
  }

  //--------------------------------------------------------
  //  set_thermostat_rhs:
  //    sets up the right-hand side including boundary
  //    fluxes (coupling & prescribed), heat sources, and
  //    fixed (uncoupled) nodes
  //--------------------------------------------------------
  void ThermostatPowerVerletFiltered::set_thermostat_rhs(DENS_MAT & rhs_nodes)
  {
    // (a) for flux based : 
    // form rhs :  \int N_I r dV
    // vs  Wagner, CMAME, 2008 eq(24) RHS_I = 2/(3kB) flux_I
    // fluxes are set in ATC transfer
//     for (int i = 0; i < nNodes_; i++) {
//       rhs_nodes(i) = heatSource_(i,0) + filterScale_*heatSourceRoc_(i,0);
//     }
    rhs_nodes = heatSource_ + filterScale_*heatSourceRoc_;

    // (b) for ess. bcs
    // form rhs : {sum_a (N_Ia * v_ia * f_ia) - 0.5*(dtheta/dt)_I}
    DENS_MAT atomicPower;
    atcTransfer_->compute_atomic_power(atomicPower, v_, f_);
    DENS_MAT rhs_nodes_prescribed(nNodes_,1);
    atcTransfer_->restrict_volumetric_quantity(atomicPower, rhs_nodes_prescribed);
//     DENS_VEC myTemperatureRoc(nNodes_);
//     for (int i = 0; i < nNodes_; i++) {
//       myTemperatureRoc(i) = nodalTemperatureRoc_(i,0) + filterScale_*nodalTemperature2Roc_(i,0);
//     }
    rhs_nodes_prescribed -= 0.5*(mdMassMatrix_*(nodalTemperatureRoc_ + filterScale_*nodalTemperature2Roc_));

    // replace rhs for prescribed nodes
    for (int i = 0; i < nNodes_; i++) {
      if (isFixedNode_(i,0)) {
        rhs_nodes(i,0) = rhs_nodes_prescribed(i,0);
      }
    }
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void ThermostatPowerVerletFiltered::output(double dt, OUTPUT_LIST & outputData)
  {
    outputData["lambda"] = &(lambda_);
    outputData["nodalLambdaPower"] = &(lambdaPowerFiltered_);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class HooverVerletFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatHooverVerletFiltered::ThermostatHooverVerletFiltered(Thermostat * thermostat) :
    ThermostatPowerVerletFiltered(thermostat),
    nodeToOverlapMap_(atcTransfer_->get_node_to_overlap_map())
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  compute_thermostat:
  //    sets up and solves the thermostat equations since
  //    they are the same at different parts of the time
  //    step
  //--------------------------------------------------------
  void ThermostatHooverVerletFiltered::compute_thermostat(double dt)
  {
    // apply prescribed/extrinsic sources and fixed nodes
    ThermostatPowerVerletFiltered::compute_thermostat(dt);
    
    // set up Hoover rhs
    DENS_MAT rhs(nNodes_,1);
    set_hoover_rhs(rhs);
    
    // solve linear system for lambda
    solve_for_lambda(rhs);
    
    // NOTE modify these sections to be for hoover
    // compute force applied by lambda
    DENS_MAT lambdaForceHoover(nLocal_,nsd_);
    compute_lambda_force(v_,lambdaForceHoover);
    lambdaForce_ += lambdaForceHoover;

    // compute nodal atomic power from Hoover coupling
    // only add in contribution to uncoupled nodes
    if (atcTransfer_->use_localized_lambda())
      add_to_lambda_power(lambdaForceHoover,dt);
  }

  //--------------------------------------------------------
  //  set_hoover_rhs:
  //    sets up the right-hand side for fixed value,
  //    i.e. Hoover coupling
  //--------------------------------------------------------
  void ThermostatHooverVerletFiltered::set_hoover_rhs(DENS_MAT & rhs_nodes)
  {
    // form rhs : sum_a (N_Ia * v_ia * f_ia) - 0.5*M_MD*(dtheta/dt)_I
    DENS_MAT atomicPower(atcTransfer_->get_nlocal(),1);
    atcTransfer_->compute_atomic_power(atomicPower, v_, f_);
    atcTransfer_->restrict_volumetric_quantity(atomicPower, rhs_nodes);
//     DENS_VEC myTemperatureRoc(nNodes_);
//     for (int i = 0; i < nNodes_; i++) {
//       myTemperatureRoc(i) = nodalTemperatureRoc_(i,0) + filterScale_*nodalTemperature2Roc_(i,0);
//     }
    rhs_nodes -= 0.5*(mdMassMatrix_*(nodalTemperatureRoc_ + filterScale_*nodalTemperature2Roc_));
  }

  //--------------------------------------------------------
  //  add_to_nodal_lambda_power:
  //    determines the power exerted by the Hoover 
  //    thermostat at each FE node
  //--------------------------------------------------------
  void ThermostatHooverVerletFiltered::add_to_lambda_power(DENS_MAT & myLambdaForce,
                                                           double dt)
  {
    DENS_VEC atomicLambdaPower;
    atcTransfer_->compute_atomic_power(atomicLambdaPower,v_,myLambdaForce);
    DENS_VEC myNodalLambdaPower(nNodes_);
    atcTransfer_->restrict_volumetric_quantity(atomicLambdaPower,myNodalLambdaPower);
    myNodalLambdaPower *= 2.;
    for (int i = 0; i < nNodes_; ++i) {
      if (nodeToOverlapMap_(i)==-1)
        nodalAtomicLambdaPower_(i) += myNodalLambdaPower(i);
      else
        myNodalLambdaPower(i) = 0.;
    }
    timeFilter_->apply_post_step1(lambdaPowerFiltered_,myNodalLambdaPower,dt);
  }

};                   
