#include "Thermostat.h"
#include "ATC_Coupling.h"
#include "ATC_Error.h"
#include "PrescribedDataManager.h"
#include "ThermalTimeIntegrator.h"
#include "TransferOperator.h"

using namespace std;

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class Thermostat
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  Thermostat::Thermostat(ATC_Coupling * atc,
                         const string & regulatorPrefix) :
    AtomicRegulator(atc,regulatorPrefix),
    lambdaMaxIterations_(myLambdaMaxIterations)
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

    int argIndex = 0;
    if (strcmp(arg[argIndex],"thermal")==0) {
      argIndex++;
    
      // thermostat type
      /*! \page man_control_thermal fix_modify AtC control thermal 
        \section syntax
        fix_modify AtC control thermal <control_type> <optional_args>
        - control_type (string) = none | rescale | hoover | flux\n
      
        fix_modify AtC control thermal rescale <frequency> \n
        - frequency (int) = time step frequency for applying velocity rescaling \n

        fix_modify AtC control thermal hoover \n

        fix_modify AtC control thermal flux <boundary_integration_type(optional)> <face_set_id(optional)>\n
        - boundary_integration_type (string) = faceset | interpolate\n
        - face_set_id (string), optional = id of boundary face set, if not specified
        (or not possible when the atomic domain does not line up with 
        mesh boundaries) defaults to an atomic-quadrature approximate 
        evaulation, does not work with interpolate\n
        \section examples
        <TT> fix_modify AtC control thermal none </TT> \n
        <TT> fix_modify AtC control thermal rescale 10 </TT> \n
        <TT> fix_modify AtC control thermal hoover </TT> \n
        <TT> fix_modify AtC control thermal flux </TT> \n
        <TT> fix_modify AtC control thermal flux faceset bndy_faces </TT> \n
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
      if (strcmp(arg[argIndex],"none")==0) { // restore defaults
        regulatorTarget_ = NONE;
        couplingMode_ = UNCOUPLED;
        howOften_ = 1;
        boundaryIntegrationType_ = NO_QUADRATURE;
        foundMatch = true;
      }
      else if (strcmp(arg[argIndex],"rescale")==0) {
        argIndex++;
        howOften_ = atoi(arg[argIndex]);
        if (howOften_ < 1) {
          throw ATC_Error("Bad rescaling thermostat frequency");
        }
        else {
          regulatorTarget_ = FIELD;
          couplingMode_ = UNCOUPLED;
          boundaryIntegrationType_ = NO_QUADRATURE;
          foundMatch = true;
        }
      }
      else if (strcmp(arg[argIndex],"hoover")==0) {
        regulatorTarget_ = DYNAMICS;
        couplingMode_ = FIXED;
        howOften_ = 1;
        boundaryIntegrationType_ = NO_QUADRATURE;
        foundMatch = true;
      }
      else if (strcmp(arg[argIndex],"flux")==0) {
        regulatorTarget_ = DYNAMICS;
        couplingMode_ = FLUX;
        howOften_ = 1;
        argIndex++;
        
        boundaryIntegrationType_ = atc_->parse_boundary_integration(narg-argIndex,&arg[argIndex],boundaryFaceSet_);
        foundMatch = true;
      }
      // set parameters for numerical matrix solutions unique to this thermostat
      /*! \page man_control_thermal_correction_max_iterations fix_modify AtC control thermal correction_max_iterations
        \section syntax
        fix_modify AtC control thermal correction_max_iterations <max_iterations>
        - max_iterations (int) = maximum number of iterations that will be used by iterative matrix solvers\n

        \section examples
        <TT> fix_modify AtC control thermal correction_max_iterations 10 </TT> \n
        \section description
        Sets the maximum number of iterations to compute the 2nd order in time correction term for lambda with the fractional step method.  The method uses the same tolerance as the controller's matrix solver.
        \section restrictions
        only for use with thermal physics using the fractional step method.
        \section related
        \section default
        correction_max_iterations is 20
      */
      else if (strcmp(arg[argIndex],"correction_max_iterations")==0) {
        argIndex++;
        lambdaMaxIterations_ = atoi(arg[argIndex]);
        if (lambdaMaxIterations_ < 1) {
          throw ATC_Error("Bad correction maximum iteration count");
        }
        foundMatch = true;
      }
    }

    if (!foundMatch)
      foundMatch = AtomicRegulator::modify(narg,arg);
    if (foundMatch)
      needReset_ = true;
    return foundMatch;
  }

  //--------------------------------------------------------
  //  reset_lambda_contribution:
  //    resets the thermostat generated power to a
  //    prescribed value
  //--------------------------------------------------------
  void Thermostat::reset_lambda_contribution(const DENS_MAT & target)
  {
    DENS_MAN * lambdaPowerFiltered = regulator_data("LambdaPowerFiltered",1);
    *lambdaPowerFiltered = target;
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
  void Thermostat::construct_methods()
  {
    // get data associated with stages 1 & 2 of ATC_Method::initialize
    AtomicRegulator::construct_methods();

    if (atc_->reset_methods()) {
      // eliminate existing methods
      delete_method();

      // update time filter
      TimeIntegrator::TimeIntegrationType myIntegrationType = (atc_->time_integrator(TEMPERATURE))->time_integration_type();
      TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
      if (timeFilterManager->need_reset() ) {
        if (myIntegrationType == TimeIntegrator::GEAR) {
          timeFilter_ = timeFilterManager->construct(TimeFilterManager::EXPLICIT);
        }
        else if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP) {
          timeFilter_ = timeFilterManager->construct(TimeFilterManager::EXPLICIT_IMPLICIT);
        }
      }
      
      if (timeFilterManager->filter_dynamics()) {
        switch (regulatorTarget_) {
        case NONE: {
          regulatorMethod_ = new RegulatorMethod(this);
          break;
        }
        case FIELD: { // error check, rescale and filtering not supported together
          throw ATC_Error("Cannot use rescaling thermostat with time filtering");
          break;
        }
        case DYNAMICS: {
          switch (couplingMode_) {
          case FIXED: {
            if (use_lumped_lambda_solve()) {
              throw ATC_Error("Thermostat:construct_methods - lumped lambda solve cannot be used with Hoover thermostats");
            }
            if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP) {
              if (md_flux_nodes(TEMPERATURE)) {
                if (!md_fixed_nodes(TEMPERATURE) && (boundaryIntegrationType_ == NO_QUADRATURE)) {
                  // there are fluxes but no fixed or coupled nodes
                  regulatorMethod_ = new ThermostatIntegratorFluxFiltered(this,lambdaMaxIterations_);
                }
                else {
                  // there are both fixed and flux nodes
                  regulatorMethod_ = new ThermostatFluxFixedFiltered(this,lambdaMaxIterations_);
                }
              }
              else {
                // there are only fixed nodes
                regulatorMethod_ = new ThermostatIntegratorFixedFiltered(this,lambdaMaxIterations_);
              }
            }
            else {
              regulatorMethod_ = new ThermostatHooverVerletFiltered(this);
            }
            break;
          }
          case FLUX: {
            if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP) {
              if (use_lumped_lambda_solve()) {
                throw ATC_Error("Thermostat:construct_methods - lumped lambda solve has been depricated for fractional step thermostats");
              }
              if (md_fixed_nodes(TEMPERATURE)) {
                if (!md_flux_nodes(TEMPERATURE) && (boundaryIntegrationType_ == NO_QUADRATURE)) {
                // there are fixed nodes but no fluxes
                regulatorMethod_ = new ThermostatIntegratorFixedFiltered(this,lambdaMaxIterations_);
                }
                else {
                  // there are both fixed and flux nodes
                  regulatorMethod_ = new ThermostatFluxFixedFiltered(this,lambdaMaxIterations_);
                }
              }
              else {
                // there are only flux nodes
                regulatorMethod_ = new ThermostatIntegratorFluxFiltered(this,lambdaMaxIterations_);
              }
            }
            else {
              if (use_localized_lambda()) {
                if (!((atc_->prescribed_data_manager())->no_fluxes(TEMPERATURE)) &&
                    atc_->boundary_integration_type() != NO_QUADRATURE) {
                  throw ATC_Error("Cannot use flux coupling with localized lambda");
                }
              }
              regulatorMethod_ = new ThermostatPowerVerletFiltered(this);
            }
            break;
          }
          default:
            throw ATC_Error("Unknown coupling mode in Thermostat::initialize");
          }
          break;
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
          if (atc_->temperature_def()==KINETIC)
            regulatorMethod_ = new ThermostatRescale(this);
          else if (atc_->temperature_def()==TOTAL)
            regulatorMethod_ = new ThermostatRescaleMixedKePe(this);
          else
            throw ATC_Error("Unknown temperature definition");
          break;
        }
        case DYNAMICS: {
          switch (couplingMode_) {
          case FIXED: {
            if (use_lumped_lambda_solve()) {
              throw ATC_Error("Thermostat:construct_methods - lumped lambda solve cannot be used with Hoover thermostats");
            }
            if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP) {
              if (md_flux_nodes(TEMPERATURE)) {
                if (!md_fixed_nodes(TEMPERATURE) && (boundaryIntegrationType_ == NO_QUADRATURE)) {
                  // there are fluxes but no fixed or coupled nodes
                  regulatorMethod_ = new ThermostatIntegratorFlux(this,lambdaMaxIterations_);
                }
                else {
                  // there are both fixed and flux nodes
                  regulatorMethod_ = new ThermostatFluxFixed(this,lambdaMaxIterations_);
                }
              }
              else {
                // there are only fixed nodes
                regulatorMethod_ = new ThermostatIntegratorFixed(this,lambdaMaxIterations_);
              }
            }
            else {
              regulatorMethod_ = new ThermostatHooverVerlet(this);
            }
            break;
          }
          case FLUX: {
            if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP) {
              if (use_lumped_lambda_solve()) {
                throw ATC_Error("Thermostat:construct_methods - lumped lambda solve has been depricated for fractional step thermostats");
              }
              if (md_fixed_nodes(TEMPERATURE)) {
                if (!md_flux_nodes(TEMPERATURE) && (boundaryIntegrationType_ == NO_QUADRATURE)) {
                  // there are fixed nodes but no fluxes
                  regulatorMethod_ = new ThermostatIntegratorFixed(this,lambdaMaxIterations_);
                }
                else {
                  // there are both fixed and flux nodes
                  regulatorMethod_ = new ThermostatFluxFixed(this,lambdaMaxIterations_);
                }
              }
              else {
                // there are only flux nodes
                regulatorMethod_ = new ThermostatIntegratorFlux(this,lambdaMaxIterations_);
              }
            }
            else {
              if (use_localized_lambda()) {
                if (!((atc_->prescribed_data_manager())->no_fluxes(TEMPERATURE)) &&
                    atc_->boundary_integration_type() != NO_QUADRATURE) {
                  throw ATC_Error("Cannot use flux coupling with localized lambda");
                }
              }
              regulatorMethod_ = new ThermostatPowerVerlet(this);
            }
            break;
          }
          default:
            throw ATC_Error("Unknown coupling mode in Thermostat::initialize");
          }
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
  //  Class ThermostatShapeFunction
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ThermostatShapeFunction::ThermostatShapeFunction(AtomicRegulator * thermostat,
                                                   const string & regulatorPrefix) :
    RegulatorShapeFunction(thermostat,regulatorPrefix),
    mdMassMatrix_(atc_->set_mass_mat_md(TEMPERATURE)),
    atomVelocities_(NULL)
  {
    fieldMask_(TEMPERATURE,FLUX) = true;
    lambda_ = atomicRegulator_->regulator_data(regulatorPrefix_+"LambdaEnergy",1); // data associated with stage 3 in ATC_Method::initialize
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void ThermostatShapeFunction::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());
    RegulatorShapeFunction::construct_transfers();

    // get atom velocity data from manager
    atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY);

    // construct lambda evaluated at atom locations
    atomLambdas_ = new FtaShapeFunctionProlongation(atc_,
                                                    lambda_,
                                                    interscaleManager.per_atom_sparse_matrix("Interpolant"));
    interscaleManager.add_per_atom_quantity(atomLambdas_,regulatorPrefix_+"AtomLambdaEnergy");
  }

  //---------------------------------------------------------
  //  set_weights:
  //    set the diagonal weighting matrix to be the atomic
  //    temperatures
  //---------------------------------------------------------
  void ThermostatShapeFunction::set_weights()
  {
    if (this->use_local_shape_functions()) {
      VelocitySquaredMapped * myWeights = new VelocitySquaredMapped(atc_,lambdaAtomMap_);
      weights_ = myWeights;
      (atc_->interscale_manager()).add_per_atom_quantity(myWeights,
                                                         regulatorPrefix_+"AtomVelocitySquaredMapped");
    }
    else {
      VelocitySquared * myWeights = new VelocitySquared(atc_);
      weights_ = myWeights;
      (atc_->interscale_manager()).add_per_atom_quantity(myWeights,
                                                         regulatorPrefix_+"AtomVelocitySquared");
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
  ThermostatRescale::ThermostatRescale(AtomicRegulator * thermostat) :
    ThermostatShapeFunction(thermostat),
    nodalTemperature_(atc_->field(TEMPERATURE)),
    atomVelocityRescalings_(NULL)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void ThermostatRescale::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // set up node mappings
    create_node_maps();

    // set up data for linear solver
    shapeFunctionMatrix_ = new LambdaCouplingMatrix(atc_,nodeToOverlapMap_);
    (atc_->interscale_manager()).add_per_atom_sparse_matrix(shapeFunctionMatrix_,
                                                            regulatorPrefix_+"LambdaCouplingMatrixEnergy");
    linearSolverType_ = AtomicRegulator::CG_SOLVE;

    // base class transfers
    ThermostatShapeFunction::construct_transfers();

    // velocity rescaling factor
    atomVelocityRescalings_ = new AtomicVelocityRescaleFactor(atc_,atomLambdas_);
    interscaleManager.add_per_atom_quantity(atomVelocityRescalings_,
                                            regulatorPrefix_+"AtomVelocityRescaling");
  }

  //---------------------------------------------------------
  //  set_weights:
  //    set the diagonal weighting matrix to be the atomic
  //    temperatures
  //---------------------------------------------------------
  void ThermostatRescale::set_weights()
  {
    weights_ = (atc_->interscale_manager()).per_atom_quantity("AtomicEnergyForTemperature");
  }

  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the thermostat in the post corrector phase
  //--------------------------------------------------------
  void ThermostatRescale::apply_post_corrector(double dt)
  {
    compute_thermostat(dt);

    // application of rescaling lambda due
    apply_to_atoms(atomVelocities_);
  }

  //--------------------------------------------------------
  //  compute_thermostat
  //            manages the solution of the
  //            thermostat equations and variables
  //--------------------------------------------------------
  void ThermostatRescale::compute_thermostat(double dt)
  {
    // compute right-hand side
    this->set_rhs(_rhs_);
    
    // solve equations
    solve_for_lambda(_rhs_,lambda_->set_quantity());
  }

  //--------------------------------------------------------
  //  set_rhs:
  //    constructs the RHS vector with the target
  //    temperature
  //--------------------------------------------------------
  void ThermostatRescale::set_rhs(DENS_MAT & rhs)
  {
    rhs = mdMassMatrix_.quantity()*nodalTemperature_.quantity();
  }

  //--------------------------------------------------------
  //  apply_lambda_to_atoms:
  //    applies the velocity rescale with an existing lambda
  //    note oldAtomicQuantity and dt are not used
  //--------------------------------------------------------
  void ThermostatRescale::apply_to_atoms(PerAtomQuantity<double> * atomVelocities)
  {
    *atomVelocities *= atomVelocityRescalings_->quantity();
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void ThermostatRescale::output(OUTPUT_LIST & outputData)
  {
    DENS_MAT & lambda(lambda_->set_quantity());
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData["LambdaEnergy"] = &lambda;
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatRescaleMixedKePe
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ThermostatRescaleMixedKePe::ThermostatRescaleMixedKePe(AtomicRegulator * thermostat) :
    ThermostatRescale(thermostat),
    nodalAtomicFluctuatingPotentialEnergy_(NULL)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void ThermostatRescaleMixedKePe::construct_transfers()
  {
    ThermostatRescale::construct_transfers();
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // get fluctuating PE at nodes
    nodalAtomicFluctuatingPotentialEnergy_ = 
      interscaleManager.dense_matrix("NodalAtomicFluctuatingPotentialEnergy");
  }

  //---------------------------------------------------------
  //  set_weights:
  //    set the diagonal weighting matrix to be the atomic
  //    temperatures
  //---------------------------------------------------------
  void ThermostatRescaleMixedKePe::set_weights()
  {
    weights_ = (atc_->interscale_manager()).per_atom_quantity("AtomicTwiceKineticEnergy");
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void ThermostatRescaleMixedKePe::initialize()
  {
    ThermostatRescale::initialize();
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // multipliers for KE and PE
    AtomicEnergyForTemperature * atomEnergyForTemperature =
      static_cast<AtomicEnergyForTemperature * >(interscaleManager.per_atom_quantity("AtomicEnergyForTemperature"));
    keMultiplier_ = atomEnergyForTemperature->kinetic_energy_multiplier();
    peMultiplier_ = 2. - keMultiplier_;
    keMultiplier_ /= 2.; // account for use of 2 X KE in matrix equation
  }

  //--------------------------------------------------------
  //  set_rhs:
  //    accounts for potential energy contribution to
  //    definition of atomic temperature
  //--------------------------------------------------------
  void ThermostatRescaleMixedKePe::set_rhs(DENS_MAT & rhs)
  {
    ThermostatRescale::set_rhs(rhs);
    rhs -=  peMultiplier_*(nodalAtomicFluctuatingPotentialEnergy_->quantity());
    rhs /= keMultiplier_;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatFsSolver
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ThermostatFsSolver::ThermostatFsSolver(AtomicRegulator * thermostat,
                                         int lambdaMaxIterations,
                                         const string & regulatorPrefix) :
    RegulatorShapeFunction(thermostat,regulatorPrefix),
    lambdaMaxIterations_(lambdaMaxIterations),
    rhsLambdaSquared_(NULL),
    dtFactor_(1.)
  {
    fieldMask_(TEMPERATURE,FLUX) = true;
    lambda_ = atomicRegulator_->regulator_data(regulatorPrefix_+"LambdaEnergy",1); // data associated with stage 3 in ATC_Method::initialize
  }

  //--------------------------------------------------------
  //  initialize
  //    creates mapping from all nodes to those to which
  //    the thermostat applies
  //--------------------------------------------------------
  void ThermostatFsSolver::initialize()
  {
    RegulatorShapeFunction::initialize();

    rhsMap_.resize(overlapToNodeMap_->nRows(),1);
    DENS_MAT rhsMapGlobal(nNodes_,1);
    const set<int> & applicationNodes(applicationNodes_->quantity());
    for (int i = 0; i < nNodes_; i++) {
      if (applicationNodes.find(i) != applicationNodes.end()) {
        rhsMapGlobal(i,0) = 1.;
      }
      else {
        rhsMapGlobal(i,0) = 0.;
      }
    }
    map_unique_to_overlap(rhsMapGlobal,rhsMap_);
  }

  //---------------------------------------------------------
  //  set_weights:
  //    set the diagonal weighting matrix to be the atomic
  //    temperatures
  //---------------------------------------------------------
  void ThermostatFsSolver::set_weights()
  {
    if (this->use_local_shape_functions()) {
      VelocitySquaredMapped * myWeights = new VelocitySquaredMapped(atc_,lambdaAtomMap_);
      weights_ = myWeights;
      (atc_->interscale_manager()).add_per_atom_quantity(myWeights,
                                                         regulatorPrefix_+"AtomVelocitySquaredMapped");
    }
    else {
      VelocitySquared * myWeights = new VelocitySquared(atc_);
      weights_ = myWeights;
      (atc_->interscale_manager()).add_per_atom_quantity(myWeights,
                                                         regulatorPrefix_+"AtomVelocitySquared");
    }
  }

  //--------------------------------------------------------
  //  compute_lambda:
  //   solves linear system for lambda, if the
  //   bool is true it iterators to a non-linear solution
  //--------------------------------------------------------
  void ThermostatFsSolver::compute_lambda(const DENS_MAT & rhs,
                                          bool iterateSolution)
  {
    // solve linear system for lambda guess
    DENS_MAT & lambda(lambda_->set_quantity());
    solve_for_lambda(rhs,lambda);

    // iterate to solution
    if (iterateSolution) {
      iterate_lambda(rhs);
    }
  }

  //--------------------------------------------------------
  //  iterate_lambda:
  //    iteratively solves the equations for lambda
  //    for the higher order dt corrections, assuming
  //    an initial guess for lambda
  //--------------------------------------------------------
  void ThermostatFsSolver::iterate_lambda(const MATRIX & rhs)
  {
    int nNodeOverlap = overlapToNodeMap_->nRows();
    DENS_VEC _lambdaOverlap_(nNodeOverlap);
    DENS_MAT & lambda(lambda_->set_quantity());
    map_unique_to_overlap(lambda,_lambdaOverlap_);
    double factor = 0.5*dtFactor_*atc_->dt();

    _lambdaOld_.resize(nNodes_,1);
    _rhsOverlap_.resize(nNodeOverlap,1);
    map_unique_to_overlap(rhs,_rhsOverlap_);
    _rhsTotal_.resize(nNodeOverlap);

    // solve assuming we get initial guess for lambda
    double error(-1.);
    for (int i = 0; i < lambdaMaxIterations_; ++i) {
      _lambdaOld_ = lambda;

      // solve the system with the new rhs
      const DENS_MAT & rhsLambdaSquared(rhsLambdaSquared_->quantity());
      for (int i = 0; i < nNodeOverlap; i++) {
        if (rhsMap_(i,0) == 1.) {
          _rhsTotal_(i) = _rhsOverlap_(i,0) + factor*rhsLambdaSquared(i,0);
        }
        else {
          _rhsTotal_(i) = 0.;
        }
      }
      matrixSolver_->execute(_rhsTotal_,_lambdaOverlap_);

      // check convergence
      map_overlap_to_unique(_lambdaOverlap_,lambda);
      lambda_->force_reset();
      DENS_MAT difference = lambda-_lambdaOld_;
      error = difference.col_norm()/_lambdaOld_.col_norm();
      if (error < tolerance_)
        break;
    }
    
    if (error >= tolerance_) {
      stringstream message;
      message << "WARNING: Iterative solve for lambda failed to converge after " << lambdaMaxIterations_ << " iterations, final tolerance was " << error << "\n";
      ATC::LammpsInterface::instance()->print_msg(message.str());
    }
  }
 
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatGlcFs
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ThermostatGlcFs::ThermostatGlcFs(AtomicRegulator * thermostat,
                                   int lambdaMaxIterations,
                                   const string & regulatorPrefix) :
    RegulatorMethod(thermostat,regulatorPrefix),
    lambdaSolver_(NULL),
    mdMassMatrix_(atc_->set_mass_mat_md(TEMPERATURE)),
    atomVelocities_(NULL),
    temperature_(atc_->field(TEMPERATURE)),
    timeFilter_(atomicRegulator_->time_filter()),
    nodalAtomicLambdaPower_(NULL),
    lambdaPowerFiltered_(NULL),
    atomLambdas_(NULL),
    atomThermostatForces_(NULL),
    atomMasses_(NULL),
    isFirstTimestep_(true),
    nodalAtomicEnergy_(NULL),
    atomPredictedVelocities_(NULL),
    nodalAtomicPredictedEnergy_(NULL),
    firstHalfAtomForces_(NULL)
  {
    // construct/obtain data corresponding to stage 3 of ATC_Method::initialize
    nodalAtomicLambdaPower_ = thermostat->regulator_data(regulatorPrefix_+"NodalAtomicLambdaPower",1);
    lambdaPowerFiltered_ = atomicRegulator_->regulator_data(regulatorPrefix_+"LambdaPowerFiltered",1);
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void ThermostatGlcFs::construct_transfers()
  {
    lambdaSolver_->construct_transfers();
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // get atom velocity data from manager
    atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY);

    // construct lambda evaluated at atom locations
    atomLambdas_ = interscaleManager.per_atom_quantity(regulatorPrefix_+"AtomLambdaEnergy");
    if (!atomLambdas_) {
      DENS_MAN * lambda = atomicRegulator_->regulator_data(regulatorPrefix_+"LambdaEnergy",1);
      atomLambdas_ = new FtaShapeFunctionProlongation(atc_,
                                                      lambda,
                                                      interscaleManager.per_atom_sparse_matrix("Interpolant"));
      interscaleManager.add_per_atom_quantity(atomLambdas_,regulatorPrefix_+"AtomLambdaEnergy");
    }

    // get data from manager
    atomMasses_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS),
    nodalAtomicEnergy_ = interscaleManager.dense_matrix("NodalAtomicEnergy");
    
    // thermostat forces based on lambda and the atomic velocities
    atomThermostatForces_ = new AtomicThermostatForce(atc_,atomLambdas_);
    interscaleManager.add_per_atom_quantity(atomThermostatForces_,
                                            regulatorPrefix_+"AtomThermostatForce");

    // predicted temperature quantities:  atom velocities, atom energies, and restricted atom energies
    atomPredictedVelocities_ = new AtcAtomQuantity<double>(atc_,atc_->nsd());
    // MAKE THINGS WORK WITH ONLY ONE PREDICTED VELOCITY, CHECK IT EXISTS
    interscaleManager.add_per_atom_quantity(atomPredictedVelocities_,
                                            regulatorPrefix_+"AtomicPredictedVelocities");
    AtomicEnergyForTemperature * atomPredictedEnergyForTemperature = new TwiceKineticEnergy(atc_,
                                                                                            atomPredictedVelocities_);
    interscaleManager.add_per_atom_quantity(atomPredictedEnergyForTemperature,
                                            regulatorPrefix_+"AtomicPredictedTwiceKineticEnergy");
    nodalAtomicPredictedEnergy_ = new AtfShapeFunctionRestriction(atc_,
                                                                  atomPredictedEnergyForTemperature,
                                                                  interscaleManager.per_atom_sparse_matrix("Interpolant"));
    interscaleManager.add_dense_matrix(nodalAtomicPredictedEnergy_,
                                       regulatorPrefix_+"NodalAtomicPredictedEnergy");
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void ThermostatGlcFs::initialize()
  {
    RegulatorMethod::initialize();
    
    TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
    if (!timeFilterManager->end_equilibrate()) {
      // we should reset lambda and lambdaForce to zero in this case
      // implies an initial condition of 0 for the filtered nodal lambda power
      // initial conditions will always be needed when using time filtering
      // however, the fractional step scheme must assume the instantaneous
      // nodal lambda power is 0 initially because all quantities are in delta form
      lambdaSolver_->initialize(); // ensures initial lambda force is zero
      lambdaSolver_->set_lambda_to_value(0.);
      *nodalAtomicLambdaPower_ = 0.; // energy change due to thermostats
      *lambdaPowerFiltered_ = 0.; // filtered energy change due to thermostats
    }
    else {
      lambdaSolver_->initialize();
      // we can grab lambda power variables using time integrator and atc transfer in cases for equilibration
    }

    // sets up time filter for cases where variables temporally filtered
    if (timeFilterManager->need_reset()) {
      // the form of this integrator implies no time filters that require history data can be used
      timeFilter_->initialize(nodalAtomicLambdaPower_->quantity());
    }

    atomThermostatForces_->quantity(); // initialize
    atomThermostatForces_->fix_quantity();
    firstHalfAtomForces_ = atomThermostatForces_; // initialize
#ifdef OBSOLETE
    compute_rhs_map();
#endif
  }

  //--------------------------------------------------------
  //  reset_atom_materials:
  //    resets the localized atom to material map
  //--------------------------------------------------------
  void ThermostatGlcFs::reset_atom_materials(const Array<int> & elementToMaterialMap,
                                             const MatrixDependencyManager<DenseMatrix, int> * atomElement)
  {
    lambdaSolver_->reset_atom_materials(elementToMaterialMap,
                                        atomElement);
  }

  //--------------------------------------------------------
  //  apply_to_atoms:
  //     determines what if any contributions to the
  //     atomic moition is needed for
  //     consistency with the thermostat
  //     and computes the instantaneous induced power
  //--------------------------------------------------------
  void ThermostatGlcFs::apply_to_atoms(PerAtomQuantity<double> * atomicVelocity,
                                      const DENS_MAN * nodalAtomicEnergy,
                                      const DENS_MAT & lambdaForce,
                                      DENS_MAT & nodalAtomicLambdaPower,
                                      double dt)
  {
    // compute initial contributions to lambda power
    nodalAtomicLambdaPower = nodalAtomicEnergy->quantity();
    nodalAtomicLambdaPower *= -1.;

    // apply lambda force to atoms
    _velocityDelta_ = lambdaForce;
    _velocityDelta_ /= atomMasses_->quantity();
    _velocityDelta_ *= dt;
    (*atomicVelocity) += _velocityDelta_;

    // finalize lambda power
    nodalAtomicLambdaPower += nodalAtomicEnergy->quantity();
  }

  //--------------------------------------------------------
  //  full_prediction:
  //    flag to perform a full prediction calcalation
  //    for lambda rather than using the old value
  //--------------------------------------------------------
  bool ThermostatGlcFs::full_prediction()
  {
    if (isFirstTimestep_ || ((atc_->atom_to_element_map_type() == EULERIAN)
                             && (atc_->atom_to_element_map_frequency() > 1)
                             && (atc_->step() % atc_->atom_to_element_map_frequency() == 0 ))) {
      return true;
    }
    return false;
  }

  //--------------------------------------------------------
  //  apply_predictor:
  //    apply the thermostat to the atoms in the first step
  //    of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatGlcFs::apply_pre_predictor(double dt)
  {
    DENS_MAT & myLambdaPowerFiltered(lambdaPowerFiltered_->set_quantity());
    DENS_MAT & myNodalAtomicLambdaPower(nodalAtomicLambdaPower_->set_quantity());

    // update filtered power
    timeFilter_->apply_pre_step1(myLambdaPowerFiltered,myNodalAtomicLambdaPower,dt); // equivalent update to measure power as change in energy due to thermostat

    // apply lambda force to atoms and compute instantaneous lambda power for first half of time step
    this->apply_to_atoms(atomVelocities_,nodalAtomicEnergy_,
                         firstHalfAtomForces_->quantity(),
                         myNodalAtomicLambdaPower,0.5*dt);

    // update nodal variables for first half of time step
    this->add_to_energy(myNodalAtomicLambdaPower,deltaEnergy1_,0.5*dt);

    // start update of filtered lambda power
    myNodalAtomicLambdaPower = 0.; // temporary power for first part of update
    timeFilter_->apply_post_step1(myLambdaPowerFiltered,myNodalAtomicLambdaPower,dt);
  }

  //--------------------------------------------------------
  //  apply_pre_corrector:
  //    apply the thermostat to the atoms in the first part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatGlcFs::apply_pre_corrector(double dt)
  {
    (*atomPredictedVelocities_) = atomVelocities_->quantity();

    // do full prediction if we just redid the shape functions
    if (full_prediction()) {
      this->compute_lambda(dt);
      
      atomThermostatForces_->unfix_quantity();  // allow update of atomic force applied by lambda
    }

    // apply lambda force to atoms and compute instantaneous lambda power to predict second half of time step
    DENS_MAT & myNodalAtomicLambdaPower(nodalAtomicLambdaPower_->set_quantity());
    apply_to_atoms(atomPredictedVelocities_,
                   nodalAtomicPredictedEnergy_,
                   firstHalfAtomForces_->quantity(),
                   myNodalAtomicLambdaPower,0.5*dt);

    if (full_prediction())
      atomThermostatForces_->fix_quantity();
    
    // update predicted nodal variables for second half of time step
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
  void ThermostatGlcFs::apply_post_corrector(double dt)
  {
    // remove predicted power effects
    DENS_MAT & myTemperature(temperature_.set_quantity());
    atc_->apply_inverse_mass_matrix(deltaEnergy2_,TEMPERATURE);
    myTemperature -= deltaEnergy2_;

    // set up equation and update lambda
    this->compute_lambda(dt);

    // apply lambda force to atoms and compute instantaneous lambda power for second half of time step
    DENS_MAT & myNodalAtomicLambdaPower(nodalAtomicLambdaPower_->set_quantity());
    // allow computation of force applied by lambda using current velocities
    atomThermostatForces_->unfix_quantity(); 
    atomThermostatForces_->quantity();
    atomThermostatForces_->fix_quantity();
    apply_to_atoms(atomVelocities_,nodalAtomicEnergy_,
                   atomThermostatForces_->quantity(),
                   myNodalAtomicLambdaPower,0.5*dt);

    // finalize filtered lambda power by adding latest contribution
    timeFilter_->apply_post_step2(lambdaPowerFiltered_->set_quantity(),
                                  myNodalAtomicLambdaPower,dt);

    // update nodal variables for second half of time step
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
  void ThermostatGlcFs::compute_lambda(double dt,
                                       bool iterateSolution)
  {
    // set up rhs for lambda equation
    this->set_thermostat_rhs(rhs_,0.5*dt);

    // solve system
    lambdaSolver_->compute_lambda(rhs_,iterateSolution);
#ifdef OBSOLETE
    // solve linear system for lambda guess
    DENS_MAT & lambda(lambda_->set_quantity());
    solve_for_lambda(rhs_,lambda);
    
    // iterate to solution
    if (iterateSolution) {
      iterate_lambda(rhs_);
    }
#endif
  }

  //--------------------------------------------------------
  //  compute_boundary_flux
  //    default computation of boundary flux based on
  //    finite
  //--------------------------------------------------------
  void ThermostatGlcFs::compute_boundary_flux(FIELDS & fields)
  {
    
    lambdaSolver_->compute_boundary_flux(fields);
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void ThermostatGlcFs::output(OUTPUT_LIST & outputData)
  {
    _lambdaPowerOutput_ = nodalAtomicLambdaPower_->quantity();
    // approximate value for lambda power
    double dt =  LammpsInterface::instance()->dt();
    _lambdaPowerOutput_ *= (2./dt);
    DENS_MAT & lambda((atomicRegulator_->regulator_data(regulatorPrefix_+"LambdaEnergy",1))->set_quantity());
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData[regulatorPrefix_+"LambdaEnergy"] = &lambda;
      outputData[regulatorPrefix_+"NodalLambdaPower"] = &(_lambdaPowerOutput_);
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatSolverFlux
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatSolverFlux::ThermostatSolverFlux(AtomicRegulator * thermostat,
                                             int lambdaMaxIterations,
                                             const string & regulatorPrefix) :
    ThermostatFsSolver(thermostat,lambdaMaxIterations,regulatorPrefix)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void ThermostatSolverFlux::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // set up node mappings
    create_node_maps();

    // set up data for linear solver
    shapeFunctionMatrix_ = new LambdaCouplingMatrix(atc_,nodeToOverlapMap_);
    interscaleManager.add_per_atom_sparse_matrix(shapeFunctionMatrix_,
                                                 regulatorPrefix_+"LambdaCouplingMatrixEnergy");
    if (elementMask_) {
      lambdaAtomMap_ = new AtomToElementset(atc_,elementMask_);
      interscaleManager.add_per_atom_int_quantity(lambdaAtomMap_,
                                                  regulatorPrefix_+"LambdaAtomMap");
    }
    if (atomicRegulator_->use_localized_lambda()) {
      linearSolverType_ = AtomicRegulator::RSL_SOLVE;
    }
    else {
      linearSolverType_ = AtomicRegulator::CG_SOLVE;
    }

    // base class transfers
    ThermostatFsSolver::construct_transfers();

    // add transfers for computation of extra RHS term accounting of O(lambda^2)
    // lambda squared followed by fractional step RHS contribution
    atomLambdas_ = (interscaleManager.per_atom_quantity(regulatorPrefix_+"AtomLambdaEnergy"));
    if (!atomLambdas_) {
      atomLambdas_ = new FtaShapeFunctionProlongation(atc_,
                                                      lambda_,
                                                      interscaleManager.per_atom_sparse_matrix("Interpolant"));
      interscaleManager.add_per_atom_quantity(atomLambdas_,regulatorPrefix_+"AtomLambdaEnergy");
    }
    LambdaSquared * lambdaSquared = new LambdaSquared(atc_,
                                                      interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS),
                                                      weights_,
                                                      atomLambdas_);
    interscaleManager.add_per_atom_quantity(lambdaSquared,
                                            regulatorPrefix_+"LambdaSquaredMapped");
    rhsLambdaSquared_ = new AtfShapeFunctionRestriction(atc_,lambdaSquared,shapeFunctionMatrix_);
    interscaleManager.add_dense_matrix(rhsLambdaSquared_,
                                            regulatorPrefix_+"RhsLambdaSquared");
  }

  //--------------------------------------------------------
  //  construct_regulated_nodes:
  //    constructs the set of nodes being regulated
  //--------------------------------------------------------
  void ThermostatSolverFlux::construct_regulated_nodes()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // matrix requires all entries even if localized for correct lumping
    regulatedNodes_ = interscaleManager.set_int(regulatorPrefix_+"ThermostatRegulatedNodes");
    if (!regulatedNodes_) {
      regulatedNodes_ = new RegulatedNodes(atc_);
      interscaleManager.add_set_int(regulatedNodes_,
                                    regulatorPrefix_+"ThermostatRegulatedNodes");
    }
    
    // if localized monitor nodes with applied fluxes
    if (atomicRegulator_->use_localized_lambda()) {
      if ((atomicRegulator_->coupling_mode() == Thermostat::FLUX) && (atomicRegulator_->boundary_integration_type() != NO_QUADRATURE)) {
        // include boundary nodes
        applicationNodes_ = new FluxBoundaryNodes(atc_);
        
        boundaryNodes_ = new BoundaryNodes(atc_);
        interscaleManager.add_set_int(boundaryNodes_,
                                      regulatorPrefix_+"ThermostatBoundaryNodes");
      }
      else {
        // fluxed nodes only
        applicationNodes_ = new FluxNodes(atc_);
      }
      interscaleManager.add_set_int(applicationNodes_,
                                    regulatorPrefix_+"ThermostatApplicationNodes");
    }
    else {
      applicationNodes_ = regulatedNodes_;
    }

    // special set of boundary elements for boundary flux quadrature  
    if ((atomicRegulator_->boundary_integration_type() == FE_INTERPOLATION)
        && (atomicRegulator_->use_localized_lambda())) {
      elementMask_ = interscaleManager.dense_matrix_bool(regulatorPrefix_+"BoundaryElementMask");
      if (!elementMask_) {
        elementMask_ = new ElementMaskNodeSet(atc_,applicationNodes_);
        interscaleManager.add_dense_matrix_bool(elementMask_,
                                                regulatorPrefix_+"BoundaryElementMask");
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatIntegratorFlux
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatIntegratorFlux::ThermostatIntegratorFlux(AtomicRegulator * thermostat,
                                                     int lambdaMaxIterations,
                                                     const string & regulatorPrefix) :
    ThermostatGlcFs(thermostat,lambdaMaxIterations,regulatorPrefix),
    heatSource_(atc_->atomic_source(TEMPERATURE))
  {
    lambdaSolver_ = new ThermostatSolverFlux(thermostat,
                                             lambdaMaxIterations,
                                             regulatorPrefix);
  }

  //--------------------------------------------------------
  //  add_to_temperature
  //    add in contributions from lambda power and boundary
  //    flux to the FE temperature
  //--------------------------------------------------------
  void ThermostatIntegratorFlux::add_to_energy(const DENS_MAT & nodalLambdaPower,
                                               DENS_MAT & deltaEnergy,
                                               double dt)
  {
    deltaEnergy.resize(nNodes_,1);
    const DENS_MAT & myBoundaryFlux(boundaryFlux_[TEMPERATURE].quantity());
    for (int i = 0; i < nNodes_; i++) {
      deltaEnergy(i,0) = nodalLambdaPower(i,0) + dt*myBoundaryFlux(i,0);
    }
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void ThermostatIntegratorFlux::initialize()
  {
    ThermostatGlcFs::initialize();

    // timestep factor
    lambdaSolver_->set_timestep_factor(1.);
  }

  //--------------------------------------------------------
  //  set_thermostat_rhs:
  //    sets up the right-hand side including boundary
  //    fluxes (coupling & prescribed), heat sources, and
  //    fixed (uncoupled) nodes
  //--------------------------------------------------------
  void ThermostatIntegratorFlux::set_thermostat_rhs(DENS_MAT & rhs,
                                                    double dt)
  {
    
    // only tested with flux != 0 + ess bc = 0

    // (a) for flux based : 
    // form rhs :  2/3kB * W_I^-1 * \int N_I r dV
    // vs  Wagner, CMAME, 2008 eq(24) RHS_I = 2/(3kB) flux_I
    // fluxes are set in ATC transfer
    const DENS_MAT & heatSource(heatSource_.quantity());
#if true
    const set<int> & applicationNodes((lambdaSolver_->application_nodes())->quantity());
    rhs.resize(nNodes_,1);
    for (int i = 0; i < nNodes_; i++) {
      if (applicationNodes.find(i) != applicationNodes.end()) {
        rhs(i,0) = heatSource(i,0);
      }
      else {
        rhs(i,0) = 0.;
      }
    }
#else
    rhs.resize(nNodes_,1);
    for (int i = 0; i < nNodes_; i++) {
      rhs(i,0) = heatSource(i,0);
    }
#endif
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatSolverFixed
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatSolverFixed::ThermostatSolverFixed(AtomicRegulator * thermostat,
                                               int lambdaMaxIterations,
                                               const string & regulatorPrefix) :
    ThermostatFsSolver(thermostat,lambdaMaxIterations,regulatorPrefix)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void ThermostatSolverFixed::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // set up node mappings
    create_node_maps();

    // determine if map is needed and set up if so
    if (this->use_local_shape_functions()) {
      lambdaAtomMap_ = new AtomToElementset(atc_,elementMask_);
      interscaleManager.add_per_atom_int_quantity(lambdaAtomMap_,
                                                  regulatorPrefix_+"LambdaAtomMap");
      shapeFunctionMatrix_ = new LocalLambdaCouplingMatrix(atc_,
                                                           lambdaAtomMap_,
                                                           nodeToOverlapMap_);
    }
    else {
      shapeFunctionMatrix_ = new LambdaCouplingMatrix(atc_,nodeToOverlapMap_);
    }
    interscaleManager.add_per_atom_sparse_matrix(shapeFunctionMatrix_,
                                                 regulatorPrefix_+"LambdaCouplingMatrixEnergy");
    linearSolverType_ = AtomicRegulator::CG_SOLVE;

    // base class transfers, e.g. weights
    ThermostatFsSolver::construct_transfers();

    // add transfers for computation of extra RHS term accounting of O(lambda^2)
    // lambda squared followed by fractional step RHS contribution
    atomLambdas_ = interscaleManager.per_atom_quantity(regulatorPrefix_+"AtomLambdaEnergy");
    if (!atomLambdas_) {
      atomLambdas_ = new FtaShapeFunctionProlongation(atc_,
                                                      lambda_,
                                                      interscaleManager.per_atom_sparse_matrix("Interpolant"));
      interscaleManager.add_per_atom_quantity(atomLambdas_,regulatorPrefix_+"AtomLambdaEnergy");
    }
    if (lambdaAtomMap_) {
      LambdaSquaredMapped * lambdaSquared = new LambdaSquaredMapped(atc_,
                                                                    lambdaAtomMap_,
                                                                    interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS),
                                                                    weights_,
                                                                    atomLambdas_);
      interscaleManager.add_per_atom_quantity(lambdaSquared,
                                              regulatorPrefix_+"LambdaSquared");
      rhsLambdaSquared_ = new AtfShapeFunctionRestriction(atc_,lambdaSquared,shapeFunctionMatrix_);
    }
    else {
      LambdaSquared * lambdaSquared = new LambdaSquared(atc_,
                                                       interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS),
                                                        weights_,
                                                        atomLambdas_);
      interscaleManager.add_per_atom_quantity(lambdaSquared,
                                              regulatorPrefix_+"LambdaSquaredMapped");
      rhsLambdaSquared_ = new AtfShapeFunctionRestriction(atc_,lambdaSquared,shapeFunctionMatrix_);
    }
    interscaleManager.add_dense_matrix(rhsLambdaSquared_,
                                       regulatorPrefix_+"RhsLambdaSquared");
  }

  //--------------------------------------------------------
  //  construct_regulated_nodes:
  //    constructs the set of nodes being regulated
  //--------------------------------------------------------
  void ThermostatSolverFixed::construct_regulated_nodes()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());
    regulatedNodes_ = interscaleManager.set_int(regulatorPrefix_+"ThermostatRegulatedNodes");

    if (!regulatedNodes_) {
      if (!atomicRegulator_->use_localized_lambda()) {
        regulatedNodes_ = new RegulatedNodes(atc_);
      }
      else if (atomicRegulator_->coupling_mode() == AtomicRegulator::FLUX) {
        regulatedNodes_ = new FixedNodes(atc_);
      }
      else if (atomicRegulator_->coupling_mode() == AtomicRegulator::FIXED) {
          // include boundary nodes
          regulatedNodes_ = new FixedBoundaryNodes(atc_);
      }
      else {
        throw ATC_Error("ThermostatSolverFixed::construct_regulated_nodes - couldn't determine set of regulated nodes");
      }
     
      interscaleManager.add_set_int(regulatedNodes_,
                                    regulatorPrefix_+"ThermostatRegulatedNodes");
    }

    applicationNodes_ = regulatedNodes_;

    // special set of boundary elements for defining regulated atoms 
    if (atomicRegulator_->use_localized_lambda()) {
      elementMask_ = interscaleManager.dense_matrix_bool(regulatorPrefix_+"BoundaryElementMask");
      if (!elementMask_) {
        elementMask_ = new ElementMaskNodeSet(atc_,applicationNodes_);
        interscaleManager.add_dense_matrix_bool(elementMask_,
                                                regulatorPrefix_+"BoundaryElementMask");
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatIntegratorFixed
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatIntegratorFixed::ThermostatIntegratorFixed(AtomicRegulator * thermostat,
                                                       int lambdaMaxIterations,
                                                       const string & regulatorPrefix) :
    ThermostatGlcFs(thermostat,lambdaMaxIterations,regulatorPrefix),
    atomThermostatForcesPredVel_(NULL),
    filterCoefficient_(1.)
  {
    lambdaSolver_ = new ThermostatSolverFixed(thermostat,
                                              lambdaMaxIterations,
                                              regulatorPrefix);
  }

  //--------------------------------------------------------
  //  construct_regulated_nodes:
  //    constructs the set of nodes being regulated
  //--------------------------------------------------------
  void ThermostatIntegratorFixed::construct_transfers()
  {
    ThermostatGlcFs::construct_transfers();

    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // predicted forces for halving update
    atomThermostatForcesPredVel_ = new AtomicThermostatForce(atc_,atomLambdas_,atomPredictedVelocities_);
    interscaleManager.add_per_atom_quantity(atomThermostatForcesPredVel_,
                                            regulatorPrefix_+"AtomThermostatForcePredictedVelocity");
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void ThermostatIntegratorFixed::initialize()
  {
    ThermostatGlcFs::initialize();
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // set KE multiplier
    AtomicEnergyForTemperature * atomEnergyForTemperature =
      static_cast<AtomicEnergyForTemperature * >(interscaleManager.per_atom_quantity("AtomicEnergyForTemperature"));
    keMultiplier_ = atomEnergyForTemperature->kinetic_energy_multiplier();

    // reset data to zero
    deltaFeEnergy_.reset(nNodes_,1);
    deltaNodalAtomicEnergy_.reset(nNodes_,1);

    // initialize filtered energy
    TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
    if (!timeFilterManager->end_equilibrate()) {
      nodalAtomicEnergyFiltered_ = nodalAtomicEnergy_->quantity();
    }

    // timestep factor
    lambdaSolver_->set_timestep_factor(0.5);
  }

  //--------------------------------------------------------
  //  halve_force:
  //    flag to halve the lambda force for improved
  //    accuracy
  //--------------------------------------------------------
  bool ThermostatIntegratorFixed::halve_force()
  {
    if (isFirstTimestep_ || ((atc_->atom_to_element_map_type() == EULERIAN)
                             && (atc_->atom_to_element_map_frequency() > 1)
                             && (atc_->step() % atc_->atom_to_element_map_frequency() == 1))) {
      return true;
    }
    return false;
  }

  //--------------------------------------------------------
  //  initialize_delta_nodal_atomic_energy:
  //    initializes storage for the variable tracking
  //    the change in the nodal atomic energy
  //    that has occured over the past timestep
  //--------------------------------------------------------
  void ThermostatIntegratorFixed::initialize_delta_nodal_atomic_energy(double dt)
  {
    // initialize delta energy
    const DENS_MAT & myNodalAtomicEnergy(nodalAtomicEnergy_->quantity());
    initialNodalAtomicEnergy_ = myNodalAtomicEnergy;
    initialNodalAtomicEnergy_ *= -1.; // initially stored as negative for efficiency
    timeFilter_->apply_pre_step1(nodalAtomicEnergyFiltered_.set_quantity(),
                                 myNodalAtomicEnergy,dt);
  }

  //--------------------------------------------------------
  //  compute_delta_nodal_atomic_energy:
  //    computes the change in the nodal atomic energy
  //    that has occured over the past timestep
  //--------------------------------------------------------
  void ThermostatIntegratorFixed::compute_delta_nodal_atomic_energy(double dt)
  {
    // set delta energy based on predicted atomic velocities
    const DENS_MAT & myNodalAtomicEnergy(nodalAtomicEnergy_->quantity());
    timeFilter_->apply_post_step1(nodalAtomicEnergyFiltered_.set_quantity(),
                                  myNodalAtomicEnergy,dt);
    deltaNodalAtomicEnergy_ = initialNodalAtomicEnergy_;
    deltaNodalAtomicEnergy_ += myNodalAtomicEnergy;
  }

  //--------------------------------------------------------
  //  compute_lambda:
  //   sets up and solves linear system for lambda
  //--------------------------------------------------------
  void ThermostatIntegratorFixed::compute_lambda(double dt,
                                                 bool iterateSolution)
  {
    // compute predicted changes in nodal atomic energy
    compute_delta_nodal_atomic_energy(dt);

    // change in finite element energy
    deltaFeEnergy_ = initialFeEnergy_;
    deltaFeEnergy_ += (mdMassMatrix_.quantity())*(temperature_.quantity());

    ThermostatGlcFs::compute_lambda(dt,iterateSolution);
  }

  //--------------------------------------------------------
  //  apply_predictor:
  //    apply the thermostat to the atoms in the first step
  //    of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatIntegratorFixed::apply_pre_predictor(double dt)
  {
    // initialize values to be track change in finite element energy over the timestep
    initialize_delta_nodal_atomic_energy(dt);
    initialFeEnergy_ = -1.*((mdMassMatrix_.quantity())*(temperature_.quantity())); // initially stored as negative for efficiency

    ThermostatGlcFs::apply_pre_predictor(dt);
  }

  //--------------------------------------------------------
  //  apply_pre_corrector:
  //    apply the thermostat to the atoms in the first part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatIntegratorFixed::apply_pre_corrector(double dt)
  {
    // do full prediction if we just redid the shape functions
    if (full_prediction()) {
      firstHalfAtomForces_ = atomThermostatForces_; // reset in case this time step needed special treatment
      _tempNodalAtomicEnergyFiltered_ = nodalAtomicEnergyFiltered_.quantity();
    }

    ThermostatGlcFs::apply_pre_corrector(dt);

    if (full_prediction()) {
      // reset temporary variables
      nodalAtomicEnergyFiltered_ = _tempNodalAtomicEnergyFiltered_;
    }

    if (halve_force()) {
      // save old velocities if we are doing halving calculation of lambda force
      // copy velocities over into temporary storage
      (*atomPredictedVelocities_) = atomVelocities_->quantity();
    }
  }

  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the thermostat to the atoms in the second part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatIntegratorFixed::apply_post_corrector(double dt)
  {
    
    bool halveForce = halve_force();

    ThermostatGlcFs::apply_post_corrector(dt);
    
    // update filtered energy with lambda power
    DENS_MAT & myNodalAtomicLambdaPower(nodalAtomicLambdaPower_->set_quantity());
    timeFilter_->apply_post_step2(nodalAtomicEnergyFiltered_.set_quantity(),
                                  myNodalAtomicLambdaPower,dt);

    if (halveForce) {
      // Halve lambda force due to fixed temperature constraints
      // 1) makes up for poor initial condition
      // 2) accounts for possibly large value of lambda when atomic shape function values change
      //    from eulerian mapping after more than 1 timestep
      //    avoids unstable oscillations arising from 
      //    thermostat having to correct for error introduced in lambda changing the 
      //    shape function matrices
      lambdaSolver_->scale_lambda(0.5);
      firstHalfAtomForces_ = atomThermostatForcesPredVel_;
      atomThermostatForcesPredVel_->unfix_quantity();
    }
    else {
      firstHalfAtomForces_ = atomThermostatForces_;
    }
  }

  //--------------------------------------------------------
  //  add_to_temperature
  //    add in contributions from lambda power and boundary
  //    flux to the FE temperature
  //--------------------------------------------------------
  void ThermostatIntegratorFixed::add_to_energy(const DENS_MAT & nodalLambdaPower,
                                                DENS_MAT & deltaEnergy,
                                                double dt)
  {
    deltaEnergy.resize(nNodes_,1);
    
    SetDependencyManager<int> * myRegulatedNodes =
      (atc_->interscale_manager()).set_int(regulatorPrefix_+"ThermostatRegulatedNodes");
    const set<int> & regulatedNodes(myRegulatedNodes->quantity());
    
    for (int i = 0; i < nNodes_; i++) {
      if (regulatedNodes.find(i) != regulatedNodes.end()) {
        deltaEnergy(i,0) = 0.;
      }
      else {
        deltaEnergy(i,0) = nodalLambdaPower(i,0);
      }
    }
  }

  //--------------------------------------------------------
  //  set_thermostat_rhs:
  //    sets up the right-hand side including boundary
  //    fluxes (coupling & prescribed), heat sources, and
  //    fixed (uncoupled) nodes
  //--------------------------------------------------------
  void ThermostatIntegratorFixed::set_thermostat_rhs(DENS_MAT & rhs,
                                                     double dt)
  {
    // for essential bcs (fixed nodes) :
    // form rhs : (delThetaV - delTheta)/dt
    SetDependencyManager<int> * myRegulatedNodes =
      (atc_->interscale_manager()).set_int(regulatorPrefix_+"ThermostatRegulatedNodes");
    const set<int> & regulatedNodes(myRegulatedNodes->quantity());
    double factor = (1./dt)/keMultiplier_;
    rhs.resize(nNodes_,1);
    
    for (int i = 0; i < nNodes_; i++) {
      if (regulatedNodes.find(i) != regulatedNodes.end()) {
        rhs(i,0) = factor*(deltaNodalAtomicEnergy_(i,0) - deltaFeEnergy_(i,0));
      }
      else {
        rhs(i,0) = 0.;
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatIntegratorFluxFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatIntegratorFluxFiltered::ThermostatIntegratorFluxFiltered(AtomicRegulator * thermostat,
                                                                     int lambdaMaxIterations,
                                                                     const string & regulatorPrefix) :
    ThermostatIntegratorFlux(thermostat,lambdaMaxIterations,regulatorPrefix)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void ThermostatIntegratorFluxFiltered::initialize()
  {
    ThermostatIntegratorFlux::initialize();

    TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
    if (!timeFilterManager->end_equilibrate()) {
      // always must start as zero because of filtering scheme
      heatSourceOld_.reset(nNodes_,1);
      instantHeatSource_.reset(nNodes_,1);
      timeStepSource_.reset(nNodes_,1);
    }
  }

  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the thermostat to the atoms in the second part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatIntegratorFluxFiltered::apply_post_corrector(double dt)
  {
    // compute lambda
    ThermostatIntegratorFlux::apply_post_corrector(dt);

    // store data needed for filter inversion of heat flux for thermostat rhs
    instantHeatSource_ = rhs_;
    heatSourceOld_ = heatSource_.quantity();
  }

  //--------------------------------------------------------
  //  add_to_temperature
  //    add in contributions from lambda power and boundary
  //    flux to the FE temperature
  //--------------------------------------------------------
  void ThermostatIntegratorFluxFiltered::add_to_energy(const DENS_MAT & nodalLambdaPower,
                                                       DENS_MAT & deltaEnergy,
                                                       double dt)
  {
    deltaEnergy.reset(nNodes_,1);
    double coef = timeFilter_->unfiltered_coefficient_post_s1(2.*dt);
    const DENS_MAT & myBoundaryFlux(boundaryFlux_[TEMPERATURE].quantity());
    for (int i = 0; i < nNodes_; i++) {
      deltaEnergy(i,0) = coef*nodalLambdaPower(i,0) + dt*myBoundaryFlux(i,0);
    }
  }

  //--------------------------------------------------------
  //  set_thermostat_rhs:
  //    sets up the right-hand side including boundary
  //    fluxes (coupling & prescribed), heat sources, and
  //    fixed (uncoupled) nodes
  //--------------------------------------------------------
  void ThermostatIntegratorFluxFiltered::set_thermostat_rhs(DENS_MAT & rhs,
                                                            double dt)
  {
    
    // only tested with flux != 0 + ess bc = 0

    // (a) for flux based : 
    // form rhs :  2/3kB * W_I^-1 * \int N_I r dV
    // vs  Wagner, CMAME, 2008 eq(24) RHS_I = 2/(3kB) flux_I
    // fluxes are set in ATC transfer

    // invert heatSource_ to get unfiltered source
    // relevant coefficients from time filter
    
    double coefF1 = timeFilter_->filtered_coefficient_pre_s1(2.*dt);
    double coefF2 = timeFilter_->filtered_coefficient_post_s1(2.*dt);
    double coefU1 = timeFilter_->unfiltered_coefficient_pre_s1(2.*dt);
    double coefU2 = timeFilter_->unfiltered_coefficient_post_s1(2.*dt);

    const DENS_MAT & heatSource(heatSource_.quantity());
    SetDependencyManager<int> * myApplicationNodes = 
      (atc_->interscale_manager()).set_int(regulatorPrefix_+"ThermostatApplicationNodes");
    const set<int> & applicationNodes(myApplicationNodes->quantity());
    rhs.resize(nNodes_,1);
    for (int i = 0; i < nNodes_; i++) {
      if (applicationNodes.find(i) != applicationNodes.end()) {
        rhs(i,0) = heatSource(i,0) - coefF1*coefF2*heatSourceOld_(i,0) - coefU1*coefF2*instantHeatSource_(i,0);
        rhs(i,0) /= coefU2;
      }
      else {
        rhs(i,0) = 0.;
      }
    }
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void ThermostatIntegratorFluxFiltered::output(OUTPUT_LIST & outputData)
  {
    _lambdaPowerOutput_ = lambdaPowerFiltered_->quantity();
    // approximate value for lambda power
    double dt =  LammpsInterface::instance()->dt();
    _lambdaPowerOutput_ *= (2./dt);
    DENS_MAT & lambda((atomicRegulator_->regulator_data(regulatorPrefix_+"LambdaEnergy",1))->set_quantity());
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData[regulatorPrefix_+"Lambda"] = &lambda;
      outputData[regulatorPrefix_+"NodalLambdaPower"] = &(_lambdaPowerOutput_);
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatIntegratorFixedFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatIntegratorFixedFiltered::ThermostatIntegratorFixedFiltered(AtomicRegulator * thermostat,
                                                                       int lambdaMaxIterations,
                                                                       const string & regulatorPrefix) :
    ThermostatIntegratorFixed(thermostat,lambdaMaxIterations,regulatorPrefix)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  initialize_delta_nodal_atomic_energy:
  //    initializes storage for the variable tracking
  //    the change in the nodal atomic energy
  //    that has occured over the past timestep
  //--------------------------------------------------------
  
  
  void ThermostatIntegratorFixedFiltered::initialize_delta_nodal_atomic_energy(double dt)
  {
    // initialize delta energy
    DENS_MAT & myNodalAtomicEnergyFiltered(nodalAtomicEnergyFiltered_.set_quantity());
    initialNodalAtomicEnergy_ = myNodalAtomicEnergyFiltered;
    initialNodalAtomicEnergy_ *= -1.; // initially stored as negative for efficiency
    timeFilter_->apply_pre_step1(myNodalAtomicEnergyFiltered,
                                 nodalAtomicEnergy_->quantity(),dt);
  }

  //--------------------------------------------------------
  //  compute_delta_nodal_atomic_energy:
  //    computes the change in the nodal atomic energy
  //    that has occured over the past timestep
  //--------------------------------------------------------
  void ThermostatIntegratorFixedFiltered::compute_delta_nodal_atomic_energy(double dt)
  {
    // set delta energy based on predicted atomic velocities
    DENS_MAT & myNodalAtomicEnergyFiltered(nodalAtomicEnergyFiltered_.set_quantity());
    timeFilter_->apply_post_step1(myNodalAtomicEnergyFiltered,
                                  nodalAtomicEnergy_->quantity(),dt);
    deltaNodalAtomicEnergy_ = initialNodalAtomicEnergy_;
    deltaNodalAtomicEnergy_ += myNodalAtomicEnergyFiltered;
  }

  //--------------------------------------------------------
  //  add_to_temperature
  //    add in contributions from lambda power and boundary
  //    flux to the FE temperature
  //--------------------------------------------------------
  void ThermostatIntegratorFixedFiltered::add_to_energy(const DENS_MAT & nodalLambdaPower,
                                                        DENS_MAT & deltaEnergy,
                                                        double dt)
  {
    deltaEnergy.resize(nNodes_,1);
    SetDependencyManager<int> * myRegulatedNodes =
      (atc_->interscale_manager()).set_int(regulatorPrefix_+"ThermostatRegulatedNodes");
    const set<int> & regulatedNodes(myRegulatedNodes->quantity());
    double coef = timeFilter_->unfiltered_coefficient_post_s1(2.*dt);
    for (int i = 0; i < nNodes_; i++) {
      if (regulatedNodes.find(i) != regulatedNodes.end()) {
        deltaEnergy(i,0) = 0.;
      }
      else {
        deltaEnergy(i,0) = coef*nodalLambdaPower(i,0);
      }
    }
  }
  //--------------------------------------------------------
  //  set_thermostat_rhs:
  //    sets up the right-hand side for fixed
  //    (coupling & prescribed) temperature values
  //--------------------------------------------------------
  void ThermostatIntegratorFixedFiltered::set_thermostat_rhs(DENS_MAT & rhs,
                                                             double dt)
  {
    // (b) for essential bcs (fixed nodes):
    // form rhs : (delThetaV - delTheta)/dt
    SetDependencyManager<int> * myRegulatedNodes =
      (atc_->interscale_manager()).set_int(regulatorPrefix_+"ThermostatRegulatedNodes");
    const set<int> & regulatedNodes(myRegulatedNodes->quantity());
    double factor = (1./dt)/keMultiplier_;
    factor /= timeFilter_->unfiltered_coefficient_post_s1(2.*dt);
    rhs.resize(nNodes_,1);
    
    for (int i = 0; i < nNodes_; i++) {
      if (regulatedNodes.find(i) != regulatedNodes.end()) {
        rhs(i,0) = factor*(deltaNodalAtomicEnergy_(i,0) - deltaFeEnergy_(i,0));
      }
      else {
        rhs(i,0) = 0.;
      }
    }
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void ThermostatIntegratorFixedFiltered::output(OUTPUT_LIST & outputData)
  {
     _lambdaPowerOutput_ = lambdaPowerFiltered_->quantity();
    // approximate value for lambda power
    double dt =  LammpsInterface::instance()->dt();
    _lambdaPowerOutput_ *= (2./dt);
    DENS_MAT & lambda((atomicRegulator_->regulator_data(regulatorPrefix_+"LambdaEnergy",1))->set_quantity());
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData[regulatorPrefix_+"Lambda"] = &lambda;
      outputData[regulatorPrefix_+"NodalLambdaPower"] = &(_lambdaPowerOutput_);
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatFluxFixed
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ThermostatFluxFixed::ThermostatFluxFixed(AtomicRegulator * thermostat,
                                           int lambdaMaxIterations,
                                           bool constructThermostats) :
    RegulatorMethod(thermostat),
    thermostatFlux_(NULL),
    thermostatFixed_(NULL),
    thermostatBcs_(NULL)
  {
    if (constructThermostats) {
      thermostatFlux_ = new ThermostatIntegratorFlux(thermostat,lambdaMaxIterations,regulatorPrefix_+"Flux");
      thermostatFixed_ = new ThermostatIntegratorFixed(thermostat,lambdaMaxIterations,regulatorPrefix_+"Fixed");

      // need to choose BC type based on coupling mode
      if (thermostat->coupling_mode() == AtomicRegulator::FLUX) {
        thermostatBcs_ = thermostatFlux_;
      }
      else if (thermostat->coupling_mode() == AtomicRegulator::FIXED) {
        thermostatBcs_ = thermostatFixed_;
      }
      else {
        throw ATC_Error("ThermostatFluxFixed:create_thermostats - invalid thermostat type provided");
      }
    }
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ThermostatFluxFixed::~ThermostatFluxFixed()
  {
    if (thermostatFlux_) delete thermostatFlux_;
    if (thermostatFixed_) delete thermostatFixed_;
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void ThermostatFluxFixed::construct_transfers()
  {
    thermostatFlux_->construct_transfers();
    thermostatFixed_->construct_transfers();
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void ThermostatFluxFixed::initialize()
  {
    thermostatFixed_->initialize();
    thermostatFlux_->initialize();
  }

  //--------------------------------------------------------
  //  apply_predictor:
  //    apply the thermostat to the atoms in the first step
  //    of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatFluxFixed::apply_pre_predictor(double dt)
  {
    thermostatFixed_->apply_pre_predictor(dt);
    thermostatFlux_->apply_pre_predictor(dt);
  }

  //--------------------------------------------------------
  //  apply_pre_corrector:
  //    apply the thermostat to the atoms in the first part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatFluxFixed::apply_pre_corrector(double dt)
  {
    thermostatFlux_->apply_pre_corrector(dt);
    if (thermostatFixed_->full_prediction()) {
      atc_->set_fixed_nodes();
    }
    thermostatFixed_->apply_pre_corrector(dt);
  }

  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the thermostat to the atoms in the second part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatFluxFixed::apply_post_corrector(double dt)
  {
    thermostatFlux_->apply_post_corrector(dt);
    atc_->set_fixed_nodes();
    thermostatFixed_->apply_post_corrector(dt);
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void ThermostatFluxFixed::output(OUTPUT_LIST & outputData)
  {
    thermostatFlux_->output(outputData);
    thermostatFixed_->output(outputData);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatFluxFixedFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ThermostatFluxFixedFiltered::ThermostatFluxFixedFiltered(AtomicRegulator * thermostat,
                                                           int lambdaMaxIterations) :
    ThermostatFluxFixed(thermostat,lambdaMaxIterations,false)
  {
    thermostatFlux_ = new ThermostatIntegratorFluxFiltered(thermostat,lambdaMaxIterations,regulatorPrefix_+"Flux");
    thermostatFixed_ = new ThermostatIntegratorFixedFiltered(thermostat,lambdaMaxIterations,regulatorPrefix_+"Fixed");

    // need to choose BC type based on coupling mode
    if (thermostat->coupling_mode() == AtomicRegulator::FLUX) {
      thermostatBcs_ = thermostatFlux_;
    }
    else if (thermostat->coupling_mode() == AtomicRegulator::FIXED) {
      thermostatBcs_ = thermostatFixed_;
    }
    else {
      throw ATC_Error("ThermostatFluxFixed:create_thermostats - invalid thermostat type provided");
    }
  }
  //--------------------------------------------------------
  //  Class ThermostatGlc
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ThermostatGlc::ThermostatGlc(AtomicRegulator * thermostat) :
    ThermostatShapeFunction(thermostat),
    timeFilter_(atomicRegulator_->time_filter()),
    lambdaPowerFiltered_(NULL),
    atomThermostatForces_(NULL),
    prescribedDataMgr_(atc_->prescribed_data_manager()),
    atomMasses_(NULL)
  {
    // consistent with stage 3 of ATC_Method::initialize
    lambdaPowerFiltered_= atomicRegulator_->regulator_data(regulatorPrefix_+"LambdaPowerFiltered",1);
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void ThermostatGlc::construct_transfers()
  {
    ThermostatShapeFunction::construct_transfers();
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // get data from manager
    atomMasses_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS);

    // thermostat forces based on lambda and the atomic velocities
    AtomicThermostatForce * atomThermostatForces = new AtomicThermostatForce(atc_);
    interscaleManager.add_per_atom_quantity(atomThermostatForces,
                                            regulatorPrefix_+"AtomThermostatForce");
    atomThermostatForces_ = atomThermostatForces;
  }

  //--------------------------------------------------------
  //  apply_to_atoms:
  //            determines what if any contributions to the
  //            atomic moition is needed for
  //            consistency with the thermostat
  //--------------------------------------------------------
  void ThermostatGlc::apply_to_atoms(PerAtomQuantity<double> * atomVelocities,
                                     const DENS_MAT & lambdaForce,
                                     double dt)
  {
    _velocityDelta_ = lambdaForce;
    _velocityDelta_ /= atomMasses_->quantity();
    _velocityDelta_ *= dt;
    (*atomVelocities) += _velocityDelta_;
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
  ThermostatPowerVerlet::ThermostatPowerVerlet(AtomicRegulator * thermostat) :
    ThermostatGlc(thermostat),
    nodalTemperatureRoc_(atc_->field_roc(TEMPERATURE)),
    heatSource_(atc_->atomic_source(TEMPERATURE)),
    nodalAtomicPower_(NULL),
    nodalAtomicLambdaPower_(NULL)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void ThermostatPowerVerlet::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // set up node mappings
    create_node_maps();

    // determine if mapping is needed and set up if so
    if (atomicRegulator_->use_localized_lambda()) {
        lambdaAtomMap_ = new AtomToElementset(atc_,elementMask_);
        interscaleManager.add_per_atom_int_quantity(lambdaAtomMap_,
                                                    regulatorPrefix_+"LambdaAtomMap");
    }

    // set up linear solver
    if (atomicRegulator_->use_lumped_lambda_solve()) {
      shapeFunctionMatrix_ = new LambdaCouplingMatrix(atc_,nodeToOverlapMap_);
      linearSolverType_ = AtomicRegulator::RSL_SOLVE;
    }
    else {
      if (lambdaAtomMap_) {
        shapeFunctionMatrix_ = new LocalLambdaCouplingMatrix(atc_,
                                                             lambdaAtomMap_,
                                                             nodeToOverlapMap_);
      }
      else {
        shapeFunctionMatrix_ = new LambdaCouplingMatrix(atc_,nodeToOverlapMap_);
      }
      linearSolverType_ = AtomicRegulator::CG_SOLVE;
    }
    interscaleManager.add_per_atom_sparse_matrix(shapeFunctionMatrix_,
                                                 regulatorPrefix_+"LambdaCouplingMatrixEnergy");

    // base class transfers, e.g. weights
    ThermostatGlc::construct_transfers();

    // get managed data
    nodalAtomicPower_ = interscaleManager.dense_matrix("NodalAtomicPower");
    
    // power induced by lambda
    DotTwiceKineticEnergy * atomicLambdaPower = 
      new DotTwiceKineticEnergy(atc_,atomThermostatForces_);
    interscaleManager.add_per_atom_quantity(atomicLambdaPower,
                                            regulatorPrefix_+"AtomicLambdaPower");
    
    // restriction to nodes of power induced by lambda
    nodalAtomicLambdaPower_ = new AtfShapeFunctionRestriction(atc_,
                                                              atomicLambdaPower,
                                                              interscaleManager.per_atom_sparse_matrix("Interpolant"));
    interscaleManager.add_dense_matrix(nodalAtomicLambdaPower_,
                                            regulatorPrefix_+"NodalAtomicLambdaPower");
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void ThermostatPowerVerlet::initialize()
  {
    ThermostatGlc::initialize();

    // sets up time filter for cases where variables temporally filtered
    TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
    if (!timeFilterManager->end_equilibrate()) {
      _nodalAtomicLambdaPowerOut_ = 0.;
      *lambdaPowerFiltered_ = 0.;
      timeFilter_->initialize(lambdaPowerFiltered_->quantity());
    }
  }

  //--------------------------------------------------------
  //  apply_predictor:
  //    apply the thermostat to the atoms in the first step
  //    of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatPowerVerlet::apply_pre_predictor(double dt)
  {
    atomThermostatForces_->unfix_quantity();
    compute_thermostat(0.5*dt);

    // apply lambda force to atoms
    const DENS_MAT & thermostatForces(atomThermostatForces_->quantity());
    atomThermostatForces_->fix_quantity();
    apply_to_atoms(atomVelocities_,thermostatForces,0.5*dt);
  }

  //--------------------------------------------------------
  //  apply_pre_corrector:
  //    apply the thermostat to the atoms in the first part
  //    of the corrector step of the Verlet algorithm
  //--------------------------------------------------------
  void ThermostatPowerVerlet::apply_pre_corrector(double dt)
  {
    atomThermostatForces_->unfix_quantity();
    compute_thermostat(0.5*dt);

    // apply lambda force to atoms
    const DENS_MAT & thermostatForces(atomThermostatForces_->quantity());
    atomThermostatForces_->fix_quantity();
    apply_to_atoms(atomVelocities_,thermostatForces,0.5*dt);
  }

  //--------------------------------------------------------
  //  add_to_rhs:
  //            determines what if any contributions to the
  //            finite element equations are needed for
  //            consistency with the thermostat
  //--------------------------------------------------------
  void ThermostatPowerVerlet::add_to_rhs(FIELDS & rhs)
  {
    rhs[TEMPERATURE] += nodalAtomicLambdaPower_->quantity() + boundaryFlux_[TEMPERATURE].quantity();
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
    rhs_nodes = heatSource_.quantity();
    
    // (b) for ess. bcs
    // form rhs : {sum_a (2 * N_Ia * v_ia * f_ia) - (dtheta/dt)_I}
    
    // replace rhs for prescribed nodes
    const DENS_MAT & myNodalAtomicPower(nodalAtomicPower_->quantity());
    const DIAG_MAT & myMdMassMatrix(mdMassMatrix_.quantity());
    const DENS_MAT & myNodalTemperatureRoc(nodalTemperatureRoc_.quantity());
    for (int i = 0; i < nNodes_; i++) {
      if (prescribedDataMgr_->is_fixed(i,TEMPERATURE,0)) {
        rhs_nodes(i,0) = 0.5*(myNodalAtomicPower(i,0) - myMdMassMatrix(i,i)*myNodalTemperatureRoc(i,0));
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
    set_thermostat_rhs(_rhs_);
    
    // solve linear system for lambda
    DENS_MAT & myLambda(lambda_->set_quantity());
    solve_for_lambda(_rhs_,myLambda);

    nodalAtomicLambdaPower_->unfix_quantity(); // enable computation of force applied by lambda
    timeFilter_->apply_pre_step1(lambdaPowerFiltered_->set_quantity(),
                                 nodalAtomicLambdaPower_->quantity(),dt);
    nodalAtomicLambdaPower_->fix_quantity();
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void ThermostatPowerVerlet::output(OUTPUT_LIST & outputData)
  {
    _nodalAtomicLambdaPowerOut_ = nodalAtomicLambdaPower_->quantity();
    DENS_MAT & lambda(lambda_->set_quantity());
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData["lambda"] = &lambda;
      outputData["nodalLambdaPower"] = &(_nodalAtomicLambdaPowerOut_);
    }
  }

  //--------------------------------------------------------
  //  finish:
  //    final tasks after a run
  //--------------------------------------------------------
  void ThermostatPowerVerlet::finish()
  {
    _nodalAtomicLambdaPowerOut_ = nodalAtomicLambdaPower_->quantity();
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
  ThermostatHooverVerlet::ThermostatHooverVerlet(AtomicRegulator * thermostat) :
    ThermostatPowerVerlet(thermostat),
    lambdaHoover_(NULL),
    nodalAtomicHooverLambdaPower_(NULL)
  {
    // set up data consistent with stage 3 of ATC_Method::initialize
    lambdaHoover_ = atomicRegulator_->regulator_data(regulatorPrefix_+"LambdaHoover",1);
    
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void ThermostatHooverVerlet::construct_transfers()
  {
    ThermostatPowerVerlet::construct_transfers();
    InterscaleManager & interscaleManager(atc_->interscale_manager());


    FtaShapeFunctionProlongation * atomHooverLambdas = new FtaShapeFunctionProlongation(atc_,
                                                                                        lambdaHoover_,
                                                                                        interscaleManager.per_atom_sparse_matrix("Interpolant"));
    interscaleManager.add_per_atom_quantity(atomHooverLambdas,
                                            regulatorPrefix_+"AtomHooverLambda");
    AtomicThermostatForce * atomHooverThermostatForces = new AtomicThermostatForce(atc_,atomHooverLambdas);
    interscaleManager.add_per_atom_quantity(atomHooverThermostatForces,
                                            regulatorPrefix_+"AtomHooverThermostatForce");
    SummedAtomicQuantity<double> * atomTotalThermostatForces = 
      new SummedAtomicQuantity<double>(atc_,atomThermostatForces_,atomHooverThermostatForces);
    interscaleManager.add_per_atom_quantity(atomTotalThermostatForces,
                                            regulatorPrefix_+"AtomTotalThermostatForce");
    atomThermostatForces_ = atomTotalThermostatForces;
          
    // transfers dependent on time integration method
    DotTwiceKineticEnergy * atomicHooverLambdaPower = 
      new DotTwiceKineticEnergy(atc_,atomHooverThermostatForces);
    interscaleManager.add_per_atom_quantity(atomicHooverLambdaPower,
                                            regulatorPrefix_+"AtomicHooverLambdaPower");
    
     nodalAtomicHooverLambdaPower_ = new AtfShapeFunctionRestriction(atc_,
                                                                     atomicHooverLambdaPower,
                                                                     interscaleManager.per_atom_sparse_matrix("Interpolant"));
     interscaleManager.add_dense_matrix(nodalAtomicHooverLambdaPower_,
                                             regulatorPrefix_+"NodalAtomicHooverLambdaPower");
  }

  //--------------------------------------------------------
  //  add_to_rhs:
  //            determines what if any contributions to the
  //            finite element equations are needed for
  //            consistency with the thermostat
  //--------------------------------------------------------
  void ThermostatHooverVerlet::add_to_rhs(FIELDS & rhs)
  {
    rhs[TEMPERATURE] += _nodalAtomicLambdaPowerOut_;
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
    ThermostatPowerVerlet::compute_thermostat(0.5*dt);
    _nodalAtomicLambdaPowerOut_ = nodalAtomicLambdaPower_->quantity(); // save power from lambda in power-based thermostat
    
    // set up Hoover rhs
    set_hoover_rhs(_rhs_);
    
    // solve linear system for lambda
    DENS_MAT & myLambda(lambdaHoover_->set_quantity());
    solve_for_lambda(_rhs_,myLambda);
    
    // compute force applied by lambda
    // compute nodal atomic power from Hoover coupling
    // only add in contribution to uncoupled nodes
    if (atomicRegulator_->use_localized_lambda())
      add_to_lambda_power(atomThermostatForces_->quantity(),0.5*dt);
  }

  //--------------------------------------------------------
  //  set_hoover_rhs:
  //    sets up the right-hand side for fixed value,
  //    i.e. Hoover coupling
  //--------------------------------------------------------
  void ThermostatHooverVerlet::set_hoover_rhs(DENS_MAT & rhs)
  {
    // form rhs : sum_a ( N_Ia * v_ia * f_ia) - 0.5*M_MD*(dtheta/dt)_I
    rhs = nodalAtomicPower_->quantity();
    rhs -= mdMassMatrix_.quantity()*nodalTemperatureRoc_.quantity();
    rhs /= 2.;
  }

  //--------------------------------------------------------
  //  add_to_nodal_lambda_power:
  //    determines the power exerted by the Hoover 
  //    thermostat at each FE node
  //--------------------------------------------------------
  void ThermostatHooverVerlet::add_to_lambda_power(const DENS_MAT & myLambdaForce,
                                                   double dt)
  {
    _myNodalLambdaPower_ = nodalAtomicHooverLambdaPower_->quantity();
    const INT_ARRAY & nodeToOverlapMap(nodeToOverlapMap_->quantity());
    for (int i = 0; i < nNodes_; ++i) {
      if (nodeToOverlapMap(i,0)==-1)
        _nodalAtomicLambdaPowerOut_(i,0) += _myNodalLambdaPower_(i,0);
      else
        _myNodalLambdaPower_(i,0) = 0.;
    }
    timeFilter_->apply_post_step1(lambdaPowerFiltered_->set_quantity(),_myNodalLambdaPower_,dt);
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
  ThermostatPowerVerletFiltered::ThermostatPowerVerletFiltered(AtomicRegulator * thermostat) :
    ThermostatPowerVerlet(thermostat),
    nodalTemperature2Roc_(atc_->field_2roc(TEMPERATURE)),
    fieldsRoc_(atc_->fields_roc()),
    filterScale_((atc_->time_filter_manager())->filter_scale())
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
    atc_->compute_boundary_flux(fieldMask_,
                                fieldsRoc_,
                                fluxRoc_,
                                atomMaterialGroups_,
                                shpFcnDerivs_);

    

    // compute extrinsic model rate of change
    (atc_->extrinsic_model_manager()).set_sources(fieldsRoc_,fluxRoc_);
    heatSourceRoc_ = fluxRoc_[TEMPERATURE].quantity(); 
  }

  //--------------------------------------------------------
  //  add_to_rhs:
  //            determines what if any contributions to the
  //            finite element equations are needed for
  //            consistency with the thermostat
  //--------------------------------------------------------
  void ThermostatPowerVerletFiltered::add_to_rhs(FIELDS & rhs)
  {
    rhs[TEMPERATURE] += lambdaPowerFiltered_->quantity() + boundaryFlux_[TEMPERATURE].quantity();
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
    rhs_nodes = heatSource_.quantity() + filterScale_*heatSourceRoc_.quantity();

    // (b) for ess. bcs
    // form rhs : {sum_a (N_Ia * v_ia * f_ia) - 0.5*(dtheta/dt)_I}
    const DENS_MAT & myNodalAtomicPower(nodalAtomicPower_->quantity());
    const DIAG_MAT & myMdMassMatrix(mdMassMatrix_.quantity());
    const DENS_MAT & myNodalTemperatureRoc(nodalTemperatureRoc_.quantity());
    const DENS_MAT & myNodalTemperature2Roc(nodalTemperature2Roc_.quantity());
    for (int i = 0; i < nNodes_; i++) {
      if (prescribedDataMgr_->is_fixed(i,TEMPERATURE,0)) {
        rhs_nodes(i,0) = 0.5*(myNodalAtomicPower(i,0) - myMdMassMatrix(i,i)*(myNodalTemperatureRoc(i,0)+ filterScale_*myNodalTemperature2Roc(i,0)));
      }
    }
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void ThermostatPowerVerletFiltered::output(OUTPUT_LIST & outputData)
  {
    outputData["lambda"] = &(lambda_->set_quantity());
    outputData["nodalLambdaPower"] = &(lambdaPowerFiltered_->set_quantity());
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatHooverVerletFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and thermostat data
  //--------------------------------------------------------
  ThermostatHooverVerletFiltered::ThermostatHooverVerletFiltered(AtomicRegulator * thermostat) :
    ThermostatPowerVerletFiltered(thermostat),
    lambdaHoover_(NULL),
    nodalAtomicHooverLambdaPower_(NULL)
  {
    // consistent with stage 3 of ATC_Method::initialize
    lambdaHoover_ = atomicRegulator_->regulator_data("LambdaHoover",1);
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void ThermostatHooverVerletFiltered::construct_transfers()
  {
    ThermostatPowerVerletFiltered::construct_transfers();
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    FtaShapeFunctionProlongation * atomHooverLambdas = new FtaShapeFunctionProlongation(atc_,
                                                                                        lambdaHoover_,
                                                                                        interscaleManager.per_atom_sparse_matrix("Interpolant"));
    interscaleManager.add_per_atom_quantity(atomHooverLambdas,
                                            regulatorPrefix_+"AtomHooverLambda");
    AtomicThermostatForce * atomHooverThermostatForces = new AtomicThermostatForce(atc_,atomHooverLambdas);
    interscaleManager.add_per_atom_quantity(atomHooverThermostatForces,
                                            regulatorPrefix_+"AtomHooverThermostatForce");
    SummedAtomicQuantity<double> * atomTotalThermostatForces = 
      new SummedAtomicQuantity<double>(atc_,atomThermostatForces_,atomHooverThermostatForces);
    interscaleManager.add_per_atom_quantity(atomTotalThermostatForces,
                                            regulatorPrefix_+"AtomTotalThermostatForce");
    atomThermostatForces_ = atomTotalThermostatForces;
          
    // transfers dependent on time integration method
    DotTwiceKineticEnergy * atomicHooverLambdaPower = 
      new DotTwiceKineticEnergy(atc_,atomHooverThermostatForces);
    interscaleManager.add_per_atom_quantity(atomicHooverLambdaPower,
                                            regulatorPrefix_+"AtomicHooverLambdaPower");
    
     nodalAtomicHooverLambdaPower_ = new AtfShapeFunctionRestriction(atc_,
                                                                     atomicHooverLambdaPower,
                                                                     interscaleManager.per_atom_sparse_matrix("Interpolant"));
     interscaleManager.add_dense_matrix(nodalAtomicHooverLambdaPower_,
                                             regulatorPrefix_+"NodalAtomicHooverLambdaPower");
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
    ThermostatPowerVerletFiltered::compute_thermostat(0.5*dt);
    _nodalAtomicLambdaPowerOut_ = nodalAtomicLambdaPower_->quantity(); // save power from lambda in power-based thermostat
    
    // set up Hoover rhs
    set_hoover_rhs(_rhs_);
    
    // solve linear system for lambda
    DENS_MAT & myLambda(lambdaHoover_->set_quantity());
    solve_for_lambda(_rhs_,myLambda);
    
    
    // compute force applied by lambda
    // compute nodal atomic power from Hoover coupling
    // only add in contribution to uncoupled nodes
    if (atomicRegulator_->use_localized_lambda())
      add_to_lambda_power(atomThermostatForces_->quantity(),0.5*dt);
  }

  //--------------------------------------------------------
  //  set_hoover_rhs:
  //    sets up the right-hand side for fixed value,
  //    i.e. Hoover coupling
  //--------------------------------------------------------
  void ThermostatHooverVerletFiltered::set_hoover_rhs(DENS_MAT & rhs)
  {
    // form rhs : sum_a (N_Ia * v_ia * f_ia) - 0.5*M_MD*(dtheta/dt)_I
    rhs = nodalAtomicPower_->quantity();
    rhs -= mdMassMatrix_.quantity()*(nodalTemperatureRoc_.quantity() + filterScale_*nodalTemperature2Roc_.quantity());
    rhs /= 2.;
  }

  //--------------------------------------------------------
  //  add_to_nodal_lambda_power:
  //    determines the power exerted by the Hoover 
  //    thermostat at each FE node
  //--------------------------------------------------------
  void ThermostatHooverVerletFiltered::add_to_lambda_power(const DENS_MAT & myLambdaForce,
                                                           double dt)
  {
    _myNodalLambdaPower_ = nodalAtomicHooverLambdaPower_->quantity();
    const INT_ARRAY nodeToOverlapMap(nodeToOverlapMap_->quantity());
    for (int i = 0; i < nNodes_; ++i) {
      if (nodeToOverlapMap(i,0)==-1)
        _nodalAtomicLambdaPowerOut_(i,0) += _myNodalLambdaPower_(i,0);
      else
        _myNodalLambdaPower_(i,0) = 0.;
    }
    timeFilter_->apply_post_step1(lambdaPowerFiltered_->set_quantity(),_myNodalLambdaPower_,dt);
  }

  //--------------------------------------------------------
  //  add_to_rhs:
  //            determines what if any contributions to the
  //            finite element equations are needed for
  //            consistency with the thermostat
  //--------------------------------------------------------
  void ThermostatHooverVerletFiltered::add_to_rhs(FIELDS & rhs)
  {
    rhs[TEMPERATURE] += lambdaPowerFiltered_->quantity();
  }

};
