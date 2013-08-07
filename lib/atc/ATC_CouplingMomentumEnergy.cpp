// ATC headers
#include "ATC_CouplingMomentumEnergy.h"
#include "KinetoThermostat.h"
#include "ATC_Error.h"
#include "PrescribedDataManager.h"

// Other Headers
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <typeinfo>

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ATC_CouplingMomentumEnergy
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ATC_CouplingMomentumEnergy::ATC_CouplingMomentumEnergy(string groupName,
                                                         double ** & perAtomArray,
                                                         LAMMPS_NS::Fix * thisFix,
                                                         string matParamFile,
                                                         ExtrinsicModelType extrinsicModel)
    : ATC_Coupling(groupName,perAtomArray,thisFix),
      nodalAtomicMass_(NULL),
      nodalAtomicCount_(NULL),
      nodalAtomicHeatCapacity_(NULL),
      nodalAtomicKineticTemperature_(NULL),
      nodalAtomicConfigurationalTemperature_(NULL),
      boundaryDynamics_(PRESCRIBED),
      gamma_(0),mu_(1),kappa_(1),
      refPE_(0)
  {
    // Allocate PhysicsModel 
    create_physics_model(THERMO_ELASTIC, matParamFile);

    // create extrinsic physics model
    if (extrinsicModel != NO_MODEL) {
      extrinsicModelManager_.create_model(extrinsicModel,matParamFile);  
    }

    // Defaults
    set_time();
    bndyIntType_ = FE_INTERPOLATION;
    trackDisplacement_ = true;
  
    // set up field data based on physicsModel
    physicsModel_->num_fields(fieldSizes_,fieldMask_);
    fieldSizes_[DISPLACEMENT] = fieldSizes_[VELOCITY];

    // set up atomic regulator
    atomicRegulator_ = new KinetoThermostat(this);

    // default to not track charge
    trackCharge_ = false;

    // set up physics specific time integrator and thermostat
    timeIntegrators_[VELOCITY] = new MomentumTimeIntegrator(this,TimeIntegrator::FRACTIONAL_STEP);
    timeIntegrators_[TEMPERATURE] = new ThermalTimeIntegrator(this,TimeIntegrator::FRACTIONAL_STEP);

    // default physics
    temperatureDef_ = KINETIC;

    // output variable vector info:
    // output[1] = total coarse scale mechanical kinetic energy
    // output[2] = total coarse scale mechanical potential energy
    // output[3] = total coarse scale mechanical energy
    // output[1] = total coarse scale thermal energy
    // output[2] = average temperature
    scalarFlag_ = 1;
    vectorFlag_ = 1;
    sizeVector_ = 5;
    scalarVectorFreq_ = 1;
    extVector_ = 1;
    if (extrinsicModel != NO_MODEL)
      sizeVector_ += extrinsicModelManager_.size_vector(sizeVector_);

    // create PE per atom ccompute
    lammpsInterface_->create_compute_pe_peratom();
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ATC_CouplingMomentumEnergy::~ATC_CouplingMomentumEnergy()
  {
    // clear out all managed memory to avoid conflicts with dependencies on class member data
    interscaleManager_.clear();
  }

  //--------------------------------------------------------
  //  initialize
  //    sets up all the necessary data
  //--------------------------------------------------------
  void ATC_CouplingMomentumEnergy::initialize()
  {
    // clear displacement entries if requested
    if (!trackDisplacement_) {
      fieldSizes_.erase(DISPLACEMENT);
      for (int i = 0; i < NUM_FLUX; i++)
        fieldMask_(DISPLACEMENT,i) = false;
    }

    // Base class initalizations
    ATC_Coupling::initialize();
    
    // resetting precedence:
    // time integrator -> kinetostat/thermostat -> time filter
    // init_filter uses fieldRateNdFiltered which comes from the time integrator,
    // which is why the time integrator is initialized first

    // set the reference potential, if necessary, because the nodal energy is needed to initialize the time integrator
    if (!initialized_) {
      if (temperatureDef_==TOTAL) {
        PerAtomQuantity<double> * atomicReferencePotential = interscaleManager_.per_atom_quantity("AtomicReferencePotential");
        PerAtomQuantity<double> * atomicPotentialEnergy = interscaleManager_.per_atom_quantity("AtomicPotentialEnergy");
        atomicReferencePotential->set_quantity() = atomicPotentialEnergy->quantity();
      }
    }

    // other initializations
    
    if (reset_methods()) {
      for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
        (_tiIt_->second)->force_reset();
      }
      atomicRegulator_->force_reset();
    }
    if (reset_methods()) {
      for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
        (_tiIt_->second)->initialize();
      }
      atomicRegulator_->initialize();
      extrinsicModelManager_.initialize();
    }
    if (timeFilterManager_.need_reset()) // reset thermostat power
      init_filter();
    timeFilterManager_.initialize(); // clears need for reset

    if (!initialized_) {
      // initialize sources based on initial FE temperature
      double dt = lammpsInterface_->dt();
      prescribedDataMgr_->set_sources(time()+0.5*dt,sources_);
      extrinsicModelManager_.set_sources(fields_,extrinsicSources_);
      atomicRegulator_->compute_boundary_flux(fields_);
      compute_atomic_sources(fieldMask_,fields_,atomicSources_);

      // read in field data if necessary
      if (useRestart_) {
        RESTART_LIST data;
        read_restart_data(restartFileName_,data);
        useRestart_ = false;
      }

      // set consistent initial conditions, if requested
      if (!timeFilterManager_.filter_dynamics()) {
        if (consistentInitialization_) {
          
          DENS_MAT & velocity(fields_[VELOCITY].set_quantity());
          DENS_MAN * nodalAtomicVelocity(interscaleManager_.dense_matrix("NodalAtomicVelocity"));
          const DENS_MAT & atomicVelocity(nodalAtomicVelocity->quantity());
          DENS_MAT & temperature(fields_[TEMPERATURE].set_quantity());
          DENS_MAN * nodalAtomicTemperature(interscaleManager_.dense_matrix("NodalAtomicTemperature"));
          const DENS_MAT & atomicTemperature(nodalAtomicTemperature->quantity());
          const INT_ARRAY & nodeType(nodalGeometryType_->quantity());
          for (int i = 0; i<nNodes_; ++i) {
            
            if (nodeType(i,0)==MD_ONLY) {
              for (int j = 0; j < nsd_; j++) {
                velocity(i,j) = atomicVelocity(i,j);
              }
              temperature(i,0) = atomicTemperature(i,0);
            }
          }
          if (trackDisplacement_) {
            DENS_MAT & displacement(fields_[DISPLACEMENT].set_quantity());
            DENS_MAN * nodalAtomicDisplacement(interscaleManager_.dense_matrix("NodalAtomicDisplacement"));
            const DENS_MAT & atomicDisplacement(nodalAtomicDisplacement->quantity());
            for (int i = 0; i<nNodes_; ++i) {
              
              if (nodeType(i,0)==MD_ONLY) {
                for (int j = 0; j < nsd_; j++) {
                  displacement(i,j) = atomicDisplacement(i,j);
                }
              }
            }
          }
        }
      }
      initialized_ = true;
    }

    // reset integration field mask
    velocityMask_.reset(NUM_FIELDS,NUM_FLUX);
    velocityMask_ = false;
    for (int i = 0; i < NUM_FLUX; i++)
      velocityMask_(VELOCITY,i) = fieldMask_(VELOCITY,i);
    temperatureMask_.reset(NUM_FIELDS,NUM_FLUX);
    temperatureMask_ = false;
    for (int i = 0; i < NUM_FLUX; i++)
      temperatureMask_(TEMPERATURE,i) = fieldMask_(TEMPERATURE,i);

    refPE_=0;
    refPE_=potential_energy();
  }

  //--------------------------------------------------------
  //  construct_methods
  //    have managers instantiate requested algorithms
  //    and methods
  //--------------------------------------------------------
  void ATC_CouplingMomentumEnergy::construct_methods()
  {
    ATC_Coupling::construct_methods();

    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->construct_methods();
    }
    atomicRegulator_->construct_methods();
  }

  //--------------------------------------------------------
  //  construct_transfers
  //    constructs needed transfer operators
  //--------------------------------------------------------
  void ATC_CouplingMomentumEnergy::construct_transfers()
  {
    ATC_Coupling::construct_transfers();

    // momentum of each atom
    AtomicMomentum * atomicMomentum = new AtomicMomentum(this);
    interscaleManager_.add_per_atom_quantity(atomicMomentum,
                                             "AtomicMomentum");
    
    // nodal momentum for RHS
    AtfShapeFunctionRestriction * nodalAtomicMomentum = new AtfShapeFunctionRestriction(this,
                                                                                        atomicMomentum,
                                                                                        shpFcn_);
    interscaleManager_.add_dense_matrix(nodalAtomicMomentum,
                                        "NodalAtomicMomentum");
    
    // nodal forces
    FundamentalAtomQuantity * atomicForce = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_FORCE);
    AtfShapeFunctionRestriction * nodalAtomicForce = new AtfShapeFunctionRestriction(this,
                                                                                     atomicForce,
                                                                                     shpFcn_);
    interscaleManager_.add_dense_matrix(nodalAtomicForce,
                                        "NodalAtomicForce");
    
    // nodal velocity derived only from atoms
    AtfShapeFunctionMdProjection * nodalAtomicVelocity = new AtfShapeFunctionMdProjection(this,
                                                                                          nodalAtomicMomentum,
                                                                                          VELOCITY);
    interscaleManager_.add_dense_matrix(nodalAtomicVelocity,
                                        "NodalAtomicVelocity");
    
    if (trackDisplacement_) {
      // mass-weighted (center-of-mass) displacement of each atom
      AtomicMassWeightedDisplacement * atomicMassWeightedDisplacement;
      if (needXrefProcessorGhosts_ || groupbitGhost_) { // explicit construction on internal group
        PerAtomQuantity<double> * atomReferencePositions = interscaleManager_.per_atom_quantity("AtomicInternalReferencePositions");
        atomicMassWeightedDisplacement = new AtomicMassWeightedDisplacement(this,atomPositions_,
                                                                            atomMasses_,
                                                                            atomReferencePositions,
                                                                            INTERNAL);
      }
      else
        atomicMassWeightedDisplacement = new AtomicMassWeightedDisplacement(this);
      interscaleManager_.add_per_atom_quantity(atomicMassWeightedDisplacement,
                                               "AtomicMassWeightedDisplacement");
      
      // nodal (RHS) mass-weighted displacement
      AtfShapeFunctionRestriction * nodalAtomicMassWeightedDisplacement = new AtfShapeFunctionRestriction(this,
                                                                                                          atomicMassWeightedDisplacement,
                                                                                                          shpFcn_);
      interscaleManager_.add_dense_matrix(nodalAtomicMassWeightedDisplacement,
                                          "NodalAtomicMassWeightedDisplacement");
      
      // nodal displacement derived only from atoms
      AtfShapeFunctionMdProjection * nodalAtomicDisplacement = new AtfShapeFunctionMdProjection(this,
                                                                                                nodalAtomicMassWeightedDisplacement,
                                                                                                VELOCITY);
      interscaleManager_.add_dense_matrix(nodalAtomicDisplacement,
                                          "NodalAtomicDisplacement");
    }

    // always need kinetic energy
    FtaShapeFunctionProlongation * atomicMeanVelocity = new FtaShapeFunctionProlongation(this,&fields_[VELOCITY],shpFcn_);
    interscaleManager_.add_per_atom_quantity(atomicMeanVelocity,
                                             "AtomicMeanVelocity");
    AtomicEnergyForTemperature * atomicTwiceKineticEnergy = new TwiceFluctuatingKineticEnergy(this);
    AtomicEnergyForTemperature * atomEnergyForTemperature = NULL;

    // Appropriate per-atom quantity based on desired temperature definition
    if (temperatureDef_==KINETIC) {
      atomEnergyForTemperature = atomicTwiceKineticEnergy;
    }
    else if (temperatureDef_==TOTAL) {
      if (timeIntegrators_[TEMPERATURE]->time_integration_type() != TimeIntegrator::FRACTIONAL_STEP)
        throw ATC_Error("ATC_CouplingMomentumEnergy:construct_transfers()  on the fractional step time integrator can be used with non-kinetic defitions of the temperature");

      // kinetic energy
      interscaleManager_.add_per_atom_quantity(atomicTwiceKineticEnergy,
                                               "AtomicTwiceKineticEnergy");

      // atomic potential energy
      ComputedAtomQuantity * atomicPotentialEnergy = new ComputedAtomQuantity(this,lammpsInterface_->compute_pe_name(),
                                                                              1./(lammpsInterface_->mvv2e()));
      interscaleManager_.add_per_atom_quantity(atomicPotentialEnergy,
                                               "AtomicPotentialEnergy");

      // reference potential energy
      AtcAtomQuantity<double> * atomicReferencePotential;
      if (!initialized_) {
        atomicReferencePotential = new AtcAtomQuantity<double>(this);
        interscaleManager_.add_per_atom_quantity(atomicReferencePotential,
                                                 "AtomicReferencePotential");
        atomicReferencePotential->set_memory_type(PERSISTENT);
      }
      else {
        atomicReferencePotential = static_cast<AtcAtomQuantity<double> * >(interscaleManager_.per_atom_quantity("AtomicReferencePotential"));
      }
      AtfShapeFunctionRestriction * nodalAtomicReferencePotential = new AtfShapeFunctionRestriction(this,
                                                                                                    atomicReferencePotential,
                                                                                                    shpFcn_);
      interscaleManager_.add_dense_matrix(nodalAtomicReferencePotential,
                                               "NodalAtomicReferencePotential");

      // fluctuating potential energy
      AtomicEnergyForTemperature * atomicFluctuatingPotentialEnergy = new FluctuatingPotentialEnergy(this,
                                                                                                     atomicPotentialEnergy,
                                                                                                     atomicReferencePotential);
      interscaleManager_.add_per_atom_quantity(atomicFluctuatingPotentialEnergy,
                                               "AtomicFluctuatingPotentialEnergy");

      // atomic total energy
      atomEnergyForTemperature = new MixedKePeEnergy(this,1,1);

      // kinetic temperature measure for post-processing
      // nodal restriction of the atomic energy quantity for the temperature definition
      AtfShapeFunctionRestriction * nodalAtomicTwiceKineticEnergy = new AtfShapeFunctionRestriction(this,
                                                                                                    atomicTwiceKineticEnergy,
                                                                                                    shpFcn_);
      interscaleManager_.add_dense_matrix(nodalAtomicTwiceKineticEnergy,
                                          "NodalAtomicTwiceKineticEnergy");
      nodalAtomicKineticTemperature_ = new AtfShapeFunctionMdProjection(this,
                                                                        nodalAtomicTwiceKineticEnergy,
                                                                        TEMPERATURE);
      interscaleManager_.add_dense_matrix(nodalAtomicKineticTemperature_,
                                          "NodalAtomicKineticTemperature");

      // potential temperature measure for post-processing (must multiply by 2 for configurational temperature
      // nodal restriction of the atomic energy quantity for the temperature definition
      AtfShapeFunctionRestriction * nodalAtomicFluctuatingPotentialEnergy = new AtfShapeFunctionRestriction(this,
                                                                                                            atomicFluctuatingPotentialEnergy,
                                                                                                            shpFcn_);
      interscaleManager_.add_dense_matrix(nodalAtomicFluctuatingPotentialEnergy,
                                          "NodalAtomicFluctuatingPotentialEnergy");
      nodalAtomicConfigurationalTemperature_ = new AtfShapeFunctionMdProjection(this,
                                                                                nodalAtomicFluctuatingPotentialEnergy,
                                                                                TEMPERATURE);
      interscaleManager_.add_dense_matrix(nodalAtomicConfigurationalTemperature_,
                                          "NodalAtomicConfigurationalTemperature");
    }

    // register the per-atom quantity for the temperature definition
    interscaleManager_.add_per_atom_quantity(atomEnergyForTemperature,
                                             "AtomicEnergyForTemperature");
    
    // nodal restriction of the atomic energy quantity for the temperature definition
    AtfShapeFunctionRestriction * nodalAtomicEnergy = new AtfShapeFunctionRestriction(this,
                                                                                      atomEnergyForTemperature,
                                                                                      shpFcn_);
    interscaleManager_.add_dense_matrix(nodalAtomicEnergy,
                                        "NodalAtomicEnergy");
    
    // nodal atomic temperature field
    
    AtfShapeFunctionMdProjection * nodalAtomicTemperature = new AtfShapeFunctionMdProjection(this,
                                                                                             nodalAtomicEnergy,
                                                                                             TEMPERATURE);
    interscaleManager_.add_dense_matrix(nodalAtomicTemperature,
                                        "NodalAtomicTemperature");
    
    if (!useFeMdMassMatrix_) {
      // atomic momentum mass matrix
      FundamentalAtomQuantity * atomicMass = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_MASS);
      nodalAtomicMass_ = new AtfShapeFunctionRestriction(this,
                                                         atomicMass,
                                                         shpFcn_);
      interscaleManager_.add_dense_matrix(nodalAtomicMass_,
                                          "AtomicMomentumMassMat");

      // atomic dimensionless mass matrix
      ConstantQuantity<double> * atomicOnes = new ConstantQuantity<double>(this,1);
      interscaleManager_.add_per_atom_quantity(atomicOnes,"AtomicOnes");
      nodalAtomicCount_ = new AtfShapeFunctionRestriction(this,
                                                          atomicOnes,
                                                          shpFcn_);
      interscaleManager_.add_dense_matrix(nodalAtomicCount_,
                                          "AtomicDimensionlessMassMat");

      // classical thermodynamic heat capacity of the atoms
      HeatCapacity * heatCapacity = new HeatCapacity(this);
      interscaleManager_.add_per_atom_quantity(heatCapacity,
                                               "AtomicHeatCapacity");

      // atomic thermal mass matrix
      nodalAtomicHeatCapacity_ = new AtfShapeFunctionRestriction(this,
                                                                 heatCapacity,
                                                                 shpFcn_);
      interscaleManager_.add_dense_matrix(nodalAtomicHeatCapacity_,
                                          "NodalAtomicHeatCapacity");
    }

    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->construct_transfers();
    }
    atomicRegulator_->construct_transfers();
  }

  //---------------------------------------------------------
  //  init_filter
  //    sets up the time filtering operations in all objects
  //---------------------------------------------------------
  void ATC_CouplingMomentumEnergy::init_filter()
  {
    if (timeIntegrators_[TEMPERATURE]->time_integration_type() != TimeIntegrator::FRACTIONAL_STEP) {
      throw ATC_Error("ATC_CouplingMomentumEnergy::initialize - method only valid with fractional step time integration");
    }

    
    ATC_Coupling::init_filter();
    
    
    
    
    if (timeFilterManager_.end_equilibrate() && equilibriumStart_) {
      if (atomicRegulator_->coupling_mode(VELOCITY)==AtomicRegulator::FLUX || atomicRegulator_->coupling_mode(VELOCITY)==AtomicRegulator::GHOST_FLUX)
        // nothing needed in other cases since kinetostat force is balanced by boundary flux in FE equations
       atomicRegulator_->reset_lambda_contribution(nodalAtomicFieldsRoc_[VELOCITY].quantity(),VELOCITY);
 
      DENS_MAT powerMat(-1.*(nodalAtomicFields_[TEMPERATURE].quantity()));
      atomicRegulator_->reset_lambda_contribution(powerMat,TEMPERATURE);
    }
  }

  //---------------------------------------------------------
  //  compute_md_mass_matrix
  //    compute the mass matrix arising from only atomistic
  //    quadrature and contributions as a summation
  //---------------------------------------------------------
  void ATC_CouplingMomentumEnergy::compute_md_mass_matrix(FieldName thisField,
                                                  DIAG_MAT & massMat)
  {
    
    
    if (thisField == DISPLACEMENT || thisField == VELOCITY) {
      massMat.reset(nodalAtomicMass_->quantity());
    }
    else if (thisField == MASS_DENSITY) { // dimensionless mass matrix
      massMat.reset(nodalAtomicCount_->quantity());
    }
    else if (thisField == TEMPERATURE) {
      massMat.reset(nodalAtomicHeatCapacity_->quantity());
    }
  }     

  //--------------------------------------------------------
  //  finish
  //    final clean up after a run
  //--------------------------------------------------------
  void ATC_CouplingMomentumEnergy::finish()
  {
    // base class
    ATC_Coupling::finish();

    atomicRegulator_->finish();
  }

  //--------------------------------------------------------
  //  modify
  //    parses inputs and modifies state of the filter
  //--------------------------------------------------------
  bool ATC_CouplingMomentumEnergy::modify(int narg, char **arg)
  {
    bool foundMatch = false;
    int argIndex = 0;
    return foundMatch;
  }

  //--------------------------------------------------
  // pack_fields
  //   bundle all allocated field matrices into a list
  //   for output needs
  //--------------------------------------------------
  void ATC_CouplingMomentumEnergy::pack_quantity_fields(RESTART_LIST & data)
  {
    atomicRegulator_->pack_fields(data);
  }
  
  //--------------------------------------------------
  // write_restart_file
  //   bundle matrices that need to be saved and call
  //   fe_engine to write the file
  //--------------------------------------------------
  void ATC_CouplingMomentumEnergy::write_restart_data(string fileName, RESTART_LIST & data)
  {
    pack_quantity_fields(data);
    ATC_Method::write_restart_data(fileName,data);
  }
    
  //--------------------------------------------------
  // write_restart_file
  //   bundle matrices that need to be saved and call
  //   fe_engine to write the file
  //--------------------------------------------------
  void ATC_CouplingMomentumEnergy::read_restart_data(string fileName, RESTART_LIST & data)
  {
    pack_quantity_fields(data);
    ATC_Method::read_restart_data(fileName,data);
  }

  //--------------------------------------------------
  void ATC_CouplingMomentumEnergy::reset_nlocal()
  {
    ATC_Coupling::reset_nlocal();
    atomicRegulator_->reset_nlocal();
  }

  //--------------------------------------------------
  // reset_atom_materials
  //   update the atom materials map 
  //--------------------------------------------------
  void ATC_CouplingMomentumEnergy::reset_atom_materials()
  {
    ATC_Coupling::reset_atom_materials();
    atomicRegulator_->reset_atom_materials(elementToMaterialMap_,
                                           atomElement_);
  }

  //--------------------------------------------------
  void ATC_CouplingMomentumEnergy::pre_init_integrate()
  {
    ATC_Coupling::pre_init_integrate();
    double dt = lammpsInterface_->dt();

    // Perform any initialization, no actual integration
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->pre_initial_integrate1(dt);
    }
    
    // Apply controllers to atom velocities, if needed
    atomicRegulator_->apply_pre_predictor(dt,lammpsInterface_->ntimestep());

    // Predict nodal temperatures and time derivatives based on FE data
    // predict nodal velocities
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->pre_initial_integrate2(dt);
    }
    extrinsicModelManager_.pre_init_integrate();
  }

  //--------------------------------------------------------
  //  mid_init_integrate
  //    time integration between the velocity update and
  //    the position lammps update of Verlet step 1
  //--------------------------------------------------------
  void ATC_CouplingMomentumEnergy::mid_init_integrate()
  {
    // CONTINUOUS VELOCITY UPDATE
    
    ATC_Coupling::mid_init_integrate();
    double dt = lammpsInterface_->dt();

    // Compute nodal velocity at n+1/2
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->mid_initial_integrate1(dt);
    }
    atomicRegulator_->apply_mid_predictor(dt,lammpsInterface_->ntimestep());
    extrinsicModelManager_.mid_init_integrate();
  }

  //--------------------------------------------------------
  //  post_init_integrate
  //    time integration after the lammps atomic updates of
  //    Verlet step 1
  //--------------------------------------------------------
  void ATC_CouplingMomentumEnergy::post_init_integrate()
  {

    // CONTINUOUS DISPLACEMENT UPDATE
  
    double dt = lammpsInterface_->dt();

    // Compute nodal velocity at n+1
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->post_initial_integrate1(dt);
    }

    // Update kinetostat quantities if displacement is being regulated
    atomicRegulator_->apply_post_predictor(dt,lammpsInterface_->ntimestep());

    // Update extrisic model
    extrinsicModelManager_.post_init_integrate();

    // fixed values, non-group bcs handled through FE
    set_fixed_nodes();
 
    
    // enforce atomic boundary conditions
    if      (boundaryDynamics_==PRESCRIBED) set_ghost_atoms();
    else if (boundaryDynamics_==DAMPED_HARMONIC) initial_integrate_ghost();
    else if (boundaryDynamics_==COUPLED)         initial_integrate_ghost();
        
    // update time by a half dt
    update_time(0.5);

    ATC_Coupling::post_init_integrate();
  }

  //--------------------------------------------------------
  //  pre_final_integrate
  //    integration before the second stage lammps atomic 
  //    update of Verlet step 2
  //--------------------------------------------------------
  void ATC_CouplingMomentumEnergy::pre_final_integrate()
  {
    ATC_Coupling::pre_final_integrate();

    if      (boundaryDynamics_==DAMPED_HARMONIC) {
      apply_ghost_forces();
      final_integrate_ghost();
    }
    else if (boundaryDynamics_==COUPLED) {
      add_ghost_forces();
      final_integrate_ghost();
    }
  }

  //--------------------------------------------------
  void ATC_CouplingMomentumEnergy::post_final_integrate()
  {
    // CONTINUOUS VELOCITY RHS UPDATE

    double dt = lammpsInterface_->dt();

    // update of atomic contributions for fractional step methods
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->pre_final_integrate1(dt);
    }

    // Set sources
    prescribedDataMgr_->set_sources(time()+0.5*dt,sources_);
    extrinsicModelManager_.pre_final_integrate();
    if (timeIntegrators_[TEMPERATURE]->has_final_predictor() || timeIntegrators_[MOMENTUM]->has_final_predictor()) {
      extrinsicModelManager_.set_sources(fields_,extrinsicSources_);
      atomicRegulator_->compute_boundary_flux(fields_);
      compute_atomic_sources(velocityMask_,fields_,atomicSources_);
    }
    atomicRegulator_->apply_pre_corrector(dt,lammpsInterface_->ntimestep());

    // Compute atom-integrated rhs
    // parallel communication happens within FE_Engine
    compute_rhs_vector(velocityMask_,fields_,rhs_,FE_DOMAIN);
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->add_to_rhs();
    }
    atomicRegulator_->add_to_rhs(rhs_);
    
    // Compute and add atomic contributions to FE equations
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->post_final_integrate1(dt);
    }
    
   // fix nodes, non-group bcs applied through FE
    set_fixed_nodes();
        
    // corrector step extrinsic model
    extrinsicModelManager_.post_final_integrate();
    if (timeIntegrators_[TEMPERATURE]->has_final_corrector() || timeIntegrators_[MOMENTUM]->has_final_corrector()) {
      // set state-based sources
      extrinsicModelManager_.set_sources(fields_,extrinsicSources_);
      atomicRegulator_->compute_boundary_flux(fields_);
      compute_atomic_sources(velocityMask_,fields_,atomicSources_);
    }

    // Finish update of FE velocity
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->post_final_integrate2(dt);
    }
    
    // apply corrector phase of thermostat
    atomicRegulator_->apply_post_corrector(dt,lammpsInterface_->ntimestep());

    // final phase of time integration
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->post_final_integrate3(dt);
    }

    // Fix nodes, non-group bcs applied through FE
    set_fixed_nodes();
    
    update_time(0.5);
    
    output();
    ATC_Coupling::post_final_integrate(); // adds next step to computes
  }

  //--------------------------------------------------------------------
  //     compute_scalar : added energy
  //--------------------------------------------------------------------
  double ATC_CouplingMomentumEnergy::compute_scalar(void)
  {
    double energy = 0.0;
    energy += extrinsicModelManager_.compute_scalar();
    return energy;
  }

    //--------------------------------------------------------------------
  //     total kinetic energy
  //--------------------------------------------------------------------
  double ATC_CouplingMomentumEnergy::kinetic_energy(void)
  {
    const MATRIX & M = massMats_[VELOCITY].quantity();

    const DENS_MAT & velocity(fields_[VELOCITY].quantity());
    double mvv2e = lammpsInterface_->mvv2e();
    double kineticEnergy = 0;
    DENS_VEC velocitySquared(nNodes_);
    for (int i = 0; i < nNodes_; i++)
      for (int j = 0; j < nsd_; j++)
        velocitySquared(i) += velocity(i,j)*velocity(i,j);
    kineticEnergy = (M*velocitySquared).sum();
    kineticEnergy *= mvv2e; // convert to LAMMPS units
    return kineticEnergy;
  }
  //--------------------------------------------------------------------
  //     total potential energy
  //--------------------------------------------------------------------
  double ATC_CouplingMomentumEnergy::potential_energy(void)
  {
    Array<FieldName> mask(1);
    mask(0) = VELOCITY;
    FIELD_MATS energy;
    feEngine_->compute_energy(mask, 
                                fields_,
                                physicsModel_,
                                elementToMaterialMap_,
                                energy,
                                &(elementMask_->quantity()));
    double potentialEnergy = energy[VELOCITY].col_sum();
    double mvv2e = lammpsInterface_->mvv2e();
    potentialEnergy *= mvv2e; // convert to LAMMPS units
    return potentialEnergy-refPE_;
  }
  
  //--------------------------------------------------------------------
  //     compute_vector
  //--------------------------------------------------------------------
  // this is for direct output to lammps thermo
  double ATC_CouplingMomentumEnergy::compute_vector(int n)
  {
    // output[1] = total coarse scale kinetic energy
    // output[2] = total coarse scale potential energy
    // output[3] = total coarse scale energy
    // output[4] = total coarse scale thermal energy
    // output[5] = average temperature

    double mvv2e = lammpsInterface_->mvv2e(); // convert to lammps energy units

    if (n == 0) {
      return kinetic_energy();
    }
    else if (n == 1) {
      return potential_energy();
    }
    else if (n == 2) {
      return kinetic_energy()+potential_energy();
    }
    else if (n == 4) {
      Array<FieldName> mask(1);
      FIELD_MATS energy;
      mask(0) = TEMPERATURE;
      
      feEngine_->compute_energy(mask, 
                                fields_,
                                physicsModel_,
                                elementToMaterialMap_,
                                energy,
                                &(elementMask_->quantity()));
      
      double phononEnergy = mvv2e * energy[TEMPERATURE].col_sum();
      return phononEnergy;
    }
    else if (n == 5) {
      double aveT = (fields_[TEMPERATURE].quantity()).col_sum()/nNodes_;
      return aveT;
    }
    else if (n > 5) {
      double extrinsicValue = extrinsicModelManager_.compute_vector(n);
      return extrinsicValue;
    }

    return 0.;

  }

  //--------------------------------------------------------------------
  //     output
  //--------------------------------------------------------------------
  void ATC_CouplingMomentumEnergy::output()
  {
    if (output_now()) {
      feEngine_->departition_mesh();
      // avoid possible mpi calls
      if (nodalAtomicKineticTemperature_)
        _keTemp_ = nodalAtomicKineticTemperature_->quantity();
      if (nodalAtomicConfigurationalTemperature_)
        _peTemp_ = nodalAtomicConfigurationalTemperature_->quantity();
      
      OUTPUT_LIST outputData;
      
      // base class output
      ATC_Method::output();

      // push atc fields time integrator modifies into output arrays
      for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
        (_tiIt_->second)->post_process();
      }

      // auxilliary data
      
      for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
        (_tiIt_->second)->output(outputData);
      }
      atomicRegulator_->output(outputData);
      extrinsicModelManager_.output(outputData);
        
      DENS_MAT & velocity(nodalAtomicFields_[VELOCITY].set_quantity());
      DENS_MAT & rhs(rhs_[VELOCITY].set_quantity());
      DENS_MAT & temperature(nodalAtomicFields_[TEMPERATURE].set_quantity());
      DENS_MAT & dotTemperature(dot_fields_[TEMPERATURE].set_quantity());
      DENS_MAT & ddotTemperature(ddot_fields_[TEMPERATURE].set_quantity());
      DENS_MAT & rocTemperature(nodalAtomicFieldsRoc_[TEMPERATURE].set_quantity());
      DENS_MAT & fePower(rhs_[TEMPERATURE].set_quantity());
      if (lammpsInterface_->rank_zero()) {
        // global data
        double T_mean   = (fields_[TEMPERATURE].quantity()).col_sum(0)/nNodes_;
        feEngine_->add_global("temperature_mean",  T_mean);
        double T_stddev   = (fields_[TEMPERATURE].quantity()).col_stdev(0);
        feEngine_->add_global("temperature_std_dev",  T_stddev);
        double Ta_mean =  (nodalAtomicFields_[TEMPERATURE].quantity()).col_sum(0)/nNodes_;
        feEngine_->add_global("atomic_temperature_mean",  Ta_mean);
        double Ta_stddev =  (nodalAtomicFields_[TEMPERATURE].quantity()).col_stdev(0); 
        feEngine_->add_global("atomic_temperature_std_dev",  Ta_stddev);

        // different temperature measures, if appropriate
        if (nodalAtomicKineticTemperature_)
          outputData["kinetic_temperature"] = & _keTemp_;
        
        if (nodalAtomicConfigurationalTemperature_) {
          _peTemp_ *= 2; // account for full temperature
          outputData["configurational_temperature"] = & _peTemp_;
        }
        
        // mesh data
        outputData["NodalAtomicVelocity"] = &velocity;
        outputData["FE_Force"] = &rhs;
        if (trackDisplacement_)
          outputData["NodalAtomicDisplacement"] = & nodalAtomicFields_[DISPLACEMENT].set_quantity();
        outputData["NodalAtomicTemperature"] = &temperature;
        outputData["dot_temperature"] = &dotTemperature;
        outputData["ddot_temperature"] = &ddotTemperature;
        outputData["NodalAtomicPower"] = &rocTemperature;
        outputData["fePower"] = &fePower;
      
        feEngine_->write_data(output_index(), fields_, & outputData);
      }
      
      //      hence propagation is performed on proc 0 but not others.
      //      The real fix is to have const data in the output list
      // force optional variables to reset to keep in sync
      if (trackDisplacement_) {
        nodalAtomicFields_[DISPLACEMENT].force_reset();
      }
      fields_[VELOCITY].propagate_reset();

      feEngine_->partition_mesh();
    }
  }

    //--------------------------------------------------------
  //  set_ghost_atoms
  //    sets ghost atom positions to finite element
  //    displacements based on shape functions
  //--------------------------------------------------------
  void ATC_CouplingMomentumEnergy::set_ghost_atoms()
  {
    // set atomic displacements based on FE displacements
    double ** x = lammpsInterface_->xatom();
    // prolong
    DenseMatrix<double> ghostAtomData(nLocalGhost_,nsd_);
    if (nLocalGhost_>0)
      ghostAtomData = (shpFcnGhost_->quantity())*(fields_[DISPLACEMENT].quantity());

    for (int i = 0; i < nLocalGhost_; ++i)
      for (int j = 0; j < nsd_; ++j)
        x[ghostToAtom_(i)][j] = ghostAtomData(i,j)+xref_[ghostToAtom_(i)][j];

    
  }

  //--------------------------------------------------------
  //  add_ghost_forces
  //    add forces to dynamic ghosts
  //--------------------------------------------------------
  void ATC_CouplingMomentumEnergy::add_ghost_forces()
  {
    double **x = lammpsInterface_->xatom();
    double **v = lammpsInterface_->vatom();
    double **f = lammpsInterface_->fatom();

    // add forces
    DENS_MAT coarseDisp(nLocalGhost_,nsd_);
    DENS_MAT coarseVel(nLocalGhost_,nsd_);
    if (nLocalGhost_>0) {
      coarseDisp = (shpFcnGhost_->quantity())*(fields_[DISPLACEMENT].quantity());
      coarseVel  = (shpFcnGhost_->quantity())*(fields_[VELOCITY].quantity());
    }
    // dynamics one-way coupled to real atoms in a well tied to coarse scale
    for (int i = 0; i < nLocalGhost_; ++i) {
      for (int j = 0; j < nsd_; ++j) {
        double du = coarseDisp(i,j)+xref_[ghostToAtom_(i)][j]-x[ghostToAtom_(i)][j];
        double dv = coarseVel(i,j)-v[ghostToAtom_(i)][j];
        f[ghostToAtom_(i)][j] += mu_*du + gamma_*dv;
          
      }
    }
  }

  void ATC_CouplingMomentumEnergy::apply_ghost_forces()
  {
    double **x = lammpsInterface_->xatom();
    double **v = lammpsInterface_->vatom();
    double **f = lammpsInterface_->fatom();

    // add forces
    DENS_MAT coarseDisp(nLocalGhost_,nsd_);
    DENS_MAT coarseVel(nLocalGhost_,nsd_);
    if (nLocalGhost_>0) {
      coarseDisp = (shpFcnGhost_->quantity())*(fields_[DISPLACEMENT].quantity());
      coarseVel  = (shpFcnGhost_->quantity())*(fields_[VELOCITY].quantity());
    }
    // dynamics one-way coupled to real atoms in a well tied to coarse scale
    for (int i = 0; i < nLocalGhost_; ++i) {
      for (int j = 0; j < nsd_; ++j) {
        double du = coarseDisp(i,j)+xref_[ghostToAtom_(i)][j]-x[ghostToAtom_(i)][j];
        double dv = coarseVel(i,j)-v[ghostToAtom_(i)][j];
        f[ghostToAtom_(i)][j] = mu_*du + gamma_*dv;
          
      }
    }
  }

  //--------------------------------------------------------
  //  initial_integrate_ghost
  //    does the first step of the Verlet integration for
  //    ghost atoms, to be used with non-reflecting BCs
  //--------------------------------------------------------
  void ATC_CouplingMomentumEnergy::initial_integrate_ghost()
  {
    double dtfm;

    double **x = lammpsInterface_->xatom();
    double **v = lammpsInterface_->vatom();
    double **f = lammpsInterface_->fatom();
    int *type =  lammpsInterface_->atom_type();
    const int *mask =  lammpsInterface_->atom_mask();
    int nlocal = lammpsInterface_->nlocal();
    double dtv = lammpsInterface_->dt();
    double dtf = 0.5 * lammpsInterface_->dt() * lammpsInterface_->ftm2v();
    double *mass = lammpsInterface_->atom_mass();

    if (mass) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbitGhost_) {
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
        }
      }

    } else {
      double *rmass = lammpsInterface_->atom_rmass();
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbitGhost_) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
        }
      }
    }

  }

  //--------------------------------------------------------
  //  final_integrate_ghost
  //    does the second step of the Verlet integration for
  //    ghost atoms, to be used with non-reflecting BCs
  //--------------------------------------------------------
  void ATC_CouplingMomentumEnergy::final_integrate_ghost()
  {
    double dtfm;

    double **v = lammpsInterface_->vatom();
    double **f = lammpsInterface_->fatom();
    int *type =  lammpsInterface_->atom_type();
    const int *mask =  lammpsInterface_->atom_mask();
    int nlocal = lammpsInterface_->nlocal();
    double dtf = 0.5 * lammpsInterface_->dt() * lammpsInterface_->ftm2v();

    double *mass = lammpsInterface_->atom_mass();
    if (mass) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbitGhost_) {
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
        }
      }

    } else {
      double *rmass = lammpsInterface_->atom_rmass();
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbitGhost_) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
        }
      }
    }

  }
};
