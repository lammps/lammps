// ATC headers
#include "ATC_CouplingMomentum.h"
#include "ATC_Error.h"
#include "LammpsInterface.h"
#include "PrescribedDataManager.h"
#include "PerAtomQuantity.h"
#include "TransferOperator.h"

// Other Headers
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <iostream>

using std::string;

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ATC_CouplingMomentum
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ATC_CouplingMomentum::ATC_CouplingMomentum(string groupName, 
                                             double **& perAtomArray,
                                             LAMMPS_NS::Fix * thisFix,
                                             string matParamFile,
                                             PhysicsType intrinsicModel,
                                             ExtrinsicModelType extrinsicModel)
    : ATC_Coupling(groupName,perAtomArray,thisFix),
      nodalAtomicMass_(NULL),
      nodalAtomicCount_(NULL),
      refPE_(0)
  {
    // Allocate PhysicsModel 
    create_physics_model(intrinsicModel, matParamFile);

    // create extrinsic physics model
    if (extrinsicModel != NO_MODEL) {
      extrinsicModelManager_.create_model(extrinsicModel,matParamFile);  
    }
  
    // set up field data based on physicsModel
    physicsModel_->num_fields(fieldSizes_,fieldMask_);

    // Defaults
    set_time();
    bndyIntType_ = FE_INTERPOLATION;
    trackCharge_ = false;

    // use a kinetostat
    atomicRegulator_ = new Kinetostat(this);

    // set time integrator and change any defaults based on model type
    if (intrinsicModel == ELASTIC) {
      trackDisplacement_ = true;
      fieldSizes_[DISPLACEMENT] = fieldSizes_[VELOCITY];
      timeIntegrators_[VELOCITY] = new MomentumTimeIntegrator(this,TimeIntegrator::VERLET);
      ghostManager_.set_boundary_dynamics(GhostManager::PRESCRIBED);
    }
    else if (intrinsicModel == SHEAR) {
      atomToElementMapType_ = EULERIAN;
      atomToElementMapFrequency_ = 1;
      timeIntegrators_[VELOCITY] = new MomentumTimeIntegrator(this,TimeIntegrator::GEAR);
      ghostManager_.set_boundary_dynamics(GhostManager::NO_BOUNDARY_DYNAMICS);
    }

    // output variable vector info:
    // output[1] = total coarse scale kinetic energy
    // output[2] = total coarse scale potential energy
    // output[3] = total coarse scale energy
    scalarFlag_ = 1;
    vectorFlag_ = 1;
    sizeVector_ = 5;
    scalarVectorFreq_ = 1;
    extVector_ = 1;
    thermoEnergyFlag_ = 1;
    if (extrinsicModel != NO_MODEL)
      sizeVector_ += extrinsicModelManager_.size_vector(sizeVector_);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ATC_CouplingMomentum::~ATC_CouplingMomentum()
  {
    interscaleManager_.clear();
  }

  //--------------------------------------------------------
  //  initialize
  //    sets up all the necessary data
  //--------------------------------------------------------
  void ATC_CouplingMomentum::initialize()
  {
    // clear displacement entries if requested
    if (!trackDisplacement_) {
      fieldSizes_.erase(DISPLACEMENT);
      for (int i = 0; i < NUM_FLUX; i++)
        fieldMask_(DISPLACEMENT,i) = false;
    }

    // Base class initalizations
    ATC_Coupling::initialize();

    // check resetting precedence:
    // time integrator -> kinetostat -> time filter

    // other initializations
    if (reset_methods()) {
      for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
        (_tiIt_->second)->initialize();
      }
      atomicRegulator_->initialize();
    }
    extrinsicModelManager_.initialize();
    if (timeFilterManager_.need_reset()) { // reset kinetostat power
      init_filter();
    }
    // clears need for reset
    timeFilterManager_.initialize();
    ghostManager_.initialize();
    
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
          const INT_ARRAY & nodeType(nodalGeometryType_->quantity());
          for (int i = 0; i<nNodes_; ++i) {
            
            if (nodeType(i,0)==MD_ONLY) {
              for (int j = 0; j < nsd_; j++) {
                velocity(i,j) = atomicVelocity(i,j);
              }
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


    refPE_=0;
    refPE_=potential_energy();
  }

  //--------------------------------------------------------
  //  construct_methods
  //    have managers instantiate requested algorithms
  //    and methods
  //--------------------------------------------------------
  void ATC_CouplingMomentum::construct_methods()
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
  void ATC_CouplingMomentum::construct_transfers()
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

    // atomic mass matrix data
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
  void ATC_CouplingMomentum::init_filter()
  {
    
    ATC_Coupling::init_filter();

    if (timeFilterManager_.end_equilibrate() && equilibriumStart_) // set up correct initial lambda forces to enforce initial accerlation
      if (atomicRegulator_->coupling_mode()==AtomicRegulator::FLUX || atomicRegulator_->coupling_mode()==AtomicRegulator::GHOST_FLUX)
        // nothing needed in other cases since kinetostat force is balanced by boundary flux in FE equations
        atomicRegulator_->reset_lambda_contribution(nodalAtomicFieldsRoc_[VELOCITY].quantity());
  }

  //---------------------------------------------------------
  //  compute_md_mass_matrix
  //    compute the mass matrix arising from only atomistic
  //    quadrature and contributions as a summation
  //---------------------------------------------------------
  void ATC_CouplingMomentum::compute_md_mass_matrix(FieldName thisField,
                                                    DIAG_MAT & massMat)
  {
    
    
    if (thisField == DISPLACEMENT || thisField == VELOCITY)
      massMat.reset(nodalAtomicMass_->quantity());
    else if (thisField == MASS_DENSITY) { // dimensionless mass matrix
      massMat.reset(nodalAtomicCount_->quantity());
    }
  }

  //--------------------------------------------------------
  //  finish
  //    final clean up after a run
  //--------------------------------------------------------
  void ATC_CouplingMomentum::finish()
  {
    // base class
    ATC_Coupling::finish();

    atomicRegulator_->finish();
  }

  //--------------------------------------------------------
  //  modify
  //    parses inputs and modifies state of the filter
  //--------------------------------------------------------
  bool ATC_CouplingMomentum::modify(int narg, char **arg)
  {
    bool foundMatch = false;
    int argIndex = 0;

    // check to see if it is a transfer class command
    // check derived class before base class

    // pass-through to kinetostat
    if (strcmp(arg[argIndex],"control")==0) {
      argIndex++;
      foundMatch = atomicRegulator_->modify(narg-argIndex,&arg[argIndex]);
    }

    // pass-through to timeIntegrator class
    else if (strcmp(arg[argIndex],"time_integration")==0) {
      argIndex++;
      foundMatch = timeIntegrators_[VELOCITY]->modify(narg-argIndex,&arg[argIndex]);
    }

    // switch for if displacement is tracked or not
    /*! \page man_track_displacement fix_modify AtC track_displacement
      \section syntax
      fix_modify AtC track_displacement <on/off> \n
      \section examples
      <TT> fix_modify atc track_displacement on </TT> \n
      \section description
      Determines whether displacement is tracked or not.  For solids problems this is a useful quantity, but for fluids it is not relevant.
      \section restrictions
      Some constitutive models require the displacement field
      \section default
      on
    */
    else if (strcmp(arg[argIndex],"track_displacement")==0) {
      argIndex++;
      if (strcmp(arg[argIndex],"on")==0) {
        trackDisplacement_ = true;
        foundMatch = true;
      }
      else if (strcmp(arg[argIndex],"off")==0) {
        trackDisplacement_ = false;
        foundMatch = true;
      }
      if (foundMatch) {
        needReset_ = true;
      }
    }

    else if (strcmp(arg[argIndex],"boundary_dynamics")==0) {
      argIndex++;
      foundMatch = ghostManager_.modify(narg-argIndex,&arg[argIndex]);
    }

    // no match, call base class parser
    if (!foundMatch) {
      foundMatch = ATC_Coupling::modify(narg, arg);
    }

    return foundMatch;

  }

  //--------------------------------------------------
  // pack_fields
  //   bundle all allocated field matrices into a list
  //   for output needs
  //--------------------------------------------------
  void ATC_CouplingMomentum::pack_elastic_fields(RESTART_LIST & data)
  {
    atomicRegulator_->pack_fields(data);
  }
  
  //--------------------------------------------------
  // write_restart_file
  //   bundle matrices that need to be saved and call
  //   fe_engine to write the file
  //--------------------------------------------------
  void ATC_CouplingMomentum::write_restart_data(string fileName, RESTART_LIST & data)
  {
    pack_elastic_fields(data);
    ATC_Method::write_restart_data(fileName,data);
  }
    
  //--------------------------------------------------
  // write_restart_file
  //   bundle matrices that need to be saved and call
  //   fe_engine to write the file
  //--------------------------------------------------
  void ATC_CouplingMomentum::read_restart_data(string fileName, RESTART_LIST & data)
  {
    pack_elastic_fields(data);
    ATC_Method::read_restart_data(fileName,data);
  }

  //--------------------------------------------------------
  void ATC_CouplingMomentum::reset_nlocal()
  {
    ATC_Coupling::reset_nlocal();
    atomicRegulator_->reset_nlocal();
  }

  //--------------------------------------------------
  // reset_atom_materials
  //   update the atom materials map 
  //--------------------------------------------------
  void ATC_CouplingMomentum::reset_atom_materials()
  {
    ATC_Coupling::reset_atom_materials();
    atomicRegulator_->reset_atom_materials(elementToMaterialMap_,
                                           atomElement_);
  }

#ifdef OBSOLETE
  //--------------------------------------------------------
  //  mid_init_integrate
  //    time integration between the velocity update and
  //    the position lammps update of Verlet step 1
  //--------------------------------------------------------
  void ATC_CouplingMomentum::mid_init_integrate()
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
  void ATC_CouplingMomentum::post_init_integrate()
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
        
    // update time by a half dt
    update_time(0.5);

    ATC_Coupling::post_init_integrate();
  }
#endif
  //--------------------------------------------------------
  //  pre_final_integrate
  //    integration before the second stage lammps atomic 
  //    update of Verlet step 2
  //--------------------------------------------------------
  void ATC_CouplingMomentum::pre_final_integrate()
  {
    ATC_Coupling::pre_final_integrate();
  }

  //--------------------------------------------------------
  //  post_final_integrate
  //    integration after the second stage lammps atomic 
  //    update of Verlet step 2
  //--------------------------------------------------------
  void ATC_CouplingMomentum::post_final_integrate()
  {
    // COMPUTE FORCES FOR FE VELOCITY RHS

    double dt = lammpsInterface_->dt();

    // updating of data based on atomic forces
   for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->pre_final_integrate1(dt);
    }

    // Set prescribed sources for current time
    prescribedDataMgr_->set_sources(time()+0.5*dt,sources_);

    // predictor step in extrinsic model
    extrinsicModelManager_.pre_final_integrate();

    
    if (timeIntegrators_[VELOCITY]->has_final_predictor()) {
      // set state-based sources
      extrinsicModelManager_.set_sources(fields_,extrinsicSources_);
      atomicRegulator_->compute_boundary_flux(fields_);
      compute_atomic_sources(velocityMask_,fields_,atomicSources_);
    }

    // Compute kinetostat forces and add kinetostat contributions to FE equations
    
    atomicRegulator_->apply_pre_corrector(dt,lammpsInterface_->ntimestep());  // computes but does not apply kstat, and only for StressFlux

    // set state-based RHS
    // Determine FE contributions to dv/dt-----------------------
    // Compute atom-integrated rhs
    // parallel communication happens within FE_Engine
    compute_rhs_vector(velocityMask_,fields_,rhs_,FE_DOMAIN);
    // Compute and add atomic contributions to FE equations
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->add_to_rhs();
    }
    // add in kinetostat contributions to FE equations
    atomicRegulator_->add_to_rhs(rhs_);

    // final phase predictor step
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->post_final_integrate1(dt);
    }

    // fix nodes, non-group bcs applied through FE
    set_fixed_nodes();

    // CONTINUOUS VELOCITY RHS UPDATE

    // corrector step extrinsic model
    extrinsicModelManager_.post_final_integrate();

    if (timeIntegrators_[VELOCITY]->has_final_corrector()) {
      // set state-based sources
      extrinsicModelManager_.set_sources(fields_,extrinsicSources_);
      atomicRegulator_->compute_boundary_flux(fields_);
      compute_atomic_sources(velocityMask_,fields_,atomicSources_);
    }
        
    // Finish update of FE velocity
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->post_final_integrate2(dt);
    }

    // Apply kinetostat to atoms
    atomicRegulator_->apply_post_corrector(dt,lammpsInterface_->ntimestep());

    // finalize time integration
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->post_final_integrate3(dt);
    }

    // Fix nodes, non-group bcs applied through FE
    set_fixed_nodes();

    // update time by a half dt
    update_time(0.5);

    output();
    ATC_Coupling::post_final_integrate(); // addstep for computes
  }

  //--------------------------------------------------------
  //  min_pre_force
  //    add to interatomic forces for minimize
  //--------------------------------------------------------
  void ATC_CouplingMomentum::min_pre_force()
  {
  }

  //--------------------------------------------------------
  //  min_post_force
  //    add to interatomic forces for minimize
  //    this determines the search direction
  //--------------------------------------------------------
  void ATC_CouplingMomentum::min_post_force()
  {
    // reset positions and shape functions
    ATC_Method::min_post_force();

    // Set sources
    
    prescribedDataMgr_->set_sources(time(),sources_);
    extrinsicModelManager_.set_sources(fields_,extrinsicSources_);
    extrinsicModelManager_.pre_final_integrate();



    
    if (outputNow_) {
      update_time(1.0);
      update_step();
      output();
      outputNow_ = false;
    }
    
    
    localStep_ += 1;
  }
  
  //--------------------------------------------------------
  //  output
  //    does post-processing steps and outputs data
  //--------------------------------------------------------
  void ATC_CouplingMomentum::output() 
  {
    if (output_now()) {
      feEngine_->departition_mesh();
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
      if (lammpsInterface_->rank_zero()) {
        // mesh data
        outputData["NodalAtomicVelocity"] = &velocity;
        outputData["FE_Force"] = &rhs;
        if (trackDisplacement_) {
          outputData["NodalAtomicDisplacement"] = & nodalAtomicFields_[DISPLACEMENT].set_quantity();
        }
        
        feEngine_->write_data(output_index(), fields_, & outputData);
      }
      // force optional variables to reset to keep in sync
      if (trackDisplacement_) {
        nodalAtomicFields_[DISPLACEMENT].force_reset();
      }

      feEngine_->partition_mesh();
    }
  }

  //--------------------------------------------------------------------
  //     compute_scalar : added energy
  //        this is used in the line search
  //--------------------------------------------------------------------
  double ATC_CouplingMomentum::compute_scalar(void)
  {
    double energy = extrinsicModelManager_.compute_scalar();
    return energy;
  }

  //--------------------------------------------------------------------
  //     kinetic energy
  //--------------------------------------------------------------------
  double ATC_CouplingMomentum::kinetic_energy(const IntegrationDomainType domain) // const
  {
    const MATRIX & M = massMats_[VELOCITY].quantity();
    const DENS_MAT & velocity(fields_[VELOCITY].quantity());
    double kineticEnergy = 0;
    for (int j = 0; j < nsd_; j++) {
      CLON_VEC v = column(velocity,j);
      kineticEnergy += v.dot(M*v);
    }
    if (domain == FE_DOMAIN) {
      
      Array<FieldName> massMask(1);
      massMask(0) = VELOCITY;
      feEngine_->compute_lumped_mass_matrix(massMask,fields_,physicsModel_,atomMaterialGroups_,
                                            atomVolume_->quantity(),shpFcn_->quantity(),
                                            Ma_);
      const MATRIX & Ma = Ma_[VELOCITY].quantity();
      for (int j = 0; j < nsd_; j++) {
        CLON_VEC v = column(velocity,j);
        kineticEnergy -= v.dot(Ma*v);
      }
    }
    double mvv2e = lammpsInterface_->mvv2e(); 
    kineticEnergy *= 0.5*mvv2e; // convert to LAMMPS units

    return kineticEnergy;
  }
  //--------------------------------------------------------------------
  //     potential/strain energy 
  //--------------------------------------------------------------------
  double ATC_CouplingMomentum::potential_energy(const IntegrationDomainType domain) const
  {
    Array<FieldName> mask(1);
    mask(0) = VELOCITY;
    FIELD_MATS energy;
    feEngine_->compute_energy(mask, 
                                fields_,
                                physicsModel_,
                                elementToMaterialMap_,
                                energy,
                                &(elementMask_->quantity()),
                                domain);
    double potentialEnergy = energy[VELOCITY].col_sum();
    double mvv2e = lammpsInterface_->mvv2e();
    potentialEnergy *= mvv2e; // convert to LAMMPS units
    return potentialEnergy-refPE_;
  }
  //--------------------------------------------------------------------
  //     compute_vector
  //--------------------------------------------------------------------
  // this is for direct output to lammps thermo
  double ATC_CouplingMomentum::compute_vector(int n)
  {
    // output[1] = total coarse scale kinetic energy
    // output[2] = total coarse scale potential energy
    // output[3] = total coarse scale energy
    // output[4] = fe-only coarse scale kinetic energy
    // output[5] = fe-only coarse scale potential energy
  


    if (n == 0) {
      return kinetic_energy();
    }
    else if (n == 1) {
      return potential_energy();
    }
    else if (n == 2) {
      return kinetic_energy()+potential_energy();
    }
    else if (n == 3) {
      return kinetic_energy(FE_DOMAIN);
    }
    else if (n == 4) {
      return potential_energy(FE_DOMAIN);
    }
    else if (n > 4) {
      double extrinsicValue = extrinsicModelManager_.compute_vector(n);
      return extrinsicValue;
    }
    return 0.;

  }

};
