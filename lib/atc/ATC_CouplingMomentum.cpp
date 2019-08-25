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

    // reset integration field mask
    intrinsicMask_.reset(NUM_FIELDS,NUM_FLUX);
    intrinsicMask_ = false;
    for (int i = 0; i < NUM_FLUX; i++)
      intrinsicMask_(VELOCITY,i) = fieldMask_(VELOCITY,i);


    refPE_=0;
    refPE_=potential_energy();
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
