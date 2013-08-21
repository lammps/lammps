// ATC_Transfer headers
#include "ATC_CouplingMass.h"
#include "ATC_Error.h"
#include "FE_Engine.h"
#include "SpeciesTimeIntegrator.h"
#include "PrescribedDataManager.h"
#include "ExtrinsicModelElectrostatic.h"
#include "PoissonSolver.h"
#include "ChargeRegulator.h"
#include "ConcentrationRegulator.h"
#include "PerAtomQuantityLibrary.h"
#include "TransferOperator.h"
#include "AtomToMoleculeTransfer.h"
#include "MoleculeSet.h"
#include "FieldManager.h"

// Other Headers
#include <vector>
#include <set>
#include <utility>

using ATC_Utility::to_string;
using std::map;
using std::string;
using std::pair;

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ATC_CouplingMass
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ATC_CouplingMass::ATC_CouplingMass(string groupName, 
                                     double **& perAtomArray,
                                     LAMMPS_NS::Fix * thisFix,
                                     string matParamFile,
                                     ExtrinsicModelType extrinsicModel)
    : ATC_Coupling(groupName,perAtomArray,thisFix),
      resetNlocal_(false)
  {
    // Allocate PhysicsModel 
    create_physics_model(SPECIES, matParamFile); 

    // create extrinsic physics model
    if (extrinsicModel != NO_MODEL) {
      extrinsicModelManager_.create_model(extrinsicModel,matParamFile);  
    }

    // Defaults
    set_time();
    bndyIntType_ = NO_QUADRATURE;
    
  
    // set up field data based on physicsModel
    physicsModel_->num_fields(fieldSizes_,fieldMask_);

    // regulator
    atomicRegulator_ = new ConcentrationRegulator(this);

    // set up physics specific time integrator
    //WIP_JAT should be species concentration
    timeIntegrators_[MASS_DENSITY] = new SpeciesTimeIntegrator(this,TimeIntegrator::FRACTIONAL_STEP);

    // output variable vector info:
    // output[1] = system mass density
    vectorFlag_ = 1;
    sizeVector_ = 0;
    scalarVectorFreq_ = 1;
    extVector_ = 1;
    if (extrinsicModel != NO_MODEL)
      sizeVector_ += extrinsicModelManager_.size_vector(sizeVector_);
    sizeVector_ += atomicRegulator_->size_vector(sizeVector_);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ATC_CouplingMass::~ATC_CouplingMass()
  {
    interscaleManager_.clear();
  }

  //--------------------------------------------------------
  //  modify
  //    parses inputs and modifies state 
  //--------------------------------------------------------
  bool ATC_CouplingMass::modify(int narg, char **arg)
  {
    bool match = false;
    // check to see if it is a transfer class command
    
    // check derived class before base class
    int argIndex = 0;
    // pass-through to concentration regulator
    if (strcmp(arg[argIndex],"control")==0) {
      argIndex++;
      if (strcmp(arg[argIndex],"concentration")==0) {
        argIndex++;
        match = atomicRegulator_->modify(narg-argIndex,&arg[argIndex]);
      }
    }
    // no match, call base class parser
    if (!match) {
      match = ATC_Coupling::modify(narg, arg);
    }
    return match;
  }

  //--------------------------------------------------------
  //  initialize
  //    sets up all the necessary data
  //--------------------------------------------------------
  void ATC_CouplingMass::initialize()
  {
    
    fieldSizes_[SPECIES_CONCENTRATION] = ntracked();

    // Base class initalizations
    ATC_Coupling::initialize();

    // check that only all atoms

    if (bndyIntType_ != NO_QUADRATURE) throw ATC_Error("ATC_CouplingMass: only all atoms simulations are supported");

    // set consistent initial conditions, if requested
    if (!timeFilterManager_.filter_dynamics()) {
      if (consistentInitialization_) {
        
        DENS_MAT & massDensity(fields_[MASS_DENSITY].set_quantity());
        const DENS_MAT & atomicMassDensity(atomicFields_[MASS_DENSITY]->quantity());
        
        DENS_MAT & speciesConcentration(fields_[SPECIES_CONCENTRATION].set_quantity());
        const DENS_MAT & atomicSpeciesConcentration(atomicFields_[SPECIES_CONCENTRATION]->quantity());

        const INT_ARRAY & nodeType(nodalGeometryType_->quantity());
        for (int i = 0; i<nNodes_; ++i) {
          
          if (nodeType(i,0)==MD_ONLY) {
            massDensity(i,0) = atomicMassDensity(i,0);
            for (int j = 0; j < atomicSpeciesConcentration.nCols(); ++j) {
              speciesConcentration(i,j) = atomicSpeciesConcentration(i,j);
            }
          }
        }
      }
    }


    // other initializatifields_[SPECIES_CONCENTRATION].quantity()ons
    if (reset_methods()) { 
      for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
        (_tiIt_->second)->initialize();
      }
    }
    extrinsicModelManager_.initialize();  // always needed to construct new Poisson solver
    if (timeFilterManager_.need_reset()) {
      init_filter();
    }
    // clears need for reset
    timeFilterManager_.initialize();
    atomicRegulator_->initialize();
    ghostManager_.initialize();
    
    if (!initialized_) {
      // initialize sources based on initial FE temperature
      double dt = lammpsInterface_->dt(); 
      // set sources

      prescribedDataMgr_->set_sources(time()+0.5*dt,sources_);
      extrinsicModelManager_.set_sources(fields_,extrinsicSources_);
      compute_atomic_sources(fieldMask_,fields_,atomicSources_);

      // read in field data if necessary
      if (useRestart_) {
        RESTART_LIST data;
        read_restart_data(restartFileName_,data);
        useRestart_ = false;
      }
      
      initialized_ = true;
    }

    // reset integration field mask
    speciesMask_.reset(NUM_FIELDS,NUM_FLUX);
    speciesMask_ = false;

  }

  //--------------------------------------------------------
  //  construct_methods
  //    have managers instantiate requested algorithms
  //    and methods
  //--------------------------------------------------------
  void ATC_CouplingMass::construct_methods()
  {
    ATC_Coupling::construct_methods();
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->construct_methods();
    }
    atomicRegulator_->construct_methods();
  }

  void ATC_CouplingMass::construct_transfers()
  {
    ATC_Coupling::construct_transfers();
    FieldManager fmgr(this);
    atomicFields_[MASS_DENSITY]  = fmgr.nodal_atomic_field(MASS_DENSITY, field_to_intrinsic_name(MASS_DENSITY));

    if (has_tracked_species()) { 
      atomicFields_[SPECIES_CONCENTRATION]  = fmgr.nodal_atomic_field(SPECIES_CONCENTRATION, field_to_intrinsic_name(SPECIES_CONCENTRATION));
      
      //if (atomicRegulator_->needs_temperature()) {

        atomicFields_[TEMPERATURE]  = fmgr.nodal_atomic_field(KINETIC_TEMPERATURE, field_to_intrinsic_name(TEMPERATURE));
        //atomicFields_[TEMPERATURE]  = fmgr.nodal_atomic_field(TEMPERATURE, field_to_intrinsic_name(TEMPERATURE));
        field(TEMPERATURE) = atomicFields_[TEMPERATURE]->quantity();
        //}
    }
    else {
      
      throw ATC_Error("ATC_CouplingMass: no tracked species");
    }

    //==========================================================================
    // add molecule mass density transfer operators
    //==========================================================================
    map<string,pair<MolSize,int> >::const_iterator molecule;
    FundamentalAtomQuantity * mass = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_MASS,
                                                                                  PROC_GHOST);

    for (molecule = moleculeIds_.begin(); molecule != moleculeIds_.end(); molecule++) {
      const string moleculeName = molecule->first;
      SmallMoleculeSet * smallMoleculeSet = interscaleManager_.small_molecule_set(moleculeName);
      SPAR_MAN * shpFcnMol = interscaleManager_.sparse_matrix("ShapeFunction"+moleculeName);
      AtomToSmallMoleculeTransfer<double> * moleculeMass = 
        new AtomToSmallMoleculeTransfer<double>(this,mass,smallMoleculeSet);
      interscaleManager_.add_dense_matrix(moleculeMass,"MoleculeMass"+moleculeName);
      MotfShapeFunctionRestriction * nodalAtomicMoleculeMass = 
        new MotfShapeFunctionRestriction(moleculeMass,shpFcnMol);
      interscaleManager_.add_dense_matrix(nodalAtomicMoleculeMass,"NodalMoleculeMass"+moleculeName);

      
      AtfShapeFunctionMdProjection * nodalAtomicMoleculeMassDensity =
        new AtfShapeFunctionMdProjection(this,nodalAtomicMoleculeMass,MASS_DENSITY);
      interscaleManager_.add_dense_matrix(nodalAtomicMoleculeMassDensity,"NodalMoleculeMassDensity"+moleculeName);
    }
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->construct_transfers();
    }
  }

  void ATC_CouplingMass::init_filter()
  {
    
    ATC_Coupling::init_filter();
  }

  void ATC_CouplingMass::compute_md_mass_matrix(FieldName thisField,
                                                DIAG_MAT & massMat)
  {

    if (thisField == MASS_DENSITY ||
        thisField == SPECIES_CONCENTRATION) {
      massMat.reset(nodalAtomicVolume_->quantity());
    }
  }

  //--------------------------------------------------
  // write_restart_file
  //   bundle matrices that need to be saved and call
  //   fe_engine to write the file
  //--------------------------------------------------
  void ATC_CouplingMass::write_restart_data(string fileName, RESTART_LIST & data)
  {
    ATC_Method::write_restart_data(fileName,data);
  }
    
  //--------------------------------------------------
  // write_restart_file
  //   bundle matrices that need to be saved and call
  //   fe_engine to write the file
  //--------------------------------------------------
  void ATC_CouplingMass::read_restart_data(string fileName, RESTART_LIST & data)
  {
    ATC_Method::read_restart_data(fileName,data);
  }

  //--------------------------------------------------------
  //  pre_force
  //    prior to calculation of forces
  //--------------------------------------------------------
  void ATC_CouplingMass::pre_force()
  {
    ATC_Coupling::pre_force();
    atomicRegulator_->pre_force(); 
  }

  //--------------------------------------------------------
  //  pre_exchange
  //    prior to exchange of atoms
  //--------------------------------------------------------
  void ATC_CouplingMass::pre_exchange()
  {
    ATC_Coupling::pre_exchange();
    
    //if (atomicRegulator_->needs_temperature()) {
      field(TEMPERATURE) = atomicFields_[TEMPERATURE]->quantity(); 
///}
    atomicRegulator_->pre_exchange(); 
    if (resetNlocal_) {
      this->reset_nlocal();
      resetNlocal_ = false;
    }
  }

#ifdef OBSOLETE
  //--------------------------------------------------------
  //  mid_init_integrate
  //    time integration between the velocity update and
  //    the position lammps update of Verlet step 1
  //--------------------------------------------------------
  void ATC_CouplingMass::mid_init_integrate()
  {
    ATC_Coupling::mid_init_integrate();
    double dt = lammpsInterface_->dt();

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
  void ATC_CouplingMass::post_init_integrate()
  {
    double dt = lammpsInterface_->dt();
  
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->post_initial_integrate1(dt);
    }

    atomicRegulator_->apply_post_predictor(dt,lammpsInterface_->ntimestep());

    extrinsicModelManager_.post_init_integrate();

    set_fixed_nodes();
    update_time(0.5); // half step
    ATC_Coupling::post_init_integrate();
  }
#endif
  //--------------------------------------------------------
  //  post_final_integrate
  //    integration after the second stage lammps atomic 
  //    update of Verlet step 2
  //--------------------------------------------------------
  void ATC_CouplingMass::post_final_integrate()
  {
    double dt = lammpsInterface_->dt();

    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->pre_final_integrate1(dt);
    }

    prescribedDataMgr_->set_sources(time()+0.5*dt,sources_);
    extrinsicModelManager_.pre_final_integrate();

    
    if (timeIntegrators_[MASS_DENSITY]->has_final_predictor()) {
      // set state-based sources
      extrinsicModelManager_.set_sources(fields_,extrinsicSources_);
      
      compute_atomic_sources(speciesMask_,fields_,atomicSources_);
    }

    

    // set state-based RHS
    // Determine FE contributions to dv/dt-----------------------
    // Compute atom-integrated rhs
    // parallel communication happens within FE_Engine
    compute_rhs_vector(speciesMask_,fields_,rhs_,FE_DOMAIN);
    // Compute and add atomic contributions to FE equations
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->add_to_rhs();
    }
    atomicRegulator_->add_to_rhs(rhs_);

    // final phase predictor step
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->post_final_integrate1(dt);
    }

    set_fixed_nodes();

    // corrector step extrinsic model
    extrinsicModelManager_.post_final_integrate();

    if (timeIntegrators_[MASS_DENSITY]->has_final_corrector()) {
      // set state-based sources
      extrinsicModelManager_.set_sources(fields_,extrinsicSources_);
      
      compute_atomic_sources(speciesMask_,fields_,atomicSources_);
    }

    // finish FE temperature update
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->post_final_integrate2(dt);
    }

    // apply corrector phase of thermostat
    atomicRegulator_->apply_post_corrector(dt,lammpsInterface_->ntimestep());

    // finalize time integration
    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->post_final_integrate3(dt);
    }

    // Fix nodes, non-group bcs applied through FE
    set_fixed_nodes();

    update_time(0.5);
    output();
    ATC_Coupling::post_final_integrate(); // addstep for computes
  }
  
  //--------------------------------------------------------
  //  output
  //    does post-processing steps and outputs data
  //--------------------------------------------------------
  void ATC_CouplingMass::output()
  { 
    if (output_now()) {
      feEngine_->departition_mesh();
      OUTPUT_LIST outputData;

      // base class output
      ATC_Coupling::output();

      // push atc fields time integrator modifies into output arrays
      for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
        (_tiIt_->second)->post_process();
      }

      // auxilliary data
      for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
        (_tiIt_->second)->output(outputData);
      }
      extrinsicModelManager_.output(outputData);
      atomicRegulator_->output(outputData); 

      FIELD_POINTERS::iterator itr;
      for (itr=atomicFields_.begin(); itr!=atomicFields_.end();itr++) { 
        FieldName name = itr->first;
        const DENS_MAT & data = (itr->second)->quantity();
        outputData[field_to_intrinsic_name(name)] = & data;
      }
      // compute partial forces
      int * type =lammpsInterface_->atom_type();
      double ** f =lammpsInterface_->fatom();
      for (unsigned int j = 0; j < typeList_.size(); j++) {
        string speciesName = typeNames_[j];
        int sType = typeList_[j];
        double localF[3] = {0,0,0}, F[3] = {0,0,0};
        for (int i = 0; i < nLocal_; i++) {
          int a = internalToAtom_(i);
          if (sType == type[a]) {
            double * fa = f[a];
            localF[0] += fa[0];
            localF[1] += fa[1];
            localF[2] += fa[2];
          }
         }
         lammpsInterface_->allsum(localF,F,3);
         if (lammpsInterface_->rank_zero()) {
           for (int i = 0; i < 3; ++i)  {
            feEngine_->add_global(speciesName+"_F"+to_string(i+1), F[i]);
           }
         }
      }
      if (lammpsInterface_->rank_zero()) {
        // tagged data --only for molecule
        map<string,DENS_MAN>::iterator densMan;
        for (densMan = taggedDensMan_.begin(); densMan != taggedDensMan_.end(); densMan++) {
          outputData[densMan->first] = & (densMan->second).set_quantity();
        }

        feEngine_->write_data(output_index(), fields_, & outputData);
      }
      // force reset of tagged data to keep in sync
      map<string,DENS_MAN>::iterator densMan;
      for (densMan = taggedDensMan_.begin(); densMan != taggedDensMan_.end(); densMan++)
        (densMan->second).force_reset();
      feEngine_->partition_mesh();
    }
  }

  //--------------------------------------------------------------------
  //     compute_vector
  //--------------------------------------------------------------------
  // this is for direct output to lammps thermo
  double ATC_CouplingMass::compute_vector(int n)
  {
    return atomicRegulator_->compute_vector(n);
  }
};
