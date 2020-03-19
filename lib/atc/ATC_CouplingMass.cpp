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

    // reset integration field mask
    intrinsicMask_.reset(NUM_FIELDS,NUM_FLUX);
    intrinsicMask_ = false;

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

  //WIP_JAT consolidate to coupling when we handle the temperature correctly
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

      // auxiliary data
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
