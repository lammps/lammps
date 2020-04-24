// ATC_Transfer headers
#include "ATC_CouplingEnergy.h"
#include "Thermostat.h"
#include "ATC_Error.h"
#include "PrescribedDataManager.h"
#include "FieldManager.h"

// Other Headers
#include <vector>
#include <set>
#include <utility>
#include <typeinfo>

using std::string;

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ATC_CouplingEnergy
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ATC_CouplingEnergy::ATC_CouplingEnergy(string groupName,
                                         double ** & perAtomArray,
                                         LAMMPS_NS::Fix * thisFix,
                                         string matParamFile,
                                         ExtrinsicModelType extrinsicModel)
    : ATC_Coupling(groupName,perAtomArray,thisFix),
      nodalAtomicKineticTemperature_(NULL),
      nodalAtomicConfigurationalTemperature_(NULL)
  {
    // Allocate PhysicsModel 
    create_physics_model(THERMAL, matParamFile);

    // create extrinsic physics model
    if (extrinsicModel != NO_MODEL) {
      extrinsicModelManager_.create_model(extrinsicModel,matParamFile);  
    }

    // Defaults
    set_time();
    bndyIntType_ = FE_INTERPOLATION;
  
    // set up field data based on physicsModel
    physicsModel_->num_fields(fieldSizes_,fieldMask_);

    // set up atomic regulator
    atomicRegulator_ = new Thermostat(this);

    // set up physics specific time integrator and thermostat
    timeIntegrators_[TEMPERATURE] = new ThermalTimeIntegrator(this,TimeIntegrator::GEAR);

    // default physics
    temperatureDef_ = KINETIC;

    // output variable vector info:
    // output[1] = total coarse scale thermal energy
    // output[2] = average temperature
    vectorFlag_ = 1;
    sizeVector_ = 2;
    scalarVectorFreq_ = 1;
    extVector_ = 1;
    if (extrinsicModel != NO_MODEL)
      sizeVector_ += extrinsicModelManager_.size_vector(sizeVector_);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ATC_CouplingEnergy::~ATC_CouplingEnergy()
  {
    // clear out all managed memory to avoid conflicts with dependencies on class member data
    interscaleManager_.clear();
  }

  //--------------------------------------------------------
  //  initialize
  //    sets up all the necessary data
  //--------------------------------------------------------
  void ATC_CouplingEnergy::initialize()
  {
    // Base class initalizations
    ATC_Coupling::initialize();

    // reset integration field mask
    intrinsicMask_.reset(NUM_FIELDS,NUM_FLUX);
    intrinsicMask_ = false;
    for (int i = 0; i < NUM_FLUX; i++)
      intrinsicMask_(TEMPERATURE,i) = fieldMask_(TEMPERATURE,i);
  }

  //--------------------------------------------------------
  //  construct_transfers
  //    constructs needed transfer operators
  //--------------------------------------------------------
  void ATC_CouplingEnergy::construct_transfers()
  {
    ATC_Coupling::construct_transfers();

    // always need kinetic energy
    AtomicEnergyForTemperature * atomicTwiceKineticEnergy = new TwiceKineticEnergy(this);
    AtomicEnergyForTemperature * atomEnergyForTemperature = NULL;

    // Appropriate per-atom quantity based on desired temperature definition
    if (temperatureDef_==KINETIC) {
      atomEnergyForTemperature = atomicTwiceKineticEnergy;
    }
    else if (temperatureDef_==TOTAL) {
      if (timeIntegrators_[TEMPERATURE]->time_integration_type() != TimeIntegrator::FRACTIONAL_STEP)
        throw ATC_Error("ATC_CouplingEnergy:construct_transfers()  on the fractional step time integrator can be used with non-kinetic defitions of the temperature");

      // kinetic energy
      interscaleManager_.add_per_atom_quantity(atomicTwiceKineticEnergy,
                                               "AtomicTwiceKineticEnergy");

      // atomic potential energy
      ComputedAtomQuantity * atomicPotentialEnergy = new ComputedAtomQuantity(this,
                                                                              lammpsInterface_->compute_pe_name(),
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
      nodalRefPotentialEnergy_ = new AtfShapeFunctionRestriction(this,
                                                                 atomicReferencePotential,
                                                                 shpFcn_);
      interscaleManager_.add_dense_matrix(nodalRefPotentialEnergy_,
                                          "NodalAtomicReferencePotential");

      // fluctuating potential energy
      AtomicEnergyForTemperature * atomicFluctuatingPotentialEnergy =
        new FluctuatingPotentialEnergy(this,
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

    for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
      (_tiIt_->second)->construct_transfers();
    }
    atomicRegulator_->construct_transfers();
  }

  //---------------------------------------------------------
  //  init_filter
  //    sets up the time filtering operations in all objects
  //---------------------------------------------------------
  void ATC_CouplingEnergy::init_filter()
  {
    
    TimeIntegrator::TimeIntegrationType timeIntegrationType = timeIntegrators_[TEMPERATURE]->time_integration_type();
    
    
    
    
    if (timeFilterManager_.end_equilibrate()) { 
      if (timeIntegrationType==TimeIntegrator::GEAR) {
        if (equilibriumStart_) {
          
          
          
          if (atomicRegulator_->regulator_target()==AtomicRegulator::DYNAMICS) { // based on FE equation
            DENS_MAT vdotflamMat(-2.*(nodalAtomicFields_[TEMPERATURE].quantity())); // note 2 is for 1/2 vdotflam addition
            atomicRegulator_->reset_lambda_contribution(vdotflamMat);
          }
          else { // based on MD temperature equation
            DENS_MAT vdotflamMat(-1.*(nodalAtomicFields_[TEMPERATURE].quantity()));
            atomicRegulator_->reset_lambda_contribution(vdotflamMat);
          }
        }
      }
      else if (timeIntegrationType==TimeIntegrator::FRACTIONAL_STEP) {
        if (equilibriumStart_) {
          DENS_MAT powerMat(-1.*(nodalAtomicFields_[TEMPERATURE].quantity()));
          atomicRegulator_->reset_lambda_contribution(powerMat);
        }
      }
    }
  }

  //--------------------------------------------------------
  //  modify
  //    parses inputs and modifies state of the filter
  //--------------------------------------------------------
  bool ATC_CouplingEnergy::modify(int narg, char **arg)
  {
    bool foundMatch = false;
    int argIndx = 0;
    
    // check to see if input is a transfer class command
    // check derived class before base class

    // pass-through to thermostat
    if (strcmp(arg[argIndx],"control")==0) {
      argIndx++;
      foundMatch = atomicRegulator_->modify(narg-argIndx,&arg[argIndx]);
    }

    // pass-through to timeIntegrator class
    else if (strcmp(arg[argIndx],"time_integration")==0) {
      argIndx++;
      foundMatch = timeIntegrators_[TEMPERATURE]->modify(narg-argIndx,&arg[argIndx]);
    }

    // switch for the kind of temperature being used
    /*! \page man_temperature_definition fix_modify AtC temperature_definition
      \section syntax
      fix_modify AtC temperature_definition <kinetic|total>

      \section examples
      <TT> fix_modify atc temperature_definition kinetic </TT> \n

      \section description
      Change the definition for the atomic temperature used to create the finite element temperature.  The kinetic option is based only on the kinetic energy of the atoms while the total option uses the total energy (kinetic + potential) of an atom.

      \section restrictions
      This command is only valid when using thermal coupling.  Also, while not a formal restriction, the user should ensure that associating a potential energy with each atom makes physical sense for the total option to be meaningful.

        \section default
        kinetic
      */
    else if (strcmp(arg[argIndx],"temperature_definition")==0) {
      argIndx++;
      string_to_temperature_def(arg[argIndx],temperatureDef_);
      if (temperatureDef_ == TOTAL) {
        setRefPE_ = true;
      }
      foundMatch = true;
      needReset_ = true;
    }

    // no match, call base class parser
    if (!foundMatch) {
      foundMatch = ATC_Coupling::modify(narg, arg);
    }

    return foundMatch;

  }

  //--------------------------------------------------------------------
  //     compute_vector
  //--------------------------------------------------------------------
  // this is for direct output to lammps thermo
  double ATC_CouplingEnergy::compute_vector(int n)
  {
    // output[1] = total coarse scale thermal energy
    // output[2] = average temperature

    double mvv2e = lammpsInterface_->mvv2e(); // convert to lammps energy units
  
    if (n == 0) {
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
    else if (n == 1) {
      double aveT = (fields_[TEMPERATURE].quantity()).col_sum()/nNodes_;
      return aveT;
    }
    else if (n > 1) {
      double extrinsicValue = extrinsicModelManager_.compute_vector(n);
      return extrinsicValue;
    }

    return 0.;

  }

  //--------------------------------------------------------------------
  //     output
  //--------------------------------------------------------------------
  void ATC_CouplingEnergy::output()
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

      // auxiliary data
      for (_tiIt_ = timeIntegrators_.begin(); _tiIt_ != timeIntegrators_.end(); ++_tiIt_) {
        (_tiIt_->second)->output(outputData);
      }
      atomicRegulator_->output(outputData);
      extrinsicModelManager_.output(outputData);
      
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
        outputData["NodalAtomicTemperature"] = &temperature;
        outputData["dot_temperature"] = &dotTemperature;
        outputData["ddot_temperature"] = &ddotTemperature;
        outputData["NodalAtomicPower"] = &rocTemperature;
        outputData["fePower"] = &fePower;

        // write data
        feEngine_->write_data(output_index(), fields_, & outputData);
      }
      feEngine_->partition_mesh();
    }
  }
};
