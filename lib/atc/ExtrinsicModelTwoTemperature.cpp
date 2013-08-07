// ATC Headers
#include "ExtrinsicModelTwoTemperature.h"
#include "ATC_Error.h"
#include "FieldEulerIntegrator.h"
#include "ATC_Coupling.h"
#include "LammpsInterface.h"
#include "PrescribedDataManager.h"
#include "PhysicsModel.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModelTwoTemperature
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ExtrinsicModelTwoTemperature::ExtrinsicModelTwoTemperature
                                (ExtrinsicModelManager * modelManager,
                                 ExtrinsicModelType modelType,
                                 string matFileName) :
    ExtrinsicModel(modelManager,modelType,matFileName),
    electronTimeIntegration_(TimeIntegrator::IMPLICIT),
    temperatureIntegrator_(NULL),
    nsubcycle_(1),
    exchangeFlag_(true), 
    baseSize_(0)
  {
     physicsModel_ = new PhysicsModelTwoTemperature(matFileName);

     // set up correct masks for coupling
     rhsMaskIntrinsic_.reset(NUM_FIELDS,NUM_FLUX);
     rhsMaskIntrinsic_ = false;
     rhsMaskIntrinsic_(TEMPERATURE,SOURCE) = true;
     atc_->fieldMask_(TEMPERATURE,EXTRINSIC_SOURCE) = true;
  }
  
  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ExtrinsicModelTwoTemperature::~ExtrinsicModelTwoTemperature()
  {
    if (temperatureIntegrator_) delete temperatureIntegrator_;
  }

  //--------------------------------------------------------
  // modify 
  //--------------------------------------------------------
  bool ExtrinsicModelTwoTemperature::modify(int narg, char **arg)
  {
    bool match = false;
    int argIndx = 0;

    // energy exchange switch
    /*! \page man_extrinsic_exchange fix_modify AtC extrinsic exchange
      \section syntax
      fix_modify AtC extrinsic exchange <on|off>

      \section examples
      <TT> fix_modify AtC extrinsic exchange on </TT> \n

      \section description
      Switches energy exchange between the MD system and electron system on and off

      \section restrictions
      Only valid for use with two_temperature type of AtC fix.

      \section related
      see \ref man_fix_atc

      \section default
      on
     */
    if (strcmp(arg[argIndx],"exchange")==0) {
      argIndx++;
      if (strcmp(arg[argIndx],"off")==0) {
        exchangeFlag_ = false;
        rhsMaskIntrinsic_(TEMPERATURE,SOURCE) = false;
        atc_->fieldMask_(ELECTRON_TEMPERATURE,SOURCE) = false;
        atc_->fieldMask_(TEMPERATURE,EXTRINSIC_SOURCE) = false;
      }
      else {
        exchangeFlag_ = true;
        rhsMaskIntrinsic_(TEMPERATURE,SOURCE) = true;
        atc_->fieldMask_(ELECTRON_TEMPERATURE,SOURCE) = true;
        atc_->fieldMask_(TEMPERATURE,EXTRINSIC_SOURCE) = true;
      }
      match = true;
    } // end "exchange"

    // electron integration type
    /*! \page man_electron_integration fix_modify AtC extrinsic electron_integration
      \section syntax
      fix_modify AtC extrinsic electron_integration <integration_type> <num_subcyle_steps(optional)> \n
        - integration_type (string) = explicit | implicit | steady  \n
        - num_subcycle_steps (int), optional = number of subcycle steps for the electron time integration

      \section examples
      <TT> fix_modify AtC extrinsic electron_integration implicit </TT> \n
      <TT> fix_modify AtC extrinsic electron_integration explicit 100 </TT> \n

      \section description
      Switches between integration scheme for the electron temperature.  The number of subcyling steps used to integrate the electron temperature 1 LAMMPS timestep can be manually adjusted to capture fast electron dynamics.

      \section restrictions
      For use only with two_temperature type of AtC fix ( see \ref man_fix_atc ) \n
      \section default
      implicit\n
      subcycle_steps = 1
     */
    else if (strcmp(arg[argIndx],"electron_integration")==0) {
      argIndx++;
      nsubcycle_ = 1;
      if (strcmp(arg[argIndx],"explicit")==0) {
        electronTimeIntegration_ = TimeIntegrator::EXPLICIT;
        match = true;
      }
      else if (strcmp(arg[argIndx],"implicit")==0) {
        electronTimeIntegration_ = TimeIntegrator::IMPLICIT;
        match = true;
      }
      else if (strcmp(arg[argIndx],"steady")==0) {
        electronTimeIntegration_ = TimeIntegrator::STEADY;
        match = true;
      }
      if (narg > ++argIndx) nsubcycle_ = atoi(arg[argIndx]);
    } // end "electron_integration"

    if (!match) {
      match = ExtrinsicModel::modify(narg, arg);
    }

    return match;
  }

  //--------------------------------------------------------
  //  initialize
  //--------------------------------------------------------
  void ExtrinsicModelTwoTemperature::initialize()
  {
    ExtrinsicModel::initialize();

    int nNodes = atc_->num_nodes();
    rhs_[TEMPERATURE].reset(nNodes,1);
    rhs_[ELECTRON_TEMPERATURE].reset(nNodes,1);

    // set up electron temperature integrator
    Array2D <bool> rhsMask(NUM_FIELDS,NUM_FLUX);
    rhsMask = false;
    for (int i = 0; i < NUM_FLUX; i++) {
      rhsMask(ELECTRON_TEMPERATURE,i) = atc_->fieldMask_(ELECTRON_TEMPERATURE,i);
    }
    if (temperatureIntegrator_) delete temperatureIntegrator_;
    if (electronTimeIntegration_ == TimeIntegrator::STEADY) {
      throw ATC_Error("not implemented");
    }
    else if (electronTimeIntegration_ == TimeIntegrator::IMPLICIT) {
      double alpha = 1; // backwards Euler
      temperatureIntegrator_ = new FieldImplicitEulerIntegrator(
        ELECTRON_TEMPERATURE, physicsModel_, atc_->feEngine_, atc_, 
        rhsMask, alpha);
    }
    else {
      temperatureIntegrator_ = new FieldExplicitEulerIntegrator(
        ELECTRON_TEMPERATURE, physicsModel_, atc_->feEngine_, atc_, 
        rhsMask);
    }


    // set up mass matrix
    Array<FieldName> massMask(1);
    massMask = ELECTRON_TEMPERATURE;
    (atc_->feEngine_)->compute_lumped_mass_matrix(massMask,atc_->fields_,physicsModel_,atc_->elementToMaterialMap_,atc_->massMats_);
    atc_->massMatsInv_[ELECTRON_TEMPERATURE] = inv(atc_->massMats_[ELECTRON_TEMPERATURE].quantity());
  }

  //--------------------------------------------------------
  //  pre initial integration
  //--------------------------------------------------------
  void ExtrinsicModelTwoTemperature::pre_init_integrate()
  {
    double dt = atc_->lammpsInterface_->dt();
    double time = atc_->time();

    // integrate fast electron variable/s
    // note: atc calls set_sources in pre_final_integrate
    atc_->set_fixed_nodes();
    double idt = dt/nsubcycle_;
    for (int i = 0; i < nsubcycle_ ; ++i) {
      temperatureIntegrator_->update(idt,time,atc_->fields_,rhs_);
    }

  }

  //--------------------------------------------------------
  //  set coupling source terms
  //--------------------------------------------------------
  void ExtrinsicModelTwoTemperature::set_sources(FIELDS & fields, FIELDS & sources)
  {
    // compute source term with appropriate masking and physics model
    atc_->evaluate_rhs_integral(rhsMaskIntrinsic_, fields,
                            sources,
                            atc_->source_integration(), physicsModel_);
  }

  //--------------------------------------------------------
  //  output
  //--------------------------------------------------------
  void ExtrinsicModelTwoTemperature::output(OUTPUT_LIST & outputData)
  {
    // nodal data
    outputData["dot_electron_temperature"] = & rhs_[ELECTRON_TEMPERATURE].set_quantity();

    // global data
    if (atc_->lammpsInterface_->rank_zero()) {
      double T_mean   = ((atc_->field(ELECTRON_TEMPERATURE)).quantity()).col_sum(0)/atc_->nNodes_;
      atc_->feEngine_->add_global("electron_temperature_mean",  T_mean);
      double T_stddev = ((atc_->field(ELECTRON_TEMPERATURE)).quantity()).col_stdev(0);
      atc_->feEngine_->add_global("electron_temperature_std_dev",  T_stddev);
    }
  }

  //--------------------------------------------------------
  //  size_vector
  //--------------------------------------------------------
  int ExtrinsicModelTwoTemperature::size_vector(int intrinsicSize)
  {
    baseSize_ = intrinsicSize;
    return 2;
  }

  //--------------------------------------------------------
  //  compute_vector
  //--------------------------------------------------------
  bool ExtrinsicModelTwoTemperature::compute_vector(int n, double & value)
  {
    // output[1] = total electron energy
    // output[2] = average electron temperature

    if (n == baseSize_) { 
      Array<FieldName> mask(1);
      FIELD_MATS energy;
      mask(0) = ELECTRON_TEMPERATURE;
      
      (atc_->feEngine_)->compute_energy(mask, 
                                          atc_->fields(),
                                          physicsModel_,
                                          atc_->elementToMaterialMap_,
                                          energy);
      // convert to lammps energy units
      double mvv2e = (atc_->lammps_interface())->mvv2e(); 
      double electronEnergy = mvv2e * energy[ELECTRON_TEMPERATURE].col_sum();
      value = electronEnergy;
      return true;
    }
    else if (n == baseSize_+1) {
      double electronTemperature = ((atc_->field(ELECTRON_TEMPERATURE)).quantity()).col_sum()/(atc_->nNodes_);
      value = electronTemperature;
      return true;
    }
    
    return false;
  }

};
