// ATC_Transfer headers
#include "ATC_TransferThermal.h"
#include "ATC_Error.h"
#include "LammpsInterface.h"
#include "PrescribedDataManager.h"

// Other Headers
#include <vector>
#include <set>
#include <utility>

namespace ATC {

  ATC_TransferThermal::ATC_TransferThermal(string groupName,
                                           string matParamFile,
                                           ExtrinsicModelType extrinsicModel)
    : ATC_Transfer(),
      thermostat_(this),
      pmfcOn_(false)
  {
    // assign default "internal" group
    int igroup = lammpsInterface_->find_group(groupName.c_str());
    groupbit_ |= lammpsInterface_->group_bit(igroup);
    igroups_.insert(igroup);

    // Allocate PhysicsModel 
    create_physics_model(THERMAL, matParamFile);

    // create extrinsic physics model
    if (extrinsicModel != NO_MODEL) {
      extrinsicModelManager_.create_model(extrinsicModel,matParamFile);  
    }

    simTime_ = 0.;
    integrationType_ = TimeIntegrator::VERLET;
    atomicOutputMask_.insert("temperature");
  
    // set up field data based on physicsModel
    physicsModel_->get_num_fields(fieldSizes_,fieldMask_);

    // output variable vector info:
    // output[1] = total coarse scale thermal energy
    // output[2] = average temperature
    vectorFlag_ = 1;
    sizeVector_ = 2;
    globalFreq_ = 1;
    extVector_ = 1;
    if (extrinsicModel != NO_MODEL)
      sizeVector_ += extrinsicModelManager_.size_vector(sizeVector_);
  }

  ATC_TransferThermal::~ATC_TransferThermal()
  {
  }

  void ATC_TransferThermal::initialize()
  {
    // Base class initalizations
    ATC_Transfer::initialize();
    
    if (!timeFilterManager_.filter_dynamics()) {
      DENS_VEC atomicKineticEnergy(nLocal_);
      compute_atomic_kinetic_energy(atomicKineticEnergy, lammpsInterface_->vatom());
      project_volumetric_quantity(atomicKineticEnergy,fieldNdFiltered_[TEMPERATURE],TEMPERATURE);
      fieldNdFiltered_[TEMPERATURE] *= 2.;
    }
   
    thermostat_.initialize();
    extrinsicModelManager_.initialize();

    if (!initialized_) {
      // initialize sources based on initial FE temperature
      double dt = lammpsInterface_->dt();
      prescribedDataMgr_->set_sources(simTime_+.5*dt,sources_);
      extrinsicModelManager_.set_sources(fields_,extrinsicSources_);
      thermostat_.compute_boundary_flux(fields_);
      compute_atomic_sources(fieldMask_,fields_,atomicSources_);
    }

    if (timeFilterManager_.need_reset()) {
      init_filter();
      timeFilterManager_.initialize();
    }

    // reset integration field mask
    temperatureMask_.reset(NUM_FIELDS,NUM_FLUX);
    temperatureMask_ = false;
    for (int i = 0; i < NUM_FLUX; i++)
      temperatureMask_(TEMPERATURE,i) = fieldMask_(TEMPERATURE,i);
    initialized_ = true;
    

    if (pmfcOn_) {
      oldFieldTemp_.reset(nNodes_,1);
    }

    // read in field data if necessary
    if (useRestart_) {
      OUTPUT_LIST data;
      read_restart_data(restartFileName_,data);
      useRestart_ = false;
    }
  }

  void ATC_TransferThermal::init_filter()
  {
    // NOTE assume total filtered time derivatives are zero
    ATC_Transfer::init_filter();

    if (integrationType_==TimeIntegrator::VERLET) {
      // initialize restricted fields
      // NOTE: comment in to initialize fieldRateNdOld_ to current time deriviative
      // current assumption is that it is zero
      //DenseMatrix<double> atomicPower(nLocal_,1);
      //compute_atomic_power(atomicPower,
      //                   lammpsInterface_->vatom(),
      //                   lammpsInterface_->fatom());
      //restrict(atomicPower,fieldRateNdOld_[TEMPERATURE]);
      //fieldRateNdOld *= 2.;
      
      DENS_MAT atomicKineticEnergy(nLocal_,1);
      DENS_MAT nodalAtomicTemperature(nNodes_,1);
      compute_atomic_temperature(atomicKineticEnergy,
                                 lammpsInterface_->vatom());
      restrict(atomicKineticEnergy, nodalAtomicTemperature);
      nodalAtomicTemperature *= 2.;
      fieldNdOld_[TEMPERATURE] = nodalAtomicTemperature;
    }
    
    if (timeFilterManager_.end_equilibrate()) { // set up correct initial lambda power to enforce initial temperature rate of 0
      if (integrationType_==TimeIntegrator::VERLET) {
        if (equilibriumStart_) {
          if (thermostat_.get_thermostat_type()==Thermostat::FLUX) { // based on FE equation
            DENS_MAT vdotflamMat(-2.*fieldRateNdFiltered_[TEMPERATURE]); // note 2 is for 1/2 vdotflam addition
            thermostat_.reset_lambda_power(vdotflamMat);
          }
          else { // based on MD temperature equation
            DENS_MAT vdotflamMat(-1.*fieldRateNdFiltered_[TEMPERATURE]);
            thermostat_.reset_lambda_power(vdotflamMat);
          }
        }
      }
    }
    else { // set up space for filtered atomic power
      // should set up all data structures common to equilibration and filtering,
      // specifically filtered atomic power
      if (integrationType_==TimeIntegrator::FRACTIONAL_STEP) {
        dot_atomicTemp.reset(nNodes_,1);
        dot_atomicTempOld.reset(nNodes_,1);
        dot_dot_atomicTemp.reset(nNodes_,1);
        dot_dot_atomicTempOld.reset(nNodes_,1);
      }
    }
  }

  void ATC_TransferThermal::compute_md_mass_matrix(FieldName thisField,
                                                   map<FieldName,DIAG_MAT> & massMats)
  {
    if (thisField == TEMPERATURE) {
      DENS_VEC atomicCapacity(nLocal_);
      atomicCapacity = nsd_*(LammpsInterface::instance()->kBoltzmann());
      DENS_VEC nodalAtomicCapacity(nNodes_);
      restrict_volumetric_quantity(atomicCapacity,nodalAtomicCapacity);
      massMats[thisField].reset(nodalAtomicCapacity);
    }
  }
        

  void ATC_TransferThermal::finish()
  {
    // base class
    ATC_Transfer::finish();

    // nothing specific to thermal
  }

  bool ATC_TransferThermal::modify(int narg, char **arg)
  {
    bool foundMatch = false;
    int argIndx = 0;
    
    // check to see if input is a transfer class command
    // check derived class before base class
    if (strcmp(arg[argIndx],"transfer")==0) {
      argIndx++;

      // pass-through to thermostat
      if (strcmp(arg[argIndx],"thermal")==0) {
        argIndx++;
        if (strcmp(arg[argIndx],"control")==0) {
          argIndx++;
          foundMatch = thermostat_.modify(narg-argIndx,&arg[argIndx]);
        }
      }

      // switch for equilibrium filtering start
      /*! \page man_equilibrium_start fix_modify AtC transfer equilibrium_start
        \section syntax
	fix_modify AtC transfer equilibrium_start <on|off>

	\section examples
	<TT> fix_modify atc transfer equilibrium_start on </TT> \n

	\section description
	Starts filtered calculations assuming they start in equilibrium, i.e. perfect finite element force balance.

	\section restrictions
	only needed before filtering is begun

	\section related
        see \ref man_time_filter

        \section default
        on
      */
      else if (strcmp(arg[argIndx],"equilibrium_start")==0) {
        argIndx++;
        if (strcmp(arg[argIndx],"on")==0) {
          equilibriumStart_ = true;
          foundMatch = true;
        }
        else if (strcmp(arg[argIndx],"off")==0) {
          equilibriumStart_ = false;
          foundMatch = true;
        }
      }

      // time integration scheme
      /*! \page man_time_integration fix_modify AtC transfer pmfc
	\section syntax
	fix_modify AtC transfer pmfc <on|off>

	\section examples
        <TT> fix_modify atc transfer pmfc on </TT>

	\section description
	Switches the poor man's fractional step algorithm on where the finite element data lags the exact atomic data by one time step for overlap nodes

	\section restrictions
        \section related
        \section default
	off
       */
      else if (strcmp(arg[argIndx],"pmfc")==0) {
	if (integrationType_!=TimeIntegrator::VERLET || timeFilterManager_.filter_dynamics())
	  throw ATC_Error(0,"Can only use poor man's fractional step with Verlet integration without filtering");
	pmfcOn_ = !pmfcOn_;
	foundMatch = true;
      }
    }

    // no match, call base class parser
    if (!foundMatch) {
      foundMatch = ATC_Transfer::modify(narg, arg);
    }

    return foundMatch;

  }

  //--------------------------------------------------
  // pack_fields
  //   bundle all allocated field matrices into a list
  //   for output needs
  //--------------------------------------------------
  void ATC_TransferThermal::pack_thermal_fields(OUTPUT_LIST & data)
  {
    data["ddFieldTemp"] = & oldFieldTemp_;
    data["lambdaPowerFiltered"] = &(thermostat_.get_filtered_lambda_power());
  }
  
  //--------------------------------------------------
  // write_restart_file
  //   bundle matrices that need to be saved and call
  //   fe_engine to write the file
  //--------------------------------------------------
  void ATC_TransferThermal::write_restart_data(string fileName, OUTPUT_LIST & data)
  {
    pack_thermal_fields(data);
    ATC_Transfer::write_restart_data(fileName,data);
  }
    
  //--------------------------------------------------
  // write_restart_file
  //   bundle matrices that need to be saved and call
  //   fe_engine to write the file
  //--------------------------------------------------
  void ATC_TransferThermal::read_restart_data(string fileName, OUTPUT_LIST & data)
  {
    pack_thermal_fields(data);
    ATC_Transfer::read_restart_data(fileName,data);
  }

  void ATC_TransferThermal::reset_nlocal()
  {
    ATC_Transfer::reset_nlocal();
    thermostat_.reset_nlocal();
  }

  void ATC_TransferThermal::pre_init_integrate()
  {
    //ATC_Transfer::pre_init_integrate();
    double dt = lammpsInterface_->dt();
    double dtLambda = 0.5*dt;

    if (pmfcOn_) {
      oldFieldTemp_ = fields_[TEMPERATURE];
    }
    
    // Apply thermostat to atom velocities
    thermostat_.apply_pre_predictor(dtLambda,lammpsInterface_->ntimestep());

    
    // Predict nodal temperatures and time derivatives based on FE data
    // 4th order Gear
    // switch call based on filter
    if (timeFilterManager_.filter_dynamics())
      gear1_4_predict(fields_[TEMPERATURE],dot_fields_[TEMPERATURE],
                      ddot_fields_[TEMPERATURE],dddot_fields_[TEMPERATURE],dt);
    else
      gear1_3_predict(fields_[TEMPERATURE],dot_fields_[TEMPERATURE],
                      ddot_fields_[TEMPERATURE],dt);

    extrinsicModelManager_.pre_init_integrate();
        
    // fixed values, non-group bcs handled through FE
    set_fixed_nodes();
      
    simTime_ += .5*dt;
  }

  void ATC_TransferThermal::post_init_integrate()
  {
    // Nothing to do here
    //ATC_Transfer::post_init_integrate();
  }

  void ATC_TransferThermal::pre_final_integrate()
  {
    ATC_Transfer::pre_final_integrate();
  }


  void ATC_TransferThermal::post_final_integrate()
  {
    //ATC_Transfer::post_final_integrate();
    double dt = lammpsInterface_->dt();
    double dtLambda = dt*0.5;
    // predict thermostat contributions
    // compute sources based on predicted FE temperature
    prescribedDataMgr_->set_sources(simTime_+.5*dt,sources_);
    extrinsicModelManager_.set_sources(fields_,extrinsicSources_);
    thermostat_.compute_boundary_flux(fields_);
    compute_atomic_sources(temperatureMask_,fields_,atomicSources_);
    thermostat_.apply_pre_corrector(dtLambda,lammpsInterface_->ntimestep());

    // Determine FE contributions to d theta/dt    
    // Compute atom-integrated rhs
    // parallel communication happens within FE_Engine
    compute_rhs_vector(temperatureMask_,fields_,rhs_,FE_DOMAIN);

    // For flux matching, add 1/2 of "drag" power
    thermostat_.add_to_rhs(rhs_);

    extrinsicModelManager_.post_final_integrate();
    // add in atomic FE contributions for verlet method
    if (integrationType_==TimeIntegrator::VERLET) {
      DENS_MAT atomicPower(nLocal_,1);
      compute_atomic_power(atomicPower,
                           lammpsInterface_->vatom(),
                           lammpsInterface_->fatom());
      DENS_MAT nodalAtomicPower(nNodes_,1);
      restrict_volumetric_quantity(atomicPower,nodalAtomicPower);
      nodalAtomicPower *= 2.;

      if (timeFilterManager_.filter_variables())
        update_filter_implicit(fieldRateNdFiltered_[TEMPERATURE],nodalAtomicPower,dt);
      else
        fieldRateNdFiltered_[TEMPERATURE] = nodalAtomicPower;
      if (!timeFilterManager_.filter_variables() || timeFilterManager_.filter_dynamics())
        rhs_[TEMPERATURE] += fieldRateNdFiltered_[TEMPERATURE];
      else
        rhs_[TEMPERATURE] += nodalAtomicPower;
    }
 
    // Finish updating temperature
    DENS_MAT R_theta(nNodes_,1);
    apply_inverse_mass_matrix(rhs_[TEMPERATURE],TEMPERATURE);
    R_theta = (rhs_[TEMPERATURE] - dot_fields_[TEMPERATURE])*dt;
    if (timeFilterManager_.filter_dynamics())
      gear1_4_correct(fields_[TEMPERATURE],dot_fields_[TEMPERATURE],
                      ddot_fields_[TEMPERATURE],dddot_fields_[TEMPERATURE],
                      R_theta,dt);
    else
      gear1_3_correct(fields_[TEMPERATURE],dot_fields_[TEMPERATURE],
                      ddot_fields_[TEMPERATURE],R_theta,dt);

    // only requirecd for poor man's fractional step
    if (pmfcOn_) {
      DENS_MAT atomicKineticEnergy(nLocal_,1);
      compute_atomic_kinetic_energy(atomicKineticEnergy,
				    lammpsInterface_->vatom());
      DENS_MAT temperatureNd(nNodes_,1);
      project_md_volumetric_quantity(atomicKineticEnergy,
				     temperatureNd,
				     TEMPERATURE);
      temperatureNd *= 2.;
      dot_fields_[TEMPERATURE] = 1.0/dt * ( temperatureNd - oldFieldTemp_); 
      fields_[TEMPERATURE] = temperatureNd;
    }

    // fix nodes, non-group bcs applied through FE
    set_fixed_nodes();

    // compute sources based on final FE updates
    prescribedDataMgr_->set_sources(simTime_+.5*dt,sources_);
    extrinsicModelManager_.set_sources(fields_,extrinsicSources_);
    thermostat_.compute_boundary_flux(fields_);
    compute_atomic_sources(temperatureMask_,fields_,atomicSources_);

    // apply corrector phase of thermostat
    thermostat_.apply_post_corrector(dtLambda,lammpsInterface_->ntimestep());

    // add in MD contributions to time derivative
    // update filtered temperature
    DENS_MAT atomicKineticEnergy(nLocal_,1);
    compute_atomic_kinetic_energy(atomicKineticEnergy, 
                                  lammpsInterface_->vatom());
    DENS_MAT temperatureNd(nNodes_,1);
    project_md_volumetric_quantity(atomicKineticEnergy,temperatureNd,TEMPERATURE);
    temperatureNd *= 2.;
    if (!timeFilterManager_.filter_dynamics()) 
      fieldNdFiltered_[TEMPERATURE] = temperatureNd;
    else if (integrationType_==TimeIntegrator::VERLET)
      update_filter_implicit(fieldNdFiltered_[TEMPERATURE],temperatureNd,dt);

    simTime_ += .5*dt;
    output();
  }
  
  //--------------------------------------------------------------------
  //     compute_vector
  //--------------------------------------------------------------------
  // this is for direct output to lammps thermo
  double ATC_TransferThermal::compute_vector(int n)
  {
    // output[1] = total coarse scale thermal energy
    // output[2] = average temperature

    double mvv2e = lammpsInterface_->mvv2e(); // convert to lammps energy units
  
    if (n == 0) {
      Array<FieldName> mask(1);
      FIELDS energy;
      mask(0) = TEMPERATURE;
      feEngine_->compute_energy(mask, // NOTE make this a base class function
                                fields_,
                                physicsModel_,
                                elementToMaterialMap_,
                                energy,
                                &elementMask_);
      
      double phononEnergy = mvv2e * energy[TEMPERATURE].col_sum();
      return phononEnergy;
    }
    else if (n == 1) {
      double aveT = fields_[TEMPERATURE].col_sum()/nNodes_;
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
  void ATC_TransferThermal::output()
  {
    double dt = lammpsInterface_->dt();

    double step = (double) lammpsInterface_->ntimestep();
    ++stepCounter_;

    if ((outputFrequency_ > 0) 
        && ((stepCounter_ ==1) || (stepCounter_ % outputFrequency_ == 0)) ) {

      OUTPUT_LIST output_data;

      // Add some outputs only on Proc 0
      if (lammpsInterface_->comm_rank() == 0) {
        // base class output
        ATC_Transfer::output();
        
        // global data
        double T_mean   = fields_[TEMPERATURE].col_sum(0)/nNodes_;
        feEngine_->add_global("temperature_mean",  T_mean);
        double T_stddev   = fields_[TEMPERATURE].col_stdev(0);
        feEngine_->add_global("temperature_std_dev",  T_stddev);
        double Ta_mean =  fieldNdFiltered_[TEMPERATURE].col_sum(0)/nNodes_; 
        feEngine_->add_global("atomic_temperature_mean",  Ta_mean);
        double Ta_stddev =  fieldNdFiltered_[TEMPERATURE].col_stdev(0); 
        feEngine_->add_global("atomic_temperature_std_dev",  Ta_stddev);
        
        // mesh data
        output_data["nodalAtomicTemperature"] = & fieldNdFiltered_[TEMPERATURE];
        output_data["dot_temperature"] = & dot_fields_[TEMPERATURE];
        output_data["ddot_temperature"] = & ddot_fields_[TEMPERATURE];
        output_data["nodalAtomicPower"] = & fieldRateNdFiltered_[TEMPERATURE];

        // auxilliary data
        thermostat_.output(dt,output_data);
        extrinsicModelManager_.output(dt,output_data);
      }


      // Output fe data on proc 0
      if (lammpsInterface_->comm_rank() == 0) {
        feEngine_->write_data(simTime_, fields_, & output_data); 
      }
      
    }

    if ((outputFrequencyAtom_ > 0)
        && ((stepCounter_ ==1) || (stepCounter_ % outputFrequencyAtom_ == 0)) )
      atomic_output();
  } 


  void ATC_TransferThermal::set_ghost_atoms()
  {
    // Nothing to do here
  }

};
