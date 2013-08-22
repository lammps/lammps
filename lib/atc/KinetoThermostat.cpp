#include "KinetoThermostat.h"
#include "ATC_Coupling.h"
#include "ATC_Error.h"
#include "PrescribedDataManager.h"
#include "ElasticTimeIntegrator.h"
#include "ThermalTimeIntegrator.h"
#include "TransferOperator.h"

using namespace std;

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class KinetoThermostat
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  KinetoThermostat::KinetoThermostat(ATC_Coupling * atc,
                                     const string & regulatorPrefix) :
    AtomicRegulator(atc,regulatorPrefix),
    kinetostat_(atc,"Momentum"),
    thermostat_(atc,"Energy")
  {
    // nothing to do
  }

  //--------------------------------------------------------
  //  modify:
  //    parses and adjusts thermostat state based on
  //    user input, in the style of LAMMPS user input
  //--------------------------------------------------------
  bool KinetoThermostat::modify(int narg, char **arg)
  {
    bool foundMatch = false;
    return foundMatch;
  }

  //--------------------------------------------------------
  //  reset_lambda_contribution:
  //    resets the thermostat generated power to a
  //    prescribed value
  //--------------------------------------------------------
  void KinetoThermostat::reset_lambda_contribution(const DENS_MAT & target,
                                                   const FieldName field)
  {
    //TEMP_JAT
    if (field==VELOCITY) {
      kinetostat_.reset_lambda_contribution(target);
    }
    else if (field == TEMPERATURE) {
      thermostat_.reset_lambda_contribution(target);
    }
    else {
      throw ATC_Error("KinetoThermostat::reset_lambda_contribution - invalid field given");
    }
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
  void KinetoThermostat::construct_methods()
  {
// #ifdef WIP_JAT
//     // get data associated with stages 1 & 2 of ATC_Method::initialize
//     AtomicRegulator::construct_methods();

//     if (atc_->reset_methods()) {
//       // eliminate existing methods
//       delete_method();

//       // update time filter
//       TimeIntegrator::TimeIntegrationType myIntegrationType = (atc_->time_integrator(TEMPERATURE))->time_integration_type();
//       TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
//       if (timeFilterManager->need_reset() ) {
//         if (myIntegrationType == TimeIntegrator::GEAR)
//          timeFilter_ = timeFilterManager->construct(TimeFilterManager::EXPLICIT);
//         else if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP)
//          timeFilter_ = timeFilterManager->construct(TimeFilterManager::EXPLICIT_IMPLICIT);
//       }
      
//       if (timeFilterManager->filter_dynamics()) {
//         switch (regulatorTarget_) {
//         case NONE: {
//           regulatorMethod_ = new RegulatorMethod(this);
//           break;
//         }
//         case FIELD: { // error check, rescale and filtering not supported together
//           throw ATC_Error("Cannot use rescaling thermostat with time filtering");
//           break;
//         }
//         case DYNAMICS: {
//           switch (couplingMode_) {
//           case FIXED: {
//             if (use_lumped_lambda_solve()) {
//               throw ATC_Error("Thermostat:construct_methods - lumped lambda solve cannot be used with Hoover thermostats");
//             }
//             if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP) {
//               if (md_flux_nodes(TEMPERATURE)) {
//                 if (!md_fixed_nodes(TEMPERATURE) && (boundaryIntegrationType_ == NO_QUADRATURE)) {
//                   // there are fluxes but no fixed or coupled nodes
//                   regulatorMethod_ = new ThermostatFluxFiltered(this);
//                 }
//                 else {
//                   // there are both fixed and flux nodes
//                   regulatorMethod_ = new ThermostatFluxFixedFiltered(this);
//                 }
//               }
//               else {
//                 // there are only fixed nodes
//                 regulatorMethod_ = new ThermostatFixedFiltered(this);
//               }
//             }
//             else {
//               regulatorMethod_ = new ThermostatHooverVerletFiltered(this);
//             }
//             break;
//           }
//           case FLUX: {
//             if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP) {
//               if (use_lumped_lambda_solve()) {
//                 throw ATC_Error("Thermostat:construct_methods - lumped lambda solve has been depricated for fractional step thermostats");
//               }
//               if (md_fixed_nodes(TEMPERATURE)) {
//                 if (!md_flux_nodes(TEMPERATURE) && (boundaryIntegrationType_ == NO_QUADRATURE)) {
//                 // there are fixed nodes but no fluxes
//                 regulatorMethod_ = new ThermostatFixedFiltered(this);
//                 }
//                 else {
//                   // there are both fixed and flux nodes
//                   regulatorMethod_ = new ThermostatFluxFixedFiltered(this);
//                 }
//               }
//               else {
//                 // there are only flux nodes
//                 regulatorMethod_ = new ThermostatFluxFiltered(this);
//               }
//             }
//             else {
//               if (use_localized_lambda()) {
//                 if (!((atc_->prescribed_data_manager())->no_fluxes(TEMPERATURE)) &&
//                     atc_->boundary_integration_type() != NO_QUADRATURE) {
//                   throw ATC_Error("Cannot use flux coupling with localized lambda");
//                 }
//               }
//               regulatorMethod_ = new ThermostatPowerVerletFiltered(this);
//             }
//             break;
//           }
//           default:
//             throw ATC_Error("Unknown coupling mode in Thermostat::initialize");
//           }
//           break;
//         }
//         default:
//           throw ATC_Error("Unknown thermostat type in Thermostat::initialize");
//         }
//       }
//       else {
//         switch (regulatorTarget_) {
//         case NONE: {
//           regulatorMethod_ = new RegulatorMethod(this);
//           break;
//         }
//         case FIELD: {
//           if (atc_->temperature_def()==KINETIC)
//             regulatorMethod_ = new ThermostatRescale(this);
//           else if (atc_->temperature_def()==TOTAL)
//             regulatorMethod_ = new ThermostatRescaleMixedKePe(this);
//           else
//             throw ATC_Error("Unknown temperature definition");
//           break;
//         }
//         case DYNAMICS: {
//           switch (couplingMode_) {
//           case FIXED: {
//             if (use_lumped_lambda_solve()) {
//               throw ATC_Error("Thermostat:construct_methods - lumped lambda solve cannot be used with Hoover thermostats");
//             }
//             if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP) {
//               if (md_flux_nodes(TEMPERATURE)) {
//                 if (!md_fixed_nodes(TEMPERATURE) && (boundaryIntegrationType_ == NO_QUADRATURE)) {
//                   // there are fluxes but no fixed or coupled nodes
//                   regulatorMethod_ = new ThermostatFlux(this);
//                 }
//                 else {
//                   // there are both fixed and flux nodes
//                   regulatorMethod_ = new ThermostatFluxFixed(this);
//                 }
//               }
//               else {
//                 // there are only fixed nodes
//                 regulatorMethod_ = new ThermostatFixed(this);
//               }
//             }
//             else {
//               regulatorMethod_ = new ThermostatHooverVerlet(this);
//             }
//             break;
//           }
//           case FLUX: {
//             if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP) {
//               if (use_lumped_lambda_solve()) {
//                 throw ATC_Error("Thermostat:construct_methods - lumped lambda solve has been depricated for fractional step thermostats");
//               }
//               if (md_fixed_nodes(TEMPERATURE)) {
//                 if (!md_flux_nodes(TEMPERATURE) && (boundaryIntegrationType_ == NO_QUADRATURE)) {
//                   // there are fixed nodes but no fluxes
//                   regulatorMethod_ = new ThermostatFixed(this);
//                 }
//                 else {
//                   // there are both fixed and flux nodes
//                   regulatorMethod_ = new ThermostatFluxFixed(this);
//                 }
//               }
//               else {
//                 // there are only flux nodes
//                 regulatorMethod_ = new ThermostatFlux(this);
//               }
//             }
//             else {
//               if (use_localized_lambda()) {
//                 if (!((atc_->prescribed_data_manager())->no_fluxes(TEMPERATURE)) &&
//                     atc_->boundary_integration_type() != NO_QUADRATURE) {
//                   throw ATC_Error("Cannot use flux coupling with localized lambda");
//                 }
//               }
//               regulatorMethod_ = new ThermostatPowerVerlet(this);
//             }
//             break;
//           }
//           default:
//             throw ATC_Error("Unknown coupling mode in Thermostat::initialize");
//           }
//           break;
//         }
//         default:
//           throw ATC_Error("Unknown thermostat target in Thermostat::initialize");
//         }
//       }
      
//       AtomicRegulator::reset_method();
//     }
//     else {
//       set_all_data_to_used();
//     }
// #endif
    //TEMP_JAT
    kinetostat_.construct_methods();
    thermostat_.construct_methods();
  }

  //TEMP_JAT all functions which have explicit thermo and kinetostats need to be removed
    //--------------------------------------------------------
  //  reset_nlocal:
  //    resizes lambda force if necessary
  //--------------------------------------------------------
  void KinetoThermostat::reset_nlocal()
  {
    kinetostat_.reset_nlocal();
    thermostat_.reset_nlocal();
  }
  //--------------------------------------------------------
  //  reset_atom_materials:
  //    resets the localized atom to material map
  //--------------------------------------------------------
  void KinetoThermostat::reset_atom_materials(const Array<int> & elementToMaterialMap,
                                             const MatrixDependencyManager<DenseMatrix, int> * atomElement)
  {
    kinetostat_.reset_atom_materials(elementToMaterialMap,
                                     atomElement);
    thermostat_.reset_atom_materials(elementToMaterialMap,
                                     atomElement);
  }
  //--------------------------------------------------------
  //  construct_transfers:
  //    pass through to appropriate transfer constuctors
  //--------------------------------------------------------
  void KinetoThermostat::construct_transfers()
  {
    kinetostat_.construct_transfers();
    thermostat_.construct_transfers();
  }
  //--------------------------------------------------------
  //  initialize:
  //    sets up methods before a run
  //--------------------------------------------------------
  void KinetoThermostat::initialize()
  {
    kinetostat_.initialize();
    thermostat_.initialize();
    needReset_ = false;
  }
  //--------------------------------------------------------
  //  output:
  //    pass through to appropriate output methods
  //--------------------------------------------------------
  void KinetoThermostat::output(OUTPUT_LIST & outputData) const
  {
    kinetostat_.output(outputData);
    thermostat_.output(outputData);
  }
  //--------------------------------------------------------
  //  finish:
  //    pass through to appropriate end-of-run methods
  //--------------------------------------------------------
  void KinetoThermostat::finish()
  {
    kinetostat_.finish();
    thermostat_.finish();
  }
    //--------------------------------------------------
  // pack_fields
  //   bundle all allocated field matrices into a list
  //   for output needs
  //--------------------------------------------------
  void KinetoThermostat::pack_fields(RESTART_LIST & data)
  {
    kinetostat_.pack_fields(data);
    thermostat_.pack_fields(data);
  }
  //--------------------------------------------------------
  //  compute_boundary_flux:
  //    computes the boundary flux to be consistent with
  //    the controller
  //--------------------------------------------------------
  void KinetoThermostat::compute_boundary_flux(FIELDS & fields)
  {
    kinetostat_.compute_boundary_flux(fields);
    thermostat_.compute_boundary_flux(fields);
  }
  //--------------------------------------------------------
  //  add_to_rhs:
  //    adds any controller contributions to the FE rhs 
  //--------------------------------------------------------
  void KinetoThermostat::add_to_rhs(FIELDS & rhs)
  {
    thermostat_.add_to_rhs(rhs);
    kinetostat_.add_to_rhs(rhs);
  }
  //--------------------------------------------------------
  //  apply_pre_predictor:
  //    applies the controller in the pre-predictor
  //    phase of the time integrator  
  //--------------------------------------------------------
  void KinetoThermostat::apply_pre_predictor(double dt, int timeStep)
  {
    thermostat_.apply_pre_predictor(dt,timeStep);
    kinetostat_.apply_pre_predictor(dt,timeStep);
  }
  //--------------------------------------------------------
  //  apply_mid_predictor:
  //    applies the controller in the mid-predictor
  //    phase of the time integrator  
  //--------------------------------------------------------
  void KinetoThermostat::apply_mid_predictor(double dt, int timeStep)
  {
    thermostat_.apply_mid_predictor(dt,timeStep);
    kinetostat_.apply_mid_predictor(dt,timeStep);
  }
  //--------------------------------------------------------
  //  apply_post_predictor:
  //    applies the controller in the post-predictor
  //    phase of the time integrator  
  //--------------------------------------------------------
  void KinetoThermostat::apply_post_predictor(double dt, int timeStep)
  {
    thermostat_.apply_post_predictor(dt,timeStep);
    kinetostat_.apply_post_predictor(dt,timeStep);
  }
  //--------------------------------------------------------
  //  apply_pre_corrector:
  //    applies the controller in the pre-corrector phase
  //    of the time integrator
  //--------------------------------------------------------
  void KinetoThermostat::apply_pre_corrector(double dt, int timeStep)
  {
    thermostat_.apply_pre_corrector(dt,timeStep);
    kinetostat_.apply_pre_corrector(dt,timeStep);
  }
  //--------------------------------------------------------
  //  apply_post_corrector:
  //    applies the controller in the post-corrector phase
  //    of the time integrator
  //--------------------------------------------------------
  void KinetoThermostat::apply_post_corrector(double dt, int timeStep)
  {
    thermostat_.apply_post_corrector(dt,timeStep);
    kinetostat_.apply_post_corrector(dt,timeStep);
  }
  //--------------------------------------------------------
  //  pre_exchange
  //--------------------------------------------------------
  void KinetoThermostat::pre_exchange()
  {
    thermostat_.pre_exchange();
    kinetostat_.pre_exchange();
  }
  //--------------------------------------------------------
  //  pre_force
  //--------------------------------------------------------
  AtomicRegulator::RegulatorCouplingType KinetoThermostat::coupling_mode(const FieldName field) const
  {
    //TEMP_JAT
    if (field==VELOCITY) {
      return kinetostat_.coupling_mode();
    }
    else if (field == TEMPERATURE) {
      return thermostat_.coupling_mode();
    }
    else {
      throw ATC_Error("KinetoThermostat::coupling_mode - invalid field given");
    }
  }

};
