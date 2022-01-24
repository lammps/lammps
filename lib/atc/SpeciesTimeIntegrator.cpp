// ATC transfer headers
#include "SpeciesTimeIntegrator.h"
#include "TransferOperator.h"
#include "ATC_CouplingMass.h"
#include "TimeFilter.h"
#include "LammpsInterface.h"
#include "ATC_Error.h"
#include "PerAtomQuantityLibrary.h"
#include "AtomToMoleculeTransfer.h"
#include "MoleculeSet.h"

using std::pair;
using std::map;
using std::string;

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class SpeciesTimeIntegrator
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //        Grab data from ATC
  //--------------------------------------------------------
  SpeciesTimeIntegrator::SpeciesTimeIntegrator(ATC_CouplingMass * atc,
                                               TimeIntegrationType timeIntegrationType) :
    TimeIntegrator(atc, timeIntegrationType),
    moleculeIds_(atc->molecule_ids())
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  modify
  //    parses inputs and modifies state of the filter
  //--------------------------------------------------------
  bool SpeciesTimeIntegrator::modify(int /* narg */, char ** /* arg */)
  {
    bool match = false;

    // no parsing needed

    return match;
  }

  //--------------------------------------------------------
  //  construct_methods
  //    creates algorithm objects
  //--------------------------------------------------------
  void SpeciesTimeIntegrator::construct_methods()
  {
    if (atc_->reset_methods()) {
      if (timeIntegrationMethod_) delete timeIntegrationMethod_;
      if (timeFilterManager_->need_reset()) {
        switch (timeIntegrationType_) {
          case FRACTIONAL_STEP: {
            timeFilter_ = timeFilterManager_->construct(TimeFilterManager::EXPLICIT_IMPLICIT);
            atc_->set_mass_mat_time_filter(SPECIES_CONCENTRATION,TimeFilterManager::EXPLICIT);
            break;
          }
          default:
            throw ATC_Error("Unknown time integration type in SpeciesTimeIntegrator::Initialize()");
        }
      }
      if (timeFilterManager_->filter_dynamics()) {
        switch (timeIntegrationType_) {
          case FRACTIONAL_STEP: {
            timeIntegrationMethod_ = new SpeciesTimeIntegratorFractionalStepFiltered(this,
                                                                                     moleculeIds_);
           }
          default:
            throw ATC_Error("Unknown time integration type in SpeciesTimeIntegrator::Initialize()");
        }
      }
      else {
        timeIntegrationMethod_ = new SpeciesTimeIntegratorFractionalStep(this,
                                                                         moleculeIds_);
      }
    }
  }

  //--------------------------------------------------------
  //  pack_fields
  //    add persistent variables to data list
  //--------------------------------------------------------
  void SpeciesTimeIntegrator::pack_fields(RESTART_LIST & data)
  {

    TimeIntegrator::pack_fields(data);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class SpeciesIntegrationMethod
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //        Grab data from ATC
  //--------------------------------------------------------
  SpeciesIntegrationMethod::SpeciesIntegrationMethod(SpeciesTimeIntegrator * speciesTimeIntegrator,
    const map<string,pair<MolSize,int> > & moleculeIds) :
    TimeIntegrationMethod(speciesTimeIntegrator),
    timeFilter_(speciesTimeIntegrator->time_filter()),
    massDensity_(atc_->field(MASS_DENSITY)),
    nodalAtomicMassDensityOut_(atc_->nodal_atomic_field(MASS_DENSITY)),
    nodalAtomicMassDensity_(nullptr),
    speciesConcentration_(atc_->field(SPECIES_CONCENTRATION)),
    nodalAtomicSpeciesConcentration_(nullptr),
    nodalAtomicSpeciesConcentrationFiltered_(speciesTimeIntegrator->nodal_atomic_species_concentration_filtered()),
    moleculeIds_(moleculeIds)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  construct_transfers
  //    sets up all the necessary transfer operators
  //--------------------------------------------------------
  void SpeciesIntegrationMethod::construct_transfers()
  {
    InterscaleManager & interscaleManager = atc_->interscale_manager();

    // get existing data
    nodalAtomicMassDensity_ = interscaleManager.dense_matrix(field_to_intrinsic_name(MASS_DENSITY));
    if (atc_->has_tracked_species())
      nodalAtomicSpeciesConcentration_ = interscaleManager.dense_matrix(field_to_intrinsic_name(SPECIES_CONCENTRATION));
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class SpeciesTimeIntegratorFractionalStep
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  SpeciesTimeIntegratorFractionalStep::SpeciesTimeIntegratorFractionalStep(SpeciesTimeIntegrator * speciesTimeIntegrator,
    const map<string,pair<MolSize,int> > & moleculeIds) :
    SpeciesIntegrationMethod(speciesTimeIntegrator,moleculeIds)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  initialize
  //        initialize all data
  //--------------------------------------------------------
  void SpeciesTimeIntegratorFractionalStep::initialize()
  {
    SpeciesIntegrationMethod::initialize();

    TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
    if (timeFilterManager->need_reset()) {
      timeFilter_->initialize();
    }

    if (!timeFilterManager->end_equilibrate()) {
      nodalAtomicSpeciesConcentrationFiltered_ = nodalAtomicSpeciesConcentration_->quantity();
    }


    pre_final_integrate1(0.);
  }

  //--------------------------------------------------------
  //  pre_initial_integrate1
  //--------------------------------------------------------
  void SpeciesTimeIntegratorFractionalStep::pre_initial_integrate1(double dt)
  {
    const DENS_MAT & my(nodalAtomicSpeciesConcentration_->quantity());
    // updated filtered energy using explicit-implicit scheme
    timeFilter_->apply_pre_step1(nodalAtomicSpeciesConcentrationFiltered_.set_quantity(),
                                 my,dt);
  }


  //--------------------------------------------------------
  //  pre_final_integrate1
  //    first time integration computations
  //    before FractionalStep step 2
  //--------------------------------------------------------
  void SpeciesTimeIntegratorFractionalStep::pre_final_integrate1(double /* dt */)
  {
    // Compute MD contribution to FEM equation


    massDensity_ = nodalAtomicMassDensity_->quantity();
    speciesConcentration_ = nodalAtomicSpeciesConcentration_->quantity();
    atc_->set_fixed_nodes();
  }

  //--------------------------------------------------------
  //  post_final_integrate2
  //--------------------------------------------------------
  void SpeciesTimeIntegratorFractionalStep::post_final_integrate2(double dt)
  {

    timeFilter_->apply_post_step1(
      nodalAtomicSpeciesConcentrationFiltered_.set_quantity(),
      nodalAtomicSpeciesConcentration_->quantity(),dt);
    speciesConcentration_ = nodalAtomicSpeciesConcentrationFiltered_.quantity();
  }

  //--------------------------------------------------------
  //  post_process
  //    post processing of variables before output
  //--------------------------------------------------------
  void SpeciesTimeIntegratorFractionalStep::post_process()
  {

    map<string,pair<MolSize,int> >::const_iterator molecule;
    for (molecule = moleculeIds_.begin(); molecule != moleculeIds_.end(); molecule++) {
      DENS_MAN & nodalMoleculeMassDensityOut(atc_->tagged_dens_man(molecule->first));
      DENS_MAN * nodalMoleculeMassDensity((atc_->interscale_manager()).dense_matrix("NodalMoleculeMassDensity"+molecule->first));
      nodalMoleculeMassDensityOut = nodalMoleculeMassDensity->quantity();
    }
  }

 //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class SpeciesTimeIntegratorFractionalStepFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //        Grab data from ATC
  //--------------------------------------------------------

  SpeciesTimeIntegratorFractionalStepFiltered::SpeciesTimeIntegratorFractionalStepFiltered(
    SpeciesTimeIntegrator * speciesTimeIntegrator,
    const map<string,pair<MolSize,int> > & moleculeIds) :
    SpeciesTimeIntegratorFractionalStep(speciesTimeIntegrator,moleculeIds)
  {
    throw ATC_Error("SpeciesTimeIntegratorFractionalStepFiltered work in progress");
    // do nothing
  }
  //--------------------------------------------------------
  //  pre_initial_integrate1
  //--------------------------------------------------------
  void SpeciesTimeIntegratorFractionalStepFiltered::pre_final_integrate1(double /* dt */)
  {
  }

};
