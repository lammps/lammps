// ATC_Transfer Headers
#include "ATC_Error.h"
#include "Kinetostat.h"
#include "CG.h"
#include "ATC_Transfer.h"
#include "LammpsInterface.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class Kinetostat
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  Kinetostat::Kinetostat(ATC_Transfer * atcTransfer) :
    AtomicRegulator(atcTransfer),
    kinetostatType_(NONE),
    couplingMode_(UNCOUPLED)
  {
    // do nothing
  }
  
  //--------------------------------------------------------
  //  modify:
  //    parses and adjusts kinetostat state based on
  //    user input, in the style of LAMMPS user input
  //--------------------------------------------------------
  bool Kinetostat::modify(int narg, char **arg)
  {
    bool foundMatch = false;

    // fluxstat type
    /*! \page man_disp_control fix_modify AtC transfer momentum control
      \section syntax
      fix_modify AtC transfer momentum control none \n
      
      fix_modify AtC transfer momentum control rescale <frequency>\n
      - frequency (int) = time step frequency for applying displacement and velocity rescaling \n
      
      fix_modify AtC transfer momentum control glc_displacement \n
  
      fix_modify AtC transfer momentum control glc_velocity \n

      fix_modify AtC transfer momentum control flux [faceset face_set_id, interpolate] 
      - face_set_id (string) = id of boundary face set, if not specified
      (or not possible when the atomic domain does not line up with 
      mesh boundaries) defaults to an atomic-quadrature approximate 
      evaulation\n
      \section examples
      fix_modify AtC transfer momentum control glc_velocity \n
      fix_modify AtC transfer momentum control stress_flux faceset bndy_faces \n
      \section description
      \section restrictions
      only for be used with specific transfers :
      elastic \n
      rescale not valid with time filtering activated
      \section related
      \section default
      none
    */
    
    int argIndex = 0;
    if (strcmp(arg[argIndex],"none")==0) { // restore defaults
      kinetostatType_ = NONE;
      couplingMode_ = UNCOUPLED;
      howOften_ = 1;
      boundaryIntegrationType_ = ATC_Transfer::NO_QUADRATURE;
      foundMatch = true;
    }
    else if (strcmp(arg[argIndex],"glc_displacement")==0) {
      kinetostatType_ = GLC_DISPLACEMENT;
      couplingMode_ = FIXED;
      howOften_ = 1;
      boundaryIntegrationType_ = ATC_Transfer::NO_QUADRATURE;
      foundMatch = true;
    }
    else if (strcmp(arg[argIndex],"glc_velocity")==0) {
      kinetostatType_ = GLC_VELOCITY;
      couplingMode_ = FIXED;
      howOften_ = 1;
      boundaryIntegrationType_ = ATC_Transfer::NO_QUADRATURE;
      foundMatch = true;
    }
    else if (strcmp(arg[argIndex],"hoover")==0) {
      kinetostatType_ = FORCE;
      couplingMode_ = FIXED;
      howOften_ = 1;
      boundaryIntegrationType_ = ATC_Transfer::NO_QUADRATURE;
      foundMatch = true;
    }
    else if (strcmp(arg[argIndex],"flux")==0) {
      kinetostatType_ = FORCE;
      couplingMode_ = FLUX;
      howOften_ = 1;
      argIndex++;
      boundaryIntegrationType_ = atcTransfer_->parse_boundary_integration(narg-argIndex,&arg[argIndex],boundaryFaceSet_);
      foundMatch = true;
    }
  
    if (!foundMatch)
      foundMatch = AtomicRegulator::modify(narg,arg);
    if (foundMatch)
      needReset_ = true;
    return foundMatch;
  }

  //--------------------------------------------------------
  //  reset_nodal_lambda_force:
  //    resets the kinetostat generated force to a
  //    prescribed value
  //--------------------------------------------------------
  void Kinetostat::reset_lambda_force(DENS_MAT & target)
  {
    for (int i = 0; i < nNodes_; ++i)
      for (int j = 0; j < nsd_; ++j)
        lambdaForceFiltered_(i,j) = target(i,j);
  }
 
  //--------------------------------------------------------
  //  initialize:
  //    sets up methods before a run
  //    NOTE we currently only include time filter
  //    dependence, but in general there is also a
  //    time integrator dependence.  In general the 
  //    precedence order is:
  //    time filter -> time integrator -> kinetostat
  //    In the future this may need to be added if
  //    different types of time integrators can be
  //    specified.
  //--------------------------------------------------------
  void Kinetostat::initialize()
  {
    // NOTE: only support basic time integration
    TimeFilterManager * timeFilterManager = atcTransfer_->get_time_filter_manager();

    // reset data if needed and perform any error/conflict checking
    if (resetData_) {
      AtomicRegulator::reset_data();

      // set up storage
      lambda_.reset(nNodes_,nsd_);
      nodalAtomicLambdaForce_.reset(nNodes_,nsd_);
      lambdaForceFiltered_.reset(nNodes_,nsd_);
    }
    
    if (needReset_ || timeFilterManager->need_reset() || timeFilterManager->end_equilibrate()) {
      // eliminate existing methods
      destroy();
      
      DENS_MAT nodalGhostForceFiltered;
      if (timeFilterManager->end_equilibrate() && kinetostatType_==FORCE) {
        StressFlux * myMethod;
        myMethod = dynamic_cast<StressFlux *>(regulatorMethod_);
        nodalGhostForceFiltered = myMethod->get_filtered_ghost_force();
      }

      if (timeFilterManager->need_reset()) {
        if (timeFilter_)
          delete timeFilter_;
        timeFilter_ = timeFilterManager->construct(TimeFilterManager::IMPLICIT_UPDATE);
      }
                         
      if (timeFilterManager->filter_dynamics()) {
        switch (kinetostatType_) {
        case NONE: {
          break;
        }
        case GLC_DISPLACEMENT: {
          regulatorMethod_ = new DisplacementGlcFiltered(this);
          break;
        }
        case GLC_VELOCITY: {
          regulatorMethod_ = new VelocityGlcFiltered(this);
          break;
        }
        case FORCE: {
          regulatorMethod_ = new StressFluxFiltered(this);
          if (timeFilterManager->end_equilibrate()) {
            StressFlux * myMethod;
            myMethod = dynamic_cast<StressFlux *>(regulatorMethod_);
            myMethod->reset_filtered_ghost_force(nodalGhostForceFiltered);
          }
          break;
        }
        default:
          throw ATC_Error(0,"Unknown kinetostat type in Kinetostat::initialize");
        }
      }
      else {
        switch (kinetostatType_) {
        case NONE: {
          break;
        }
        case GLC_DISPLACEMENT: {
          regulatorMethod_ = new DisplacementGlc(this);
          break;
        }
        case GLC_VELOCITY: {
          regulatorMethod_ = new VelocityGlc(this);
          break;
        }
        case FORCE: {
          regulatorMethod_ = new StressFlux(this);
          break;
        }
        default:
          throw ATC_Error(0,"Unknown kinetostat type in Kinetostat::initialize");
        }
      }

      AtomicRegulator::reset_method();
    }

     AtomicRegulator::initialize();
  }
 
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class GlcKinetostat
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and kinetostat data
  //--------------------------------------------------------
  GlcKinetostat::GlcKinetostat(Kinetostat *kinetostat) 
    :
    RegulatorShapeFunction(kinetostat),
    kinetostat_(kinetostat),
    timeFilter_(atomicRegulator_->get_time_filter()),
    nodalAtomicLambdaForce_(kinetostat_->get_nodal_atomic_lambda_force()),
    lambdaForceFiltered_(kinetostat_->get_lambda_force_filtered()),
    lambdaForce_(atomicRegulator_->get_lambda_force()),
    mdMassMatrix_(atcTransfer_->get_mass_mat_md(VELOCITY)),
    nodeToOverlapMap_(atcTransfer_->get_node_to_overlap_map())
  {
    // set up list of nodes using Hoover coupling
    // (a) nodes with prescribed values
    Array2D<bool> isFixedNode = atcTransfer_->get_fixed_node_flags(VELOCITY);
    for (int i = 0; i < nNodes_; ++i)
      for (int j = 0; j < isFixedNode.nCols(); ++j)
        if (isFixedNode(i,j))
          hooverNodes_.insert(pair<int,int>(i,j));

    // (b) AtC coupling nodes
    if (kinetostat_->get_coupling_mode()==Kinetostat::FIXED) {
      Array<int> & nodeType(atcTransfer_->get_node_type());
      if (atcTransfer_->use_localized_lambda()) {
        for (int i = 0; i < nNodes_; ++i)
          if (nodeType(i)==ATC_Transfer::BOUNDARY)
            for (int j = 0; j < nsd_; ++j)
              hooverNodes_.insert(pair<int,int>(i,j));
      }
      else {
        for (int i = 0; i < nNodes_; ++i)
          if (nodeType(i)==ATC_Transfer::BOUNDARY || nodeType(i)==ATC_Transfer::MD_ONLY)
            for (int j = 0; j < nsd_; ++j)
              hooverNodes_.insert(pair<int,int>(i,j));
      }
    }
  }

  //--------------------------------------------------------
  //  apply_lambda_to_atoms
  //            uses existing lambda to modify given
  //            atomic quantity
  //--------------------------------------------------------
  void GlcKinetostat::apply_to_atoms(double ** atomicQuantity,
                                     const DENS_MAT & lambdaAtom,
                                     double dt)
  {
    if (nLocal_>0) {
      for (int i = 0; i < nLocal_; i++) {
        for (int j = 0; j < nsd_; j++) {
          atomicQuantity[internalToAtom_(i)][j] -= lambdaAtom(i,j);
        }
      }
    }
  }

  //--------------------------------------------------------
  //  set_weights
  //            sets diagonal weighting matrix used in
  //            solve_for_lambda
  //--------------------------------------------------------
  void GlcKinetostat::set_weights(DIAG_MAT & weights)
  {
    if (nLocalLambda_>0) {
      DENS_VEC weightVector(nLocal_);
      weightVector = 1.;
      DENS_VEC maskedWeightVector(nLocalLambda_);
      maskedWeightVector = internalToOverlapMap_*weightVector;
      weights.reset(maskedWeightVector);
    }
  }

  //--------------------------------------------------------
  //  reset_nlocal:
  //    resets data dependent on local atom count
  //--------------------------------------------------------
  void GlcKinetostat::reset_nlocal()
  {
    RegulatorShapeFunction::reset_nlocal();
    if (nLocal_ > 0) {
      atcTransfer_->compute_atomic_mass(atomicMass_);
    }
  }
 
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class DisplacementGlc
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and kinetostat data
  //--------------------------------------------------------
  DisplacementGlc::DisplacementGlc(Kinetostat * kinetostat) :
    GlcKinetostat(kinetostat),
    nodalDisplacements_(atcTransfer_->get_field(DISPLACEMENT)),
    x_(atcTransfer_->get_x())           
  {
    // sets up time filter for cases where variables temporally filtered
    TimeFilterManager * timeFilterManager = atcTransfer_->get_time_filter_manager();
    if (!timeFilterManager->end_equilibrate()) {
      nodalAtomicLambdaForce_ = 0.;
      lambdaForceFiltered_ = 0.;
      timeFilter_->initialize(lambdaForceFiltered_);
    }
  }
 
  //--------------------------------------------------------
  //  apply:
  //    apply the kinetostat to the atoms
  //--------------------------------------------------------
  void DisplacementGlc::apply_post_predictor(double dt)
  {
    compute_kinetostat(dt);
  }

  //--------------------------------------------------------
  //  compute_kinetostat
  //            manages the solution and application of the
  //            kinetostat equations and variables
  //--------------------------------------------------------
  void DisplacementGlc::compute_kinetostat(double dt)
  {
    // initial filtering update
    apply_pre_filtering(dt);

    // set up rhs
    DENS_MAT rhs(nNodes_,nsd_);
    set_kinetostat_rhs(rhs,dt);

    // solve linear system for lambda
    solve_for_lambda(rhs);

    // compute force applied by lambda
    DENS_MAT lambdaAtom;
    compute_lambda_force(lambdaForce_,lambdaAtom,dt);

    // compute nodal atomic power
    compute_nodal_lambda_force(dt);

    // apply kinetostat to atoms
    apply_to_atoms(x_,lambdaAtom);
  }

  //--------------------------------------------------------
  //  set_kinetostat_rhs
  //            sets up the right-hand side of the
  //            kinetostat equations
  //--------------------------------------------------------
  void DisplacementGlc::set_kinetostat_rhs(DENS_MAT & rhs, double dt)
  {
    // form rhs : sum_a (hatN_Ia * x_ai) - (Upsilon)_Ii
    DENS_MAT atomicDisplacement;
    atcTransfer_->compute_atomic_centerOfMass_displacement(atomicDisplacement,x_);
    rhs.reset(nNodes_,nsd_);
    atcTransfer_->restrict_volumetric_quantity(atomicDisplacement,rhs);
    rhs -= (atcTransfer_->get_mass_mat_md(VELOCITY))*nodalDisplacements_;
  }

  //--------------------------------------------------------
  //  compute_lambda_force
  //            compute the equivalent force on the atoms
  //            induced by lambda
  //--------------------------------------------------------
  void DisplacementGlc::compute_lambda_force(DENS_MAT & lambdaForce,
                                             DENS_MAT & lambdaAtom,
                                             double dt)
  {
    if (nLocal_>0) {
      // prolongation to (unique) nodes
      lambdaAtom.reset(nLocal_,nsd_);
      atcTransfer_->prolong(lambda_,lambdaAtom);

      for (int i = 0; i < nLocal_; i++)
        for (int j = 0; j < nsd_; j++)
          lambdaForce(i,j) =  (-1./(dt*dt))*lambdaAtom(i,j)*atomicMass_(i);
    }
  }

  //--------------------------------------------------------
  //  compute_nodal_lambda_force
  //            compute the effective FE force applied
  //            by the kinetostat
  //--------------------------------------------------------
  void DisplacementGlc::compute_nodal_lambda_force(double dt)
  {
    atcTransfer_->restrict_volumetric_quantity(lambdaForce_,nodalAtomicLambdaForce_);
    timeFilter_->apply_post_step1(lambdaForceFiltered_,nodalAtomicLambdaForce_,dt);

    // update FE displacements for localized thermostats
    apply_localization_correction(nodalAtomicLambdaForce_,
                                  nodalDisplacements_,
                                  dt*dt);
  }

  //--------------------------------------------------------
  //  apply_pre_filtering
  //            applies first step of filtering to
  //            relevant variables
  //--------------------------------------------------------
  void DisplacementGlc::apply_pre_filtering(double dt)
  {
    // apply time filtered lambda force
    DENS_MAT lambdaZero(nNodes_,nsd_);
    timeFilter_->apply_pre_step1(lambdaForceFiltered_,(-1./dt/dt)*lambdaZero,dt);
  }

  //--------------------------------------------------------
  //  set_weights
  //            sets diagonal weighting matrix used in
  //            solve_for_lambda
  //--------------------------------------------------------
  void DisplacementGlc::set_weights(DIAG_MAT & weights)
  {
    if (nLocalLambda_>0) {
      DENS_VEC maskedWeightVector(nLocalLambda_);
      maskedWeightVector = internalToOverlapMap_*atomicMass_;
      weights.reset(maskedWeightVector);
    }
  }

  //--------------------------------------------------------
  //  apply_localization_correction
  //            corrects for localized kinetostats only
  //            solving kinetostat equations on a subset
  //            of the MD region
  //--------------------------------------------------------
  void DisplacementGlc::apply_localization_correction(const DENS_MAT & source,
                                                      DENS_MAT & nodalField,
                                                      double weight)
  {
    // NOTE can add check to see if set is empty for faster performance
    DENS_MAT nodalLambdaRoc(nNodes_,nsd_);
    atcTransfer_->apply_inverse_mass_matrix(source,
                                            nodalLambdaRoc,
                                            VELOCITY);
    set<pair<int,int> >::const_iterator iter;
    for (iter = hooverNodes_.begin(); iter != hooverNodes_.end(); ++iter) {
      nodalLambdaRoc(iter->first,iter->second) = 0.;
    }

    nodalField += weight*nodalLambdaRoc;
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void DisplacementGlc::output(double dt, OUTPUT_LIST & outputData)
  {
    outputData["lambda"] = &(lambda_);
    outputData["nodalLambdaForce"] = &(nodalAtomicLambdaForce_);
  }
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class DisplacementGlcFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and kinetostat data
  //--------------------------------------------------------
  DisplacementGlcFiltered::DisplacementGlcFiltered(Kinetostat * kinetostat) :
    DisplacementGlc(kinetostat),
    nodalAtomicDisplacements_(atcTransfer_->get_atomic_field(DISPLACEMENT))
  {
    // do nothing
  }
 
  //--------------------------------------------------------
  //  apply_pre_filtering
  //            applies first step of filtering to
  //            relevant variables
  //--------------------------------------------------------
  void DisplacementGlcFiltered::apply_pre_filtering(double dt)
  {
    // apply time filtered lambda to atomic fields
    DisplacementGlc::apply_pre_filtering(dt);
    DENS_MAT nodalAcceleration(nNodes_,nsd_);
    atcTransfer_->apply_inverse_md_mass_matrix(lambdaForceFiltered_,
                                               nodalAcceleration,
                                               VELOCITY);
    nodalAtomicDisplacements_ += dt*dt*nodalAcceleration;
  }

  //--------------------------------------------------------
  //  set_kinetostat_rhs
  //            sets up the right-hand side of the
  //            kinetostat equations
  //--------------------------------------------------------
  void DisplacementGlcFiltered::set_kinetostat_rhs(DENS_MAT & rhs, double dt)
  {
    // form rhs : sum_a (hatN_Ia * x_ai) - (Upsilon)_Ii
    double coef = 1./(timeFilter_->get_unfiltered_coefficient_pre_s1(dt));
    rhs = coef*(atcTransfer_->get_mass_mat_md(VELOCITY))*(nodalAtomicDisplacements_ - nodalDisplacements_);
  }

  //--------------------------------------------------------
  //  compute_nodal_lambda_force
  //            compute the effective FE force applied
  //            by the kinetostat
  //--------------------------------------------------------
  void DisplacementGlcFiltered::compute_nodal_lambda_force(double dt)
  {
    atcTransfer_->restrict_volumetric_quantity(lambdaForce_,nodalAtomicLambdaForce_);
    timeFilter_->apply_post_step1(lambdaForceFiltered_,nodalAtomicLambdaForce_,dt);

    // update filtered atomic displacements
    DENS_MAT nodalLambdaRoc(nodalAtomicLambdaForce_.nRows(),nodalAtomicLambdaForce_.nCols());
    atcTransfer_->apply_inverse_md_mass_matrix(nodalAtomicLambdaForce_,
                                               nodalLambdaRoc,
                                               VELOCITY);
    timeFilter_->apply_post_step1(nodalAtomicDisplacements_,dt*dt*nodalLambdaRoc,dt);

    // update FE displacements for localized thermostats
    apply_localization_correction(lambdaForceFiltered_,
                                  nodalDisplacements_,
                                  dt*dt);
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void DisplacementGlcFiltered::output(double dt, OUTPUT_LIST & outputData)
  {
    outputData["lambda"] = &(lambda_);
    outputData["nodalLambdaForce"] = &(lambdaForceFiltered_);
  }
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class VelocityGlc
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and kinetostat data
  //--------------------------------------------------------
  VelocityGlc::VelocityGlc(Kinetostat * kinetostat) :
    GlcKinetostat(kinetostat),
    nodalVelocities_(atcTransfer_->get_field(VELOCITY)),
    v_(atcTransfer_->get_v())
  {
    // sets up time filter for cases where variables temporally filtered
    TimeFilterManager * timeFilterManager = atcTransfer_->get_time_filter_manager();
    if (!timeFilterManager->end_equilibrate()) {
      nodalAtomicLambdaForce_ = 0.;
      lambdaForceFiltered_ = 0.;
      timeFilter_->initialize(lambdaForceFiltered_);
    }
  }
 
  //--------------------------------------------------------
  //  apply_mid_corrector:
  //    apply the kinetostat during the middle of the
  //    predictor phase
  //--------------------------------------------------------
  void VelocityGlc::apply_mid_predictor(double dt)
  {
    double dtLambda = 0.5*dt;
    compute_kinetostat(dtLambda);
  }

  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the kinetostat after the corrector phase
  //--------------------------------------------------------
  void VelocityGlc::apply_post_corrector(double dt)
  {
    double dtLambda = 0.5*dt;
    compute_kinetostat(dt);
  }

  //--------------------------------------------------------
  //  apply_pre_filtering
  //            applies first step of filtering to
  //            relevant variables
  //--------------------------------------------------------
  void VelocityGlc::apply_pre_filtering(double dt)
  {
    // apply time filtered lambda to atomic fields
    DENS_MAT lambdaZero(nNodes_,nsd_);
    timeFilter_->apply_pre_step1(lambdaForceFiltered_,(-1./dt)*lambdaZero,dt);
  }

  //--------------------------------------------------------
  //  compute_kinetostat
  //            manages the solution and application of the
  //            kinetostat equations and variables
  //--------------------------------------------------------
  void VelocityGlc::compute_kinetostat(double dt)
  {
    // initial filtering update
    apply_pre_filtering(dt);

    // set up rhs
    DENS_MAT rhs(nNodes_,nsd_);
    set_kinetostat_rhs(rhs,dt);
    
    // solve linear system for lambda
    solve_for_lambda(rhs);
    
    // compute force applied by lambda
    DENS_MAT lambdaAtom;
    compute_lambda_force(lambdaForce_,lambdaAtom,dt);

    // compute nodal atomic power
    compute_nodal_lambda_force(dt);

    // apply kinetostat to atoms
    apply_to_atoms(v_,lambdaAtom);
  }

  //--------------------------------------------------------
  //  set_kinetostat_rhs
  //            sets up the right-hand side of the
  //            kinetostat equations
  //--------------------------------------------------------
  void VelocityGlc::set_kinetostat_rhs(DENS_MAT & rhs, double dt)
  {
    // form rhs : sum_a (hatN_Ia * x_ai) - (\dot{Upsilon})_Ii
    DENS_MAT atomicMomentum;
    atcTransfer_->compute_atomic_momentum(atomicMomentum,v_);
    rhs.reset(nNodes_,nsd_);
    atcTransfer_->restrict_volumetric_quantity(atomicMomentum,rhs);
    rhs -= (atcTransfer_->get_mass_mat_md(VELOCITY))*nodalVelocities_;
  }

  //--------------------------------------------------------
  //  compute_lambda_force
  //            compute the equivalent force on the atoms
  //            induced by lambda
  //--------------------------------------------------------
  void VelocityGlc::compute_lambda_force(DENS_MAT & lambdaForce,
                                          DENS_MAT & lambdaAtom,
                                          double dt)
  {
    if (nLocal_>0) {
      // prolongation to (unique) nodes
      lambdaAtom.reset(nLocal_,nsd_);
      atcTransfer_->prolong(lambda_,lambdaAtom);

      for (int i = 0; i < nLocal_; i++)
        for (int j = 0; j < nsd_; j++)
          lambdaForce(i,j) =  (-1./(dt))*lambdaAtom(i,j)*atomicMass_(i);
    }
  }

  //--------------------------------------------------------
  //  compute_nodal_lambda_force
  //            compute the effective FE force applied
  //            by the kinetostat
  //--------------------------------------------------------
  void VelocityGlc::compute_nodal_lambda_force(double dt)
  {
    atcTransfer_->restrict_volumetric_quantity(lambdaForce_,nodalAtomicLambdaForce_);
    timeFilter_->apply_pre_step1(lambdaForceFiltered_,nodalAtomicLambdaForce_,dt);

    // update FE displacements for localized thermostats
    apply_localization_correction(nodalAtomicLambdaForce_,
                                  nodalVelocities_,
                                  dt);
  }

  //--------------------------------------------------------
  //  set_weights
  //            sets diagonal weighting matrix used in
  //            solve_for_lambda
  //--------------------------------------------------------
  void VelocityGlc::set_weights(DIAG_MAT & weights)
  {
    if (nLocalLambda_>0) {
      DENS_VEC maskedWeightVector = internalToOverlapMap_*atomicMass_;
      weights.reset(maskedWeightVector);
    }
  }

  //--------------------------------------------------------
  //  apply_localization_correction
  //            corrects for localized kinetostats only
  //            solving kinetostat equations on a subset
  //            of the MD region
  //--------------------------------------------------------
  void VelocityGlc::apply_localization_correction(const DENS_MAT & source,
                                                  DENS_MAT & nodalField,
                                                  double weight)
  {
    // NOTE can add check to see if set is empty for faster performance
    DENS_MAT nodalLambdaRoc(nNodes_,nsd_);
    atcTransfer_->apply_inverse_mass_matrix(source,
                                            nodalLambdaRoc,
                                            VELOCITY);
    set<pair<int,int> >::const_iterator iter;
    for (iter = hooverNodes_.begin(); iter != hooverNodes_.end(); ++iter) {
      nodalLambdaRoc(iter->first,iter->second) = 0.;
    }

    nodalField += weight*nodalLambdaRoc;
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void VelocityGlc::output(double dt, OUTPUT_LIST & outputData)
  {
    outputData["lambda"] = &(lambda_);
    outputData["nodalLambdaForce"] = &(nodalAtomicLambdaForce_);
  }
 
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class VelocityGlcFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and kinetostat data
  //--------------------------------------------------------
  VelocityGlcFiltered::VelocityGlcFiltered(Kinetostat *kinetostat) 
    : VelocityGlc(kinetostat),
      nodalAtomicVelocities_(atcTransfer_->get_atomic_field(VELOCITY))
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  apply_pre_filtering
  //            applies first step of filtering to
  //            relevant variables
  //--------------------------------------------------------
  void VelocityGlcFiltered::apply_pre_filtering(double dt)
  {
    // apply time filtered lambda to atomic fields
    VelocityGlc::apply_pre_filtering(dt);
    DENS_MAT nodalAcceleration(nNodes_,nsd_);
    atcTransfer_->apply_inverse_md_mass_matrix(lambdaForceFiltered_,
                                               nodalAcceleration,
                                               VELOCITY);
    nodalAtomicVelocities_ += dt*nodalAcceleration;
  }

  //--------------------------------------------------------
  //  set_kinetostat_rhs
  //            sets up the right-hand side of the
  //            kinetostat equations
  //--------------------------------------------------------
  void VelocityGlcFiltered::set_kinetostat_rhs(DENS_MAT & rhs, double dt)
  {
    // form rhs : sum_a (hatN_Ia * x_ai) - (Upsilon)_Ii
    double coef = 1./(timeFilter_->get_unfiltered_coefficient_pre_s1(dt));
    rhs = coef*(atcTransfer_->get_mass_mat_md(VELOCITY))*(nodalAtomicVelocities_ - nodalVelocities_);
  }

  //--------------------------------------------------------
  //  compute_nodal_lambda_force
  //            compute the effective FE force applied
  //            by the kinetostat
  //--------------------------------------------------------
  void VelocityGlcFiltered::compute_nodal_lambda_force(double dt)
  {
    atcTransfer_->restrict_volumetric_quantity(lambdaForce_,nodalAtomicLambdaForce_);
    timeFilter_->apply_post_step1(lambdaForceFiltered_,nodalAtomicLambdaForce_,dt);

    // update filtered atomic displacements
    DENS_MAT nodalLambdaRoc(nodalAtomicLambdaForce_.nRows(),nodalAtomicLambdaForce_.nCols());
    atcTransfer_->apply_inverse_md_mass_matrix(nodalAtomicLambdaForce_,
                                               nodalLambdaRoc,
                                               VELOCITY);
    timeFilter_->apply_post_step1(nodalAtomicVelocities_,dt*nodalLambdaRoc,dt);

    // update FE displacements for localized thermostats
    apply_localization_correction(lambdaForceFiltered_,
                                 nodalVelocities_,
                                 dt);
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void VelocityGlcFiltered::output(double dt, OUTPUT_LIST & outputData)
  {
    outputData["lambda"] = &(lambda_);
    outputData["nodalLambdaForce"] = &(lambdaForceFiltered_);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class StressFlux
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and kinetostat data
  //--------------------------------------------------------
  StressFlux::StressFlux(Kinetostat * kinetostat) :
    GlcKinetostat(kinetostat),
    nodalForce_(atcTransfer_->get_field_rhs(VELOCITY)),
    nodalAtomicForce_(atcTransfer_->get_fe_atomic_field_roc(VELOCITY)),
    momentumSource_(atcTransfer_->get_atomic_source(VELOCITY)),
    v_(atcTransfer_->get_v()),
    f_(atcTransfer_->get_f())
  {
    // flag for performing boundary flux calculation
    if (kinetostat_->get_coupling_mode()==Kinetostat::FLUX)
      fieldMask_(VELOCITY,FLUX) = true;

    // sets up space for ghost force related variables
    nodalGhostForce_.reset(nNodes_,nsd_);
    nodalGhostForceFiltered_.reset(nNodes_,nsd_);

    // NOTE ifdefs are for using reference lammps force as base pressure
#if false
    int nLocalTotal = atcTransfer->get_nlocal_total();
    int nsd = atcTransfer->get_nsd();
    f0_ = new double*[nLocalTotal];
    for (int i = 0; i < nLocalTotal; ++i) {
      f0_[i] = new double[nsd];
      for (int j = 0; j < nsd; ++j)
        f0_[i][j] = f_[i][j];
    }
#endif
  }

  StressFlux::~StressFlux()
  {
#if false
    if (f0_) {
      ATC_Transfer * atcTransfer = kinetostat_->get_atc_transfer();
      int nLocalTotal = atcTransfer->get_nlocal_total();
      for (int i = 0; i < nLocalTotal; ++i)
        delete [] f0_[i];
      delete [] f0_;
    }
#endif
  }

  //--------------------------------------------------------
  //  apply_pre_predictor:
  //    apply the kinetostat to the atoms in the
  //    mid-predictor integration phase
  //--------------------------------------------------------
  void StressFlux::apply_mid_predictor(double dt)
  {
    double dtLambda = 0.5*dt;
    // apply lambda force to atoms
    apply_to_atoms(v_,lambdaForce_,dtLambda);
  }

  //--------------------------------------------------------
  //  apply_pre_corrector:
  //    apply the kinetostat to the atoms in the
  //    pre-corrector integration phase
  //--------------------------------------------------------
  void StressFlux::apply_pre_corrector(double dt)
  {
    // compute the kinetostat force
    compute_kinetostat(dt);
  }
 
  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the kinetostat to the atoms in the
  //    post-corrector integration phase
  //--------------------------------------------------------
  void StressFlux::apply_post_corrector(double dt)
  {
    double dtLambda = 0.5*dt;
    // apply lambda force to atoms
    apply_to_atoms(v_,lambdaForce_,dtLambda);
  }

  //--------------------------------------------------------
  //  compute_kinetostat
  //            manages the solution and application of the
  //            kinetostat equations and variables
  //--------------------------------------------------------
  void StressFlux::compute_kinetostat(double dt)
  {
    // set up ghost force
    compute_ghost_force(nodalGhostForce_);

    // initial filtering update
    apply_pre_filtering(dt);

    // set up rhs
    DENS_MAT rhs(nNodes_,nsd_);
    set_kinetostat_rhs(rhs,dt);
    
    // solve linear system for lambda
    solve_for_lambda(rhs);
    
    // compute force applied by lambda
    compute_lambda_force(lambdaForce_);

    // compute nodal atomic power
    compute_nodal_lambda_force(dt);
  }

  //--------------------------------------------------------
  //  apply_pre_filtering
  //            applies first step of filtering to
  //            relevant variables
  //--------------------------------------------------------
  void StressFlux::apply_pre_filtering(double dt)
  {
    // apply time filtered lambda force
    DENS_MAT lambdaZero(nNodes_,nsd_);
    timeFilter_->apply_pre_step1(lambdaForceFiltered_,lambdaZero,dt);
    timeFilter_->apply_pre_step1(nodalGhostForceFiltered_,nodalGhostForce_,dt);
  }

  //--------------------------------------------------------
  //  set_kinetostat_rhs
  //            sets up the RHS of the kinetostat equations
  //            for the coupling parameter lambda
  //--------------------------------------------------------
  void StressFlux::set_kinetostat_rhs(DENS_MAT & rhs, double dt)
  {
    // (a) for flux based : 
    // form rhs :  \int N_I r dV - \sum_g N_Ig^* f_g
    // sources are set in ATC transfer
    rhs.reset(nNodes_,nsd_);
    rhs = momentumSource_ - nodalGhostForce_;

    // (b) for ess. bcs
    // NOTE can add check to see if set is empty for faster performance
    // form rhs : {sum_a (N_Ia * f_ia) - M_md * (ddupsilon/dt)_I}
    DENS_MAT rhsPrescribed = -1.*nodalForce_;
    atcTransfer_->apply_inverse_mass_matrix(rhsPrescribed,VELOCITY);
    rhsPrescribed = mdMassMatrix_*rhsPrescribed;
    rhsPrescribed += nodalAtomicForce_;

    set<pair<int,int> >::const_iterator iter;
    for (iter = hooverNodes_.begin(); iter != hooverNodes_.end(); ++iter) {
      rhs(iter->first,iter->second) = rhsPrescribed(iter->first,iter->second);
    }
  }

  //--------------------------------------------------------
  //  compute_lambda_force
  //            computes the force induced by lambda
  //            on the atoms
  //--------------------------------------------------------
  void StressFlux::compute_lambda_force(DENS_MAT & lambdaForce)
  {
    if (nLocal_>0) {
      // prolongation to (unique) nodes
      lambdaForce.reset(nLocal_,nsd_);
      atcTransfer_->prolong(lambda_,lambdaForce);
      lambdaForce *= -1.;
    }
  }

  //--------------------------------------------------------
  //  compute_nodal_lambda_force
  //            computes the force induced on the FE
  //            by applying lambdaForce on the atoms
  //--------------------------------------------------------
  void StressFlux::compute_nodal_lambda_force(double dt)
  {
    atcTransfer_->restrict_volumetric_quantity(lambdaForce_,nodalAtomicLambdaForce_);
    set<pair<int,int> >::const_iterator iter;
    for (iter = hooverNodes_.begin(); iter != hooverNodes_.end(); ++iter) {
      nodalAtomicLambdaForce_(iter->first,iter->second) = 0.;
    }

    timeFilter_->apply_post_step1(lambdaForceFiltered_,nodalAtomicLambdaForce_,dt);
  }

  //--------------------------------------------------------
  //  add_to_rhs:
  //            determines what if any contributions to the
  //            finite element equations are needed for
  //            consistency with the kinetostat
  //--------------------------------------------------------
  void StressFlux::add_to_rhs(FIELDS & rhs)
  {
    rhs[VELOCITY] += nodalAtomicLambdaForce_ + boundaryFlux_[VELOCITY];
  }

  //--------------------------------------------------------
  //  compute_ghost_force
  //            computes computes the restricted force on
  //            the ghost atoms
  //--------------------------------------------------------
  void StressFlux::compute_ghost_force(DENS_MAT & nodalGhostForce)
  { 
    DENS_MAT nodalGhostForceLocal(nNodes_,nsd_);
    nodalGhostForce.reset(nNodes_,nsd_);
    
    if (nLocalGhost_ > 0) {
      DENS_MAT ghostForce(nLocalGhost_,nsd_);
      for (int i = 0; i < nLocalGhost_; i++) {
        int atomIndex = ghostToAtom_(i);
        for (int j = 0; j < nsd_; j++) {
          ghostForce(i,j) = f_[atomIndex][j];
#if false
          ghostForce(i,j) -= f0_[atomIndex][j];
#endif
        }
      }
      nodalGhostForceLocal = shapeFunctionGhost_.transMat(ghostForce);
    }
    
    LammpsInterface::instance()->allsum(nodalGhostForceLocal.get_ptr(),
                                        nodalGhostForce.get_ptr(), nodalGhostForce.size());
    // convert from Lammps force units to ATC force units
    nodalGhostForce *= LammpsInterface::instance()->ftm2v();
  }

  //--------------------------------------------------------
  //  apply_lambda_to_atoms
  //            uses existing lambda to modify given
  //            atomic quantity
  //--------------------------------------------------------
  void StressFlux::apply_to_atoms(double ** atomicVelocity,
                                  const DENS_MAT & lambdaForce,
                                  double dt)
  {
    if (nLocal_>0) {
      // explicit update
      for (int i = 0; i < nLocal_; i++) {
        for (int j = 0; j < nsd_; j++) {
          atomicVelocity[internalToAtom_(i)][j] += dt*lambdaForce(i,j)/atomicMass_(i);
        }
      }
    }
  }

  //--------------------------------------------------------
  //  reset_filtered_ghost_force:
  //    resets the kinetostat generated ghost force to a
  //    prescribed value
  //--------------------------------------------------------
  void StressFlux::reset_filtered_ghost_force(DENS_MAT & target)
  {
    for (int i = 0; i < nNodes_; ++i)
      for (int j = 0; j < nsd_; ++j)
        nodalGhostForceFiltered_(i,j) = target(i,j);
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void StressFlux::output(double dt, OUTPUT_LIST & outputData)
  {
    outputData["lambda"] = &(lambda_);
    outputData["nodalLambdaForce"] = &(nodalAtomicLambdaForce_);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class StressFluxFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  StressFluxFiltered::StressFluxFiltered(Kinetostat * kinetostat) :
    StressFlux(kinetostat),
    nodalAtomicVelocity_(atcTransfer_->get_atomic_field(VELOCITY))
  {
    // do nothing
  }
 
#if false
  //--------------------------------------------------------
  //  apply:
  //    apply the kinetostat to the atoms
  //--------------------------------------------------------
  void StressFluxFiltered::apply_pre_corrector(double dt)
  {
    // NOTE currently not implemented, below based on GlcVelocityFiltered

    // get neccesary data
    double coef = 1./(timeFilter_->get_unfiltered_coefficient_pre_s1(dt));
   
    // form rhs : sum_a (hatN_Ia * x_ai) - (\dot{Upsilon})_Ii
    DENS_MAT rhs(nNodes_,nsd_);
    for (int i = 0; i < nNodes_; i++)
      for (int j = 0; j < nsd_; j++)
        rhs(i,j) += coef*(nodalAtomicVelocities_(i,j) - nodalVelocities_(i,j));
    DENS_MAT rhsOverlap(nNodeOverlap_,nsd_);
    atcTransfer_->map_unique_to_overlap(rhs, rhsOverlap);
    
    // solve matrix equation and map back to all nodes
    DENS_MAT lambdaOverlap(nNodeOverlap_,nsd_);
    DIAG_MAT weights;
    set_weights(weights);
    solve_for_lambda(rhsOverlap, weights, lambdaOverlap);
    atcTransfer_->map_overlap_to_unique(lambdaOverlap,lambda_);
        
    // apply lambda to the atoms and nodal atomic fields
    DENS_MAT lambdaRestricted(nNodes_,nsd_);
    apply_lambda_to_atoms(v_,lambdaRestricted);
    timeFilter_->apply_post_step1(nodalAtomicVelocities_,-1.*lambdaRestricted,dt);
    timeFilter_->apply_post_step1(nodalLambdaForce_,(-1./dt)*lambdaRestricted,dt);
    DENS_MAT nodalLambdaRoc(nodalLambdaForce_.nRows,nodalLambdaForce_.nCols());
    atcTransfer_->apply_inverse_mass_matrix(nodalLambdaForce_,
                                            nodalLambdaRoc,
                                            VELOCITY);
    atcTransfer_->apply_internal_atomic_source(nodalVelocities_,
                                               nodalLambdaRoc,
                                               dt);

  }
#endif
//   //--------------------------------------------------------
//   //  apply_pre_filtering
//   //            applies first step of filtering to
//   //            relevant variables
//   //--------------------------------------------------------
//   void StressFluxFiltered::apply_pre_filtering(double dt)
//   {
//     // apply time filtered lambda to atomic fields
//     StressFlux::apply_pre_filtering(dt);
//     nodalAtomicForce_ += lambdaForceFiltered;
//   }

  //--------------------------------------------------------
  //  set_kinetostat_rhs
  //            sets up the RHS of the kinetostat equations
  //            for the coupling parameter lambda
  //--------------------------------------------------------
  void StressFluxFiltered::set_kinetostat_rhs(DENS_MAT & rhs, double dt)
  {
    // set basic terms
    // (a) for flux based : 
    // form rhs :  \int N_I r dV - \sum_g N_Ig^* f_g
    // sources are set in ATC transfer
    rhs.reset(nNodes_,nsd_);
    rhs = momentumSource_ - nodalGhostForceFiltered_;
    
    // (b) for ess. bcs
    // NOTE can add check to see if set is empty for faster performance
    // form rhs : {sum_a (N_Ia * f_ia) - M_md * (ddupsilon/dt)_I}
    DENS_MAT rhsPrescribed = -1.*nodalForce_;
    atcTransfer_->apply_inverse_mass_matrix(rhsPrescribed,VELOCITY);
    rhsPrescribed = mdMassMatrix_*rhsPrescribed;
    rhsPrescribed += nodalAtomicForce_;

    set<pair<int,int> >::const_iterator iter;
    for (iter = hooverNodes_.begin(); iter != hooverNodes_.end(); ++iter) {
      rhs(iter->first,iter->second) = rhsPrescribed(iter->first,iter->second);
    }

    // adjust for application of current lambda force
    rhs += lambdaForceFiltered_;

    // correct for time filtering
    rhs *= 1./(timeFilter_->get_unfiltered_coefficient_pre_s1(dt));
  }

  //--------------------------------------------------------
  //  apply_lambda_to_atoms
  //            uses existing lambda to modify given
  //            atomic quantity
  //--------------------------------------------------------
  void StressFluxFiltered::apply_to_atoms(double ** atomicVelocity,
                                          const DENS_MAT & lambdaForce,
                                          double dt)
  {
    StressFlux::apply_to_atoms(atomicVelocity,lambdaForce,dt);

    // add in corrections to filtered nodal atomice velocity
    DENS_MAT velocityRoc(nNodes_,nsd_);
    atcTransfer_->apply_inverse_md_mass_matrix(lambdaForceFiltered_,
                                               velocityRoc,
                                               VELOCITY);
    nodalAtomicVelocity_ += dt*velocityRoc;
  }

  //--------------------------------------------------------
  //  add_to_rhs:
  //            determines what if any contributions to the
  //            finite element equations are needed for
  //            consistency with the kinetostat
  //--------------------------------------------------------
  void StressFluxFiltered::add_to_rhs(FIELDS & rhs)
  {
    rhs[VELOCITY] += lambdaForceFiltered_ + boundaryFlux_[VELOCITY];
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void StressFluxFiltered::output(double dt, OUTPUT_LIST & outputData)
  {
    outputData["lambda"] = &(lambda_);
    outputData["nodalLambdaForce"] = &(lambdaForceFiltered_);
  }
  
};
