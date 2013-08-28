#include "Kinetostat.h"
#include "ATC_Error.h"
#include "ATC_Coupling.h"
#include "LammpsInterface.h"
#include "PerAtomQuantityLibrary.h"
#include "PrescribedDataManager.h"
#include "ElasticTimeIntegrator.h"
#include "TransferOperator.h"

using std::set;
using std::pair;
using std::string;

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class Kinetostat
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  Kinetostat::Kinetostat(ATC_Coupling * atc,
                         const string & regulatorPrefix) :
    AtomicRegulator(atc,regulatorPrefix)
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

    int argIndex = 0;
    if (strcmp(arg[argIndex],"momentum")==0) {
      argIndex++;

      // fluxstat type
      /*! \page man_control_momentum fix_modify AtC control momentum
        \section syntax
        fix_modify AtC control momentum none \n
      
        fix_modify AtC control momentum rescale <frequency>\n
        - frequency (int) = time step frequency for applying displacement and velocity rescaling \n
      
        fix_modify AtC control momentum glc_displacement \n
  
        fix_modify AtC control momentum glc_velocity \n

        fix_modify AtC control momentum hoover \n

        fix_modify AtC control momentum flux [faceset face_set_id, interpolate] 
        - face_set_id (string) = id of boundary face set, if not specified
        (or not possible when the atomic domain does not line up with 
        mesh boundaries) defaults to an atomic-quadrature approximate 
        evaulation\n
        \section examples
        fix_modify AtC control momentum glc_velocity \n
        fix_modify AtC control momentum flux faceset bndy_faces \n
        \section description
        \section restrictions
        only to be used with specific transfers :
        elastic \n
        rescale not valid with time filtering activated
        \section related
        \section default
        none
      */
      boundaryIntegrationType_ = NO_QUADRATURE;
      howOften_ = 1;
      if (strcmp(arg[argIndex],"none")==0) { // restore defaults
        regulatorTarget_ = NONE;
        couplingMode_ = UNCOUPLED;
        foundMatch = true;
      }
      else if (strcmp(arg[argIndex],"glc_displacement")==0) {
        regulatorTarget_ = FIELD;
        couplingMode_ = FIXED;
        foundMatch = true;
      }
      else if (strcmp(arg[argIndex],"glc_velocity")==0) {
        regulatorTarget_ = DERIVATIVE;
        couplingMode_ = FIXED;
        foundMatch = true;
      }
      else if (strcmp(arg[argIndex],"hoover")==0) {
        regulatorTarget_ = DYNAMICS;
        couplingMode_ = FIXED;
        foundMatch = true;
      }
      else if (strcmp(arg[argIndex],"flux")==0) {
        regulatorTarget_ = DYNAMICS;
        couplingMode_ = FLUX;
        argIndex++;
        boundaryIntegrationType_ = atc_->parse_boundary_integration(narg-argIndex,&arg[argIndex],boundaryFaceSet_);
        foundMatch = true;
      }
      else if (strcmp(arg[argIndex],"ghost_flux")==0) {
        regulatorTarget_ = DYNAMICS;
        couplingMode_ = GHOST_FLUX;
        foundMatch = true;
      }
    }
  
    if (!foundMatch)
      foundMatch = AtomicRegulator::modify(narg,arg);
    if (foundMatch)
      needReset_ = true;
    return foundMatch;
  }

  //--------------------------------------------------------
  //  reset_lambda_contribution
  //    resets the kinetostat generated force to a
  //    prescribed value
  //--------------------------------------------------------
  void Kinetostat::reset_lambda_contribution(const DENS_MAT & target)
  {
    DENS_MAN * lambdaForceFiltered = regulator_data("LambdaForceFiltered",nsd_);
    lambdaForceFiltered->set_quantity() = target;
  }
 
  //--------------------------------------------------------
  //  initialize:
  //    sets up methods before a run
  
  //    dependence, but in general there is also a
  //    time integrator dependence.  In general the 
  //    precedence order is:
  //    time filter -> time integrator -> kinetostat
  //    In the future this may need to be added if
  //    different types of time integrators can be
  //    specified.
  //--------------------------------------------------------
  void Kinetostat::construct_methods()
  {
    // get data associated with stages 1 & 2 of ATC_Method::initialize
    AtomicRegulator::construct_methods();

    if (atc_->reset_methods()) {
      // eliminate existing methods
      delete_method();
      
      DENS_MAT nodalGhostForceFiltered;
      TimeIntegrator::TimeIntegrationType myIntegrationType = (atc_->time_integrator(VELOCITY))->time_integration_type();
      TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
      if (timeFilterManager->end_equilibrate() && regulatorTarget_==AtomicRegulator::DYNAMICS) {
        StressFlux * myMethod;
        myMethod = dynamic_cast<StressFlux *>(regulatorMethod_);
        nodalGhostForceFiltered = (myMethod->filtered_ghost_force()).quantity();
      }
      
      // update time filter
      
      if (timeFilterManager->need_reset()) {
        
        timeFilter_ = timeFilterManager->construct(TimeFilterManager::IMPLICIT_UPDATE);
      }
                         
      if (timeFilterManager->filter_dynamics()) {
        switch (regulatorTarget_) {
        case NONE: {
          regulatorMethod_ = new RegulatorMethod(this);
          break;
        }
        case FIELD: {
          regulatorMethod_ = new DisplacementGlcFiltered(this);
          break;
        }
        case DERIVATIVE: {
          regulatorMethod_ = new VelocityGlcFiltered(this);
          break;
        }
        case DYNAMICS: {
          throw ATC_Error("Kinetostat::initialize - force based kinetostats not yet implemented with time filtering");
          regulatorMethod_ = new StressFluxFiltered(this);
          if (timeFilterManager->end_equilibrate()) {
            StressFlux * myMethod;
            myMethod = dynamic_cast<StressFlux *>(regulatorMethod_);
            myMethod->reset_filtered_ghost_force(nodalGhostForceFiltered);
          }
          break;
        }
        default:
          throw ATC_Error("Unknown kinetostat type in Kinetostat::initialize");
        }
      }
      else {
        switch (regulatorTarget_) {
        case NONE: {
          regulatorMethod_ = new RegulatorMethod(this);
          break;
        }
        case FIELD: {
          regulatorMethod_ = new DisplacementGlc(this);
          break;
        }
        case DERIVATIVE: {
          regulatorMethod_ = new VelocityGlc(this);
          break;
        }
        case DYNAMICS: {
          if (myIntegrationType == TimeIntegrator::FRACTIONAL_STEP) {
            if (couplingMode_ == GHOST_FLUX) {
              regulatorMethod_ = new KinetostatFluxGhost(this);
            }
            else if (couplingMode_ == FIXED) {
              regulatorMethod_ = new KinetostatFixed(this);
            }
            else if (couplingMode_ == FLUX) {
              regulatorMethod_ = new KinetostatFlux(this);
            }
            break;
          }
          if (myIntegrationType == TimeIntegrator::GEAR) {
            if (couplingMode_ == FIXED) {
              regulatorMethod_ = new KinetostatFixed(this);
            }
            else if (couplingMode_ == FLUX) {
              regulatorMethod_ = new KinetostatFlux(this);
            }
            break;
          }
          else {
            if (couplingMode_ == GHOST_FLUX) {
              regulatorMethod_ = new StressFluxGhost(this);
            }
            else {
              regulatorMethod_ = new StressFlux(this);
            }
            break;
          }
        }
        default:
          throw ATC_Error("Unknown kinetostat type in Kinetostat::initialize");
        }
        AtomicRegulator::reset_method();
      }
    }
    else {
      set_all_data_to_used();
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class KinetostatShapeFunction
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and kinetostat data
  //--------------------------------------------------------
  KinetostatShapeFunction::KinetostatShapeFunction(Kinetostat *kinetostat,
                                                   const string & regulatorPrefix) :
    RegulatorShapeFunction(kinetostat,regulatorPrefix),
    kinetostat_(kinetostat),
    timeFilter_(atomicRegulator_->time_filter()),
    nodalAtomicLambdaForce_(NULL),
    lambdaForceFiltered_(NULL),
    atomKinetostatForce_(NULL),
    atomVelocities_(NULL),
    atomMasses_(NULL)
  {
    // data associated with stage 3 in ATC_Method::initialize
    lambda_ = kinetostat->regulator_data(regulatorPrefix_+"LambdaMomentum",nsd_);
    lambdaForceFiltered_ = kinetostat_->regulator_data("LambdaForceFiltered",nsd_);
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void KinetostatShapeFunction::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // needed fundamental quantities
    atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY);
    atomMasses_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS);

    // base class transfers
    RegulatorShapeFunction::construct_transfers();

    // lambda interpolated to the atomic coordinates
    atomLambdas_ = new FtaShapeFunctionProlongation(atc_,
                                                   lambda_,
                                                    interscaleManager.per_atom_sparse_matrix("Interpolant"));
    interscaleManager.add_per_atom_quantity(atomLambdas_,
                                            regulatorPrefix_+"AtomLambdaMomentum");
  }

  //--------------------------------------------------------
  //  set_weights
  //            sets diagonal weighting matrix used in
  //            solve_for_lambda
  //--------------------------------------------------------
  void KinetostatShapeFunction::set_weights()
  {
    if (this->use_local_shape_functions()) {
      ConstantQuantityMapped<double> * myWeights = new ConstantQuantityMapped<double>(atc_,1.,lambdaAtomMap_);
      weights_ = myWeights;
      (atc_->interscale_manager()).add_per_atom_quantity(myWeights,
                                                         "AtomOnesMapped");
    }
    else {
      weights_ = (atc_->interscale_manager()).per_atom_quantity("AtomicOnes");
      if (!weights_) {
        weights_ = new ConstantQuantity<double>(atc_,1.);
        (atc_->interscale_manager()).add_per_atom_quantity(weights_,
                                                           "AtomicOnes");
      }
    }
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
  GlcKinetostat::GlcKinetostat(Kinetostat *kinetostat) :
    KinetostatShapeFunction(kinetostat),
    mdMassMatrix_(atc_->set_mass_mat_md(VELOCITY)),
    atomPositions_(NULL)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void GlcKinetostat::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // needed fundamental quantities
    atomPositions_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_POSITION);
    
    // base class transfers
    KinetostatShapeFunction::construct_transfers();
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void GlcKinetostat::initialize()
  {
    KinetostatShapeFunction::initialize();

    
    // set up list of nodes using Hoover coupling
    // (a) nodes with prescribed values
    PrescribedDataManager * prescribedDataMgr(atc_->prescribed_data_manager());
    for (int i = 0; i < nNodes_; ++i)
      for (int j = 0; j < nsd_; ++j)
        if (prescribedDataMgr->is_fixed(i,VELOCITY,j))
          hooverNodes_.insert(pair<int,int>(i,j));

    // (b) AtC coupling nodes
    if (atomicRegulator_->coupling_mode()==AtomicRegulator::FIXED) {
      InterscaleManager & interscaleManager(atc_->interscale_manager());
      const INT_ARRAY & nodeType((interscaleManager.dense_matrix_int("NodalGeometryType"))->quantity());
      if (atomicRegulator_->use_localized_lambda()) {
        for (int i = 0; i < nNodes_; ++i) {
          if (nodeType(i,0)==BOUNDARY) {
            for (int j = 0; j < nsd_; ++j) {
              hooverNodes_.insert(pair<int,int>(i,j));
            }
          }
        }
      }
      else {
        for (int i = 0; i < nNodes_; ++i) {
          if (nodeType(i,0)==BOUNDARY || nodeType(i,0)==MD_ONLY) {
            for (int j = 0; j < nsd_; ++j) {
              hooverNodes_.insert(pair<int,int>(i,j));
            }
          }
        }
      }
    }
  }

  //--------------------------------------------------------
  //  apply_lambda_to_atoms
  //            uses existing lambda to modify given
  //            atomic quantity
  //--------------------------------------------------------
  void GlcKinetostat::apply_to_atoms(PerAtomQuantity<double> * quantity,
                                     const DENS_MAT & lambdaAtom,
                                     double dt)
  {
    *quantity -= lambdaAtom;
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
    nodalAtomicMassWeightedDisplacement_(NULL),
    nodalDisplacements_(atc_->field(DISPLACEMENT))
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void DisplacementGlc::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // set up node mappings
    create_node_maps();

    // set up shape function matrix
    if (this->use_local_shape_functions()) {
      lambdaAtomMap_ = new AtomToElementset(atc_,elementMask_);
      interscaleManager.add_per_atom_int_quantity(lambdaAtomMap_,
                                                  regulatorPrefix_+"LambdaAtomMap");
      shapeFunctionMatrix_ = new LocalLambdaCouplingMatrix(atc_,
                                                           lambdaAtomMap_,
                                                           nodeToOverlapMap_);
    }
    else {
      shapeFunctionMatrix_ = new LambdaCouplingMatrix(atc_,nodeToOverlapMap_);
    }
    interscaleManager.add_per_atom_sparse_matrix(shapeFunctionMatrix_,
                                                 regulatorPrefix_+"LambdaCouplingMatrixMomentum");

    // set linear solver strategy
    if (atomicRegulator_->use_lumped_lambda_solve()) {
      linearSolverType_ = AtomicRegulator::RSL_SOLVE;
    }
    else {
      linearSolverType_ = AtomicRegulator::CG_SOLVE;
    }
    

    // base class transfers
    GlcKinetostat::construct_transfers();

    // atomic force induced by kinetostat
    atomKinetostatForce_ = new AtomicKinetostatForceDisplacement(atc_);
    interscaleManager.add_per_atom_quantity(atomKinetostatForce_,
                                            regulatorPrefix_+"AtomKinetostatForce");
        
    // restricted force due to kinetostat
    nodalAtomicLambdaForce_ = new AtfShapeFunctionRestriction(atc_,
                                                              atomKinetostatForce_,
                                                              interscaleManager.per_atom_sparse_matrix("Interpolant"));
    interscaleManager.add_dense_matrix(nodalAtomicLambdaForce_,
                                       regulatorPrefix_+"NodalAtomicLambdaForce");
    
    // nodal displacement restricted from atoms
    nodalAtomicMassWeightedDisplacement_ = interscaleManager.dense_matrix("NodalAtomicMassWeightedDisplacement");
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void DisplacementGlc::initialize()
  {
    GlcKinetostat::initialize();

    // sets up time filter for cases where variables temporally filtered
    TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
    if (!timeFilterManager->end_equilibrate()) {
      *lambdaForceFiltered_ = 0.;
      timeFilter_->initialize(lambdaForceFiltered_->quantity());
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
    solve_for_lambda(rhs,lambda_->set_quantity());

    // compute nodal atomic power
    compute_nodal_lambda_force(dt);

    // apply kinetostat to atoms
    apply_to_atoms(atomPositions_,atomLambdas_->quantity());
  }

  //--------------------------------------------------------
  //  set_kinetostat_rhs
  //            sets up the right-hand side of the
  //            kinetostat equations
  //--------------------------------------------------------
  void DisplacementGlc::set_kinetostat_rhs(DENS_MAT & rhs, double dt)
  {
    // form rhs : sum_a (hatN_Ia * x_ai) - (Upsilon)_Ii
    rhs = nodalAtomicMassWeightedDisplacement_->quantity();
    rhs -= ((atc_->mass_mat_md(VELOCITY)).quantity())*(nodalDisplacements_.quantity());
  }

  //--------------------------------------------------------
  //  compute_nodal_lambda_force
  //            compute the effective FE force applied
  //            by the kinetostat
  //--------------------------------------------------------
  void DisplacementGlc::compute_nodal_lambda_force(double dt)
  {
    const DENS_MAT & myNodalAtomicLambdaForce(nodalAtomicLambdaForce_->quantity());
    timeFilter_->apply_post_step1(lambdaForceFiltered_->set_quantity(),
                                  myNodalAtomicLambdaForce,dt);

    // update FE displacements for localized thermostats
    apply_localization_correction(myNodalAtomicLambdaForce,
                                  nodalDisplacements_.set_quantity(),
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
    timeFilter_->apply_pre_step1(lambdaForceFiltered_->set_quantity(),(-1./dt/dt)*lambdaZero,dt);
  }

  //--------------------------------------------------------
  //  set_weights
  //            sets diagonal weighting matrix used in
  //            solve_for_lambda
  //--------------------------------------------------------
  void DisplacementGlc::set_weights()
  {
    if (lambdaAtomMap_) {
      MappedAtomQuantity * myWeights = new MappedAtomQuantity(atc_,atomMasses_,lambdaAtomMap_);
      weights_ = myWeights;
      (atc_->interscale_manager()).add_per_atom_quantity(myWeights,
                                                         "AtomMassesMapped");
    }
    else {
      weights_ = atomMasses_;
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
    
    DENS_MAT nodalLambdaRoc(nNodes_,nsd_);
    atc_->apply_inverse_mass_matrix(source,
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
  void DisplacementGlc::output(OUTPUT_LIST & outputData)
  {
    _nodalAtomicLambdaForceOut_ = nodalAtomicLambdaForce_->quantity();
    DENS_MAT & lambda(lambda_->set_quantity());
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData[regulatorPrefix_+"Lambda"] = &lambda;
      outputData[regulatorPrefix_+"NodalLambdaForce"] = &(_nodalAtomicLambdaForceOut_);
    }
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
    nodalAtomicDisplacements_(atc_->nodal_atomic_field(DISPLACEMENT))
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
    atc_->apply_inverse_md_mass_matrix(lambdaForceFiltered_->set_quantity(),
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
    double coef = 1./(timeFilter_->unfiltered_coefficient_pre_s1(dt));
    rhs = coef*((atc_->mass_mat_md(VELOCITY)).quantity())*(nodalAtomicDisplacements_.quantity() - nodalDisplacements_.quantity());
  }

  //--------------------------------------------------------
  //  compute_nodal_lambda_force
  //            compute the effective FE force applied
  //            by the kinetostat
  //--------------------------------------------------------
  void DisplacementGlcFiltered::compute_nodal_lambda_force(double dt)
  {
    const DENS_MAT & myNodalAtomicLambdaForce(nodalAtomicLambdaForce_->quantity());
    DENS_MAT & myLambdaForceFiltered(lambdaForceFiltered_->set_quantity());
    timeFilter_->apply_post_step1(myLambdaForceFiltered,
                                  myNodalAtomicLambdaForce,dt);

    // update filtered atomic displacements
    DENS_MAT nodalLambdaRoc(myNodalAtomicLambdaForce.nRows(),myNodalAtomicLambdaForce.nCols());
    atc_->apply_inverse_md_mass_matrix(myNodalAtomicLambdaForce,
                                       nodalLambdaRoc,
                                       VELOCITY);
    timeFilter_->apply_post_step1(nodalAtomicDisplacements_.set_quantity(),dt*dt*nodalLambdaRoc,dt);

    // update FE displacements for localized thermostats
    apply_localization_correction(myLambdaForceFiltered,
                                  nodalDisplacements_.set_quantity(),
                                  dt*dt);
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void DisplacementGlcFiltered::output(OUTPUT_LIST & outputData)
  {
    DENS_MAT & lambda(lambda_->set_quantity());
    DENS_MAT & lambdaForceFiltered(lambdaForceFiltered_->set_quantity());
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData[regulatorPrefix_+"Lambda"] = &lambda;
      outputData[regulatorPrefix_+"NodalLambdaForce"] = &lambdaForceFiltered;
    }
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
    nodalAtomicMomentum_(NULL),
    nodalVelocities_(atc_->field(VELOCITY))
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void VelocityGlc::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // set up node mappings
    create_node_maps();

    // set up shape function matrix
    if (this->use_local_shape_functions()) {
      lambdaAtomMap_ = new AtomToElementset(atc_,elementMask_);
      interscaleManager.add_per_atom_int_quantity(lambdaAtomMap_,
                                                  regulatorPrefix_+"LambdaAtomMap");
      shapeFunctionMatrix_ = new LocalLambdaCouplingMatrix(atc_,
                                                           lambdaAtomMap_,
                                                           nodeToOverlapMap_);
    }
    else {
      shapeFunctionMatrix_ = new LambdaCouplingMatrix(atc_,nodeToOverlapMap_);
    }
    interscaleManager.add_per_atom_sparse_matrix(shapeFunctionMatrix_,
                                                 regulatorPrefix_+"LambdaCouplingMatrixMomentum");

    // set linear solver strategy
    if (atomicRegulator_->use_lumped_lambda_solve()) {
      linearSolverType_ = AtomicRegulator::RSL_SOLVE;
    }
    else {
      linearSolverType_ = AtomicRegulator::CG_SOLVE;
    }
    
    // base class transfers
    GlcKinetostat::construct_transfers();

    // atomic force induced by kinetostat
    atomKinetostatForce_ = new AtomicKinetostatForceVelocity(atc_);
    interscaleManager.add_per_atom_quantity(atomKinetostatForce_,
                                            regulatorPrefix_+"AtomKinetostatForce");
        
    // restricted force due to kinetostat
    nodalAtomicLambdaForce_ = new AtfShapeFunctionRestriction(atc_,
                                                              atomKinetostatForce_,
                                                              interscaleManager.per_atom_sparse_matrix("Interpolant"));
    interscaleManager.add_dense_matrix(nodalAtomicLambdaForce_,
                                            regulatorPrefix_+"NodalAtomicLambdaForce");
    
    // nodal momentum restricted from atoms
    nodalAtomicMomentum_ = interscaleManager.dense_matrix("NodalAtomicMomentum");

    
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void VelocityGlc::initialize()
  {
    GlcKinetostat::initialize();

    // sets up time filter for cases where variables temporally filtered
    TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
    if (!timeFilterManager->end_equilibrate()) {
      lambdaForceFiltered_->set_quantity() = 0.;
      timeFilter_->initialize(lambdaForceFiltered_->quantity());
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
    compute_kinetostat(dtLambda);
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
    timeFilter_->apply_pre_step1(lambdaForceFiltered_->set_quantity(),(-1./dt)*lambdaZero,dt);
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
    solve_for_lambda(rhs,lambda_->set_quantity());

    // compute nodal atomic power
    compute_nodal_lambda_force(dt);

    // apply kinetostat to atoms
    apply_to_atoms(atomVelocities_,atomLambdas_->quantity());
  }

  //--------------------------------------------------------
  //  set_kinetostat_rhs
  //            sets up the right-hand side of the
  //            kinetostat equations
  //--------------------------------------------------------
  void VelocityGlc::set_kinetostat_rhs(DENS_MAT & rhs, double dt)
  {
    // form rhs : sum_a (hatN_Ia * x_ai) - (\dot{Upsilon})_Ii
    rhs = nodalAtomicMomentum_->quantity();
    rhs -= ((atc_->mass_mat_md(VELOCITY)).quantity())*(nodalVelocities_.quantity());
  }

  //--------------------------------------------------------
  //  compute_nodal_lambda_force
  //            compute the effective FE force applied
  //            by the kinetostat
  //--------------------------------------------------------
  void VelocityGlc::compute_nodal_lambda_force(double dt)
  {
    const DENS_MAT & myNodalAtomicLambdaForce(nodalAtomicLambdaForce_->quantity());

    timeFilter_->apply_pre_step1(lambdaForceFiltered_->set_quantity(),
                                 myNodalAtomicLambdaForce,dt);

    // update FE displacements for localized thermostats
    apply_localization_correction(myNodalAtomicLambdaForce,
                                  nodalVelocities_.set_quantity(),
                                  dt);
  }

  //--------------------------------------------------------
  //  set_weights
  //            sets diagonal weighting matrix used in
  //            solve_for_lambda
  //--------------------------------------------------------
  void VelocityGlc::set_weights()
  {
    if (lambdaAtomMap_) {
      MappedAtomQuantity * myWeights = new MappedAtomQuantity(atc_,atomMasses_,lambdaAtomMap_);
      weights_ = myWeights;
      (atc_->interscale_manager()).add_per_atom_quantity(myWeights,
                                                         "AtomMassesMapped");
    }
    else {
      weights_ = atomMasses_;
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
    
    DENS_MAT nodalLambdaRoc(nNodes_,nsd_);
    atc_->apply_inverse_mass_matrix(source,
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
  void VelocityGlc::output(OUTPUT_LIST & outputData)
  {
    _nodalAtomicLambdaForceOut_ = nodalAtomicLambdaForce_->quantity();
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData[regulatorPrefix_+"Lambda"] = &(lambda_->set_quantity());
      outputData[regulatorPrefix_+"NodalLambdaForce"] = &(_nodalAtomicLambdaForceOut_);
    }
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
      nodalAtomicVelocities_(atc_->nodal_atomic_field(VELOCITY))
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
    atc_->apply_inverse_md_mass_matrix(lambdaForceFiltered_->quantity(),
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
    double coef = 1./(timeFilter_->unfiltered_coefficient_pre_s1(dt));
    rhs = coef*((atc_->mass_mat_md(VELOCITY)).quantity())*(nodalAtomicVelocities_.quantity() - nodalVelocities_.quantity());
  }

  //--------------------------------------------------------
  //  compute_nodal_lambda_force
  //            compute the effective FE force applied
  //            by the kinetostat
  //--------------------------------------------------------
  void VelocityGlcFiltered::compute_nodal_lambda_force(double dt)
  {
    const DENS_MAT & myNodalAtomicLambdaForce(nodalAtomicLambdaForce_->quantity());
    DENS_MAT & myLambdaForceFiltered(lambdaForceFiltered_->set_quantity());
    timeFilter_->apply_post_step1(myLambdaForceFiltered,myNodalAtomicLambdaForce,dt);

    // update filtered atomic displacements
    DENS_MAT nodalLambdaRoc(myNodalAtomicLambdaForce.nRows(),myNodalAtomicLambdaForce.nCols());
    atc_->apply_inverse_md_mass_matrix(myNodalAtomicLambdaForce,
                                       nodalLambdaRoc,
                                       VELOCITY);
    timeFilter_->apply_post_step1(nodalAtomicVelocities_.set_quantity(),dt*nodalLambdaRoc,dt);

    // update FE displacements for localized thermostats
    apply_localization_correction(myLambdaForceFiltered,
                                  nodalVelocities_.set_quantity(),
                                  dt);
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void VelocityGlcFiltered::output(OUTPUT_LIST & outputData)
  {
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData[regulatorPrefix_+"Lambda"] = &(lambda_->set_quantity());
      outputData[regulatorPrefix_+"NodalLambdaForce"] = &(lambdaForceFiltered_->set_quantity());
    }
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
    nodalForce_(atc_->field_rhs(VELOCITY)),
    nodalAtomicForce_(NULL),
    nodalGhostForce_(NULL),
    momentumSource_(atc_->atomic_source(VELOCITY))
  {
    // flag for performing boundary flux calculation
    fieldMask_(VELOCITY,FLUX) = true;
  }

  StressFlux::~StressFlux()
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void StressFlux::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // set up node mappings
    create_node_maps();

    // set up shape function matrix
    if (this->use_local_shape_functions()) {
      lambdaAtomMap_ = new AtomToElementset(atc_,elementMask_);
      interscaleManager.add_per_atom_int_quantity(lambdaAtomMap_,
                                                  regulatorPrefix_+"LambdaAtomMap");
      shapeFunctionMatrix_ = new LocalLambdaCouplingMatrix(atc_,
                                                           lambdaAtomMap_,
                                                           nodeToOverlapMap_);
    }
    else {
      shapeFunctionMatrix_ = new LambdaCouplingMatrix(atc_,nodeToOverlapMap_);
    }
    interscaleManager.add_per_atom_sparse_matrix(shapeFunctionMatrix_,
                                                 regulatorPrefix_+"LambdaCouplingMatrixMomentum");

    // set linear solver strategy
    if (atomicRegulator_->use_lumped_lambda_solve()) {
      linearSolverType_ = AtomicRegulator::RSL_SOLVE;
    }
    else {
      linearSolverType_ = AtomicRegulator::CG_SOLVE;
    }

    // base class transfers
    GlcKinetostat::construct_transfers();

    // force at nodes due to atoms
    nodalAtomicForce_ = interscaleManager.dense_matrix("NodalAtomicForce");

    // atomic force induced by kinetostat
    atomKinetostatForce_ = new AtomicKinetostatForceStress(atc_,atomLambdas_);
    interscaleManager.add_per_atom_quantity(atomKinetostatForce_,
                                            regulatorPrefix_+"AtomKinetostatForce");
        
    // restricted force due to kinetostat
    nodalAtomicLambdaForce_ = new AtfShapeFunctionRestriction(atc_,
                                                              atomKinetostatForce_,
                                                              interscaleManager.per_atom_sparse_matrix("Interpolant"));
    interscaleManager.add_dense_matrix(nodalAtomicLambdaForce_,
                                       regulatorPrefix_+"NodalAtomicLambdaForce");

    // sets up space for ghost force related variables
    if (atc_->groupbit_ghost()) {
      GhostCouplingMatrix * shapeFunctionGhost = new GhostCouplingMatrix(atc_,interscaleManager.per_atom_sparse_matrix("InterpolantGhost"),
                                                                         regulatedNodes_,nodeToOverlapMap_);
      interscaleManager.add_sparse_matrix(shapeFunctionGhost,
                                          regulatorPrefix_+"GhostCouplingMatrix");
      FundamentalAtomQuantity * atomGhostForce = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_FORCE,
                                                                                                 GHOST);
      nodalGhostForce_ = new AtfShapeFunctionRestriction(atc_,atomGhostForce,
                                                         shapeFunctionGhost);
      interscaleManager.add_dense_matrix(nodalGhostForce_,
                                         regulatorPrefix_+"NodalGhostForce");
      nodalGhostForceFiltered_.reset(nNodes_,nsd_);
    }
  }

  //--------------------------------------------------------
  //  compute_boundary_flux:
  //    computes the boundary flux to be consistent with
  //    the controller
  //--------------------------------------------------------
  void StressFlux::compute_boundary_flux(FIELDS & fields)
  {
    GlcKinetostat::compute_boundary_flux(fields);
  }

  //--------------------------------------------------------
  //  apply_pre_predictor:
  //    apply the kinetostat to the atoms in the
  //    mid-predictor integration phase
  //--------------------------------------------------------
  void StressFlux::apply_pre_predictor(double dt)
  {
    double dtLambda = 0.5*dt;
    // apply lambda force to atoms
    apply_to_atoms(atomVelocities_,atomKinetostatForce_->quantity(),dtLambda);
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
    apply_to_atoms(atomVelocities_,atomKinetostatForce_->quantity(),dtLambda);
  }

  //--------------------------------------------------------
  //  compute_kinetostat
  //            manages the solution and application of the
  //            kinetostat equations and variables
  //--------------------------------------------------------
  void StressFlux::compute_kinetostat(double dt)
  {
    // initial filtering update
    apply_pre_filtering(dt);

    // set up rhs
    DENS_MAT rhs(nNodes_,nsd_);
    set_kinetostat_rhs(rhs,dt);

    // solve linear system for lambda
    solve_for_lambda(rhs,lambda_->set_quantity());

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
    timeFilter_->apply_pre_step1(lambdaForceFiltered_->set_quantity(),lambdaZero,dt);
    if (nodalGhostForce_) {
      timeFilter_->apply_pre_step1(nodalGhostForceFiltered_.set_quantity(),
                                   nodalGhostForce_->quantity(),dt);
    }
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
    rhs = momentumSource_.quantity();
    if (nodalGhostForce_) {
      rhs -= nodalGhostForce_->quantity();
    }
    
    // (b) for ess. bcs
    
    // form rhs : {sum_a (N_Ia * f_ia) - M_md * (ddupsilon/dt)_I}
    DENS_MAT rhsPrescribed = -1.*nodalForce_.quantity();
    atc_->apply_inverse_mass_matrix(rhsPrescribed,VELOCITY);
    rhsPrescribed = (mdMassMatrix_.quantity())*rhsPrescribed;
    rhsPrescribed += nodalAtomicForce_->quantity();

    set<pair<int,int> >::const_iterator iter;
    for (iter = hooverNodes_.begin(); iter != hooverNodes_.end(); ++iter) {
      rhs(iter->first,iter->second) = rhsPrescribed(iter->first,iter->second);
    }
  }

  //--------------------------------------------------------
  //  compute_nodal_lambda_force
  //            computes the force induced on the FE
  //            by applying lambdaForce on the atoms
  //--------------------------------------------------------
  void StressFlux::compute_nodal_lambda_force(double dt)
  {
    DENS_MAT myNodalAtomicLambdaForce = nodalAtomicLambdaForce_->quantity();
    set<pair<int,int> >::const_iterator iter;
    for (iter = hooverNodes_.begin(); iter != hooverNodes_.end(); ++iter) {
      myNodalAtomicLambdaForce(iter->first,iter->second) = 0.;
    }

    timeFilter_->apply_post_step1(lambdaForceFiltered_->set_quantity(),
                                  myNodalAtomicLambdaForce,dt);
  }

  //--------------------------------------------------------
  //  add_to_rhs:
  //            determines what if any contributions to the
  //            finite element equations are needed for
  //            consistency with the kinetostat
  //--------------------------------------------------------
  void StressFlux::add_to_rhs(FIELDS & rhs)
  {
    // compute the kinetostat force
    compute_kinetostat(atc_->dt());

    rhs[VELOCITY] += nodalAtomicLambdaForce_->quantity() + boundaryFlux_[VELOCITY].quantity();
  }

  //--------------------------------------------------------
  //  apply_lambda_to_atoms
  //            uses existing lambda to modify given
  //            atomic quantity
  //--------------------------------------------------------
  void StressFlux::apply_to_atoms(PerAtomQuantity<double> * atomVelocities,
                                  const DENS_MAT & lambdaForce,
                                  double dt)
  {
    _deltaVelocity_ = lambdaForce;
    _deltaVelocity_ /= atomMasses_->quantity();
    _deltaVelocity_ *= dt;
    *atomVelocities += _deltaVelocity_;
  }

  //--------------------------------------------------------
  //  reset_filtered_ghost_force:
  //    resets the kinetostat generated ghost force to a
  //    prescribed value
  //--------------------------------------------------------
  void StressFlux::reset_filtered_ghost_force(DENS_MAT & target)
  {
    nodalGhostForceFiltered_ = target;
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void StressFlux::output(OUTPUT_LIST & outputData)
  {
    _nodalAtomicLambdaForceOut_ = nodalAtomicLambdaForce_->quantity();
    DENS_MAT & lambda(lambda_->set_quantity());
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData[regulatorPrefix_+"Lambda"] = &lambda;
      outputData[regulatorPrefix_+"NodalLambdaForce"] = &(_nodalAtomicLambdaForceOut_);
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class StressFluxGhost
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and kinetostat data
  //--------------------------------------------------------
  StressFluxGhost::StressFluxGhost(Kinetostat * kinetostat) :
    StressFlux(kinetostat)
  {
    // flag for performing boundary flux calculation
    fieldMask_(VELOCITY,FLUX) = false;
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void StressFluxGhost::construct_transfers()
  {
    StressFlux::construct_transfers();
    if (!nodalGhostForce_) {
      throw ATC_Error("StressFluxGhost::StressFluxGhost - ghost atoms must be specified");
    }
  }

  //--------------------------------------------------------
  //  compute_boundary_flux:
  //    computes the boundary flux to be consistent with
  //    the controller
  //--------------------------------------------------------
  void StressFluxGhost::compute_boundary_flux(FIELDS & fields)
  {
    // This is only used in computation of atomic sources
    boundaryFlux_[VELOCITY] = 0.;
  }

  //--------------------------------------------------------
  //  add_to_rhs:
  //            determines what if any contributions to the
  //            finite element equations are needed for
  //            consistency with the kinetostat
  //--------------------------------------------------------
  void StressFluxGhost::add_to_rhs(FIELDS & rhs)
  {
    // compute the kinetostat force
    compute_kinetostat(atc_->dt());

    // uses ghost force as the boundary flux to add to the RHS
    rhs[VELOCITY] += nodalAtomicLambdaForce_->quantity() + nodalGhostForce_->quantity();
  }

    //--------------------------------------------------------
  //  set_kinetostat_rhs
  //            sets up the RHS of the kinetostat equations
  //            for the coupling parameter lambda
  //--------------------------------------------------------
  void StressFluxGhost::set_kinetostat_rhs(DENS_MAT & rhs, double dt)
  {
    // (a) for flux based : 
    // form rhs :  \int N_I r dV - \sum_g N_Ig^* f_g
    // sources are set in ATC transfer
    rhs.reset(nNodes_,nsd_);
    rhs = momentumSource_.quantity();
    
    // (b) for ess. bcs
    
    // form rhs : {sum_a (N_Ia * f_ia) - M_md * (ddupsilon/dt)_I}
    DENS_MAT rhsPrescribed = -1.*nodalForce_.quantity();
    atc_->apply_inverse_mass_matrix(rhsPrescribed,VELOCITY);
    rhsPrescribed = (mdMassMatrix_.quantity())*rhsPrescribed;
    rhsPrescribed += nodalAtomicForce_->quantity();

    set<pair<int,int> >::const_iterator iter;
    for (iter = hooverNodes_.begin(); iter != hooverNodes_.end(); ++iter) {
      rhs(iter->first,iter->second) = rhsPrescribed(iter->first,iter->second);
    }
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
    nodalAtomicVelocity_(atc_->nodal_atomic_field(VELOCITY))
  {
    // do nothing
  }

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
    rhs = momentumSource_.quantity() - nodalGhostForceFiltered_.quantity();
    
    // (b) for ess. bcs
    
    // form rhs : {sum_a (N_Ia * f_ia) - M_md * (ddupsilon/dt)_I}
    DENS_MAT rhsPrescribed = -1.*nodalForce_.quantity();
    atc_->apply_inverse_mass_matrix(rhsPrescribed,VELOCITY);
    rhsPrescribed = (mdMassMatrix_.quantity())*rhsPrescribed;
    rhsPrescribed += nodalAtomicForce_->quantity();

    set<pair<int,int> >::const_iterator iter;
    for (iter = hooverNodes_.begin(); iter != hooverNodes_.end(); ++iter) {
      rhs(iter->first,iter->second) = rhsPrescribed(iter->first,iter->second);
    }

    // adjust for application of current lambda force
    rhs += lambdaForceFiltered_->quantity();

    // correct for time filtering
    rhs *= 1./(timeFilter_->unfiltered_coefficient_pre_s1(dt));
  }

  //--------------------------------------------------------
  //  apply_lambda_to_atoms
  //            uses existing lambda to modify given
  //            atomic quantity
  //--------------------------------------------------------
  void StressFluxFiltered::apply_to_atoms(PerAtomQuantity<double> * atomVelocities,
                                          const DENS_MAT & lambdaForce,
                                          double dt)
  {
    StressFlux::apply_to_atoms(atomVelocities,lambdaForce,dt);

    // add in corrections to filtered nodal atomice velocity
    DENS_MAT velocityRoc(nNodes_,nsd_);
    atc_->apply_inverse_md_mass_matrix(lambdaForceFiltered_->quantity(),
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
    // compute kinetostat forces
    compute_kinetostat(atc_->dt());
    rhs[VELOCITY] += lambdaForceFiltered_->quantity() + boundaryFlux_[VELOCITY].quantity();
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void StressFluxFiltered::output(OUTPUT_LIST & outputData)
  {
    DENS_MAT & lambda(lambda_->set_quantity());
    DENS_MAT & lambdaForceFiltered(lambdaForceFiltered_->set_quantity());
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData[regulatorPrefix_+"Lambda"] = &lambda;
      outputData[regulatorPrefix_+"NodalLambdaForce"] = &lambdaForceFiltered;
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class KinetostatGlcFs
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and kinetostat data
  //--------------------------------------------------------
  KinetostatGlcFs::KinetostatGlcFs(Kinetostat * kinetostat,
                                   const string & regulatorPrefix) :
    KinetostatShapeFunction(kinetostat,regulatorPrefix),
    velocity_(atc_->field(VELOCITY))
  {
    // constuct/obtain data corresponding to stage 3 of ATC_Method::initialize
    nodalAtomicLambdaForce_ = kinetostat_->regulator_data(regulatorPrefix_+"NodalAtomicLambdaForce",nsd_);
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void KinetostatGlcFs::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // base class transfers
    KinetostatShapeFunction::construct_transfers();

    // get data from manager
    nodalAtomicMomentum_ = interscaleManager.dense_matrix("NodalAtomicMomentum");

    // atomic force induced by kinetostat
    PerAtomQuantity<double> * atomLambdas = interscaleManager.per_atom_quantity(regulatorPrefix_+"AtomLambdaMomentum");
    atomKinetostatForce_ = new AtomicKinetostatForceStress(atc_,atomLambdas);
    interscaleManager.add_per_atom_quantity(atomKinetostatForce_,
                                            regulatorPrefix_+"AtomKinetostatForce");
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void KinetostatGlcFs::initialize()
  {
    KinetostatShapeFunction::initialize();

    // set up workspaces
    _deltaMomentum_.reset(nNodes_,nsd_);
    _lambdaForceOutput_.reset(nNodes_,nsd_);
    
    TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
    if (!timeFilterManager->end_equilibrate()) {
      // we should reset lambda and lambdaForce to zero in this case
      // implies an initial condition of 0 for the filtered nodal lambda power
      // initial conditions will always be needed when using time filtering
      // however, the fractional step scheme must assume the instantaneous
      // nodal lambda power is 0 initially because all quantities are in delta form
      *lambda_ = 0.; // ensures initial lambda force is zero
      *nodalAtomicLambdaForce_ = 0.; // momentum change due to kinetostat
      *lambdaForceFiltered_ = 0.; // filtered momentum change due to kinetostats
    }
    else {
      // we can grab lambda power variables using time integrator and atc transfer in cases for equilibration
    }

    // sets up time filter for cases where variables temporally filtered
    if (timeFilterManager->need_reset()) {
      // the form of this integrator implies no time filters that require history data can be used
      timeFilter_->initialize(nodalAtomicLambdaForce_->quantity());
    }

    compute_rhs_map();
  }

  //--------------------------------------------------------
  //  compute_rhs_map
  //    creates mapping from all nodes to those to which
  //    the kinetostat applies
  //--------------------------------------------------------
  void KinetostatGlcFs::compute_rhs_map()
  {
    rhsMap_.resize(overlapToNodeMap_->nRows(),1);
    DENS_MAT rhsMapGlobal(nNodes_,1);
    const set<int> & applicationNodes(applicationNodes_->quantity());
    for (int i = 0; i < nNodes_; i++) {
      if (applicationNodes.find(i) != applicationNodes.end()) {
        rhsMapGlobal(i,0) = 1.;
      }
      else {
        rhsMapGlobal(i,0) = 0.;
      }
    }
    map_unique_to_overlap(rhsMapGlobal,rhsMap_);
  }

  //--------------------------------------------------------
  //  apply_pre_predictor:
  //    apply the kinetostat to the atoms in the
  //    pre-predictor integration phase
  //--------------------------------------------------------
  void KinetostatGlcFs::apply_pre_predictor(double dt)
  {
    DENS_MAT & lambdaForceFiltered(lambdaForceFiltered_->set_quantity());
    DENS_MAT & nodalAtomicLambdaForce(nodalAtomicLambdaForce_->set_quantity());

    // update filtered forces
    timeFilter_->apply_pre_step1(lambdaForceFiltered,nodalAtomicLambdaForce,dt);

    // apply lambda force to atoms and compute instantaneous lambda force 
    this->apply_to_atoms(atomVelocities_,nodalAtomicMomentum_,
                         atomKinetostatForce_->quantity(),
                         nodalAtomicLambdaForce,0.5*dt);

    // update nodal variables for first half of timestep
    this->add_to_momentum(nodalAtomicLambdaForce,_deltaMomentum_,0.5*dt);
    atc_->apply_inverse_mass_matrix(_deltaMomentum_,VELOCITY);
    velocity_ += _deltaMomentum_;

    // start update of filtered lambda force
    nodalAtomicLambdaForce = 0.;
    timeFilter_->apply_post_step1(lambdaForceFiltered,nodalAtomicLambdaForce,dt);
  }
 
  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the kinetostat to the atoms in the
  //    post-corrector integration phase
  //--------------------------------------------------------
  void KinetostatGlcFs::apply_post_corrector(double dt)
  {
    // compute the kinetostat equation and update lambda
    this->compute_lambda(dt);

    DENS_MAT & lambdaForceFiltered(lambdaForceFiltered_->set_quantity());
    DENS_MAT & nodalAtomicLambdaForce(nodalAtomicLambdaForce_->set_quantity());

    // update filtered force
    timeFilter_->apply_pre_step1(lambdaForceFiltered,nodalAtomicLambdaForce,dt);

    // apply lambda force to atoms and compute instantaneous lambda force 
    this->apply_to_atoms(atomVelocities_,nodalAtomicMomentum_,
                         atomKinetostatForce_->quantity(),
                         nodalAtomicLambdaForce,0.5*dt);

    // update nodal variables for first half of timestep
    this->add_to_momentum(nodalAtomicLambdaForce,_deltaMomentum_,0.5*dt);
    nodalAtomicLambdaForce *= 2./dt;
    atc_->apply_inverse_mass_matrix(_deltaMomentum_,VELOCITY);
    velocity_ += _deltaMomentum_;

    // start update of filtered lambda force
    timeFilter_->apply_post_step2(lambdaForceFiltered,nodalAtomicLambdaForce,dt);
  }

  //--------------------------------------------------------
  //  compute_kinetostat
  //            manages the solution and application of the
  //            kinetostat equations and variables
  //--------------------------------------------------------
  void KinetostatGlcFs::compute_lambda(double dt)
  {
    // set up rhs for lambda equation
    this->set_kinetostat_rhs(rhs_,0.5*dt);

    // solve linear system for lambda
    DENS_MAT & lambda(lambda_->set_quantity());
    solve_for_lambda(rhs_,lambda);
  }

  //--------------------------------------------------------
  //  apply_lambda_to_atoms
  //            uses existing lambda to modify given
  //            atomic quantity
  //--------------------------------------------------------
  void KinetostatGlcFs::apply_to_atoms(PerAtomQuantity<double> * atomVelocity,
                                      const DENS_MAN * nodalAtomicMomentum,
                                      const DENS_MAT & lambdaForce,
                                      DENS_MAT & nodalAtomicLambdaForce,
                                      double dt)
  {
    // compute initial contributions to lambda force
    nodalAtomicLambdaForce = nodalAtomicMomentum->quantity();
    nodalAtomicLambdaForce *= -1.;

    // apply lambda force to atoms
    _velocityDelta_ = lambdaForce;
    _velocityDelta_ /= atomMasses_->quantity();
    _velocityDelta_ *= dt;
    (*atomVelocity) += _velocityDelta_;

    // finalize lambda force
    nodalAtomicLambdaForce += nodalAtomicMomentum->quantity();
  }

  //--------------------------------------------------------
  //  output:
  //    adds all relevant output to outputData
  //--------------------------------------------------------
  void KinetostatGlcFs::output(OUTPUT_LIST & outputData)
  {
    _lambdaForceOutput_ = nodalAtomicLambdaForce_->quantity();
    // approximate value for lambda force
    double dt =  LammpsInterface::instance()->dt();
    _lambdaForceOutput_ *= (2./dt);
    DENS_MAT & lambda(lambda_->set_quantity());
    if ((atc_->lammps_interface())->rank_zero()) {
      outputData[regulatorPrefix_+"Lambda"] = &lambda;
      outputData[regulatorPrefix_+"NodalLambdaForce"] = &(_lambdaForceOutput_);
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class KinetostatFlux
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and kinetostat data
  //--------------------------------------------------------
  KinetostatFlux::KinetostatFlux(Kinetostat * kinetostat,
                                 const string & regulatorPrefix) :
    KinetostatGlcFs(kinetostat,regulatorPrefix),
    momentumSource_(atc_->atomic_source(VELOCITY)),
    nodalGhostForce_(NULL),
    nodalGhostForceFiltered_(NULL)
  {
    // flag for performing boundary flux calculation
    fieldMask_(VELOCITY,FLUX) = true;

    // constuct/obtain data corresponding to stage 3 of ATC_Method::initialize
    nodalGhostForceFiltered_ = kinetostat_->regulator_data(regulatorPrefix_+"NodalGhostForceFiltered",nsd_);
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void KinetostatFlux::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // set up node mappings
    create_node_maps();

    
    // set up linear solver
    // set up data for linear solver
    shapeFunctionMatrix_ = new LambdaCouplingMatrix(atc_,nodeToOverlapMap_);
    interscaleManager.add_per_atom_sparse_matrix(shapeFunctionMatrix_,
                                                 regulatorPrefix_+"LambdaCouplingMatrixMomentum");
    if (elementMask_) {
      lambdaAtomMap_ = new AtomToElementset(atc_,elementMask_);
      interscaleManager.add_per_atom_int_quantity(lambdaAtomMap_,
                                                  regulatorPrefix_+"LambdaAtomMap");
    }
    if (atomicRegulator_->use_localized_lambda()) {
      linearSolverType_ = AtomicRegulator::RSL_SOLVE;
    }
    else {
      linearSolverType_ = AtomicRegulator::CG_SOLVE;
    }

    // base class transfers
    KinetostatGlcFs::construct_transfers();

    // sets up space for ghost force related variables
    if (atc_->groupbit_ghost()) {
      MatrixDependencyManager<DenseMatrix, int> * nodeToOverlapMap = 
        interscaleManager.dense_matrix_int(regulatorPrefix_+"NodeToOverlapMap");
      GhostCouplingMatrix * shapeFunctionGhost = new GhostCouplingMatrix(atc_,interscaleManager.per_atom_sparse_matrix("InterpolantGhost"),
                                                                         regulatedNodes_,
                                                                         nodeToOverlapMap);
      interscaleManager.add_sparse_matrix(shapeFunctionGhost,
                                          regulatorPrefix_+"GhostCouplingMatrix");
      FundamentalAtomQuantity * atomGhostForce = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_FORCE,
                                                                                             GHOST);
      nodalGhostForce_ = new AtfShapeFunctionRestriction(atc_,atomGhostForce,
                                                         shapeFunctionGhost);
      interscaleManager.add_dense_matrix(nodalGhostForce_,
                                         regulatorPrefix_+"NodalGhostForce");
    }
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void KinetostatFlux::initialize()
  {
    KinetostatGlcFs::initialize();
    
    TimeFilterManager * timeFilterManager = atc_->time_filter_manager();
    if (!timeFilterManager->end_equilibrate()) {
      // we should reset lambda and lambdaForce to zero in this case
      // implies an initial condition of 0 for the filtered nodal lambda power
      // initial conditions will always be needed when using time filtering
      // however, the fractional step scheme must assume the instantaneous
      // nodal lambda power is 0 initially because all quantities are in delta form
      *nodalGhostForceFiltered_ = 0.; // filtered force from ghost atoms
    }
    else {
      // we can grab lambda power variables using time integrator and atc transfer in cases for equilibration
    }
  }

  //--------------------------------------------------------
  //  construct_regulated_nodes:
  //    constructs the set of nodes being regulated
  //--------------------------------------------------------
  void KinetostatFlux::construct_regulated_nodes()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // matrix requires all entries even if localized for correct lumping
    regulatedNodes_ = new RegulatedNodes(atc_);
    interscaleManager.add_set_int(regulatedNodes_,
                                  regulatorPrefix_+"KinetostatRegulatedNodes");

    // if localized monitor nodes with applied fluxes
    if (atomicRegulator_->use_localized_lambda()) {
      if ((kinetostat_->coupling_mode() == Kinetostat::FLUX) && (atomicRegulator_->boundary_integration_type() != NO_QUADRATURE)) {
        // include boundary nodes
        applicationNodes_ = new FluxBoundaryNodes(atc_);
        
        boundaryNodes_ = new BoundaryNodes(atc_);
        interscaleManager.add_set_int(boundaryNodes_,
                                      regulatorPrefix_+"KinetostatBoundaryNodes");
      }
      else {
        // fluxed nodes only
        applicationNodes_ = new FluxNodes(atc_);
      }
      interscaleManager.add_set_int(applicationNodes_,
                                    regulatorPrefix_+"KinetostatApplicationNodes");
    }
    else {
      applicationNodes_ = regulatedNodes_;
    }

    // special set of boundary elements for boundary flux quadrature  
    if ((atomicRegulator_->boundary_integration_type() == FE_INTERPOLATION)
        && (atomicRegulator_->use_localized_lambda())) {
      elementMask_ = new ElementMaskNodeSet(atc_,applicationNodes_);
      interscaleManager.add_dense_matrix_bool(elementMask_,
                                                   regulatorPrefix_+"BoundaryElementMask");
    }
  }

  //--------------------------------------------------------
  //  apply_pre_predictor:
  //    apply the kinetostat to the atoms in the
  //    pre-predictor integration phase
  //--------------------------------------------------------
  void KinetostatFlux::apply_pre_predictor(double dt)
  {
    // update filtered forces
    if (nodalGhostForce_) {
      timeFilter_->apply_pre_step1(nodalGhostForceFiltered_->set_quantity(),
                                   nodalGhostForce_->quantity(),dt);
    }

    KinetostatGlcFs::apply_pre_predictor(dt);
  }
 
  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the kinetostat to the atoms in the
  //    post-corrector integration phase
  //--------------------------------------------------------
  void KinetostatFlux::apply_post_corrector(double dt)
  {
    // update filtered ghost force
    if (nodalGhostForce_) {
      timeFilter_->apply_post_step1(nodalGhostForceFiltered_->set_quantity(),
                                    nodalGhostForce_->quantity(),dt);
    }
    
    // compute the kinetostat equation and update lambda
    KinetostatGlcFs::apply_post_corrector(dt);
  }

  //--------------------------------------------------------
  //  set_kinetostat_rhs
  //            sets up the RHS of the kinetostat equations
  //            for the coupling parameter lambda
  //--------------------------------------------------------
  void KinetostatFlux::set_kinetostat_rhs(DENS_MAT & rhs, double dt)
  {
    // (a) for flux based : 
    // form rhs :  \int N_I r dV - \sum_g N_Ig^* f_g
    // sources are set in ATC transfer
    rhs.reset(nNodes_,nsd_);
    const DENS_MAT & momentumSource(momentumSource_.quantity());
    const set<int> & applicationNodes(applicationNodes_->quantity());
    set<int>::const_iterator iNode;
    for (iNode = applicationNodes.begin(); iNode != applicationNodes.end(); iNode++) {
      for (int j = 0; j < nsd_; j++) {
        rhs(*iNode,j) = momentumSource(*iNode,j);
      }
    }

    // add ghost forces, if needed
    if (nodalGhostForce_) {
      const DENS_MAT & nodalGhostForce(nodalGhostForce_->quantity());
      for (iNode = applicationNodes.begin(); iNode != applicationNodes.end(); iNode++) {
        for (int j = 0; j < nsd_; j++) {
          rhs(*iNode,j) -= nodalGhostForce(*iNode,j);
        }
      }
    }
  }

  //--------------------------------------------------------
  //  add_to_momentum:
  //            determines what if any contributions to the
  //            finite element equations are needed for
  //            consistency with the kinetostat
  //--------------------------------------------------------
  void KinetostatFlux::add_to_momentum(const DENS_MAT & nodalLambdaForce,
                                       DENS_MAT & deltaMomentum,
                                       double dt)
  {
    deltaMomentum.resize(nNodes_,nsd_);
    const DENS_MAT & boundaryFlux(boundaryFlux_[VELOCITY].quantity());
    for (int i = 0; i < nNodes_; i++) {
      for (int j = 0; j < nsd_; j++) {
        deltaMomentum(i,j) = nodalLambdaForce(i,j) + dt*boundaryFlux(i,j);
      }
    }
  }

  //--------------------------------------------------------
  //  reset_filtered_ghost_force:
  //    resets the kinetostat generated ghost force to a
  //    prescribed value
  //--------------------------------------------------------
  void KinetostatFlux::reset_filtered_ghost_force(DENS_MAT & target)
  {
    (*nodalGhostForceFiltered_) = target;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class KinetostatFluxGhost
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and kinetostat data
  //--------------------------------------------------------
  KinetostatFluxGhost::KinetostatFluxGhost(Kinetostat * kinetostat,
                                           const string & regulatorPrefix) :
    KinetostatFlux(kinetostat,regulatorPrefix)
  {
    // flag for performing boundary flux calculation
    fieldMask_(VELOCITY,FLUX) = false;
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void KinetostatFluxGhost::construct_transfers()
  {
    KinetostatFlux::construct_transfers();
    if (!nodalGhostForce_) {
      throw ATC_Error("StressFluxGhost::StressFluxGhost - ghost atoms must be specified");
    }
  }

  //--------------------------------------------------------
  //  compute_boundary_flux:
  //    computes the boundary flux to be consistent with
  //    the controller
  //--------------------------------------------------------
  void KinetostatFluxGhost::compute_boundary_flux(FIELDS & fields)
  {
    // This is only used in computation of atomic sources
    boundaryFlux_[VELOCITY] = 0.;
  }

  //--------------------------------------------------------
  //  add_to_momentum:
  //            determines what if any contributions to the
  //            finite element equations are needed for
  //            consistency with the kinetostat
  //--------------------------------------------------------
  void KinetostatFluxGhost::add_to_momentum(const DENS_MAT & nodalLambdaForce,
                                       DENS_MAT & deltaMomentum,
                                       double dt)
  {
    deltaMomentum.resize(nNodes_,nsd_);
    const DENS_MAT & boundaryFlux(nodalGhostForce_->quantity());
    for (int i = 0; i < nNodes_; i++) {
      for (int j = 0; j < nsd_; j++) {
        deltaMomentum(i,j) = nodalLambdaForce(i,j) + dt*boundaryFlux(i,j);
      }
    }
  }

  //--------------------------------------------------------
  //  set_kinetostat_rhs
  //            sets up the RHS of the kinetostat equations
  //            for the coupling parameter lambda
  //--------------------------------------------------------
  void KinetostatFluxGhost::set_kinetostat_rhs(DENS_MAT & rhs, double dt)
  {
    // (a) for flux based : 
    // form rhs :  \int N_I r dV - \sum_g N_Ig^* f_g
    // sources are set in ATC transfer
    rhs.reset(nNodes_,nsd_);
    const DENS_MAT & momentumSource(momentumSource_.quantity());
    const set<int> & applicationNodes(applicationNodes_->quantity());
    set<int>::const_iterator iNode;
    for (iNode = applicationNodes.begin(); iNode != applicationNodes.end(); iNode++) {
      for (int j = 0; j < nsd_; j++) {
        rhs(*iNode,j) = momentumSource(*iNode,j);
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class KinetostatFixed
  //--------------------------------------------------------
  //--------------------------------------------------------
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and kinetostat data
  //--------------------------------------------------------
  KinetostatFixed::KinetostatFixed(Kinetostat * kinetostat,
                                   const string & regulatorPrefix) :
    KinetostatGlcFs(kinetostat,regulatorPrefix),
    mdMassMatrix_(atc_->set_mass_mat_md(VELOCITY)),
    isFirstTimestep_(true)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  constructor_transfers
  //    instantiates or obtains all dependency managed data
  //--------------------------------------------------------
  void KinetostatFixed::construct_transfers()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    // set up node mappings
    create_node_maps();

    // determine if map is needed and set up if so
    if (this->use_local_shape_functions()) {
      lambdaAtomMap_ = new AtomToElementset(atc_,elementMask_);
      interscaleManager.add_per_atom_int_quantity(lambdaAtomMap_,
                                                  regulatorPrefix_+"LambdaAtomMap");
      shapeFunctionMatrix_ = new LocalLambdaCouplingMatrix(atc_,
                                                           lambdaAtomMap_,
                                                           nodeToOverlapMap_);
    }
    else {
      shapeFunctionMatrix_ = new LambdaCouplingMatrix(atc_,nodeToOverlapMap_);
    }
    interscaleManager.add_per_atom_sparse_matrix(shapeFunctionMatrix_,
                                                 regulatorPrefix_+"LambdaCouplingMatrixMomentum");
    linearSolverType_ = AtomicRegulator::CG_SOLVE;

    // base class transfers
    KinetostatGlcFs::construct_transfers();
  }

  //--------------------------------------------------------
  //  initialize
  //    initializes all method data
  //--------------------------------------------------------
  void KinetostatFixed::initialize()
  {
    KinetostatGlcFs::initialize();

    // reset data to zero
    deltaFeMomentum_.reset(nNodes_,nsd_);
    deltaNodalAtomicMomentum_.reset(nNodes_,nsd_);
  }

  //--------------------------------------------------------
  //  halve_force:
  //    flag to halve the lambda force for improved
  //    accuracy
  //--------------------------------------------------------
  bool KinetostatFixed::halve_force()
  {
    if (isFirstTimestep_ || ((atc_->atom_to_element_map_type() == EULERIAN)
                             && (atc_->atom_to_element_map_frequency() > 1)
                             && (atc_->step() % atc_->atom_to_element_map_frequency() == 1))) {
      return true;
    }
    return false;
  }

  //--------------------------------------------------------
  //  construct_regulated_nodes:
  //    constructs the set of nodes being regulated
  //--------------------------------------------------------
  void KinetostatFixed::construct_regulated_nodes()
  {
    InterscaleManager & interscaleManager(atc_->interscale_manager());

    if (!atomicRegulator_->use_localized_lambda()) {
      regulatedNodes_ = new RegulatedNodes(atc_);
    }
    else if (kinetostat_->coupling_mode() == Kinetostat::FLUX) {
      regulatedNodes_ = new FixedNodes(atc_);
    }
    else if (kinetostat_->coupling_mode() == Kinetostat::FIXED) {
      // include boundary nodes
      regulatedNodes_ = new FixedBoundaryNodes(atc_);
    }
    else {
      throw ATC_Error("ThermostatFixed::construct_regulated_nodes - couldn't determine set of regulated nodes");
    }
     
    interscaleManager.add_set_int(regulatedNodes_,
                                  regulatorPrefix_+"RegulatedNodes");

    applicationNodes_ = regulatedNodes_;

    // special set of boundary elements for defining regulated atoms 
    if (atomicRegulator_->use_localized_lambda()) {
      elementMask_ = new ElementMaskNodeSet(atc_,applicationNodes_);
      interscaleManager.add_dense_matrix_bool(elementMask_,
                                                   regulatorPrefix_+"BoundaryElementMask");
    }
  }

  //--------------------------------------------------------
  //  initialize_delta_nodal_atomic_momentum:
  //    initializes storage for the variable tracking
  //    the change in the nodal atomic momentum
  //    that has occured over the past timestep
  //--------------------------------------------------------
  void KinetostatFixed::initialize_delta_nodal_atomic_momentum(double dt)
  {
    // initialize delta energy
    const DENS_MAT & myNodalAtomicMomentum(nodalAtomicMomentum_->quantity());
    initialNodalAtomicMomentum_ = myNodalAtomicMomentum;
    initialNodalAtomicMomentum_ *= -1.; // initially stored as negative for efficiency
    timeFilter_->apply_pre_step1(nodalAtomicMomentumFiltered_.set_quantity(),
                                 myNodalAtomicMomentum,dt);
  }

  //--------------------------------------------------------
  //  compute_delta_nodal_atomic_momentum:
  //    computes the change in the nodal atomic momentum
  //    that has occured over the past timestep
  //--------------------------------------------------------
  void KinetostatFixed::compute_delta_nodal_atomic_momentum(double dt)
  {
    // set delta energy based on predicted atomic velocities
    const DENS_MAT & myNodalAtomicMomentum(nodalAtomicMomentum_->quantity());
    timeFilter_->apply_post_step1(nodalAtomicMomentumFiltered_.set_quantity(),
                                  myNodalAtomicMomentum,dt);
    deltaNodalAtomicMomentum_ = initialNodalAtomicMomentum_;
    deltaNodalAtomicMomentum_ += myNodalAtomicMomentum;
  }

  //--------------------------------------------------------
  //  apply_pre_predictor:
  //    apply the kinetostat to the atoms in the
  //    pre-predictor integration phase
  //--------------------------------------------------------
  void KinetostatFixed::apply_pre_predictor(double dt)
  {
    // initialize values to be track change in finite element energy over the timestep
    initialize_delta_nodal_atomic_momentum(dt);
    initialFeMomentum_ = -1.*((mdMassMatrix_.quantity())*(velocity_.quantity())); // initially stored as negative for efficiency

    KinetostatGlcFs::apply_pre_predictor(dt);
  }
 
  //--------------------------------------------------------
  //  apply_post_corrector:
  //    apply the kinetostat to the atoms in the
  //    post-corrector integration phase
  //--------------------------------------------------------
  void KinetostatFixed::apply_post_corrector(double dt)
  {  
    KinetostatGlcFs::apply_post_corrector(dt);

    // update filtered momentum with lambda force
    DENS_MAT & myNodalAtomicLambdaForce(nodalAtomicLambdaForce_->set_quantity());
    timeFilter_->apply_post_step2(nodalAtomicMomentumFiltered_.set_quantity(),
                                  myNodalAtomicLambdaForce,dt);

    if (halve_force()) {
      // Halve lambda force due to fixed temperature constraints
      // 1) makes up for poor initial condition
      // 2) accounts for possibly large value of lambda when atomic shape function values change
      //    from eulerian mapping after more than 1 timestep
      //    avoids unstable oscillations arising from 
      //    thermostat having to correct for error introduced in lambda changing the 
      //    shape function matrices
      *lambda_ *= 0.5;
    }
    
    isFirstTimestep_ = false;
  }

  //--------------------------------------------------------
  //  compute_kinetostat
  //            manages the solution and application of the
  //            kinetostat equations and variables
  //--------------------------------------------------------
  void KinetostatFixed::compute_lambda(double dt)
  {
    // compute predicted changes in nodal atomic momentum
    compute_delta_nodal_atomic_momentum(dt);

    // change in finite element momentum
    deltaFeMomentum_ = initialFeMomentum_;
    deltaFeMomentum_ += (mdMassMatrix_.quantity())*(velocity_.quantity());

    // set up rhs for lambda equation
    KinetostatGlcFs::compute_lambda(dt);
  }

  //--------------------------------------------------------
  //  set_kinetostat_rhs
  //            sets up the RHS of the kinetostat equations
  //            for the coupling parameter lambda
  //--------------------------------------------------------
  void KinetostatFixed::set_kinetostat_rhs(DENS_MAT & rhs, double dt)
  {
    // for essential bcs (fixed nodes) :
    // form rhs : (delUpsV - delUps)/dt
    const set<int> & regulatedNodes(regulatedNodes_->quantity());
    double factor = (1./dt);
    for (int i = 0; i < nNodes_; i++) {
      if (regulatedNodes.find(i) != regulatedNodes.end()) {
        for (int j = 0; j < nsd_; j++) {
          rhs(i,j) = factor*(deltaNodalAtomicMomentum_(i,j) - deltaFeMomentum_(i,j));
        }
      }
      else {
        for (int j = 0; j < nsd_; j++) {
          rhs(i,j) = 0.;
        }
      }
    }
  }

  //--------------------------------------------------------
  //  add_to_momentum:
  //            determines what if any contributions to the
  //            finite element equations are needed for
  //            consistency with the kinetostat
  //--------------------------------------------------------
  void KinetostatFixed::add_to_momentum(const DENS_MAT & nodalLambdaForce,
                                       DENS_MAT & deltaMomentum,
                                       double dt)
  {
    deltaMomentum.resize(nNodes_,nsd_);
    const set<int> & regulatedNodes(regulatedNodes_->quantity());
    for (int i = 0; i < nNodes_; i++) {
      if (regulatedNodes.find(i) != regulatedNodes.end()) {
        for (int j = 0; j < nsd_; j++) {
          deltaMomentum(i,j) = 0.;
        }
      }
      else {
        for (int j = 0; j < nsd_; j++) {
          deltaMomentum(i,j) = nodalLambdaForce(i,j);
        }
      }
    }
  }
};
