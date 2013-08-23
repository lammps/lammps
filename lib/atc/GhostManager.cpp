// ATC transfer headers
#include "GhostManager.h"
#include "ATC_Method.h"
#include "LammpsInterface.h"
#include "ATC_Error.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class GhostManager
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  GhostManager::GhostManager(ATC_Method * atc) :
    ghostModifier_(NULL),
    atc_(atc),
    boundaryDynamics_(NO_BOUNDARY_DYNAMICS),
    needReset_(true),
    kappa_(0.),
    gamma_(0.),
    mu_(0.)
  {
    // do nothing
  }
  
  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  GhostManager::~GhostManager()
  {
    if (ghostModifier_) delete ghostModifier_;
  }

  //--------------------------------------------------------
  //  modify
  //    parses inputs and modifies state of the integrator
  //--------------------------------------------------------
  bool GhostManager::modify(int narg, char **arg)
  {
    int argIndex = 0;
    
    /*! \page man_boundary_dynamics fix_modify AtC boundary_dynamics
      \section syntax
      fix_modify AtC boundary_dynamics < on | damped_harmonic | prescribed | coupled | none > [args] \n
      \section description
      Sets different schemes for controlling boundary atoms.  On will integrate the boundary atoms using the velocity-verlet algorithm.  Damped harmonic uses a mass/spring/dashpot for the boundary atoms with added arguments of the damping and spring constants followed by the ratio of the boundary type mass to the desired mass.  Prescribed forces the boundary atoms to follow the finite element displacement.  Coupled does the same.
      \section restrictions
      Boundary atoms must be specified.  When using swaps between internal and boundary atoms, the initial configuration must have already correctly partitioned the two.
      \section related
      \man_boundary
      \section default
      prescribed
      on
    */
    if (strcmp(arg[argIndex],"none")==0) {
      boundaryDynamics_ = NO_BOUNDARY_DYNAMICS;
      needReset_ = true;
      return true;
    }
    else if (strcmp(arg[argIndex],"on")==0) {
      boundaryDynamics_ = VERLET;
      needReset_ = true;
      return true;
    }
    else if (strcmp(arg[argIndex],"prescribed")==0) {
      boundaryDynamics_ = PRESCRIBED;
      needReset_ = true;
      return true;
    }
    else if (strcmp(arg[argIndex],"damped_harmonic")==0) {
      argIndex++;
      kappa_ = atof(arg[argIndex++]);
      gamma_ = atof(arg[argIndex++]);
      mu_    = atof(arg[argIndex]);
      boundaryDynamics_ = DAMPED_HARMONIC;
      needReset_ = true;
      return true;
    }
    else if (strcmp(arg[argIndex],"coupled")==0) {
      boundaryDynamics_ = COUPLED;
      kappa_ = 0.;
      gamma_ = 0.;
      mu_    = 0.;
      needReset_ = true;
      return true;
    }

    return false;
  }

  //--------------------------------------------------------
  //  construct_methods
  //    constructs the specific method to modify the ghosts
  //--------------------------------------------------------
  void GhostManager::construct_methods()
  {
    if (ghostModifier_) {
      delete ghostModifier_;
      ghostModifier_ = NULL;
    }

    if (!atc_->groupbit_ghost()) {
      ghostModifier_ = new GhostModifier(this);
      return;
    }
    
    switch (boundaryDynamics_) {
    case VERLET: {
      ghostModifier_ = new GhostModifier(this);
      ghostModifier_->set_integrate_atoms(true);
      break;
    }
    case PRESCRIBED: {
      ghostModifier_ = new GhostModifierPrescribed(this);
      break;
    }
    case COUPLED:
    case DAMPED_HARMONIC: {
      ghostModifier_ = new GhostModifierDampedHarmonic(this,kappa_,gamma_,mu_);
      ghostModifier_->set_integrate_atoms(true);
      break;
    }
    case SWAP: {
      // if regions based on element sets, use verlet on ghosts and swap ghosts and internal when they move between regions
      const std::string & internalElementSet(atc_->internal_element_set());
      if (internalElementSet.size() && (atc_->atom_to_element_map_type()==EULERIAN)) {
        LammpsInterface * lammpsInterface = LammpsInterface::instance();
        if (atc_->atom_to_element_map_frequency() % lammpsInterface->reneighbor_frequency() != 0) {
          throw ATC_Error("GhostManager::construct_methods - eulerian frequency and lammsp reneighbor frequency must be consistent to swap boundary and internal atoms");
        }
        ghostModifier_ = new GhostIntegratorSwap(this);
      }
      break;
    }
    case SWAP_VERLET: {
      // if regions based on element sets, use verlet on ghosts and swap ghosts and internal when they move between regions
      const std::string & internalElementSet(atc_->internal_element_set());
      if (internalElementSet.size() && (atc_->atom_to_element_map_type()==EULERIAN)) {
        LammpsInterface * lammpsInterface = LammpsInterface::instance();
        if (atc_->atom_to_element_map_frequency() % lammpsInterface->reneighbor_frequency() != 0) {
          throw ATC_Error("GhostManager::construct_methods - eulerian frequency and lammsp reneighbor frequency must be consistent to swap boundary and internal atoms");
        }
        ghostModifier_ = new GhostIntegratorSwap(this);
        ghostModifier_->set_integrate_atoms(true);
      }
      break;
    }
    default: {
      ghostModifier_ = new GhostModifier(this);
    }
    }
  }
  
  //--------------------------------------------------------
  //  construct_transfers
  //    sets/constructs all required dependency managed data
  //--------------------------------------------------------
  void GhostManager::construct_transfers()
  {
    ghostModifier_->construct_transfers();
  }
  
  //--------------------------------------------------------
  //  initialize
  //    initialize all data and variables before a run
  //--------------------------------------------------------
  void GhostManager::initialize()
  {
    ghostModifier_->initialize();
    needReset_ = false;
  }

  //--------------------------------------------------------
  //  pre_exchange
  //    makes any updates required before lammps exchanges
  //    atoms
  //--------------------------------------------------------
  void GhostManager::pre_exchange()
  {
    ghostModifier_->pre_exchange();
  }
  
  //--------------------------------------------------------
  //  initial_integrate_velocity
  //    velocity update in first part of velocity-verlet
  //--------------------------------------------------------
  void GhostManager::init_integrate_velocity(double dt)
  {
    ghostModifier_->init_integrate_velocity(dt);
  }
      
  //--------------------------------------------------------
  //  initial_integrate_position
  //    position update in first part of velocity-verlet
  //--------------------------------------------------------
  void GhostManager::init_integrate_position(double dt)
  {
    ghostModifier_->init_integrate_position(dt);
  }

  //--------------------------------------------------------
  //  post_init_integrate
  //    makes any updates required after first integration
  //--------------------------------------------------------
  void GhostManager::post_init_integrate()
  {
    ghostModifier_->post_init_integrate();
  }

  //--------------------------------------------------------
  //  final_integrate
  //    velocity update in second part of velocity-verlet
  //--------------------------------------------------------
  void GhostManager::final_integrate(double dt)
  {
    ghostModifier_->final_integrate(dt);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class GhostModifier
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  GhostModifier::GhostModifier(GhostManager * ghostManager) :
    ghostManager_(ghostManager),
    atomTimeIntegrator_(NULL),
    integrateAtoms_(false)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  GhostModifier::~GhostModifier()
  {
    if (atomTimeIntegrator_) delete atomTimeIntegrator_;
  }

  //--------------------------------------------------------
  //  construct_transfers
  //    sets/constructs all required dependency managed data
  //--------------------------------------------------------
  void GhostModifier::construct_transfers()
  {
    if (atomTimeIntegrator_) delete atomTimeIntegrator_;
    if (integrateAtoms_) {
      atomTimeIntegrator_ = new AtomTimeIntegratorType(ghostManager_->atc(),GHOST);
      atomTimeIntegrator_->construct_transfers();
    }
    else {
      atomTimeIntegrator_ = new AtomTimeIntegrator();
    }

  }

  //--------------------------------------------------------
  //  initial_integrate_velocity
  //    velocity update in first part of velocity-verlet
  //--------------------------------------------------------
  void GhostModifier::init_integrate_velocity(double dt)
  {
    atomTimeIntegrator_->init_integrate_velocity(dt);
  }
      
  //--------------------------------------------------------
  //  initial_integrate_position
  //    position update in first part of velocity-verlet
  //--------------------------------------------------------
  void GhostModifier::init_integrate_position(double dt)
  {
    atomTimeIntegrator_->init_integrate_position(dt);
  }

  //--------------------------------------------------------
  //  final_integrate
  //    velocity update in second part of velocity-verlet
  //--------------------------------------------------------
  void GhostModifier::final_integrate(double dt)
  {
    atomTimeIntegrator_->final_integrate(dt);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class GhostModifierPrescribed
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  GhostModifierPrescribed::GhostModifierPrescribed(GhostManager * ghostManager) :
    GhostModifier(ghostManager),
    atomPositions_(NULL),
    atomFeDisplacement_(NULL),
    atomRefPositions_(NULL)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  construct_transfers
  //    sets/constructs all required dependency managed data
  //--------------------------------------------------------
  void GhostModifierPrescribed::construct_transfers()
  {
    GhostModifier::construct_transfers();

    InterscaleManager & interscaleManager((ghostManager_->atc())->interscale_manager());
    atomPositions_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_POSITION,GHOST);

    // prolongation from displacement field to atoms
    PerAtomSparseMatrix<double> * atomShapeFunctions = interscaleManager.per_atom_sparse_matrix("InterpolantGhost");
    if (!atomShapeFunctions) {
      atomShapeFunctions = new PerAtomShapeFunction(ghostManager_->atc(),
                                                    interscaleManager.per_atom_quantity("AtomicGhostCoarseGrainingPositions"),
                                                    interscaleManager.per_atom_int_quantity("AtomGhostElement"),
                                                    GHOST);
      interscaleManager.add_per_atom_sparse_matrix(atomShapeFunctions,"InterpolantGhost");
    }
    atomFeDisplacement_ = new FtaShapeFunctionProlongation(ghostManager_->atc(),
                                                           &(ghostManager_->atc())->field(DISPLACEMENT),
                                                           atomShapeFunctions,
                                                           GHOST);
    interscaleManager.add_per_atom_quantity(atomFeDisplacement_,field_to_prolongation_name(DISPLACEMENT)+"Ghost");
    
    atomRefPositions_ = interscaleManager.per_atom_quantity("AtomicGhostCoarseGrainingPositions");
  }
      
  //--------------------------------------------------------
  //  post_init_integrate
  //    after integration, fix ghost atoms' positions
  //--------------------------------------------------------
  void GhostModifierPrescribed::post_init_integrate()
  {
    const DENS_MAT & atomFeDisplacement(atomFeDisplacement_->quantity());
    const DENS_MAT & atomRefPositions(atomRefPositions_->quantity());
    *atomPositions_ = atomFeDisplacement + atomRefPositions;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class GhostModifierDampedHarmonic
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  GhostModifierDampedHarmonic::GhostModifierDampedHarmonic(GhostManager * ghostManager,
                                                           double kappa,
                                                           double gamma, double mu) :
    GhostModifierPrescribed(ghostManager),
    atomVelocities_(NULL),
    atomFeVelocity_(NULL),
    atomForces_(NULL),
    kappa_(kappa),
    gamma_(gamma),
    mu_(mu)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  construct_transfers
  //    sets/constructs all required dependency managed data
  //--------------------------------------------------------
  void GhostModifierDampedHarmonic::construct_transfers()
  {
    GhostModifierPrescribed::construct_transfers();
    
    InterscaleManager & interscaleManager((ghostManager_->atc())->interscale_manager());
    atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY,GHOST);
    atomForces_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_FORCE,GHOST);

    // prolongation from displacement field to atoms
    PerAtomSparseMatrix<double> * atomShapeFunctions = interscaleManager.per_atom_sparse_matrix("InterpolantGhost");
    atomFeVelocity_ = new FtaShapeFunctionProlongation(ghostManager_->atc(),
                                                       &(ghostManager_->atc())->field(VELOCITY),
                                                       atomShapeFunctions,
                                                       GHOST);
    interscaleManager.add_per_atom_quantity(atomFeVelocity_,field_to_prolongation_name(VELOCITY)+"Ghost");
    // calculate nominal bond stiffness
    int i = 0, j = 1;// HACk should be an atom and its neighbor in the boundary
    double rsq = 0.0;
    //k0_ = LammpsInterface_->bond_stiffness(i,j,rsq);
    k0_ = LammpsInterface::instance()->bond_stiffness(i,j,rsq);
  }

  //--------------------------------------------------------
  //  initial_integrate_velocity
  //    velocity update in first part of velocity-verlet
  //--------------------------------------------------------
  void GhostModifierDampedHarmonic::init_integrate_velocity(double dt)
  {
    atomTimeIntegrator_->init_integrate_velocity(mu_*dt);
  }
      
  //--------------------------------------------------------
  //  initial_integrate_position
  //    position update in first part of velocity-verlet
  //--------------------------------------------------------
  void GhostModifierDampedHarmonic::init_integrate_position(double dt)
  {
    atomTimeIntegrator_->init_integrate_position(dt);
  }

  //--------------------------------------------------------
  //  final_integrate
  //    velocity update in second part of velocity-verlet
  //--------------------------------------------------------
  void GhostModifierDampedHarmonic::final_integrate(double dt)
  {

    const DENS_MAT & atomPositions(atomPositions_->quantity());
    const DENS_MAT & atomVelocities(atomVelocities_->quantity());
    const DENS_MAT & atomFeDisplacement(atomFeDisplacement_->quantity());
    const DENS_MAT & atomFeVelocity(atomFeVelocity_->quantity());
    const DENS_MAT & atomRefPositions(atomRefPositions_->quantity());

    _forces_ = atomFeDisplacement;
    _forces_ += atomRefPositions;
    _forces_ -= atomPositions;
    _forces_ *= kappa_;
    _forces_ += gamma_*(atomFeVelocity - atomVelocities);
    *atomForces_ = _forces_;

    atomTimeIntegrator_->final_integrate(mu_*dt);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class GhostIntegratorVerletSwap
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  GhostIntegratorSwap::GhostIntegratorSwap(GhostManager * ghostManager) :
    GhostModifier(ghostManager),
    lammpsInterface_(LammpsInterface::instance()),
    elementSet_((((ghostManager_->atc())->fe_engine())->fe_mesh())->elementset((ghostManager_->atc())->internal_element_set())),
    atomElement_(NULL),
    atomGhostElement_(NULL),
    internalToAtom_((ghostManager_->atc())->internal_to_atom_map()),
    ghostToAtom_((ghostManager_->atc())->ghost_to_atom_map()),
    groupbit_((ghostManager_->atc())->groupbit()),
    groupbitGhost_((ghostManager_->atc())->groupbit_ghost())
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  construct_transfers
  //    sets/constructs all required dependency managed data
  //--------------------------------------------------------
  void GhostIntegratorSwap::construct_transfers()
  {
    GhostModifier::construct_transfers();

    InterscaleManager & interscaleManager((ghostManager_->atc())->interscale_manager());
    atomElement_ = interscaleManager.per_atom_int_quantity("AtomElement");
    atomGhostElement_ = interscaleManager.per_atom_int_quantity("AtomGhostElement");
  }

  //--------------------------------------------------------
  //  initialize
  //    initialize all data and variables before a run
  //--------------------------------------------------------
  void GhostIntegratorSwap::initialize()
  {
    if ((ghostManager_->atc())->atom_to_element_map_frequency() % lammpsInterface_->reneighbor_frequency() != 0) {
      throw ATC_Error("GhostIntegratorSwap::initialize - AtC Eulerian reset frequency must be a multiple of the Lammps reneighbor frequency when using internal/boundary atom swapping");
    }
  }

  //--------------------------------------------------------
  //  pre_exchange
  //    swaps atoms between types depending on region
  //--------------------------------------------------------
  void GhostIntegratorSwap::pre_exchange()
  {
    // lammps mask to change type
    int * mask = lammpsInterface_->atom_mask();

    const DenseMatrix<int> & atomElement(atomElement_->quantity());
    for (int i = 0; i < atomElement.nRows(); ++i) {
      if (elementSet_.find(atomElement(i,0)) == elementSet_.end()) {
        mask[internalToAtom_(i)] |= groupbitGhost_;
        // remove from internal
      }
    }
  }
};
