#ifndef GHOST_MANAGER_H
#define GHOST_MANAGER_H

// ATC headers
#include "MatrixLibrary.h"
#include "PerAtomQuantityLibrary.h"
#include "TimeIntegrator.h"
#include "ATC_TypeDefs.h"

namespace ATC {

  // forward declarations
  class ATC_Method;
  class GhostModifier;
  class LammpsInterface;

  /**
   *  @class  GhostManager
   *  @brief  Manages methods for modifying ghost atoms
   */

  class GhostManager {

  public:

    /** types of ghost boundary conditions in momentum */
    enum BoundaryDynamicsType {
      NO_BOUNDARY_DYNAMICS=0,
      VERLET,  // performs velocity-verlet 
      PRESCRIBED,  // forces ghost locations to conform to interpolated finite element locations
      DAMPED_HARMONIC, // turns ghost atoms into spring-mass-dashpot systems
      DAMPED_LAYERS, // oer layer DAMPED_HARMONIC
      COUPLED,  // applies a spring-dashpot force to the ghosts
      SWAP, // exchanges ghost and real atoms when they cross AtC boundaries
      SWAP_VERLET // like above, but integrates the ghosts using velocity verlet
    };

    // constructor
    GhostManager(ATC_Method * atc);

    // destructor
    virtual ~GhostManager();

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** create objects to implement requested numerical method */
    virtual void construct_methods();

    /** create and get necessary transfer operators */
    virtual void construct_transfers();

    /** pre time integration initialization of data */
    virtual void initialize();

    /** prior to lammps exchange */
    virtual void pre_exchange();

    /** Predictor phase, Verlet first step for velocity */
    virtual void init_integrate_velocity(double dt);

    /** Predictor phase, Verlet first step for position */
    virtual void init_integrate_position(double dt);

    /** set positions after integration */
    virtual void post_init_integrate();

    /** Corrector phase, Verlet second step for velocity */
    virtual void final_integrate(double dt);

    /** sets the boundary dynamics flag as desired */
    void set_boundary_dynamics(BoundaryDynamicsType boundaryDynamics) {boundaryDynamics_ = boundaryDynamics;}

    /** flag for reset */
    bool need_reset() const {return needReset_;};

    /** access to ATC method object */
    ATC_Method * atc() {return atc_;};

  protected:

    /** pointer to routines that modify ghosts */
    GhostModifier* ghostModifier_;

    /** pointer to access ATC methods */
    ATC_Method * atc_;

    /** boundary dynamics method type */
    BoundaryDynamicsType boundaryDynamics_;

    /** flag for reset */
    bool needReset_;

    /** spring constant for some models */
    std::vector<double> kappa_;

    /** damping constant for some models */
    std::vector<double> gamma_;

    /** ratio between mass of ghost types and desired mass for some models */
    std::vector<double> mu_;

  private:

    // DO NOT define this
    GhostManager();

  };

  /**
   *  @class  GhostModifier
   *  @brief  Base class for objects which modify the ghost atoms, integrates ghost atoms using velocity-verlet if requested
   */

  class GhostModifier {
  
  public:

    // constructor
    GhostModifier(GhostManager * ghostManager);

    // destructor
    virtual ~GhostModifier();

    /** create and get necessary transfer operators */
    virtual void construct_transfers();

    /** pre time integration initialization of data */
    virtual void initialize(){};

    /** Predictor phase, Verlet first step for velocity */
    virtual void init_integrate_velocity(double dt);

    /** Predictor phase, Verlet first step for position */
    virtual void init_integrate_position(double dt);

    /** set positions after integration */
    virtual void post_init_integrate(){};

    /** prior to lammps exchange */
    virtual void pre_exchange(){};

    /** Corrector phase, Verlet second step for velocity */
    virtual void final_integrate(double dt);

    /** sets the verlet integration flag as desired */
    void set_integrate_atoms(bool integrateAtoms) {integrateAtoms_ = integrateAtoms;}

  protected:

    /** owning ghost manager */
    GhostManager * ghostManager_;

    /** object which integrates atoms */
    AtomTimeIntegrator * atomTimeIntegrator_;

    /** flag to perform velocity-verlet integration of ghosts */
    bool integrateAtoms_;


  private:

    // DO NOT define this
    GhostModifier();

  };

  /**
   *  @class  GhostModifierPrescribed
   *  @brief  sets ghost atom positions based on FE displacement
   */

  class GhostModifierPrescribed : public GhostModifier {
  
  public:

    // constructor
    GhostModifierPrescribed(GhostManager * ghostManager);

    // destructor
    virtual ~GhostModifierPrescribed(){};

    /** create and get necessary transfer operators */
    virtual void construct_transfers();

    /** set positions after integration */
    virtual void post_init_integrate();

  protected:

    /** positions of atoms */
    PerAtomQuantity<double> * atomPositions_;

    /** FE displacement at ghost locations */
    PerAtomQuantity<double> * atomFeDisplacement_;

    /** atom reference positions */
    PerAtomQuantity<double> * atomRefPositions_;

  private:

    // DO NOT define this
    GhostModifierPrescribed();

  };

  /**
   *  @class  GhostModifierDampedHarmonic
   *  @brief  Integrates ghost atoms using velocity-verlet with a damped harmonic force
   */

  class GhostModifierDampedHarmonic : public GhostModifierPrescribed {
  
  public:

    // constructor
    GhostModifierDampedHarmonic(GhostManager * ghostManager,
                                const std::vector<double> & kappa,
                                const std::vector<double> & gamma,
                                const std::vector<double> & mu);

    // destructor
    virtual ~GhostModifierDampedHarmonic(){};

    /** create and get necessary transfer operators */
    virtual void construct_transfers();
#if true
    /** Predictor phase, Verlet first step for velocity */
    virtual void init_integrate_velocity(double dt);

    /** Predictor phase, Verlet first step for position */
    virtual void init_integrate_position(double dt);
#endif
    /** set positions after integration */
    virtual void post_init_integrate(){};

    /** Corrector phase, Verlet second step for velocity */
    virtual void final_integrate(double dt);

  protected:

    /** velocities of atoms */
    PerAtomQuantity<double> * atomVelocities_;

    /** FE velocity at ghost locations */
    PerAtomQuantity<double> * atomFeVelocity_;

    /** atom forces */
    PerAtomQuantity<double> * atomForces_;

    /** effective spring constant for potential */
    double k0_;

    /** spring constant */
    const std::vector<double> & kappa_;

    /** damping constant */
    const std::vector<double> & gamma_;

    /** ratio between mass of ghost types and desired mass */
    const std::vector<double> & mu_;

    // workspace
    DENS_MAT _forces_;

  private:

    // DO NOT define this
    GhostModifierDampedHarmonic();

  };

  /**
   *  @class  GhostModifierDampedHarmonicLayers
   *  @brief  Integrates ghost atoms using velocity-verlet with a damped harmonic force based on which layer the atom resides in
   */

  class GhostModifierDampedHarmonicLayers : public GhostModifierDampedHarmonic {
  
  public:

    // constructor
    GhostModifierDampedHarmonicLayers(GhostManager * ghostManager,
                                      const std::vector<double> & kappa,
                                      const std::vector<double> & gamma,
                                      const std::vector<double> & mu);

    // destructor
    virtual ~GhostModifierDampedHarmonicLayers(){};

    /** create and get necessary transfer operators */
    virtual void construct_transfers();

    /** pre time integration initialization of data */
    virtual void initialize();

    /** Corrector phase, Verlet second step for velocity */
    virtual void final_integrate(double dt);

  protected:

    // methods
    /** compute distance of ghost atom to boundary */
    void compute_distances();
    /** sorting heuristics to identify layers */
    int find_layers();

    // data
    /** distance from all ghost atoms to boundary, i.e. boundary face of containing element */
    PerAtomQuantity<double> * ghostToBoundaryDistance_;
    
    /** layer id for ghost atoms */
    PerAtomQuantity<int> * layerId_;

  private:

    // DO NOT define this
    GhostModifierDampedHarmonicLayers();

  };


  /**
   *  @class  GhostIntegratorSwap
   *  @brief  Integrates ghost atoms using velocity-verlet, and swaps atoms between ghost
   *          and internal depending on what element they are in
   */

  class GhostIntegratorSwap : public GhostModifier {
  
  public:

    // constructor
    GhostIntegratorSwap(GhostManager * ghostManager);

    // destructor
    virtual ~GhostIntegratorSwap(){};

    /** create and get necessary transfer operators */
    virtual void construct_transfers();

    /** pre time integration initialization of data */
    virtual void initialize();

    /** prior to lammps exchange */
    virtual void pre_exchange();

  protected:

    /** pointer to lammps interface */
    LammpsInterface * lammpsInterface_;

    /** internal element set */
    const std::set<int> & elementSet_;

    /** internal to element map */
    PerAtomQuantity<int> * atomElement_;
    
    /** ghost to element map */
    PerAtomQuantity<int> * atomGhostElement_;

    /** internal to atom map */
    const Array<int> & internalToAtom_;

    /** ghost to atom map */
    const Array<int> & ghostToAtom_;

    /** group bit for internal */
    int groupbit_;

    /** group bit for ghost */
    int groupbitGhost_;

  private:

    // DO NOT define this
    GhostIntegratorSwap();

  };


};
#endif
