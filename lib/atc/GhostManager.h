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
      VERLET,
      PRESCRIBED,
      DAMPED_HARMONIC,
      COUPLED,
      SWAP,
      SWAP_VERLET
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
    double kappa_;

    /** damping constant for some models */
    double gamma_;

    /** ratio between mass of ghost types and desired mass for some models */
    double mu_;

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
                                double kappa_, double gamma, double mu);

    // destructor
    virtual ~GhostModifierDampedHarmonic(){};

    /** create and get necessary transfer operators */
    virtual void construct_transfers();

    /** Predictor phase, Verlet first step for velocity */
    virtual void init_integrate_velocity(double dt);

    /** Predictor phase, Verlet first step for position */
    virtual void init_integrate_position(double dt);

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

    /** spring constant */
    double kappa_, k0_;

    /** damping constant */
    double gamma_;

    /** ratio between mass of ghost types and desired mass */
    double mu_;

    // workspace
    DENS_MAT _forces_;

  private:

    // DO NOT define this
    GhostModifierDampedHarmonic();

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
