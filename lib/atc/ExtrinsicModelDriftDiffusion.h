#ifndef EXTRINSIC_MODEL_DRIFT_DIFFUSION
#define EXTRINSIC_MODEL_DRIFT_DIFFUSION

#include <set>
#include <string>
#include <vector>

#include "ExtrinsicModelTwoTemperature.h"
#include "SchrodingerSolver.h"

namespace ATC {

  class ATC_Coupling;
  class PrescribedDataManager;
  class ExtrinsicModel;
  class PhysicsModel;
  class PoissonSolver;
  class LinearSolver;
  class SchrodingerSolver;
  class SchrodingerPoissonSolver;

  /**
   *  @class  ExtrinsicModelDriftDiffusion
   *  @brief  add electron temperature physics to phonon physics
   *          owned fields ELECTRON_DENSITY
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModelDriftDiffusion
  //--------------------------------------------------------
  //--------------------------------------------------------

  class ExtrinsicModelDriftDiffusion : public ExtrinsicModelTwoTemperature {
  
  public:

    // constructor
    ExtrinsicModelDriftDiffusion(ExtrinsicModelManager * modelManager,
                   ExtrinsicModelType modelType,
                   std::string matFileName);

    // destructor
    virtual ~ExtrinsicModelDriftDiffusion();

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** pre time integration */
    virtual void initialize();

    /** Predictor phase, executed before Verlet */
    virtual void pre_init_integrate();

    /** Add model-specific output data */
    virtual void output(OUTPUT_LIST & outputData);

    /** set up LAMMPS display variables */
    virtual int size_vector(int externalSize);

    /** get LAMMPS display variables */
    virtual bool compute_vector(int n, double & value);
 
  protected:
    /** Poisson solve */
    void poisson_solve();

    /** Schrodinger-Poisson solve */
    void schrodinger_poisson_solve(void); // wrapper
    void schrodinger_poisson_solve(int iterations);
    void slice_schrodinger_poisson_solve(int consistencyIter, int constraintIter);
    double update_fermi_energy(double target,bool first = false);

    /** time integrator for the continuity eqn */
    FieldEulerIntegrator * continuityIntegrator_;

    /** poisson solver type */
    SolverType poissonSolverType_;

    /** poisson solver */
    PoissonSolver * poissonSolver_;

    /** offset/size for LAMMPS display output */
    int baseSize_;

    /** ways to determine the electron density */
    int  electronDensityEqn_;
    enum electronDensityEqnType { ELECTRON_CONTINUITY,
                                  ELECTRON_EQUILIBRIUM,
                                  ELECTRON_SCHRODINGER};
  
    /** frequency for updating the electron state */
    int fluxUpdateFreq_;

    /** Schrodinger solver type */
    SolverType schrodingerSolverType_;

    /** poisson solver */
    SchrodingerSolver * schrodingerSolver_;

    /** schrodinger-poisson solver */
    SchrodingerPoissonManager  schrodingerPoissonMgr_;
    SchrodingerPoissonSolver * schrodingerPoissonSolver_;

    /** schrodinger-poisson data */
    int maxConsistencyIter_, maxConstraintIter_;
    double safe_dEf_, Ef_shift_;
    DENS_MAT phiTotal_;
    double Tmax_;

    Array2D<double> EfHistory_;

    /** one dimensional restriction */
    bool oneD_;
    int oneDcoor_;
    int oneDstride_;
    std::string oneDnodesetName_;
    std::set<int> oneDnodeset_;
    Array< std::set<int> > oneDslices_;
    int oneDconserve_;

    DENS_MAT JE_;

  };

  /**
   *  @class  ExtrinsicModelDriftDiffusionConvection
   *  @brief  add electron temperature physics to phonon physics, including convective transport
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModelDriftDiffusionConvection
  //--------------------------------------------------------
  //--------------------------------------------------------

  class ExtrinsicModelDriftDiffusionConvection : public ExtrinsicModelDriftDiffusion {
  
  public:

    // constructor
    ExtrinsicModelDriftDiffusionConvection(ExtrinsicModelManager * modelManager,
                   ExtrinsicModelType modelType,
                   std::string matFileName);

    // destructor
    virtual ~ExtrinsicModelDriftDiffusionConvection();

    /** pre time integration */
    virtual void initialize();

    /** Predictor phase, executed before Verlet */
    virtual void pre_init_integrate();

    /** Set sources to AtC equation */
    //virtual void set_sources(FIELDS & fields, FIELDS & sources);

    /** Add model-specific output data */
    virtual void output(OUTPUT_LIST & outputData);

    /** set up LAMMPS display variables */
    virtual int size_vector(int externalSize);

    /** get LAMMPS display variables */
    virtual bool compute_vector(int n, double & value);
 
  protected:

    /** compute the total kinetic energy of the electrons */
    void compute_nodal_kinetic_energy(DENS_MAT & kineticEnergy);

    /** Linear solver for velocity */
    std::vector<LinearSolver * > velocitySolvers_;

    
    /** Linear solver for solving the poisson equations */
    LinearSolver * cddmPoissonSolver_;

    /** offset/size for LAMMPS display output */
    int baseSize_;

    double timerStart_, timerCurrent_;

  };

};
#endif
