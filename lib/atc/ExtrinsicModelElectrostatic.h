#ifndef EXTRINSIC_MODEL_ELECTROSTATIC
#define EXTRINSIC_MODEL_ELECTROSTATIC

// WIP_REJ
#define CHARGED_SURFACE

#include "ExtrinsicModel.h"
#include <utility>
#include <string>
#include <vector>
#include <map>

namespace ATC {

  // forward declarations
  class ATC_Coupling;
  class PrescribedDataManager;
  class ExtrinsicModel;
  class PhysicsModel;
  class PoissonSolver;
  class FundamentalAtomQuantity;
  class AtfShapeFunctionRestriction;
  class ChargeRegulator;

  /**
   *  @class  ExtrinsicModelElectrostatic
   *  @brief  add self-consistent electrostatic forces
   *          owned field: ELECTRIC_POTENTIAL 
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModelElectrostatic
  //--------------------------------------------------------
  //--------------------------------------------------------

  class ExtrinsicModelElectrostatic : public ExtrinsicModel {
  
  public:

    // constructor
    ExtrinsicModelElectrostatic(ExtrinsicModelManager * modelManager,
                   ExtrinsicModelType modelType,
                   std::string matFileName);

    // destructor
    virtual ~ExtrinsicModelElectrostatic();

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** construct transfers needed for model */
    virtual void construct_transfers();

    /** pre time integration */
    virtual void initialize();

    /** Predictor phase, executed before Verlet */
    virtual void post_init_integrate();

    /** changes lammps forces to include long-range electrostatic interactions */
    virtual void post_force();

    /** Add model-specific output data */
    virtual void output(OUTPUT_LIST & outputData);

    /** set up LAMMPS display variables */
    virtual int size_vector(int externalSize);

    /** get LAMMPS display variables */
    virtual double compute_scalar(void);
    virtual bool compute_vector(int n, double & value);

    PoissonSolver * poisson_solver(void) const { return poissonSolver_;} 

 
  protected:
    /** poisson solver type */
    SolverType poissonSolverType_;
    double poissonSolverTol_;
    int poissonSolverMaxIter_;

    /** poisson solver */
    PoissonSolver * poissonSolver_;

    /** max solves per minimize */
    int maxSolves_;

    /** offset/size for LAMMPS display output */
    int baseSize_;

    /** rhs mask for Poisson solver */
    Array2D<bool> rhsMask_;
  
    /** estimate intrinsic charge density */
    void add_electrostatic_forces(MATRIX & nodalPotential);

    /** correct short range FE electric field */
    void correct_electrostatic_forces();

#ifdef CHARGED_SURFACE
    /** account for charged surfaces on charged atoms */
    void apply_charged_surfaces(MATRIX & nodalPotential);

    /** set charged surface data */
    void add_charged_surface(const std::string & facesetName,
                             const double chargeDensity);
#endif
    /** charge regulator */
    ChargeRegulator * chargeRegulator_; 

    /** local electric potential Green's function for each node */
    std::vector<SparseVector<double> > greensFunctions_;

#ifdef CHARGED_SURFACE
    /** stores surface charge data at fixed charge surfaces */
    std::map<std::string,double> surfaceCharges_;

    /** data structure to store information for applying charged surfaces */
    std::map<std::string,std::map<int,std::pair<DENS_VEC,double> > > chargedSurfaces_;
#endif

    /** data structure storing potential induced only by charges under the nodal shape function support */
    std::map<std::string,std::map<int,double> > nodalChargePotential_;

    
    /** allows electric force only applied only in z direction to complement LAMMPS slab command */
    bool useSlab_;

    
    /** enables method when short range interactions are off */
    bool includeShortRange_;

    /** atomic forces */
    FundamentalAtomQuantity * atomForces_;

    /** coarse-grained charges from internal atoms */
    DENS_MAN * nodalAtomicCharge_;

    /** coarse-grained charges from ghost atoms */
    DENS_MAN * nodalAtomicGhostCharge_;

    /** workspace */
    DENS_MAT _atomElectricalForce_;
    double  totalElectricalForce_[3];
  };


  /**
   *  @class  ExtrinsicModelElectrostaticMomentum
   *  @brief  add self-consistent electrostatic forces with elastic coupling
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModelElectrostaticMomentum
  //--------------------------------------------------------
  //--------------------------------------------------------

  class ExtrinsicModelElectrostaticMomentum : public ExtrinsicModelElectrostatic {
  
  public:

    // constructor
    ExtrinsicModelElectrostaticMomentum(ExtrinsicModelManager * modelManager,
                                       ExtrinsicModelType modelType,
                                       std::string matFileName);

    // destructor
    virtual ~ExtrinsicModelElectrostaticMomentum();

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** pre time integration */
    virtual void initialize();

    /** Set sources to AtC equation */
    virtual void set_sources(FIELDS & fields, FIELDS & sources);

    /** Add model-specific output data */
    virtual void output(OUTPUT_LIST & outputData);
 
  };

};
#endif
