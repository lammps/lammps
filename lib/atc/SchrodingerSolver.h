#ifndef SCHRODINGER_SOLVER_H
#define SCHRODINGER_SOLVER_H

// ATC includes
#include "Array2D.h"
#include "LinearSolver.h"

// other includes
#include <vector>
#include <map>
#include <set>

namespace ATC {

// Forward class declarations
class ATC_Coupling;
class FE_Engine;
class PrescribedDataManager;
class PhysicsModel;
class PoissonSolver;

/**
 *  @class SchrodingerSolver
 *  @brief a class to solve the (time-independent) Schrodinger equation 
 */
class SchrodingerSolver {

 public:

  /** Constructor */
  SchrodingerSolver(
    const FieldName fieldName,
    const PhysicsModel * physicsModel,
    const FE_Engine * feEngine,
    const PrescribedDataManager * prescribedDataMgr,
    /*const*/ ATC_Coupling * atc,
    bool parallel = true
  );

  /** Destructor */
  virtual ~SchrodingerSolver(){};

  /** parser */
  bool modify(int narg, char **arg){ return false;}

  /** initialize */
  void initialize(void);

  /** solve */
  virtual bool solve(FIELDS & fields); 


  

 protected:
  
  /** Pointer to ATC */
  ATC_Coupling * atc_;

  /** Pointer to FE_Engine */
  const FE_Engine * feEngine_;

  /** Pointer to PrescribedDataManager */
  const PrescribedDataManager * prescribedDataMgr_;

  /** Pointer to FE_Engine */
  const PhysicsModel * physicsModel_;

  /** field to solve for */
  FieldName fieldName_;

  /** number of nodes */
  int nNodes_;

  /** mass matrix */
  DENS_MAT M_;

  bool parallel_;

};

/**
 *  @class SchrodingerSolver
 *  @brief a class to solve the Schrodinger equation on slices
 */
class SliceSchrodingerSolver : public SchrodingerSolver {

 public:

  /** Constructor */
  SliceSchrodingerSolver(
    const FieldName fieldName,
    const PhysicsModel * physicsModel,
    const FE_Engine * feEngine,
    const PrescribedDataManager * prescribedDataMgr,
    /*const*/ ATC_Coupling * atc,
    const Array< std::set<int> > & oneDslices,
    const Array< double > & oneDdxs,
    bool parallel = true
  );

  /** Destructor */
  virtual ~SliceSchrodingerSolver(){};

  /** parser */
  bool modify(int narg, char **arg){return false;}

  /** initialize */
  void initialize(void);

  /** solve */
  virtual bool solve(FIELDS & fields); 

  Array< std::set<int> > & slices(void){ return oneDslices_;} 
  Array< double > & dxs(void){ return oneDdxs_;} 

 protected:

  Array< std::set<int> > oneDslices_;
  Array< double > oneDdxs_;

};

/**
 *  @class SchrodingerSolver
 *  @brief a class to solve the Schrodinger-Poisson system
 */
class SchrodingerPoissonSolver  {
  public:
    SchrodingerPoissonSolver(
      /*const*/ ATC_Coupling * atc,
      SchrodingerSolver * schrodingerSolver,
      PoissonSolver * poissonSolver,
      const PhysicsModel * physicsModel,
      int maxConsistencyIter
    );
    virtual ~SchrodingerPoissonSolver(void){};
    virtual void solve(
      FIELDS & rhs,
      GRAD_FIELD_MATS & fluxes
    );
  protected:
    ATC_Coupling * atc_;
    SchrodingerSolver * schrodingerSolver_;
    PoissonSolver * poissonSolver_;
    const PhysicsModel * physicsModel_;
    int maxConsistencyIter_;
    int maxConstraintIter_;
    int nNodes_;
    double Ef0_;
    double Ef_shift_;
    double safe_dEf_;
    double tol_;
};

/**
 *  @class SchrodingerSolver
 *  @brief a class to solve the Schrodinger-Poisson system on slices
 */
class SliceSchrodingerPoissonSolver : public SchrodingerPoissonSolver  {
  public:
    SliceSchrodingerPoissonSolver(
      /*const*/ ATC_Coupling * atc,
      SchrodingerSolver * schrodingerSolver,
      PoissonSolver * poissonSolver,
      const PhysicsModel * physicsModel,
      int maxConsistencyIter,
      int maxConstraintIter,
      int oneDconserve,
      double Ef_shift,
      double safe_dEf
    );
    virtual ~SliceSchrodingerPoissonSolver(void){};
    virtual void solve(
      FIELDS & rhs,
      GRAD_FIELD_MATS & fluxes
    );
  protected:
    int nslices_;
    double update_fermi_energy(double target, bool first, 
      GRAD_FIELD_MATS & fluxes);
    int oneDconserve_;
    int oneDcoor_;
    Array< std::set<int> > & oneDslices_;
    Array< double > & oneDdxs_;
    Array2D<double> EfHistory_;
};

/**
 *  @class SchrodingerSolver
 *  @brief a class to solve the Schrodinger-Poisson system on slices
 */
class GlobalSliceSchrodingerPoissonSolver : public SliceSchrodingerPoissonSolver  {
  public:
    GlobalSliceSchrodingerPoissonSolver(
      /*const*/ ATC_Coupling * atc,
      SchrodingerSolver * schrodingerSolver,
      PoissonSolver * poissonSolver,
      const PhysicsModel * physicsModel,
      int maxConsistencyIter,
      int maxConstraintIter,
      int oneDconserve,
      double Ef0,
      double alpha,
      double safe_dEf,
      double tol,
      double mu, double D
    );
    virtual ~GlobalSliceSchrodingerPoissonSolver(void);
    virtual void solve(
      FIELDS & rhs,
      GRAD_FIELD_MATS & fluxes
    );
  protected:
    void compute_flux(const DENS_MAT & n, const DENS_MAT & phi);
    void update_fermi_level();
    void report(int i);
    void exponential_electron_density();
    class LinearSolver * solver_;
    double alpha_;
    int sliceSize_, nNodes_, nfreeSlices_, nfixed_, nLambda_;
    SPAR_MAT K_;
    SPAR_MAT G_,G2_; // 1D & 2D grad mats = int N gradN dv
    DENS_VEC J_;
    DENS_VEC flux_;
    DENS_VEC dJ_;
    DENS_VEC lambda_;
    DENS_VEC F_;
    DENS_VEC Phi_;
    DENS_VEC n_;
    Array2D <bool> rhsMask_;
    double mobility_;
    double diffusivity_;
    double norm_, norm0_;
};

/**
 *  @class SchrodingerSolver
 *  @brief a manager class 
 */
class SchrodingerPoissonManager  {
  public:
    SchrodingerPoissonManager();
    ~SchrodingerPoissonManager(){};

    /** parser */
    bool modify(int narg, char **arg);

    /** initialize */
    SchrodingerPoissonSolver * initialize(   
      /*const*/ ATC_Coupling * atc,
      SchrodingerSolver * schrodingerSolver,
      PoissonSolver * poissonSolver,
      const PhysicsModel * physicsModel
    );
  protected:
    int maxConsistencyIter_;
    int maxConstraintIter_;
    bool oneD_;
    int oneDconserve_;
    double Ef_shift_;
    double safe_dEf_;
    double alpha_;
    double tol_;
    double mu_, D_;
};
} // namespace ATC
#endif
