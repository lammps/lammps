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
    const int solverType = ATC::LinearSolver::DIRECT_SOLVE,
    bool parallel = false
  );

  /** Destructor */
  virtual ~SchrodingerSolver();

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

  /** linear solver */
  LinearSolver * solver_;
  int solverType_;

  /** number of nodes */
  int nNodes_;

  /** stiffness matrix */

  //SPAR_MAT stiffness_;

  //SPAR_MAT massMatrix_;
  DENS_MAT M_;

  bool parallel_;

};

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
    const int solverType = ATC::LinearSolver::DIRECT_SOLVE,
    bool parallel = false
  );

  /** Destructor */
  virtual ~SliceSchrodingerSolver();

  /** parser */
  bool modify(int narg, char **arg){return false;}

  /** initialize */
  void initialize(void);

  /** solve */
  virtual bool solve(FIELDS & fields); 

  Array< std::set<int> > & slices(void){ return oneDslices_;} 

 protected:

  Array< std::set<int> > oneDslices_;

};

class SchrodingerPoissonSolver  {
  public:
    SchrodingerPoissonSolver(
      /*const*/ ATC_Coupling * atc,
      SchrodingerSolver * schrodingerSolver,
      PoissonSolver * poissonSolver,
      const PhysicsModel * physicsModel,
      int maxConsistencyIter
    );
    virtual ~SchrodingerPoissonSolver(void);
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
    int nNodes_;
};

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
    virtual ~SliceSchrodingerPoissonSolver(void);
    virtual void solve(
      FIELDS & rhs,
      GRAD_FIELD_MATS & fluxes
    );
  protected:
    double update_fermi_energy(double target, bool first, 
      GRAD_FIELD_MATS & fluxes);
    int maxConstraintIter_;
    int oneDconserve_;
    int oneDcoor_;
    double Ef_shift_;
    double safe_dEf_;
    Array< std::set<int> > & oneDslices_;
    Array2D<double> EfHistory_;
};

class SchrodingerPoissonManager  {
  public:
    SchrodingerPoissonManager();
    ~SchrodingerPoissonManager();

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
};
} // namespace ATC
#endif
