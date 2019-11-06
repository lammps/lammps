#include "PoissonSolver.h"
#include "ATC_Coupling.h"
#include "FE_Engine.h"
#include "PhysicsModel.h"
#include "PrescribedDataManager.h"
#include "LinearSolver.h"
#include <utility>
#include <iostream>

using std::pair;



namespace ATC {

// ====================================================================
//  PoissonSolver
// ====================================================================
PoissonSolver::PoissonSolver(
    const FieldName fieldName,
    const PhysicsModel * physicsModel,
    const FE_Engine * feEngine,
    const PrescribedDataManager * prescribedDataMgr,
    /*const*/ ATC_Coupling * atc,
    const Array2D<bool> & rhsMask,
    const int solverType,
    bool parallel
)
  : atc_(atc),
    feEngine_(feEngine),
    prescribedDataMgr_(prescribedDataMgr),
    physicsModel_(physicsModel),
    fieldName_(fieldName),
    rhsMask_(rhsMask),
    linear_(false),
    solver_(NULL),
    solverNL_(NULL),
    tangent_(NULL),
    solverType_(solverType),
    solverTol_(0),
    solverMaxIter_(0),
    integrationType_(FULL_DOMAIN),
    parallel_(parallel)
{
  if (physicsModel_->has_linear_rhs(fieldName)) {
    linear_ = true;
    rhsMask_(fieldName,FLUX) = false; 
  }

  else {
    rhsMask_(fieldName,FLUX)   = true; 
    rhsMask_(fieldName,SOURCE) = true; 
  }
  
  if (prescribedDataMgr_->has_robin_source(fieldName)) {
    
    
    
    rhsMask_(fieldName,ROBIN_SOURCE) = true;
  }
}
// --------------------------------------------------------------------
PoissonSolver::~PoissonSolver() 
{ 
  if (tangent_) delete tangent_;
  if (solverNL_) delete solverNL_;
  if (solver_) delete solver_;
}

// --------------------------------------------------------------------
//  Parser
// --------------------------------------------------------------------

  bool PoissonSolver::modify(int /* narg */, char **arg)
{
  bool match = false;
  /*! \page man_poisson_solver fix_modify AtC poisson_solver 
      \section syntax
      fix_modify AtC poisson_solver mesh create <nx> <ny> <nz> <region-id> 
      <f|p> <f|p> <f|p>
      - nx ny nz = number of elements in x, y, z
      - region-id = id of region that is to be meshed
      - f p p  = perioidicity flags for x, y, z
      \section examples
      <TT> fix_modify AtC poisson_solver mesh create 10 1 1 feRegion p p p </TT>
      \section description
      Creates a uniform mesh in a rectangular region
      \section restrictions
      creates only uniform rectangular grids in a rectangular region
      \section related
      \section default
      none
  */
  int argIdx = 0;
  if (strcmp(arg[argIdx],"poisson_solver")==0) {
    argIdx++;
    if (strcmp(arg[argIdx],"mesh")==0) {
      argIdx++;
      // create a FE_Engine

      //feEngine_ = new FE_Engine(this); need alternate constructor?
      // send args to new engine
      // arg[0] = "mesh";
      // arg[1] = "create";
      // feEngine_->modify(narg,arg);

    }
  } 
  return match;
}
// --------------------------------------------------------------------
//  Initialize
// --------------------------------------------------------------------
void PoissonSolver::initialize(void)
{
  nNodes_ = feEngine_->num_nodes();

  if (atc_->source_atomic_quadrature(fieldName_))  
    integrationType_ = FULL_DOMAIN_ATOMIC_QUADRATURE_SOURCE;

  // compute penalty for Dirichlet boundaries
  if (prescribedDataMgr_->none_fixed(fieldName_))  
    throw ATC_Error("Poisson solver needs Dirichlet data");

  const BC_SET & bcs = (prescribedDataMgr_->bcs(fieldName_))[0];

  if (linear_) { // constant rhs 
    if (! solver_ ) {
      pair<FieldName,FieldName> row_col(fieldName_,fieldName_);
      Array2D <bool> rhsMask(NUM_FIELDS,NUM_FLUX);
      rhsMask = false; rhsMask(fieldName_,FLUX) = true;

      if (prescribedDataMgr_->has_robin_source(fieldName_)) {
        rhsMask(fieldName_,ROBIN_SOURCE) = true;
      }
      // compute stiffness for Poisson solve
      atc_->compute_rhs_tangent(row_col, rhsMask, atc_->fields(),
        stiffness_, FULL_DOMAIN, physicsModel_);
      // create solver
      solver_ = new LinearSolver(stiffness_,bcs,solverType_,LinearSolver::AUTO_HANDLE_CONSTRAINTS,parallel_);
    }
    else {

    // re-initialize
    solver_->initialize(&bcs);
    }
    if (solverTol_) solver_->set_tolerance(solverTol_);
    if (solverMaxIter_) solver_->set_max_iterations(solverMaxIter_);

  }
  else {
//  print_mask(rhsMask_);
    if ( solverNL_ )  delete solverNL_;
    tangent_ = new PhysicsModelTangentOperator(atc_,physicsModel_, rhsMask_, integrationType_, fieldName_); 

    solverNL_ = new NonLinearSolver(tangent_,&bcs,0,parallel_);
    
    if (solverTol_) solverNL_->set_residual_tolerance(solverTol_);
    if (solverMaxIter_) solverNL_->set_max_iterations(solverMaxIter_);
  }
}

// --------------------------------------------------------------------
//  Solve
// --------------------------------------------------------------------
bool PoissonSolver::solve(FIELDS & fields, FIELDS & rhs)
{
  atc_->compute_rhs_vector(rhsMask_, fields, rhs,
    integrationType_, physicsModel_);
  CLON_VEC f = column(fields[fieldName_].set_quantity(),0);
  CLON_VEC r = column(rhs[fieldName_].quantity(),0);
  bool converged = false;
  if (linear_) {converged = solver_->solve(f,r);}
  else         {converged = solverNL_->solve(f);}

  if (atc_->source_atomic_quadrature(fieldName_) 
    && LammpsInterface::instance()->atom_charge() ) set_charges(fields);
  return converged;
}
bool PoissonSolver::solve(DENS_MAT & field, const DENS_MAT & rhs) 
{

  CLON_VEC f = column(field,0);
  CLON_VEC r = column(rhs,0);
  bool converged = false;
  if (linear_) {converged = solver_->solve(f,r);}
  else         {converged = solverNL_->solve(f);}

  if (atc_->source_atomic_quadrature(fieldName_) 
    && LammpsInterface::instance()->atom_charge() ) set_charges(atc_->fields());
  return converged;
}

// --------------------------------------------------------------------
//  set charges on atoms
// --------------------------------------------------------------------
void PoissonSolver::set_charges(FIELDS & fields)
{
  FIELD_MATS sources;
  
  atc_->compute_sources_at_atoms(rhsMask_, fields, physicsModel_,sources);
  FIELD_MATS::const_iterator nField = sources.find(fieldName_);
  if (nField != sources.end()) {
    const DENS_MAT & electronCharges = nField->second;
    double *  q = LammpsInterface::instance()->atom_charge();
    int nLocal = atc_->nlocal();
    if (nLocal > 0) {
      const Array<int> & i2a = atc_->internal_to_atom_map();
      for (int i=0; i < nLocal; i++) {

        int atomIdx = i2a(i);
        q[atomIdx] = -electronCharges(i,0);
      }
    }
  }
}

} // namespace ATC
