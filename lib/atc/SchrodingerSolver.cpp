// ATC Headers
#include "SchrodingerSolver.h"
#include "ATC_Error.h"
#include "ATC_Coupling.h"
#include "LammpsInterface.h"
#include "PrescribedDataManager.h"
#include "PhysicsModel.h"
#include "LinearSolver.h"
#include "PoissonSolver.h"


#include <utility>

using std::pair;
using std::set;

const double tol = 1.e-8; 
const double zero_tol = 1.e-12; 
const double f_tol = 1.e-8; 

namespace ATC {

enum oneDconservationEnum {ONED_DENSITY=0, ONED_FLUX}; 

double fermi_dirac(const double E, const double T)
{
  double f = 1.0;
  if      (T > 0) f = 1.0 / ( exp(E/kBeV_/T)+1.0 );
  else if (E > 0) f = 0;
  return f;
};


  //--------------------------------------------------------
  //  Schrodinger solve
  //--------------------------------------------------------
  SchrodingerSolver::SchrodingerSolver(
    const FieldName fieldName,
    const PhysicsModel * physicsModel,
    const FE_Engine * feEngine,
    const PrescribedDataManager * prescribedDataMgr,
    /*const*/ ATC_Coupling * atc,
    const int solverType,
    bool parallel
)
  : atc_(atc),
    feEngine_(feEngine),
    prescribedDataMgr_(prescribedDataMgr),
    physicsModel_(physicsModel),
    fieldName_(fieldName),
    solver_(NULL),
    solverType_(solverType),
    nNodes_(atc->num_nodes()),
    parallel_(parallel)
  {
  }
  SchrodingerSolver::~SchrodingerSolver()
  {
    if (solver_) delete solver_;
  }

  void SchrodingerSolver::initialize()
  {
    SPAR_MAT sparseM; 
    atc_->fe_engine()->compute_mass_matrix(sparseM);
    M_ = sparseM.dense_copy();
  }

  bool SchrodingerSolver::solve(FIELDS & fields)
  {

// typedef  struct{float real, imag;} COMPLEX;
    SPAR_MAT stiffness_;
    Array2D <bool> rhsMask(NUM_FIELDS,NUM_FLUX);
    rhsMask = false;
    rhsMask(ELECTRON_WAVEFUNCTION,FLUX) = true;
    rhsMask(ELECTRON_WAVEFUNCTION,SOURCE) = true;
    pair<FieldName,FieldName> row_col(ELECTRON_WAVEFUNCTION,
                                      ELECTRON_WAVEFUNCTION);
    //set_fixed_nodes();
    atc_->fe_engine()->compute_tangent_matrix(
      rhsMask, row_col, atc_->fields(), physicsModel_,
      atc_->element_to_material_map(), stiffness_);
    DENS_MAT K(stiffness_.dense_copy());
    set<int> fixedNodes = prescribedDataMgr_->fixed_nodes(ELECTRON_WAVEFUNCTION);
    const BC_SET & bcs 
      = (prescribedDataMgr_->bcs(ELECTRON_WAVEFUNCTION))[0];
    DENS_MAT & psi   = (atc_->field(ELECTRON_WAVEFUNCTION)).set_quantity();
    DENS_MAT & eVecs = (atc_->field(ELECTRON_WAVEFUNCTIONS)).set_quantity();
    DENS_MAT & eVals = (atc_->field(ELECTRON_WAVEFUNCTION_ENERGIES)).set_quantity();

    if (prescribedDataMgr_->all_fixed(ELECTRON_WAVEFUNCTION)) {
      ATC::LammpsInterface::instance()->print_msg("all wavefunctions fixed");
      psi.reset(nNodes_,1);
      eVecs.reset(nNodes_,1);
      eVals.reset(nNodes_,1);
      return true;
    }
    // (1) Helmholtz solve for inhomongeneous bcs
    
    LinearSolver helmholtzSolver_(K,bcs,LinearSolver::AUTO_SOLVE,-1,parallel_);
    
    psi.reset(nNodes_,1);
    // (2) Eigenvalue solve 
    helmholtzSolver_.eigen_system(eVals,eVecs,&M_);
    return true; 
  }

  //--------------------------------------------------------
  //  Schrodinger solve on slices
  //--------------------------------------------------------
  SliceSchrodingerSolver::SliceSchrodingerSolver(
    const FieldName fieldName,
    const PhysicsModel * physicsModel,
    const FE_Engine * feEngine,
    const PrescribedDataManager * prescribedDataMgr,
    /*const*/ ATC_Coupling * atc,
    const Array< set<int> > & oneDslices,
    const int solverType,
    bool parallel
)
    : SchrodingerSolver(fieldName, physicsModel, feEngine, prescribedDataMgr, 
        atc, solverType, parallel),
    oneDslices_(oneDslices)
  {
  }
  SliceSchrodingerSolver::~SliceSchrodingerSolver()
  {
  }
  void SliceSchrodingerSolver::initialize()
  {
    SchrodingerSolver::initialize();
  }

  bool SliceSchrodingerSolver::solve(FIELDS & fields)
  {
    // fields
    DENS_MAT & psi   = (atc_->field(ELECTRON_WAVEFUNCTION)).set_quantity();
    DENS_MAT & eVecs = (atc_->field(ELECTRON_WAVEFUNCTIONS)).set_quantity();
    DENS_MAT & eVals = (atc_->field(ELECTRON_WAVEFUNCTION_ENERGIES)).set_quantity();
    psi.reset(nNodes_,1);
    eVecs.reset(nNodes_,nNodes_);
    eVals.reset(nNodes_,1);
    DENS_MAT & Ef = (atc_->field(FERMI_ENERGY)).set_quantity();
    DENS_MAT & n  = (atc_->field(ELECTRON_DENSITY)).set_quantity();
    DENS_MAT & T  = (atc_->field(ELECTRON_TEMPERATURE)).set_quantity();

    
    // stiffness = K + V M
    SPAR_MAT stiffness_;
    Array2D <bool> rhsMask(NUM_FIELDS,NUM_FLUX);
    rhsMask = false;
    rhsMask(ELECTRON_WAVEFUNCTION,FLUX) = true;
    rhsMask(ELECTRON_WAVEFUNCTION,SOURCE) = true;
    pair<FieldName,FieldName> row_col(ELECTRON_WAVEFUNCTION,
                                      ELECTRON_WAVEFUNCTION);
    atc_->fe_engine()->compute_tangent_matrix(
      rhsMask, row_col, atc_->fields(), physicsModel_,
      atc_->element_to_material_map(), stiffness_);
    DENS_MAT K(stiffness_.dense_copy());

    // Eigenvalue solve 
    DENS_MAT K1,M1;
    int nslices = oneDslices_.size();
    DENS_MAT b ;
    DENS_MAT evals1,evecs1 ;
    DENS_MAT n1 ;
    BCS bcs;
    set <int> one;
    one.insert(0);
    set <int> eindex;
    int iEVal = 0;
    for (int islice = 0; islice < nslices ; islice++) {
      set<int> & slice = oneDslices_(islice);
      int snodes = slice.size();
      prescribedDataMgr_->bcs(ELECTRON_WAVEFUNCTION,slice,bcs,true);
      const BC_SET & bc = bcs[0];
      int nfixed = bc.size();
      if (nfixed != snodes) {
        K.map(slice,slice,K1);
        M_.map(slice,slice,M1);
        LinearSolver eigensolver(K1,bc,LinearSolver::AUTO_SOLVE,-1,parallel_);
        // wave functions
        evals1.reset(snodes,1);
        evecs1.reset(snodes,snodes);
        eigensolver.eigen_system(evals1,evecs1,&M1);
        eindex.clear();
        for (int j = 0; j < snodes; j++) eindex.insert(iEVal++);
        eVals.insert(eindex,one,  evals1); 
        eindex.clear();
        for (int j = 0; j < snodes; j++) eindex.insert(j);
        eVecs.insert(slice,eindex,evecs1);
        // electron density
        n1.reset(snodes,1);
        
        set<int>::const_iterator iset;
        double aveE_f = 0;
        for (iset = slice.begin(); iset != slice.end(); iset++) { 
          int gnode = *iset; 
          aveE_f += Ef(gnode,0);
        }
        aveE_f /= snodes;

        int node = 0;
        for (iset = slice.begin(); iset != slice.end(); iset++) { // node
          int gnode = *iset; 
          double temp =  T(gnode,0);
          //double E_f = Ef(gnode,0);
          for (int mode = 0; mode < snodes-nfixed; mode++) {
            double Ei = evals1(mode,0);
            double E = Ei-aveE_f;
            double f = fermi_dirac(E,temp); 
            if (f  <  f_tol) break; // take advantage of E ordering 
            double psi1 = evecs1(node,mode); // 2nd index corresp to evals order
            n1(node,0) += psi1*psi1*f;
          }
          node++;
        }
        n.insert(slice,one,  n1); // note not "assemble"
      }
    }
    return true; 
  }

  //--------------------------------------------------------
  //  Schrodinger-Poisson Manager
  //--------------------------------------------------------
  SchrodingerPoissonManager::SchrodingerPoissonManager() :
    maxConsistencyIter_(0),
    maxConstraintIter_(0),
    oneD_(false),
    oneDconserve_(ONED_FLUX),
    Ef_shift_(0.),
    safe_dEf_(0.)
  {
  }
  SchrodingerPoissonManager::~SchrodingerPoissonManager()
  {
  }

  bool SchrodingerPoissonManager::modify(int narg, char **arg)
  {
    bool match = false;
    int argIndx = 0;
    if (strcmp(arg[argIndx],"self_consistency")==0) {
      argIndx++;
      maxConsistencyIter_ = atoi(arg[argIndx]);
      match = true;
    }
    else if (strcmp(arg[argIndx],"conserve")==0) {
      oneD_ = true;
      argIndx++;
      if (strcmp(arg[argIndx],"density")==0)  oneDconserve_ = ONED_DENSITY;
      else                                    oneDconserve_ = ONED_FLUX;
      argIndx++;
      maxConstraintIter_ = atoi(arg[argIndx]);
      match = true;
    }
    else if (strcmp(arg[argIndx],"initial_fermi_level")==0) {
      argIndx++;
      Ef_shift_ = atof(arg[argIndx]);
      match = true;
    }
    else if (strcmp(arg[argIndx],"safe_fermi_increment")==0) {
      argIndx++;
      safe_dEf_ = atof(arg[argIndx]);
      match = true;
    }
    return match;
  }

  SchrodingerPoissonSolver * SchrodingerPoissonManager::initialize(
    /*const*/ ATC_Coupling * atc,
    SchrodingerSolver * schrodingerSolver,
    PoissonSolver * poissonSolver,
    const PhysicsModel * physicsModel
   )
  {
    SchrodingerPoissonSolver * ptr;
    if (oneD_) {
      ptr = new SliceSchrodingerPoissonSolver(atc,
        schrodingerSolver,poissonSolver,physicsModel,maxConsistencyIter_,
        maxConstraintIter_, oneDconserve_, Ef_shift_, safe_dEf_);
    }
    else {
      ptr = new SchrodingerPoissonSolver(atc,
        schrodingerSolver,poissonSolver,physicsModel,maxConsistencyIter_);
    }
    return ptr;
  }

  //-------------------------------------------------------------------
  // SchrodingerPoissonSolver
  //-------------------------------------------------------------------
  SchrodingerPoissonSolver::SchrodingerPoissonSolver( 
    /*const*/ ATC_Coupling * atc,
    SchrodingerSolver * schrodingerSolver,
    PoissonSolver * poissonSolver,
    const PhysicsModel * physicsModel,
    int maxConsistencyIter
  ) :
   atc_(atc),
   schrodingerSolver_(schrodingerSolver),
   poissonSolver_(poissonSolver),
   physicsModel_(physicsModel),
   maxConsistencyIter_(maxConsistencyIter),
   nNodes_(atc_->num_nodes())
  {
  }
  SchrodingerPoissonSolver::~SchrodingerPoissonSolver(void)
  {
  }

  void SchrodingerPoissonSolver::solve(FIELDS & rhs, GRAD_FIELD_MATS & fluxes)
  {
    if ((atc_->prescribed_data_manager()->all_fixed(ELECTRON_WAVEFUNCTION))
     && (atc_->prescribed_data_manager()->all_fixed(ELECTRIC_POTENTIAL)))  {
      return;
    }
    double norm = 1.0, norm0 = 1.0; // normPrev = 1.0;
    DENS_MAT nPrev,psiPrev,phiPrev;

    DENS_MAT & psi = (atc_->field(ELECTRON_WAVEFUNCTIONS)).set_quantity();
    DENS_MAT & phi = (atc_->field(ELECTRIC_POTENTIAL)).set_quantity();
    DENS_MAT & E_I = (atc_->field(ELECTRON_WAVEFUNCTION_ENERGIES)).set_quantity();
    DENS_MAT & Te  = (atc_->field(ELECTRON_TEMPERATURE)).set_quantity();
    atc_->set_fixed_nodes();
    DENS_MAT Te0 = Te; // save

    const double tol  = 1.e-4;
//  double Tmax = Te.max();

    int k = 0;
    double logRatio = 3; 
    int maxIter = (int) logRatio; 
    double base = 2.0;

    // temperature relaxation loop
    for (int i = 0; i < maxIter ; ++i) { 
      //double alpha = ((double) i) /( (double) maxIter-1);
      //double beta = 0.1;
      //alpha = (exp(beta*i)-1.0)/(exp(beta*(maxIter-1))-1.0);
      double alpha = pow(base,logRatio-i-1);
      // self consistency loop
      int j = 0; // for storage of last iterate
      
      for (j = 0; j < maxConsistencyIter_ ; ++j) { 
        // compute eigen-values and vectors
        atc_->set_fixed_nodes();
        Te = alpha*Te0;
        
        schrodingerSolver_->solve(atc_->fields());

        
        for (int l = 0; l < nNodes_; l++) {
          int count = 0;
          double T_e = Te(l,0);
          for (int m = 0; m < nNodes_; m++) {
            double f = fermi_dirac(E_I(m,0), T_e);
            if (f > tol) count++; 
          }
        }
        // compute charge density
        DENS_MAN & n = atc_->field(ELECTRON_DENSITY);
        //(n.quantity()).print("DENSITY");
        atc_->nodal_projection(ELECTRON_DENSITY,physicsModel_,n);
        atc_->set_fixed_nodes(); 
        
        
        // solve poisson eqn for electric potential
        atc_->set_fixed_nodes();
        Te = alpha*Te0;
        poissonSolver_->solve(atc_->fields(),rhs);
        
        //DENS_MAT dn = n;
        //DENS_MAT dpsi = psi;
        //DENS_MAT dphi = phi;
        if (i == 0 && j==0) {
          nPrev = n.quantity();
          psiPrev = psi;
          phiPrev = phi;
        }
        //dn -= nPrev;
        //dpsi -= psiPrev;
        //dphi -= phiPrev;
        
        norm = (n.quantity()-nPrev).norm();
        if (i == 0 && j==0) norm0 = (n.quantity()).norm();
        //normPrev = norm;
        //psi_normPrev = psi_norm;
        //phi_normPrev = phi_norm;
        nPrev = n.quantity();
        psiPrev = psi;
        phiPrev = phi;
        k++;
        if (j > 0 && norm <= tol*norm0) break;
      }
      //      Tmax_ *= 0.5;
    }
  }
  
  //----------------------------------------------------------------------------
  // SchrodingerPoissonSolver
  //-------------------------------------------------------------------
  SliceSchrodingerPoissonSolver::SliceSchrodingerPoissonSolver( 
    /*const*/ ATC_Coupling * atc,
    SchrodingerSolver * schrodingerSolver,
    PoissonSolver * poissonSolver,
    const PhysicsModel * physicsModel,
    int maxConsistencyIter,
    int maxConstraintIter,
    int oneDconserve,
    double Ef_shift,
    double safe_dEf
  ) :
  SchrodingerPoissonSolver(atc,schrodingerSolver,poissonSolver,physicsModel,maxConsistencyIter),
  maxConstraintIter_(maxConstraintIter),
  oneDconserve_(oneDconserve),
  oneDcoor_(0),
  Ef_shift_(Ef_shift),
  safe_dEf_(safe_dEf),
  oneDslices_(((SliceSchrodingerSolver *) schrodingerSolver_)->slices())
  {
    EfHistory_.reset(oneDslices_.size(),2);
  }
  SliceSchrodingerPoissonSolver::~SliceSchrodingerPoissonSolver(void)
  {
  }
  void SliceSchrodingerPoissonSolver::solve(FIELDS & rhs, GRAD_FIELD_MATS & fluxes)
  {
    const double tol = 1.e-4; // tolerance on consistency & constraint
    double norm = 1.0, norm0 = 1.0;
    DENS_MAT nPrev;
    DENS_MAT & n   = (atc_->field(ELECTRON_DENSITY)).set_quantity();
    DENS_MAT & phi = (atc_->field(ELECTRIC_POTENTIAL)).set_quantity();

    // fermi energy
    DENS_MAT & Ef = (atc_->field(FERMI_ENERGY)).set_quantity();
    Ef.reset(nNodes_,1);

    int nslices = oneDslices_.size();
    Array2D<double> nHistory(nslices,2);

    // target for constraint
    double target = 0.0; 
    set<int> & slice = oneDslices_(0); // note assume first slice is fixed
    if (oneDconserve_ == ONED_FLUX) atc_->set_sources(); 
    DENS_MAT & nSource = (atc_->source(ELECTRON_DENSITY)).set_quantity();
    for (set<int>::const_iterator iset = slice.begin(); iset != slice.end(); iset++) { 
      if (oneDconserve_ == ONED_FLUX) target  += nSource(*iset,0);
      else                            target  += n(*iset,0);
    }
    target /= slice.size(); 
 
    // self consistency loop between Phi and n(psi_i)
    double error = 1.0;
    for (int i = 0; i < maxConsistencyIter_ ; ++i) { 
      atc_->set_fixed_nodes();
      if (! atc_->prescribedDataMgr_->all_fixed(ELECTRIC_POTENTIAL) ) 
        poissonSolver_->solve(atc_->fields(),rhs); 
      if (! atc_->prescribedDataMgr_->all_fixed(ELECTRON_DENSITY) )  {
        // iterate on Ef
        //if (i==0) Ef = -1.0*phi;// E ~ -|e| \Phi,  charge of electron e = 1 
        Ef = -1.0*phi; // E ~ -|e| \Phi,  charge of electron e = 1 in eV
        Ef +=Ef_shift_;
        for (int j = 0; j < maxConstraintIter_ ; ++j) {  

          schrodingerSolver_->solve(atc_->fields()); 

          atc_->set_fixed_nodes();
          error = update_fermi_energy(target,(j==0),fluxes);
          // exit condition based on constraint satisfaction
          if (error < tol*target) break; 
        } // loop j : flux constraint
        // error based on change in field (Cauchy convergence)
        if (i == 0) {
          norm = norm0 = n.norm();
        }
        else {
          DENS_MAT dn = n;
          dn -= nPrev;
          norm = dn.norm();
        }
        nPrev = n;
        if (i > 0 && norm <= tol*norm0 && error < tol) break;
      }
    } // loop i : self consistency
  }

  //--------------------------------------------------------
  //  update fermi energy
  //--------------------------------------------------------

  double SliceSchrodingerPoissonSolver::update_fermi_energy
    (double target, bool first, GRAD_FIELD_MATS & fluxes)
  {
    DENS_MAT & Ef  = (atc_->field(FERMI_ENERGY)).set_quantity();
    double safe_dEf = safe_dEf_; 

    DENS_MAT & n   = (atc_->field(ELECTRON_DENSITY)).set_quantity();
    const DENS_MAT * y = &n;
    if (oneDconserve_ == ONED_FLUX) { // compute J_x
      Array2D <bool> rhsMask(NUM_FIELDS,NUM_FLUX); rhsMask = false;
      rhsMask(ELECTRON_DENSITY,FLUX) = true;
      atc_->compute_flux(rhsMask,atc_->fields_,fluxes,physicsModel_);
      y = & ( fluxes[ELECTRON_DENSITY][oneDcoor_] ); 
    }

    BCS bcs;
    double error = 0;
    // slice
    for (int islice = 0; islice < oneDslices_.size(); islice++) {
      set<int> & slice = oneDslices_(islice);
      int nSlice = slice.size();
      //atc_->prescribedDataMgr_->bcs(ELECTRON_DENSITY,slice,bcs,true);
      atc_->prescribedDataMgr_->bcs(ELECTRON_WAVEFUNCTION,slice,bcs,true);
      const BC_SET & bc = bcs[0];
      int nFixed = bc.size();
      if (nFixed == nSlice) continue;
      double Y = 0.0, X = 0.0;
double nave = 0.0;
      for (set<int>::const_iterator iset = slice.begin(); iset != slice.end(); iset++) { 
        int gnode = *iset;
        X +=   Ef(gnode,0);
        Y += (*y)(gnode,0);
nave += n(gnode,0);
      }
      X /= nSlice;
      Y /= nSlice;
nave /= nSlice;
      
      double dY = Y - EfHistory_(islice,0); 
      double dX = X - EfHistory_(islice,1);
      if (fabs(dY) < zero_tol*dX) throw ATC_Error("zero increment in conserved field on slice");
      double err = target - Y; 
      if (target*Y < -zero_tol*target) {
        //throw ATC_Error("target and quantity opposite signs");
        ATC::LammpsInterface::instance()->print_msg_once("WARNING: target and quantity opposite signs");
      }
      error += fabs(err);
      //error = max(error,err);
      double dEf = err / dY * dX;
      if (first) {
        dEf = (err < 0) ? -safe_dEf : safe_dEf;
      }
      else if (fabs(dEf) > safe_dEf) {
        dEf = safe_dEf * dEf / fabs(dEf);
      }
      for (set<int>::const_iterator iset = slice.begin(); iset != slice.end(); iset++) { 
        int gnode = *iset;
        Ef(gnode,0) += dEf; 
      }
      EfHistory_(islice,0) = Y;
      EfHistory_(islice,1) = X;
    } // loop slice 
    return error;
  }

};
