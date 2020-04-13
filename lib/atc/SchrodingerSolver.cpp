// ATC Headers
#include "SchrodingerSolver.h"
#include "ATC_Error.h"
#include "ATC_Coupling.h"
#include "LammpsInterface.h"
#include "PrescribedDataManager.h"
#include "PhysicsModel.h"
#include "LinearSolver.h"
#include "PoissonSolver.h"

#include "Utility.h"

#include <utility>
using std::pair;
using std::set;
using std::stringstream;
using std::min;
using ATC_Utility::to_string;
using ATC_Utility::sgn;

const double zero_tol = 1.e-12; 
const double f_tol = 1.e-8; 

namespace ATC {

enum oneDconservationEnum {ONED_DENSITY=0, ONED_FLUX, ONED_GLOBAL_FLUX}; 





double fermi_dirac(const double E, const double T)
{
  double f = 1.0;
  if      (T > 0) f = 1.0 / ( exp(E/(kBeV_*T))+1.0 );
  else if (E > 0) f = 0;
  return f;
};


  //========================================================
  //  Schrodinger solve
  //========================================================
  SchrodingerSolver::SchrodingerSolver(
    const FieldName fieldName,
    const PhysicsModel * physicsModel,
    const FE_Engine * feEngine,
    const PrescribedDataManager * prescribedDataMgr,
    /*const*/ ATC_Coupling * atc,
    bool parallel
)
  : atc_(atc),
    feEngine_(feEngine),
    prescribedDataMgr_(prescribedDataMgr),
    physicsModel_(physicsModel),
    fieldName_(fieldName),
    nNodes_(atc->num_nodes()),
    parallel_(parallel)
  {
  }
  //-----------------------------------------------------
  void SchrodingerSolver::initialize()
  {
    SPAR_MAT sparseM; 
    atc_->fe_engine()->compute_mass_matrix(sparseM);
    M_ = sparseM.dense_copy();
  }
  //-----------------------------------------------------
  bool SchrodingerSolver::solve(FIELDS & /* fields */)
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

  //========================================================
  //  Schrodinger solve on slices
  //========================================================
  SliceSchrodingerSolver::SliceSchrodingerSolver(
    const FieldName fieldName,
    const PhysicsModel * physicsModel,
    const FE_Engine * feEngine,
    const PrescribedDataManager * prescribedDataMgr,
    /*const*/ ATC_Coupling * atc,
    const Array< set<int> > & oneDslices,
    const Array< double > & oneDdxs,
    bool parallel
)
    : SchrodingerSolver(fieldName, physicsModel, feEngine, prescribedDataMgr, 
        atc, parallel),
    oneDslices_(oneDslices),
    oneDdxs_(oneDdxs)
  {}
  //--------------------------------------------------------
  void SliceSchrodingerSolver::initialize()
  {
    SchrodingerSolver::initialize();
  }
  //--------------------------------------------------------
  // compute charge density per slice
  //--------------------------------------------------------
  bool SliceSchrodingerSolver::solve(FIELDS & /* fields */)
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
        // A: solve for e-values and wavefunctions
        K.map(slice,slice,K1);
        M_.map(slice,slice,M1);
        LinearSolver eigensolver(K1,bc,LinearSolver::AUTO_SOLVE,-1);
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
        // slice charge density
        n1.reset(snodes,1);
        
        set<int>::const_iterator iset;
        double aveE_f = 0;
        for (iset = slice.begin(); iset != slice.end(); iset++) { 
          int gnode = *iset; 
          aveE_f += Ef(gnode,0);
        }
        aveE_f /= snodes;
//#define VERBOSE 
#ifdef VERBOSE
        stringstream ss;
        ss << "   slice "+to_string(islice+1)+" E_f "+to_string(aveE_f) << "\n" 
           << "#-----------------------------------------------\n"
           << "#          E-Ef          f        psi          n\n"
           << "#-----------------------------------------------\n";
#endif
        // B: compute charge density on slice
        int node = 0;
        for (iset = slice.begin(); iset != slice.end(); iset++) { // node
          int gnode = *iset; 
          double temp =  T(gnode,0);
          for (int mode = 0; mode < snodes-nfixed; mode++) {
            double Ei = evals1(mode,0);
            double E = Ei-aveE_f;
            double f = fermi_dirac(E,temp); 
            double psi1 = evecs1(node,mode); // 2nd index corresp to evals order
#ifdef VERBOSE
            ss << node<<":"<<mode << "  " << to_string(6,E) << " " << to_string(6,f) << " " << to_string(6,psi1) << " " << to_string(6,n1(node,0)+psi1*psi1*f) <<  "\n";
#endif
            if (f  <  f_tol) break; // take advantage of E ordering 
            n1(node,0) += psi1*psi1*f;
          }
          node++;
        }
#ifdef VERBOSE
        ATC::LammpsInterface::instance()->print_msg_once(ss.str());
#endif
        n.insert(slice,one,  n1); // note not "assemble"
      }
    }
    return true; 
  }

  //========================================================
  //  Schrodinger-Poisson Manager
  //========================================================
  SchrodingerPoissonManager::SchrodingerPoissonManager() :
    maxConsistencyIter_(0),
    maxConstraintIter_(0),
    oneD_(false),
    oneDconserve_(ONED_FLUX),
    Ef_shift_(0.),
    safe_dEf_(0.),
    tol_(1.e-10),
    mu_(1.),D_(0.)
  {
  }
  //----------------------------------------------------------
  bool SchrodingerPoissonManager::modify(int /* narg */, char **arg)
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
      if (strcmp(arg[argIndx],"density")==0)   oneDconserve_ = ONED_DENSITY;
      else if (strcmp(arg[argIndx],"flux")==0) oneDconserve_ = ONED_FLUX;
      else                                     oneDconserve_ = ONED_GLOBAL_FLUX;
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
    else if (strcmp(arg[argIndx],"relaxation")==0) {
      argIndx++;
      alpha_ = atof(arg[argIndx]);
      match = true;
    }
    else if (strcmp(arg[argIndx],"tolerance")==0) {
      argIndx++;
      tol_ = atof(arg[argIndx]);
      match = true;
    }
    else if (strcmp(arg[argIndx],"mobility")==0) {
      argIndx++;
      mu_ = atof(arg[argIndx]);
      match = true;
    }
    else if (strcmp(arg[argIndx],"diffusivity")==0) {
      argIndx++;
      D_ = atof(arg[argIndx]);
      match = true;
    }
    return match;
  }
  //----------------------------------------------------------------
  SchrodingerPoissonSolver * SchrodingerPoissonManager::initialize(
    /*const*/ ATC_Coupling * atc,
    SchrodingerSolver * schrodingerSolver,
    PoissonSolver * poissonSolver,
    const PhysicsModel * physicsModel
   )
  {
    SchrodingerPoissonSolver * ptr;
    if (oneD_) {
      if (oneDconserve_ == ONED_GLOBAL_FLUX) {
      ptr = new GlobalSliceSchrodingerPoissonSolver(atc,
        schrodingerSolver,poissonSolver,physicsModel,maxConsistencyIter_,
        maxConstraintIter_, oneDconserve_, Ef_shift_, alpha_, safe_dEf_, tol_,
        mu_,D_);
      }
      else {
      ptr = new SliceSchrodingerPoissonSolver(atc,
        schrodingerSolver,poissonSolver,physicsModel,maxConsistencyIter_,
        maxConstraintIter_, oneDconserve_, Ef_shift_, safe_dEf_);
      }
    }
    else {
      ptr = new SchrodingerPoissonSolver(atc,
        schrodingerSolver,poissonSolver,physicsModel,maxConsistencyIter_);
    }
    return ptr;
  }

  //===================================================================
  // SchrodingerPoissonSolver
  //===================================================================
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
  //----------------------------------------------------------------------

  void SchrodingerPoissonSolver::solve(FIELDS & rhs, GRAD_FIELD_MATS & /* fluxes */)
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
  
  //===================================================================
  // SliceSchrodingerPoissonSolver
  //===================================================================
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
  oneDconserve_(oneDconserve),
  oneDcoor_(0),
  oneDslices_(((SliceSchrodingerSolver *) schrodingerSolver_)->slices()),
  oneDdxs_(((SliceSchrodingerSolver *) schrodingerSolver_)->dxs())
  {
    Ef_shift_=Ef_shift;
    safe_dEf_=safe_dEf;
    maxConstraintIter_=maxConstraintIter;
    EfHistory_.reset(oneDslices_.size(),2);
  }
  //--------------------------------------------------------------------------
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
#ifdef VERBOSE
    if (oneDconserve_ == ONED_FLUX) {
      if (target > 0) ATC::LammpsInterface::instance()->print_msg_once(" influx target "+ to_string(target));
      else            ATC::LammpsInterface::instance()->print_msg_once(" efflux target "+ to_string(target));
    }
#endif
 
    // A: self consistency loop between Phi and n(psi_i)
    double error = 1.0;
    for (int i = 0; i < maxConsistencyIter_ ; ++i) { 
      atc_->set_fixed_nodes();
      if (! atc_->prescribedDataMgr_->all_fixed(ELECTRIC_POTENTIAL) ) 
        poissonSolver_->solve(atc_->fields(),rhs); 
      if (! atc_->prescribedDataMgr_->all_fixed(ELECTRON_DENSITY) )  {
        // iterate on Ef
        //if (i==0) Ef = -1.0*phi;// E ~ -|e| \Phi,  charge of electron e = 1 
        Ef = -1.0*phi; 
        
        Ef +=Ef_shift_;
        // B: conservation constraint
        for (int j = 0; j < maxConstraintIter_ ; ++j) {  
          schrodingerSolver_->solve(atc_->fields()); // n(E_f)
          atc_->set_fixed_nodes();
          error = update_fermi_energy(target,(j==0),fluxes);// root finder
#ifdef VERBOSE
          ATC::LammpsInterface::instance()->print_msg_once(to_string(i)+":"+to_string(j)+" constraint_error "+to_string(error)+" / "+to_string(tol*target)+"\n");
#endif
          // exit condition based on constraint satisfaction
          if (error < tol*fabs(target)) break; 
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
#ifdef VERBOSE
#if 0
        if (i > 0) ATC::LammpsInterface::instance()->print_msg_once(to_string(i)+" density_change: "+to_string(norm)+" / "+to_string(norm0));
        else       ATC::LammpsInterface::instance()->print_msg_once("initial norm "+to_string(norm));
#endif
#endif
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
    DENS_MAT & n   = (atc_->field(ELECTRON_DENSITY)).set_quantity();
    DENS_MAT & phi   = (atc_->field(ELECTRIC_POTENTIAL)).set_quantity();
    const DENS_MAT * y = &n;
    if (oneDconserve_ == ONED_FLUX) { // compute J_x
      Array2D <bool> rhsMask(NUM_FIELDS,NUM_FLUX); rhsMask = false;
      rhsMask(ELECTRON_DENSITY,FLUX) = true;
//#define WIP_REJ
      atc_->compute_flux(rhsMask,atc_->fields_,fluxes,physicsModel_);
      y = & ( fluxes[ELECTRON_DENSITY][oneDcoor_] ); 
    }
    BCS bcs;
    double error = 0;
    // slice
    for (int islice = 0; islice < oneDslices_.size(); islice++) {
#ifdef VERBOSE
      std::string cStr(" conserved ");
      std::string Estr(" Ef");
#endif
      set<int> & slice = oneDslices_(islice);
      int nSlice = slice.size();
      atc_->prescribedDataMgr_->bcs(ELECTRON_WAVEFUNCTION,slice,bcs,true);
      const BC_SET & bc = bcs[0];
      int nFixed = bc.size();
      if (nFixed == nSlice) continue; // skip if all fixed 
      double Y = 0.0, X = 0.0;
      double nAve = 0., phiAve = 0.;
      for (set<int>::const_iterator iset = slice.begin(); iset != slice.end(); iset++) { 
        int gnode = *iset;
        X +=   Ef(gnode,0);
        Y += (*y)(gnode,0);
        nAve += n(gnode,0);
        phiAve += phi(gnode,0);
      }
      X /= nSlice;
      Y /= nSlice;
      nAve /= nSlice;
      phiAve /= nSlice;
      // now adjust Ef for each slice 
      double dY = Y - EfHistory_(islice,0); 
      double dX = X - EfHistory_(islice,1);
      double err = target - Y; 
      if (target*Y < -zero_tol*target) {
#ifdef VERBOSE
        cStr = " opp. SIGNS";
#else    
        ATC::LammpsInterface::instance()->print_msg_once("WARNING: slice "+to_string(islice)+" target and quantity opposite signs "+to_string(Y));
#endif
      }
      error += fabs(err);
      double dEf = 0.;
      if (first) {
        dEf = (err < 0) ? -safe_dEf_ : safe_dEf_;
      }
      else { 
        if (fabs(dY) < zero_tol*dX) throw ATC_Error("zero increment in conserved field on slice:"+to_string(islice));
        dEf = err / dY * dX;
        if (fabs(dEf) > safe_dEf_) {
          dEf = safe_dEf_* dEf / fabs(dEf);
#ifdef VERBOSE
          Estr = " !!";
#else 
          ATC::LammpsInterface::instance()->print_msg_once("WARNING: slice "+to_string(islice)+ " large Delta E_f "+to_string(dEf));
#endif
        }
      }
      for (set<int>::const_iterator iset = slice.begin(); iset != slice.end(); iset++) { 
        int gnode = *iset;
        Ef(gnode,0) += dEf;  
      }
      EfHistory_(islice,0) = Y;
      EfHistory_(islice,1) = X;
      if ( std::isnan(Y) ) throw ATC_Error("target on slice is not a number");
#ifdef VERBOSE
      ATC::LammpsInterface::instance()->print_msg_once("   slice"+to_string(islice,2) +cStr+to_string(4,Y/target) +Estr+to_string(4,X)+" n"+to_string(5,nAve)+" phi"+to_string(4,phiAve));
      //ATC::LammpsInterface::instance()->print_msg_once("   slice "+to_string(islice) +cStr+to_string(4,Y/target) +" E_f"+to_string(4,X)+dEstr+to_string(4,X-EfHistory_(std::max(0,islice-1),1))+" n"+to_string(4,nAve)+" phi"+to_string(4,phiAve)+"  "+to_string(nFixed)+" dn "+to_string(4,dnAve)+" dphi "+to_string(4,dphiAve));
#endif
    } // loop slice 
    return error;
  }

  //===================================================================
  // GlobalSliceSchrodingerPoissonSolver
  //===================================================================
  GlobalSliceSchrodingerPoissonSolver::GlobalSliceSchrodingerPoissonSolver( 
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
  ) :
  SliceSchrodingerPoissonSolver(atc,schrodingerSolver,poissonSolver,physicsModel,maxConsistencyIter,maxConstraintIter,oneDconserve,0,0),
  solver_(NULL), 
  mobility_(mu),diffusivity_(D)
  {
    Ef0_ = Ef0;
    alpha_ = alpha;
    safe_dEf_ = safe_dEf;
    if (safe_dEf_ < 1.e-20) throw ATC_Error("safe dE_f must be positive");
    ATC::LammpsInterface::instance()->print_msg("mobility:"+to_string(mobility_)+" diffusivity:"+to_string(diffusivity_));
    tol_ = tol;
    nslices_ = oneDslices_.size();
    sliceSize_ = (oneDslices_(0)).size();
    nNodes_ = nslices_*sliceSize_;
    flux_.reset(nNodes_);
    J_.reset(nslices_);
    //nfixed_ = 2;
    nfixed_ = 1;
    nfreeSlices_ = nslices_-nfixed_;
    nLambda_ = nslices_-1;
    lambda_.reset(nLambda_);
    dJ_.reset(nLambda_);
    F_.reset(nslices_);
    Phi_.reset(nslices_);
    n_.reset(nslices_);
    // form stiffness, lhs dirichlet bc, rhs homogeneous neumann bc
    //int m = nfreeSlices_;
    int m = nLambda_;
    DENS_MAT A(m,m);
    for (int i = 1; i < m; ++i) {
      A(i,i) = -2;
      if (i>0)   A(i,i-1) = 1;
      if (i<m-1) A(i,i+1) = 1;
    }
    A(0,0) = -2;
    A(0,1) =  1;
    A(m-1,m-1) = -2; 
    A(m-1,m-2) =  1; 
    //if (nfixed_ == 1) { A(m-1,m-1) = -1;  }
    double dx = oneDdxs_(0); 
    A *=  1./dx;
    A.print("stiffness",4);
    SPAR_MAT K(A);
    K_ = K;
    // form gradient  (account for lhs bc)
    int n = nslices_;
    DENS_MAT B(m,n);
    //for (int i = 0; i < m-1; ++i) {
    for (int i = 0; i < m; ++i) {
      B(i,i)   =-1;
      B(i,i+1) = 1; //B(i,i+2) = 1;
    }
    if (nfixed_ == 1) {
      B(m-1,n-2) = -1;
      B(m-1,n-1) =  1;
    }
    B.print("gradient",4);
    SPAR_MAT G(B);
    G_ = G;
    
    DENS_MAT C(nNodes_,nNodes_);
    // local to ATC nodemap: k --> gnode = *iset
    int k = 0;
    set<int>::const_iterator iset;
    for (int islice = 0; islice < nslices_; islice++) {
      set<int> & slice = oneDslices_(islice);
      for (iset = slice.begin(); iset != slice.end(); iset++) { 
        double v = 0.5/dx;
        if ( k < sliceSize_ || k+1 > (nslices_-1)*sliceSize_ ) v *=2.0;
        if (islice > 0) {          C(k,k-sliceSize_) += v; }
        else            {          C(k,k)            += v; }
        if (islice < nslices_-1) { C(k,k+sliceSize_) -= v; }
        else            {          C(k,k)            -= v; }
        k++;
      }
    }
    //C.print("2D gradient",4);
    SPAR_MAT G2(C);
    G2_ = G2;
    solver_ = new LinearSolver(K_); // for lambda
    rhsMask_.reset(NUM_FIELDS,NUM_FLUX); rhsMask_ = false;
    rhsMask_(ELECTRON_DENSITY,FLUX) = true;

    // report
    if (nfixed_ ==2)
      ATC::LammpsInterface::instance()->print_msg_once("schrodinger-poisson solver: Dirichlet INLET, Dirichlet; OUTLET");
    else if (nfixed_ ==1)
      ATC::LammpsInterface::instance()->print_msg_once("schrodinger-poisson solver: Dirichlet INLET, Neumann; OUTLET");
    else 
      ATC_Error("schrodinger-poisson solver:too many fixed");
  }
  GlobalSliceSchrodingerPoissonSolver::~GlobalSliceSchrodingerPoissonSolver(void)  {
    if (solver_) delete solver_;
  }
  //--------------------------------------------------------------------------
  void GlobalSliceSchrodingerPoissonSolver::solve(FIELDS & rhs, GRAD_FIELD_MATS & /* fluxes */)
  {
    const DENS_MAT & phi = (atc_->fields_[ELECTRIC_POTENTIAL]).quantity();
    const DENS_MAT & n   = (atc_->fields_[ELECTRON_DENSITY]  ).quantity();
    DENS_MAT       & Ef  = (atc_->field(FERMI_ENERGY)).set_quantity();
    Ef.reset(phi.nRows(),1);
    norm_ = norm0_ = 1.0;
    for (int i = 0; i < maxConstraintIter_ ; ++i) { 
      atc_->set_fixed_nodes();
      if (! atc_->prescribedDataMgr_->all_fixed(ELECTRIC_POTENTIAL) ) {
        poissonSolver_->solve(atc_->fields(),rhs); 
      }
      else {
        ATC::LammpsInterface::instance()->print_msg_once("WARNING: phi is fixed");
      }
      if (i == 0) { report(0); }
      if (! atc_->prescribedDataMgr_->all_fixed(ELECTRON_DENSITY) )  {
        update_fermi_level(); // update Ef = Ef0 +lambda
        schrodingerSolver_->solve(atc_->fields()); // updates n(E_f)
        //exponential_electron_density(); // surrogate
        compute_flux(n,phi); // compute J(n,phi) & dJ_
        solver_->solve(lambda_,dJ_); // conservation constraint
        //lambda_.print("lambda");
        //lambda_.print("[[J}}");
      }
      else {
        ATC::LammpsInterface::instance()->print_msg_once("WARNING: rho is fixed");
      }
      norm_ = dJ_.norm();
      report(i+1);
      if (i == 0 && norm_ > tol_) norm0_ = norm_;
      else { if (norm_ < tol_*norm0_) break; }
    } 
  }
  //--------------------------------------------------------------------------
  void GlobalSliceSchrodingerPoissonSolver::exponential_electron_density() 
  {
    std::cout << "******************HACK******************\n";

    DENS_MAT & n   = (atc_->fields_[ELECTRON_DENSITY]  ).set_quantity();
    DENS_MAT       & Ef  = (atc_->field(FERMI_ENERGY)).set_quantity();
    double T = 300;
    double n0 = 1.e-2;
    set<int>::const_iterator iset;
    for (int islice = 0; islice < nslices_; islice++) {
      set<int> & slice = oneDslices_(islice);
      double aveE_f = 0.0;
      for (iset = slice.begin(); iset != slice.end(); iset++) { 
        int gnode = *iset;
        aveE_f += Ef(gnode,0);
      }
      aveE_f /= slice.size();
      for (iset = slice.begin(); iset != slice.end(); iset++) { 
        int gnode = *iset;
        //std::cout << phi(gnode,0)+aveE_f << "\n";
        //n(gnode,0) = -n0*exp(-(phi(gnode,0)+aveE_f)/(kBeV_*T));
        //n(gnode,0) = -n0*exp((-phi(gnode,0))/(kBeV_*T));
        //n(gnode,0) = -n0*exp(aveE_f/(kBeV_*T));
        //n(gnode,0) =  aveE_f+0.01;
        //n(gnode,0) = aveE_f;
        //n(gnode,0) = phi(gnode,0);
        //n(gnode,0) = -n0*(phi(gnode,0)+aveE_f)/(kBeV_*T);
        n(gnode,0) = -n0*(aveE_f)/(kBeV_*T);
      }
    }
  }
  //--------------------------------------------------------------------------
  void GlobalSliceSchrodingerPoissonSolver::report(int i) 
  {
    const DENS_MAT & phi = (atc_->fields_[ELECTRIC_POTENTIAL]).quantity();
    const DENS_MAT & n   = (atc_->fields_[ELECTRON_DENSITY]  ).quantity();
    const DENS_MAT & Ef  = (atc_->field(FERMI_ENERGY)).quantity();
    set<int>::const_iterator iset;
    for (int islice = 0; islice < nslices_; islice++) {
      set<int> & slice = oneDslices_(islice);
      double Phi = 0.0;
      double N   = 0.0;
      double EF  = 0.0;
      for (iset = slice.begin(); iset != slice.end(); iset++) {
        int gnode = *iset;
        Phi += phi(gnode,0);
        N   +=   n(gnode,0);
        EF  +=  Ef(gnode,0);
      }
      Phi /= slice.size();
      Phi_(islice) = Phi; // average potential
      N /= slice.size();
      n_(islice) = N; // average electron density
      EF /= slice.size();
      F_(islice) = EF; // average Fermi level
    }
    stringstream header;
    header << "\n";
    header << "#----------------------------------------------------------------------\n";
    header << "#          [[J]]     lambda        E_f        phi          n          J\n";
    header << "#----------------------------------------------------------------------\n";
    if (i == 0) {
      ATC::LammpsInterface::instance()->write_file("slice.dat",header.str());
    }
    stringstream ss;
    ss << "\n";
    // first slice (fixed E_F)
    double dJ0 = J_(1)-J_(0);
    ss << to_string(1,2) << "*" << to_string(6,dJ0) << " " << to_string(6,0.) << " " << to_string(6,F_(0)) << " " << to_string(6,Phi_(0)) << " " << to_string(6,n_(0)) << " " << to_string(6,J_(0)) <<  "\n";
    // interior
    for (int j = 1; j < nslices_-1; ++j) {
       ss << to_string(j+1,2) << " " << to_string(6,dJ_(j-1)) << " " << to_string(6,lambda_(j-1)) << " " << to_string(6,F_(j)) << " " << to_string(6,Phi_(j)) << " " << to_string(6,n_(j)) << " " << to_string(6,J_(j)) << "\n";
    }
    // last slice (fixed E_F)
    double dJn = J_(nslices_-1)-J_(nslices_-2);
    int j = nslices_-1;
    double lambdaN = 0.;
    std::string space = "*";
    if (nfixed_ == 1) { 
      lambdaN = lambda_(nslices_-2); 
      space = " ";
    }
    ss << to_string(nslices_,2) << space << to_string(6,dJn) << " " << to_string(6,lambdaN) << " " << to_string(6,F_(j)) << " " << to_string(6,Phi_(j)) << " " << to_string(6,n_(j)) << " " << to_string(6,J_(j)) <<  "\n";
    stringstream is;
    is << "\n# iteration: " << to_string(i)+"/ "+to_string(maxConstraintIter_)+" constraint norm:"+to_string(6,norm_/norm0_) << " " << nslices_ << " slices";
    ATC::LammpsInterface::instance()->print_msg(is.str()+header.str()+ss.str());
    ATC::LammpsInterface::instance()->write_file("slice.dat",ss.str()+is.str()+"\n",std::ofstream::app);
  }
  //--------------------------------------------------------------------------
  void GlobalSliceSchrodingerPoissonSolver::compute_flux(
     const DENS_MAT & n, const DENS_MAT & phi)
  {
    DENS_VEC f(nNodes_);
    DENS_VEC gradphi(nNodes_);
    DENS_VEC gradn(nNodes_);
    int k = 0;
    set<int>::const_iterator iset;
    // grad phi
    for (int islice = 0; islice < nslices_; islice++) {
      set<int> & slice = oneDslices_(islice);
      for (iset = slice.begin(); iset != slice.end(); iset++) { 
        int gnode = *iset;
        f(k) = phi(gnode,0);
        k++;
      }
    }
    //f.print("phi");
    gradphi = G2_*f;
    //gradphi.print("grad phi");
    k = 0;
    // grad n
    for (int islice = 0; islice < nslices_; islice++) {
      set<int> & slice = oneDslices_(islice);
      for (iset = slice.begin(); iset != slice.end(); iset++) { 
        int gnode = *iset;
        f(k) = n(gnode,0);
        k++;
      }
    }
    //f.print("n");
    gradn = G2_*f;
    ////gradn.print("grad n");
    flux_.reset(nNodes_);
    for (k = 0; k < nNodes_; k++) {
      flux_(k) = -mobility_*f(k)*gradphi(k)-diffusivity_*gradn(k);
    }
    //flux_.print("flux");
    // per slice flux and diference
    k = 0;
    for (int islice = 0; islice < nslices_; islice++) {
      set<int> & slice = oneDslices_(islice);
      J_(islice) = 0;
      for (iset = slice.begin(); iset != slice.end(); iset++) { 
        J_(islice) += flux_(k);
        k++;
      }
      J_(islice) /= slice.size();
      //std::cout << islice << "  J " << J_(islice) << "\n";
    }
    //J_.print("J");
    dJ_ = G_*J_; 
  }
  //--------------------------------------------------------------------------
  void GlobalSliceSchrodingerPoissonSolver::update_fermi_level() 
  {

    DENS_MAT & Ef  = (atc_->field(FERMI_ENERGY)      ).set_quantity();
    DENS_MAT & phi = (atc_->field(ELECTRIC_POTENTIAL)).set_quantity();
    DENS_MAT & n   = (atc_->field(ELECTRON_DENSITY)  ).set_quantity();
    set<int>::const_iterator iset;
    for (int islice = 0; islice < nslices_; islice++) {
      set<int> & slice = oneDslices_(islice);
      double Phi = 0.;
      double N = 0.;
      //F_(islice) = Ef0_;
      if (islice > 0 && islice < nslices_-1) {
        F_(islice) += alpha_*lambda_(islice-1);
      }
      for (iset = slice.begin(); iset != slice.end(); iset++) {
        int gnode = *iset;
        Phi += phi(gnode,0);
        N += n(gnode,0);
      }
      Phi /= slice.size();
      Phi_(islice) = Phi; // average potential
      N /= slice.size();
      n_(islice) = N; // average electron density
      
      //F_(j) +=  min(fabs(alpha_*lambda),safe_dEf_)*sgn(lambda);
      for (iset = slice.begin(); iset != slice.end(); iset++) {
        int gnode = *iset;
        Ef(gnode,0) = F_(islice);
      }
    }
    //Ef.print("Ef");
  }
};
