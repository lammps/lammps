// ATC Headers
#include "ExtrinsicModelElectrostatic.h"
#include "PhysicsModel.h"
#include "ATC_Error.h"
#include "FieldEulerIntegrator.h"
#include "ATC_Coupling.h"
#include "LammpsInterface.h"
#include "PrescribedDataManager.h"
#include "PoissonSolver.h"
#include "PerAtomQuantityLibrary.h"
#include "AtomToMoleculeTransfer.h"
#include "MoleculeSet.h"
#include "ChargeRegulator.h"
#include <set>

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::set;

static const double kTol_ = 1.0e-8; 
static const double tol_sparse = 1.e-30;//tolerance for compaction from dense

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModelElectrostatic
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ExtrinsicModelElectrostatic::ExtrinsicModelElectrostatic
                              (ExtrinsicModelManager * modelManager,
                               ExtrinsicModelType modelType,
                               string matFileName) :
    ExtrinsicModel(modelManager,modelType,matFileName),
    poissonSolverType_(DIRECT), // ITERATIVE | DIRECT
    poissonSolverTol_(0),
    poissonSolverMaxIter_(0),
    poissonSolver_(NULL),
    maxSolves_(0),
    baseSize_(0),
    chargeRegulator_(NULL),
    useSlab_(false),
    includeShortRange_(true),
    atomForces_(NULL),
    nodalAtomicCharge_(NULL),
    nodalAtomicGhostCharge_(NULL)
  {
     physicsModel_ = new PhysicsModelSpeciesElectrostatic(matFileName);
     // set up correct masks for coupling
     rhsMaskIntrinsic_.reset(NUM_FIELDS,NUM_FLUX);
     rhsMaskIntrinsic_ = false;
     if (atc_->track_charge()) {
        
       if (! chargeRegulator_) chargeRegulator_ = new ChargeRegulator(atc_);
     }
  }
  
  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ExtrinsicModelElectrostatic::~ExtrinsicModelElectrostatic()
  {
    if (poissonSolver_) delete poissonSolver_;
    if (chargeRegulator_) delete chargeRegulator_;
  }

  //--------------------------------------------------------
  // modify 
  //--------------------------------------------------------
  bool ExtrinsicModelElectrostatic::modify(int narg, char **arg)
  {
    bool match = false;
    int argIndx = 0;

    /** */

    if (strcmp(arg[argIndx],"poisson_solver")==0) {
      argIndx++;
      if (strcmp(arg[argIndx],"max_solves")==0) {
        argIndx++;
         maxSolves_ = atoi(arg[argIndx]) ; }
      
      else if (strcmp(arg[argIndx],"tolerance")==0) {
        argIndx++;
        poissonSolverTol_ = atof(arg[argIndx]);
      }
      else if (strcmp(arg[argIndx],"max_iterations")==0) {
        argIndx++;
        poissonSolverMaxIter_ = atoi(arg[argIndx]);
      }
      else if (strcmp(arg[argIndx],"iterative")==0) {
        poissonSolverType_ = ITERATIVE; }
      else {
        poissonSolverType_ = DIRECT; }
      match = true;
    } // end "poisson_solver"
    


    /** creates fixed charge on faceset
        units on surface charge density are lammps charge units / lammps length units ^ 2
        fix_modify ATC extrinsic fix_charge faceset_id value
    */
#ifdef CHARGED_SURFACE
    else if (strcmp(arg[argIndx],"fix_charge")==0) {
      argIndx++;
      string facesetName(arg[argIndx]);
      argIndx++;
      double chargeDensity = atof(arg[argIndx]);
      surfaceCharges_[facesetName] = chargeDensity;
      match = true;
    }

    /** */
    else if (strcmp(arg[argIndx],"unfix_charge")==0) {
      argIndx++;
      string fsetName(arg[argIndx]);
      throw ATC_Error("Ability to unfix charge not yet implemented");
      match = true;
    }
#endif
    else if (strcmp(arg[argIndx],"control")==0) {
      argIndx++;
      if (strcmp(arg[argIndx],"charge")==0) {
        argIndx++;
        if (!atc_->track_charge()) throw ATC_Error("must have charges to regulate");
        match = chargeRegulator_->modify(narg-argIndx,&arg[argIndx]);
      }
    }

    /** switch to use slabbing */
    else if (strcmp(arg[argIndx],"slab")==0) {
      argIndx++;
      if (strcmp(arg[argIndx],"on")==0) {
        useSlab_ = true;
        match = true;
      }
      else if (strcmp(arg[argIndx],"off")==0) {
        useSlab_ = false;
        match = true;
      }
    }

    /** switch to account for short range interaces */
    else if (strcmp(arg[argIndx],"short_range")==0) {
      argIndx++;
      if (strcmp(arg[argIndx],"on")==0) {
        includeShortRange_ = true;
        match = true;
      }
      else if (strcmp(arg[argIndx],"off")==0) {
        includeShortRange_ = false;
        match = true;
      }
    }

    return match;
  }

  //--------------------------------------------------------
  //  initialize
  //--------------------------------------------------------
  void ExtrinsicModelElectrostatic::construct_transfers()
  {
    // add charge density transfer operator
    if (atc_->track_charge()) {
      InterscaleManager & interscaleManager(atc_->interscale_manager());

      // make sure we have gradients at atoms
      VectorDependencyManager<SPAR_MAT * > * interpolantGradient = interscaleManager.vector_sparse_matrix("InterpolantGradient");
      if (!interpolantGradient) {
        interpolantGradient = new PerAtomShapeFunctionGradient(atc_);
        interscaleManager.add_vector_sparse_matrix(interpolantGradient,
                                                   "InterpolantGradient");
      }

      FundamentalAtomQuantity * atomicCharge = 
        interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_CHARGE);
      AtfShapeFunctionRestriction * nodalAtomicCharge = 
        new AtfShapeFunctionRestriction(atc_,atomicCharge,atc_->accumulant());
      interscaleManager.add_dense_matrix(nodalAtomicCharge,"NodalAtomicCharge");
      AtfShapeFunctionMdProjection * nodalAtomicChargeDensity =
        new AtfShapeFunctionMdProjection(atc_,nodalAtomicCharge,MASS_DENSITY);
      interscaleManager.add_dense_matrix(nodalAtomicChargeDensity,"NodalAtomicChargeDensity");


      // get the total charge and dipole moment at the node per molecule
      // small molecules require per atom quantities with ghosts
      const map<string,pair<MolSize,int> > & moleculeIds(atc_->molecule_ids());
      map<string,pair<MolSize,int> >::const_iterator molecule;
      PerAtomQuantity<double> * atomProcGhostCoarseGrainingPositions = interscaleManager.per_atom_quantity("AtomicProcGhostCoarseGrainingPositions");
      FundamentalAtomQuantity * charge = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_CHARGE,
                                                                                      PROC_GHOST);
      for (molecule = moleculeIds.begin(); molecule != moleculeIds.end(); molecule++) {
        const string moleculeName = molecule->first;
        SmallMoleculeSet * smallMoleculeSet = interscaleManager.small_molecule_set(moleculeName);
        // calculate nodal charge from the molecules
        AtomToSmallMoleculeTransfer<double> * moleculeCharge = 
          new AtomToSmallMoleculeTransfer<double>(atc_,charge,smallMoleculeSet);
        interscaleManager.add_dense_matrix(moleculeCharge,"MoleculeCharge"+moleculeName);
        MotfShapeFunctionRestriction * nodalAtomicMoleculeCharge = 
          new MotfShapeFunctionRestriction(moleculeCharge,
                                           interscaleManager.sparse_matrix("ShapeFunction"+moleculeName));
        interscaleManager.add_dense_matrix(nodalAtomicMoleculeCharge,"NodalMoleculeCharge"+moleculeName);
        AtfShapeFunctionMdProjection * nodalAtomicMoleculeChargeDensity =
          new AtfShapeFunctionMdProjection(atc_,nodalAtomicMoleculeCharge,MASS_DENSITY);
        interscaleManager.add_dense_matrix(nodalAtomicMoleculeChargeDensity,"NodalMoleculeChargeDensity"+moleculeName);

        // dipole moment density
        // calculate the dipole moment of the molecules
        SmallMoleculeCentroid * moleculeCentroid = static_cast<SmallMoleculeCentroid*>(interscaleManager.dense_matrix("MoleculeCentroid"+moleculeName));
        SmallMoleculeDipoleMoment * dipoleMoment = 
          new SmallMoleculeDipoleMoment(atc_,charge,smallMoleculeSet,atomProcGhostCoarseGrainingPositions,moleculeCentroid);
        interscaleManager.add_dense_matrix(dipoleMoment,"DipoleMoment"+moleculeName);
        MotfShapeFunctionRestriction * nodalAtomicMoleculeDipole = 
          new MotfShapeFunctionRestriction(dipoleMoment,
                                           interscaleManager.sparse_matrix("ShapeFunction"+moleculeName));
        interscaleManager.add_dense_matrix(nodalAtomicMoleculeDipole,"NodalMoleculeDipole"+moleculeName);
        AtfShapeFunctionMdProjection * nodalAtomicMoleculeDipoleDensity =
          new AtfShapeFunctionMdProjection(atc_,nodalAtomicMoleculeDipole,MASS_DENSITY);
        interscaleManager.add_dense_matrix(nodalAtomicMoleculeDipoleDensity,"NodalMoleculeDipoleDensity"+moleculeName);
      }
    }
    
  }

  //--------------------------------------------------------
  //  initialize
  //--------------------------------------------------------
  void ExtrinsicModelElectrostatic::initialize()
  {
    ExtrinsicModel::initialize();
    InterscaleManager & interscaleManager = atc_->interscale_manager();

    int nNodes = atc_->num_nodes();
    atomForces_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_FORCE);
    rhs_[ELECTRIC_POTENTIAL].reset(nNodes,1);

#ifdef CHARGED_SURFACE

    // set fixed potential surfaces form charged surfaces
    map<string,double>::const_iterator isurface;
    for (isurface = surfaceCharges_.begin(); isurface != surfaceCharges_.end(); isurface++)
      add_charged_surface(isurface->first,isurface->second);
#endif

    // set up poisson solver
    rhsMask_.reset(NUM_FIELDS,NUM_FLUX);
    rhsMask_ = false;
    for (int i = 0; i < NUM_FLUX; i++) {
      rhsMask_(ELECTRIC_POTENTIAL,i) = atc_->fieldMask_(ELECTRIC_POTENTIAL,i);
    }
    rhsMask_(ELECTRIC_POTENTIAL,FLUX) = false;// for poisson solve & rhs compute
    // need to create the bcs for the solver to configure properly
    atc_->set_fixed_nodes(); 
    if (poissonSolver_) delete poissonSolver_;


        int type = ATC::LinearSolver::ITERATIVE_SOLVE_SYMMETRIC;
    if (poissonSolverType_ == DIRECT) {
      type = ATC::LinearSolver::DIRECT_SOLVE;
    }
    poissonSolver_ = new PoissonSolver(ELECTRIC_POTENTIAL,
                                       physicsModel_, atc_->feEngine_,
                                       atc_->prescribedDataMgr_, atc_,
                                       rhsMask_,type, true);
    if (poissonSolverTol_) poissonSolver_->set_tolerance(poissonSolverTol_);
    if (poissonSolverMaxIter_) poissonSolver_->set_max_iterations(poissonSolverMaxIter_);
    poissonSolver_->initialize();

    // initialize localized Green's function for FE electric field correction
    if (atc_->track_charge() && includeShortRange_) {
      greensFunctions_.reserve(nNodes);

      // set up Green's function per node
      for (int i = 0; i < nNodes; i++) {
        set<int> localNodes;
        
        for (int j = 0; j < nNodes; j++)
          localNodes.insert(j);

        // call Poisson solver to get Green's function for node i
        DENS_VEC globalGreensFunction;
        poissonSolver_->greens_function(i,globalGreensFunction);
      
        // store green's functions as sparse vectors only on local nodes
        set<int>::const_iterator thisNode;
        SparseVector<double> sparseGreensFunction(nNodes);
        for (thisNode = localNodes.begin(); thisNode != localNodes.end(); thisNode++)
          sparseGreensFunction(*thisNode) = globalGreensFunction(*thisNode);
        greensFunctions_.push_back(sparseGreensFunction);
      }
    }


    if (atc_->track_charge()) {
      double *  q = LammpsInterface::instance()->atom_charge();
      if (!q) throw ATC_Error(" charge tracking requested but charge pointer is null");

      nodalAtomicCharge_ = interscaleManager.dense_matrix("NodalAtomicCharge");
      if (! nodalAtomicCharge_) {
        FundamentalAtomQuantity * atomCharge = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_CHARGE);
        nodalAtomicCharge_ = new AtfShapeFunctionRestriction(atc_,atomCharge,
                                                             atc_->accumulant());
        interscaleManager.add_dense_matrix(nodalAtomicCharge_,"NodalAtomicCharge");
      }
      if (atc_->groupbitGhost_) {
        nodalAtomicGhostCharge_ = interscaleManager.dense_matrix("NodalAtomicGhostCharge");
        if (! nodalAtomicGhostCharge_) {
          FundamentalAtomQuantity * ghostCharge = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_CHARGE, GHOST);
          
          PerAtomSparseMatrix<double> * ghostShapeFunctions = interscaleManager.per_atom_sparse_matrix("InterpolantGhost");
          if (!ghostShapeFunctions) {
            ghostShapeFunctions = new PerAtomShapeFunction(atc_,
                                                           interscaleManager.per_atom_quantity("AtomicGhostCoarseGrainingPositions"),
                                                           interscaleManager.per_atom_int_quantity("AtomGhostElement"),
                                                           GHOST);
            interscaleManager.add_per_atom_sparse_matrix(ghostShapeFunctions,"InterpolantGhost");
          }

          nodalAtomicGhostCharge_ = new AtfShapeFunctionRestriction(atc_,ghostCharge,
                                                                    ghostShapeFunctions);
          interscaleManager.add_dense_matrix(nodalAtomicGhostCharge_,"NodalAtomicGhostCharge");
        }
      }
    }

    if (chargeRegulator_) {
      if (! poissonSolver_) throw ATC_Error("passing of Poisson solver from ExtrinsicModelElectrostatic to ChargeRegulator failed");
      chargeRegulator_->assign_poisson_solver(poissonSolver_);
      chargeRegulator_->construct_methods();
      chargeRegulator_->initialize();
    }

    // set initial force
    post_force();
  }

  //--------------------------------------------------------
  //  pre final integration
  //--------------------------------------------------------
  void ExtrinsicModelElectrostatic::post_init_integrate()
  {
    if (chargeRegulator_) chargeRegulator_->apply_pre_force(atc_->dt());
  }
  //--------------------------------------------------------
  //  post force
  //--------------------------------------------------------
  void ExtrinsicModelElectrostatic::post_force()
  {

    if (chargeRegulator_) chargeRegulator_->apply_post_force(atc_->dt());
    // add in correction accounting for lumped mass matrix in charge density 
    // in atomistic part of domain & account for physics model fluxes,resets rhs
    
    
    // set Dirchlet data
    atc_->set_fixed_nodes();

    // set sources
    (atc_->prescribed_data_manager())->set_sources(atc_->time()+0.5*(atc_->dt()),atc_->sources());

    // compute Poisson equation RHS sources
    atc_->compute_rhs_vector(rhsMask_, atc_->fields_, rhs_, atc_->source_integration(), physicsModel_);
    
    // add atomic charges to rhs
    DENS_MAT & rhs = rhs_[ELECTRIC_POTENTIAL].set_quantity();
    if (atc_->track_charge()) {

      rhs += nodalAtomicCharge_->quantity();
      if (nodalAtomicGhostCharge_) {
        rhs += nodalAtomicGhostCharge_->quantity();
      }
    }

    

    

    // solve poisson eqn for electric potential
    // electron charge density added to Poisson RHS in solver
    DENS_MAT & potential = (atc_->field(ELECTRIC_POTENTIAL)).set_quantity();
    if ( maxSolves_ == 0 || (atc_->local_step() < maxSolves_) ) {
//potential.print("POT");
//    rhs.print("RHS");
    bool converged = poissonSolver_->solve(potential,rhs);
    if (! converged ) throw ATC_Error("Poisson solver did not converge in ExtrinsicModelElectrostatic");
    }

    // do this for intrinsic charges or effective electron charges at atoms
    if (atc_->track_charge() 
      || ( LammpsInterface::instance()->atom_charge() && atc_->source_atomic_quadrature(ELECTRIC_POTENTIAL) ) ) {
      _atomElectricalForce_.resize(atc_->nlocal(),atc_->nsd());
      add_electrostatic_forces(potential);
#ifdef CHARGED_SURFACE
      if (includeShortRange_)
        apply_charged_surfaces(potential);
#endif

      InterscaleManager & interscaleManager_ = atc_->interscale_manager(); 
      atomForces_ = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_FORCE); 
      (*atomForces_) += _atomElectricalForce_; // f_E in ours, f in lammps ultimately
    }
  }

  //--------------------------------------------------------
  //  output
  //--------------------------------------------------------
  void ExtrinsicModelElectrostatic::output(OUTPUT_LIST & outputData)
  {
    double scale = 1./(LammpsInterface::instance()->ftm2v());
    double localF[3];
    if (_atomElectricalForce_.nRows() > 0) {
      localF[0] = scale*(_atomElectricalForce_).col_sum(0);
      localF[1] = scale*(_atomElectricalForce_).col_sum(1);
      localF[2] = scale*(_atomElectricalForce_).col_sum(2);
    }
    else {
      localF[0] = 0.;localF[1] = 0.; localF[2] = 0.;
    }
    LammpsInterface::instance()->allsum(localF,totalElectricalForce_,3);  
    if (LammpsInterface::instance()->rank_zero()) {
      atc_->feEngine_->add_global("electrostatic_force_x",  totalElectricalForce_[0]);
      atc_->feEngine_->add_global("electrostatic_force_y",  totalElectricalForce_[1]);
      atc_->feEngine_->add_global("electrostatic_force_z",  totalElectricalForce_[2]);
    }

    // add in FE fields related to charge
    FIELDS & fields(atc_->fields());
    FIELDS::const_iterator rhoField = fields.find(CHARGE_DENSITY);
    if (rhoField!=fields.end()) {
      InterscaleManager & interscaleManager(atc_->interscale_manager());
      const DENS_MAN * atomicChargeDensity(interscaleManager.dense_matrix("NodalAtomicChargeDensity"));
      atc_->nodal_atomic_field(CHARGE_DENSITY) = atomicChargeDensity->quantity();
      
      fields[CHARGE_DENSITY] = atomicChargeDensity->quantity();
      
      DENS_MAT & chargeDensity(fields[CHARGE_DENSITY].set_quantity());
      DENS_MAT & nodalAtomicChargeDensity((atc_->nodal_atomic_field(CHARGE_DENSITY)).set_quantity());
      if ((atc_->lammps_interface())->rank_zero()) {
        outputData["charge_density"] = &chargeDensity;
        outputData["NodalAtomicChargeDensity"] = &nodalAtomicChargeDensity;
      }
    }

    
    if (fields.find(ELECTRON_DENSITY)==fields.end()) {
      fields[ELECTRON_DENSITY].reset(fields[CHARGE_DENSITY].nRows(),1);
      DENS_MAT & electronDensity(fields[ELECTRON_DENSITY].set_quantity());
      if ((atc_->lammps_interface())->rank_zero()) {
        outputData["electron_density"] = &electronDensity;
      }
    }
    
    const map<string,pair<MolSize,int> > & moleculeIds(atc_->molecule_ids());
    map<string,pair<MolSize,int> >::const_iterator molecule;
    for (molecule = moleculeIds.begin(); molecule != moleculeIds.end(); molecule++) {
      // net charge
      DENS_MAN & nodalMoleculeChargeDensityOut(atc_->tagged_dens_man("NodalMoleculeChargeDensity"+molecule->first));
      DENS_MAN * nodalMoleculeChargeDensity((atc_->interscale_manager()).dense_matrix("NodalMoleculeChargeDensity"+molecule->first));
      nodalMoleculeChargeDensityOut = nodalMoleculeChargeDensity->quantity();
      // dipole moment
      DENS_MAN & nodalMoleculeDipoleDensityOut(atc_->tagged_dens_man("NodalMoleculeDipoleDensity"+molecule->first));
      DENS_MAN * nodalMoleculeDipoleDensity((atc_->interscale_manager()).dense_matrix("NodalMoleculeDipoleDensity"+molecule->first));
      nodalMoleculeDipoleDensityOut = nodalMoleculeDipoleDensity->quantity();
    }

    if(chargeRegulator_) chargeRegulator_->output(outputData);
  }

  //--------------------------------------------------------
  //  size_vector
  //--------------------------------------------------------
  int ExtrinsicModelElectrostatic::size_vector(int intrinsicSize)
  {
    baseSize_ = intrinsicSize;
    return 5;
  }

  //--------------------------------------------------------
  //  compute_scalar : added energy = - f.x
  //--------------------------------------------------------
  double ExtrinsicModelElectrostatic::compute_scalar(void)
  {
    //((atc_->interscale_manager()).fundamental_atom_quantity(LammpsInterface::ATOM_POSITION))->force_reset();
    const DENS_MAT & atomPosition = ((atc_->interscale_manager()).fundamental_atom_quantity(LammpsInterface::ATOM_POSITION))->quantity();
    double local_fdotx = 0, fdotx;
    
    for (int i = 0; i < _atomElectricalForce_.nRows() ; i++) {
      for (int j = 0; j < _atomElectricalForce_.nCols() ; j++) {
        local_fdotx -= _atomElectricalForce_(i,j)*atomPosition(i,j);
      }
    }
    LammpsInterface::instance()->allsum(&local_fdotx,&fdotx,1);
    // convert
    fdotx *=  LammpsInterface::instance()->mvv2e();
    return fdotx;
  }
  //--------------------------------------------------------
  //  compute_vector
  //--------------------------------------------------------
  bool ExtrinsicModelElectrostatic::compute_vector(int n, double & value)
  {
    if (n == baseSize_) {
      
      double nSum = ((atc_->field(ELECTRON_DENSITY)).quantity()).col_sum();
      value = nSum;
      return true;
    }
    else if (n > baseSize_ && n < baseSize_+4) {
      int dof = n-baseSize_-1;
      double localF = (_atomElectricalForce_).col_sum(dof), F=0;
      LammpsInterface::instance()->allsum(&localF,&F,1);
      double ftm2v = LammpsInterface::instance()->ftm2v();
      value = F/ftm2v;
      return true;
    }
    else if (n == baseSize_+4) {
      double nSum = ((atc_->field(ELECTRIC_POTENTIAL)).quantity()).col_sum();
      value = nSum;
      return true;
    }
    return false;
  }

  //--------------------------------------------------------
  //  add_electrostatic_forces
  //--------------------------------------------------------
  void ExtrinsicModelElectrostatic::add_electrostatic_forces
    (MATRIX & potential)
  {
    
    //double qE2f = LammpsInterface::instance()->qe2f();
    double qV2e = LammpsInterface::instance()->qv2e(); // charge volts to our energy units
    //double ** f = LammpsInterface::instance()->fatom();
    double *  q = LammpsInterface::instance()->atom_charge();
    // f_ai = \sum_IJ N_Ia Bi_IJ phi_J = \sum_I N_Ia Ei_I
    int nsd = atc_->nsd();
    int nLocal = atc_->nlocal();
    DENS_MAT E(nLocal,nsd);
    const SPAR_MAT_VEC & shapeFucntionDerivatives(((atc_->interscale_manager()).vector_sparse_matrix("InterpolantGradient"))->quantity());
    if (nLocal > 0) {
      for (int i=0; i < nsd; i++) {
        CLON_VEC Ei = column(E,i);
        Ei = -1.*(*(shapeFucntionDerivatives[i])*potential);
      }
    }
   
    int dimOffset = 0;
    if (useSlab_) dimOffset = nsd - 1;
    for (int i = 0; i < nLocal; i++) {
      int atomIdx = atc_->internalToAtom_(i);
      double c = qV2e*q[atomIdx];
      for (int j = 0; j < dimOffset; j ++)
        _atomElectricalForce_(i,j) = 0.;
      for (int j = dimOffset; j < nsd; j ++)
        _atomElectricalForce_(i,j) = c*E(i,j);
    }
    
    // correct field for short range interactions
    if (includeShortRange_)
      correct_electrostatic_forces();
  }

  //--------------------------------------------------------
  //  correct_electrostatic_forces
  //--------------------------------------------------------
  void ExtrinsicModelElectrostatic::correct_electrostatic_forces()
  {
    // compute restricted sparse shape function set for each atom
    // to account for its Green's Function
    //double qE2f = LammpsInterface::instance()->qe2f();
    double qV2e = LammpsInterface::instance()->qv2e();
    double *  q = LammpsInterface::instance()->atom_charge();
    vector<SparseVector<double> > atomicFePotential;
    int nLocal = atc_->nlocal();
    
    int nGhostLammps = LammpsInterface::instance()->nghost();
    int nLocalLammps = LammpsInterface::instance()->nlocal();
    int nLocalTotal = nLocalLammps + nGhostLammps; // total number of atoms on this processor
    atomicFePotential.reserve(nLocalTotal);
    SparseVector<double> dummy(atc_->num_nodes());
    for (int i = 0; i < nLocalTotal; i++)
      atomicFePotential.push_back(dummy);

    // compute local potential contributions from atoms on this processor
    InterscaleManager & interscaleManager(atc_->interscale_manager());
    const SPAR_MAT & myShpFcn((interscaleManager.per_atom_sparse_matrix("Interpolant"))->quantity());
    for (int i = 0; i < nLocal; i++) {
      DenseVector<INDEX> nodeIndices;
      DENS_VEC nodeValues;
      myShpFcn.row(i,nodeValues,nodeIndices);
      
      int atomIdx = atc_->internalToAtom_(i);
      //double c = qE2f*q[atomIdx];
      //double c = qV2e*q[atomIdx];
      //nodeValues *= c;
      nodeValues *= q[atomIdx];
      
      for (int j = 0; j < nodeIndices.size(); j++)
        atomicFePotential[atomIdx].add_scaled(greensFunctions_[nodeIndices(j)],nodeValues(j));
    }

    // compute local potential contribtutions for lammps ghost atoms 
    // which are known to ATC,
    // this will grab both processor and periodic neighbors,
    // so we need to add in neighbor contributions using lammps indices
    // rather than atc indices or we could potentially
    // double count periodic contributions
    double ** xatom = LammpsInterface::instance()->xatom();
    const int * mask = LammpsInterface::instance()->atom_mask();
    int nodesPerElement = ((atc_->feEngine_)->fe_mesh())->num_nodes_per_element();
    int nsd = atc_->nsd();
    for (int i = nLocalLammps; i < nLocalTotal; i++) {
      if (mask[i] & atc_->groupbit_) {
        DENS_VEC coords(nsd);
        coords.copy(xatom[i],nsd);
        Array<int> nodeIndices(nodesPerElement);
        DENS_VEC nodeValues(nodesPerElement);
        (atc_->feEngine_)->shape_functions(coords,nodeValues,nodeIndices);
        
        //double c = qV2e*q[i];
        //nodeValues *= c;
        nodeValues *= q[i];

        for (int j = 0; j < nodeIndices.size(); j++) {
          atomicFePotential[i].add_scaled(greensFunctions_[nodeIndices(j)],nodeValues(j));
        }
      }
    }

    // Get sparse vectors of derivatives at each atom
    
    //      to compute this only when the shape functions change
    vector<vector<SparseVector<double> > > atomicDerivatives;
    atomicDerivatives.reserve(nLocal);
    for (int i = 0; i < nLocal; i++) {

      // determine shape function derivatives at atomic location
      // and construct sparse vectors to store derivative data
      vector<SparseVector<double> > derivativeVectors;
      derivativeVectors.reserve(nsd);
      for (int j = 0; j < nsd; j++)
        derivativeVectors.push_back(dummy);
      atomicDerivatives.push_back(derivativeVectors);
      InterscaleManager & interscaleManager(atc_->interscale_manager());
      const SPAR_MAT_VEC & shapeFucntionDerivatives((interscaleManager.vector_sparse_matrix("InterpolantGradient"))->quantity());
      for (int j = 0; j < nsd; j++) {
        DenseVector<INDEX> nodeIndices;
        DENS_VEC nodeValues;
        shapeFucntionDerivatives[j]->row(i,nodeValues,nodeIndices);
        for (int k = 0; k < nodeIndices.size(); k++)
          atomicDerivatives[i][j](nodeIndices(k)) = nodeValues(k);
      }
    }

    // loop over all atoms and correct their efield based on all their 
    // neighbor's local efield response
    
    //      need to use specific coulombic cutoff from different pairs
    //      see pair_coul_cut for an example of the data structures
    //      unfortunately don't know how to get at this data in general
    //      beyond a cast from the LAMMPS pair object (see force.h).
    //      Until this is fixed, only use this method with the coulombic force
    //      the same for all pairs and equal to the largest force cutoff.
    //      Probably the best fix is to implement our own pair style for this.

    double cutoffRadius = LammpsInterface::instance()->pair_cutoff();
    double cutoffSq = cutoffRadius*cutoffRadius;
    
    int inum = LammpsInterface::instance()->neighbor_list_inum();
    int * ilist = LammpsInterface::instance()->neighbor_list_ilist();
    int * numneigh = LammpsInterface::instance()->neighbor_list_numneigh();
    int ** firstneigh = LammpsInterface::instance()->neighbor_list_firstneigh();
    
      // loop over neighbors of my atoms
    for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii];
      if (mask[i] & atc_->groupbit_) {
        double xtmp = xatom[i][0];
        double ytmp = xatom[i][1];
        double ztmp = xatom[i][2];
        
        int * jlist = firstneigh[i];
        int jnum = numneigh[i];
        
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          if (mask[j] & atc_->groupbit_) {
            //double factor_coul = LammpsInterface::instance()->coulomb_factor(j);
            LammpsInterface::instance()->neighbor_remap(j);

            double delx = xtmp - xatom[j][0];
            double dely = ytmp - xatom[j][1];
            double delz = ztmp - xatom[j][2];
            double rsq = delx*delx + dely*dely + delz*delz;      
            if (rsq < cutoffSq) {
              DENS_VEC efield(nsd);
              efield = 0.;
              int atcIdx = atc_->atomToInternal_[i];
              for (int k = 0; k < nsd; k++)
                efield(k) = -1.*dot(atomicDerivatives[atcIdx][k],atomicFePotential[j]);

              // apply correction in atomic forces
              //double c = factor_coul*qE2f*q[i];
              //double c = factor_coul*qV2e*q[i];
              double c = qV2e*q[i];
              for (int k = 0; k < nsd; k++) {
                if ((!useSlab_) || (k==nsd)) {
                  //f[i][k] -= c*efield(k);
                  _atomElectricalForce_(atcIdx,k) -= c*efield(k);
                }
              }
            }
          }
        }
      }
    }
  }

#ifdef CHARGED_SURFACE

  //--------------------------------------------------------
  //  add_charged_surface
  //--------------------------------------------------------
  void ExtrinsicModelElectrostatic::add_charged_surface(const string & facesetName,
                                                        const double chargeDensity)
  {
    // get faceset information
    int nNodes = atc_->num_nodes();
    const FE_Mesh * feMesh = (atc_->feEngine_)->fe_mesh();
    const set< pair <int,int> > * faceset 
      = & ( feMesh->faceset(facesetName));

    // set face sources to all point at one function for use in integration
    SURFACE_SOURCE faceSources;
    XT_Function * f = XT_Function_Mgr::instance()->constant_function(1.);
    set< pair<int,int> >::const_iterator iset;
    for (iset = faceset->begin(); iset != faceset->end(); iset++) {
      pair<int,int>  face = *iset;
      // allocate
      Array < XT_Function * > & dof  = faceSources[ELECTRIC_POTENTIAL][face];
      dof.reset(1);
      dof(0) = f;
    }

    // Get associated nodeset
    set<int> nodeset;
    feMesh->faceset_to_nodeset(facesetName,nodeset);

    // Get coordinates of each node in face set
    map<int,pair<DENS_VEC,double> > & myFaceset = chargedSurfaces_[facesetName];
    set<int>::const_iterator myNode;
    for (myNode = nodeset.begin(); myNode != nodeset.end(); myNode++) {
      DENS_VEC myCoords = feMesh->nodal_coordinates(*myNode);
      pair<DENS_VEC,double> myPair(myCoords,0.);
      myFaceset[*myNode] = myPair;
    }

    // computed integrals of nodal shape functions on face
    FIELDS nodalFaceWeights;
    nodalFaceWeights[ELECTRIC_POTENTIAL].reset(nNodes,1);
    Array<bool> fieldMask(NUM_FIELDS);
    fieldMask(ELECTRIC_POTENTIAL) = true;
    
    (atc_->feEngine_)->add_fluxes(fieldMask,0.,faceSources,nodalFaceWeights);
    // set up data structure holding charged faceset information
    FIELDS sources;
    double coulombConstant = LammpsInterface::instance()->coulomb_constant();
    map<int,pair<DENS_VEC,double> >::iterator myNodeData;
    for (myNodeData = myFaceset.begin(); myNodeData != myFaceset.end(); myNodeData++) {
      // evaluate voltage at each node I
      // set up X_T function for integration: k*chargeDensity/||x_I - x_s||
      
      // integral is approximated in two parts:
      // 1) near part with all faces within r < rcrit evaluated as 2 * pi * rcrit * k sigma A/A0, A is area of this region and A0 = pi * rcrit^2, so 2 k sigma A / rcrit
      // 2) far part evaluated using Gaussian quadrature on faceset
      double rcritSq = LammpsInterface::instance()->pair_cutoff();
      rcritSq *= rcritSq;
      int nodalIndex = myNodeData->first;
      DENS_VEC myCoords((myNodeData->second).first);
      double xtArgs[8];
      xtArgs[0] = myCoords(0); xtArgs[1] = myCoords(1); xtArgs[2] = myCoords(2);
      xtArgs[3] = 1.; xtArgs[4] = 1.; xtArgs[5] = 1.;
      xtArgs[6] = coulombConstant*chargeDensity;
      xtArgs[7] = -1.;
      string radialPower = "radial_power";
      f = XT_Function_Mgr::instance()->function(radialPower,8,xtArgs);

      
      for (iset = faceset->begin(); iset != faceset->end(); iset++) {
        pair<int,int>  face = *iset;
        // allocate
        Array < XT_Function * > & dof  = faceSources[ELECTRIC_POTENTIAL][face];
        dof.reset(1);
          dof(0) = f;
      }

      // perform integration to get quantities at nodes on facesets
      // V_J' = int_S N_J k*sigma/|x_I - x_s| dS
      sources[ELECTRIC_POTENTIAL].reset(nNodes,1);
      (atc_->feEngine_)->add_fluxes(fieldMask,0.,faceSources,sources);

      double myPotential = 0.;
      // sum over all nodes in faceset to get total potential:
      // V_I = sum_J VJ'
      const DENS_MAT & myPotentialSource(sources[ELECTRIC_POTENTIAL].quantity());
      nodalChargePotential_[facesetName][nodalIndex] = myPotentialSource(nodalIndex,0);
      for (myNode = nodeset.begin(); myNode != nodeset.end(); myNode++)
        myPotential += myPotentialSource(*myNode,0);

      // assign an XT function per each node and
      // then call the prescribed data manager and fix each node individually.
      f = XT_Function_Mgr::instance()->constant_function(myPotential);
      (atc_->prescribedDataMgr_)->fix_field(nodalIndex,ELECTRIC_POTENTIAL,0,f);

      // compute effective charge at each node I
      // multiply charge density by integral of N_I over face
      (myNodeData->second).second = (nodalFaceWeights[ELECTRIC_POTENTIAL].quantity())(nodalIndex,0)*chargeDensity;
    }
  }

  //--------------------------------------------------------
  //  apply_charged_surfaces
  //--------------------------------------------------------
  void ExtrinsicModelElectrostatic::apply_charged_surfaces
    (MATRIX & potential)
  {
    //double qE2f = LammpsInterface::instance()->qe2f();
    double qV2e = LammpsInterface::instance()->qv2e();
    double qqrd2e = LammpsInterface::instance()->qqrd2e();
    //double ** fatom = LammpsInterface::instance()->fatom();
    //double *  qatom = LammpsInterface::instance()->atom_charge();
    InterscaleManager & interscaleManager(atc_->interscale_manager());
    const DENS_MAT & qatom((interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_CHARGE))->quantity());
    double cutoffRadius = LammpsInterface::instance()->pair_cutoff();
    double cutoffSq = cutoffRadius*cutoffRadius;
    int nLocal = qatom.nRows();
    int nsd = atc_->nsd();
    int nNodes = atc_->num_nodes();
    double penalty = poissonSolver_->penalty_coefficient();
    if (penalty <= 0.0) throw ATC_Error("ExtrinsicModelElectrostatic::apply_charged_surfaces expecting non zero penalty");
    SparseVector<double> dummy(atc_->num_nodes());

    map<string,map<int,pair<DENS_VEC,double> > >::const_iterator isurface;
    for (isurface = chargedSurfaces_.begin(); isurface != chargedSurfaces_.end(); isurface++) {
      string facesetName = isurface->first;
      map<int,pair<DENS_VEC,double> >::const_iterator inode;
      for (inode = (isurface->second).begin(); inode != (isurface->second).end(); inode++) {
        
        int nodeId = inode->first;
        DENS_VEC nodalCoords = (inode->second).first;
        double nodalCharge = (inode->second).second;
        double nodalPotential = nodalChargePotential_[facesetName][nodeId];
        PerAtomQuantity<double> * atomicCoords = (atc_->interscale_manager()).per_atom_quantity("AtomicCoarseGrainingPositions");
        const DENS_MAT & myAtomicCoords(atomicCoords->quantity());
        for (int i = 0; i < nLocal; i++) {
          if (abs(qatom(i,0)) > 0) { 
            double distanceSq = 0.;
            double deltaX[3];
            for (int j = 0; j < nsd; j++) {
              deltaX[j] = myAtomicCoords(i,j) - nodalCoords(j);
              distanceSq += deltaX[j]*deltaX[j];
            }
            if (distanceSq < cutoffSq) { 
              // first apply pairwise coulombic interaction
              if (!useSlab_) { 
                double coulForce = qqrd2e*nodalCharge*qatom(i,0)/(distanceSq*sqrtf(distanceSq));
                for (int j = 0; j < nsd; j++)
                  //fatom[atomIdx][j] += deltaX[j]*coulForce;
                  _atomElectricalForce_(i,j) += deltaX[j]*coulForce;
              }
              
              // second correct for FE potential induced by BCs
              // determine shape function derivatives at atomic location
              // and construct sparse vectors to store derivative data
              
              vector<SparseVector<double> > derivativeVectors;
              derivativeVectors.reserve(nsd);
              const SPAR_MAT_VEC & shapeFunctionDerivatives((interscaleManager.vector_sparse_matrix("InterpolantGradient"))->quantity());
              for (int j = 0; j < nsd; j++) {
                DenseVector<INDEX> nodeIndices;
                DENS_VEC nodeValues;
                shapeFunctionDerivatives[j]->row(i,nodeValues,nodeIndices);
                derivativeVectors.push_back(dummy);
                for (int k = 0; k < nodeIndices.size(); k++)
                  derivativeVectors[j](nodeIndices(k)) = nodeValues(k);
              }
              
              // compute greens function from charge quadrature
              
              SparseVector<double> shortFePotential(nNodes);
              shortFePotential.add_scaled(greensFunctions_[nodeId],penalty*nodalPotential);
              
              // compute electric field induced by charge
              DENS_VEC efield(nsd);
              efield = 0.;
              for (int j = 0; j < nsd; j++)
                efield(j) = -.1*dot(derivativeVectors[j],shortFePotential);
              
              // apply correction in atomic forces
              //double c = qE2f*qatom[atomIdx];
              double c = qV2e*qatom(i,0);
              for (int j = 0; j < nsd; j++)
                if ((!useSlab_) || (j==nsd)) { 
                  //fatom[atomIdx][j] -= c*efield(j);
                  _atomElectricalForce_(i,j) -= c*efield(j);
                }
            }
          }
        }
      }
    }
  }
#endif

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModelElectrostaticMomentum
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ExtrinsicModelElectrostaticMomentum::ExtrinsicModelElectrostaticMomentum
                                     (ExtrinsicModelManager * modelManager,
                                      ExtrinsicModelType modelType,
                                      string matFileName) :
    ExtrinsicModelElectrostatic(modelManager,modelType,matFileName)
  {
     if (physicsModel_) delete physicsModel_; 
     if (modelType == ELECTROSTATIC) {
       physicsModel_ = new PhysicsModelElectrostatic(matFileName);
     }
     else {
       physicsModel_ = new PhysicsModelElectrostaticEquilibrium(matFileName);
     }
     // set up correct masks for coupling
     rhsMaskIntrinsic_(VELOCITY,SOURCE) = true;
     atc_->fieldMask_(VELOCITY,EXTRINSIC_SOURCE) = true;
  }
  
  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ExtrinsicModelElectrostaticMomentum::~ExtrinsicModelElectrostaticMomentum()
  {
    // do nothing
  }

  //--------------------------------------------------------
  // modify 
  //--------------------------------------------------------
  bool ExtrinsicModelElectrostaticMomentum::modify(int narg, char **arg)
  {
    bool match = false;

    if (!match)
      match = ExtrinsicModelElectrostatic::modify(narg,arg);

    return match;
  }

  //--------------------------------------------------------
  //  initialize
  //--------------------------------------------------------
  void ExtrinsicModelElectrostaticMomentum::initialize()
  {
    ExtrinsicModelElectrostatic::initialize();

    int nNodes = atc_->num_nodes();
    int nsd    = atc_->nsd();
    rhs_[VELOCITY].reset(nNodes,nsd);
  }

  //--------------------------------------------------------
  //  set coupling source terms
  //--------------------------------------------------------
  void ExtrinsicModelElectrostaticMomentum::set_sources(FIELDS & fields, FIELDS & sources)
  {
    // compute charge density
    if (modelType_ == ELECTROSTATIC_EQUILIBRIUM) { 
      DENS_MAN & n = atc_->field(ELECTRON_DENSITY);
      atc_->nodal_projection(ELECTRON_DENSITY,physicsModel_,n);
    }
    // else {
    //   FIELDS rhs;
    //   Array2D<bool> mask;
    //   mask(ELECTRON_DENSITY,SOURCE) = true;
    //   atc_->evaluate_rhs_integral(mask,fields,rhs,FULL_DOMAIN,physicsModel_);
    //   atc_->apply_inverse_mass_matrix(rhs[ELECTRON_DENSITY].quantity(),n.set_quantity(),ELECTRON_DENSITY);
    // }

    // compute source term with appropriate masking and physics model
    atc_->evaluate_rhs_integral(rhsMaskIntrinsic_, fields, sources,
                                atc_->source_integration(), physicsModel_); 
//(sources[VELOCITY].quantity()).print("V SRC");
  }

  //--------------------------------------------------------
  //  output
  //--------------------------------------------------------
  void ExtrinsicModelElectrostaticMomentum::output(OUTPUT_LIST & outputData)
  {
    ExtrinsicModelElectrostatic::output(outputData);
  }
};
