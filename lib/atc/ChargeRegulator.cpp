#include "ChargeRegulator.h"
#include "PoissonSolver.h"
#include "LammpsInterface.h"
#include "ATC_Coupling.h"
#include "ATC_Error.h"

#include "Function.h"
#include "PrescribedDataManager.h"

   

namespace ATC {

  //========================================================
  //  Class ChargeRegulator
  //========================================================

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ChargeRegulator::ChargeRegulator(ATC_Coupling * atc) :
    AtomicRegulator(atc)
  {
    // do nothing
  }
  //--------------------------------------------------------
  // Destructor
  //--------------------------------------------------------
  ChargeRegulator::~ChargeRegulator()
  {
    map<string,ChargeRegulatorMethod *>::iterator it;
    for (it = regulators_.begin(); it != regulators_.end(); it++) {
      if (it->second) delete it->second;
    }
  }

  
  //--------------------------------------------------------
  //  modify:
  //    parses and adjusts charge regulator state based on
  //    user input, in the style of LAMMPS user input
  //--------------------------------------------------------
  bool ChargeRegulator::modify(int narg, char **arg)
  {
    bool foundMatch = false;
    return foundMatch;
  }

  //--------------------------------------------------------
  // construct methods
  //--------------------------------------------------------
  void ChargeRegulator::construct_methods()
  {
    AtomicRegulator::construct_methods();

    if (atc_->reset_methods()) {
      // eliminate existing methods
      delete_method();
      // consruct new ones
      map<string, ChargeRegulatorParameters>::iterator itr;
      for (itr = parameters_.begin();
           itr != parameters_.end(); itr++) {
        string tag = itr->first;
        if (regulators_.find(tag) != regulators_.end()) delete regulators_[tag];
        ChargeRegulatorParameters & p = itr->second;
        LammpsInterface * lammpsInterface = LammpsInterface::instance();
        p.groupBit = lammpsInterface->group_bit(tag);
        if (! p.groupBit)
          throw ATC_Error("ChargeRegulator::initialize group not found");
        switch (p.method) {
        case NONE: {
          regulators_[tag] = new ChargeRegulatorMethod(this,p);
          break;
        }
        case FEEDBACK: {
          regulators_[tag] = new ChargeRegulatorMethodFeedback(this,p);
          break;
        }
        case IMAGE_CHARGE: {
          regulators_[tag] = new ChargeRegulatorMethodImageCharge(this,p);
          break;
        }
        case EFFECTIVE_CHARGE: {
          regulators_[tag] = new ChargeRegulatorMethodEffectiveCharge(this,p);
          break;
        }
        default: 
          throw ATC_Error("ChargeRegulator::construct_method unknown charge regulator type");
        }
      }
    }
  }

  //--------------------------------------------------------
  //  initialize:
  //--------------------------------------------------------
  void ChargeRegulator::initialize()
  {

    
    map<string, ChargeRegulatorMethod *>::iterator itr;
    for (itr = regulators_.begin();
         itr != regulators_.end(); itr++) { itr->second->initialize(); }

    atc_->set_boundary_integration_type(boundaryIntegrationType_); 
    AtomicRegulator::reset_nlocal();
    AtomicRegulator::delete_unused_data();
    needReset_ = false;
  }

  //--------------------------------------------------------
  //  apply pre force
  //--------------------------------------------------------
  void ChargeRegulator::apply_pre_force(double dt)
  {
    map<string, ChargeRegulatorMethod *>::iterator itr;
    for (itr = regulators_.begin();
         itr != regulators_.end(); itr++) { itr->second->apply_pre_force(dt);}
  }
  //--------------------------------------------------------
  //  apply post force
  //--------------------------------------------------------
  void ChargeRegulator::apply_post_force(double dt)
  {
    map<string, ChargeRegulatorMethod *>::iterator itr;
    for (itr = regulators_.begin();
         itr != regulators_.end(); itr++) { itr->second->apply_post_force(dt);}
  }
  //--------------------------------------------------------
  //  output
  //--------------------------------------------------------
  void ChargeRegulator::output(OUTPUT_LIST & outputData)
  {
    map<string, ChargeRegulatorMethod *>::iterator itr;
    for (itr = regulators_.begin();
         itr != regulators_.end(); itr++) { itr->second->output(outputData);}
  }

  //========================================================
  //  Class ChargeRegulatorMethod
  //========================================================
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and ChargeRegulator 
  //--------------------------------------------------------
  ChargeRegulatorMethod::ChargeRegulatorMethod
    (ChargeRegulator *chargeRegulator, 
     ChargeRegulator::ChargeRegulatorParameters & p)
      : RegulatorShapeFunction(chargeRegulator), 
        chargeRegulator_(chargeRegulator),
        lammpsInterface_(LammpsInterface::instance()),
        rC_(0), rCsq_(0),
        targetValue_(NULL), 
        targetPhi_(p.value), 
        surface_(p.faceset),
        atomGroupBit_(p.groupBit), 
        boundary_(false), 
        depth_(p.depth),
        surfaceType_(p.surfaceType),
        permittivity_(p.permittivity),
        initialized_(false)
  {
    const FE_Mesh * feMesh = atc_->fe_engine()->fe_mesh();
    feMesh->faceset_to_nodeset(surface_,nodes_);
    // assume flat  get normal and primary coord
    PAIR face = *(surface_.begin());
    normal_.reset(nsd_);
    feMesh->face_normal(face,0,normal_);
    DENS_MAT faceCoords;
    feMesh->face_coordinates(face,faceCoords);
    point_.reset(nsd_);
    for (int i=0; i < nsd_; i++) { point_(i) = faceCoords(i,0); }
#ifdef ATC_VERBOSE
    stringstream ss; ss << "point: (" << point_(0) << "," << point_(1) << "," << point_(2) << ") normal: (" << normal_(0) << "," << normal_(1) << "," << normal_(2) << ") depth: " << depth_; 
    lammpsInterface_->print_msg_once(ss.str());
#endif
    sum_.reset(nsd_);
  }

  //--------------------------------------------------------
  //  Initialize
  
  //--------------------------------------------------------


// nomenclature might be a bit backwark: control --> nodes that exert the control, & influence --> atoms that feel the influence
  void ChargeRegulatorMethod::initialize(void)
  {
    interscaleManager_ = &(atc_->interscale_manager());

    poissonSolver_ =chargeRegulator_->poisson_solver();
    if (! poissonSolver_) throw ATC_Error("need a poisson solver to initialize charge regulator");

    // atomic vectors
    // nodal information
    nNodes_ = atc_->num_nodes();
    // constants
    rC_ = lammpsInterface_->pair_cutoff();
    rCsq_ = rC_*rC_;
    qV2e_ = lammpsInterface_->qv2e();
    qqrd2e_ = lammpsInterface_->qqrd2e();

    // note derived method set intialized to true
  }

  
  int ChargeRegulatorMethod::nlocal() { return atc_->nlocal(); }

  void ChargeRegulatorMethod::set_greens_functions(void)
  {
    // set up Green's function per node
    for (int i = 0; i < nNodes_; i++) {
      set<int> localNodes;
      for (int j = 0; j < nNodes_; j++)
        localNodes.insert(j);
      // call Poisson solver to get Green's function for node i
      DENS_VEC globalGreensFunction;
      poissonSolver_->greens_function(i,globalGreensFunction);
      // store green's functions as sparse vectors only on local nodes
      set<int>::const_iterator thisNode;
      SparseVector<double> sparseGreensFunction(nNodes_);
      for (thisNode = localNodes.begin(); thisNode != localNodes.end(); thisNode++)
        sparseGreensFunction(*thisNode) = globalGreensFunction(*thisNode);
      greensFunctions_.push_back(sparseGreensFunction);
    }
  }
  //--------------------------------------------------------
  //  output
  //--------------------------------------------------------
  void ChargeRegulatorMethod::output(OUTPUT_LIST & outputData)
  {
    //vector<double> localSum(sum_.size());
    //lammpsInteface_->allsum(localSum.pointer,sum_.pointer,sum_.size());
    DENS_VEC localSum(sum_.size());
    lammpsInterface_->allsum(localSum.ptr(),sum_.ptr(),sum_.size());
    for (int i = 0; i < sum_.size(); i++) {
      string name = "charge_regulator_influence_"+to_string(i);
//    atc_->fe_engine()->add_global(name,sum_[i]);
    }
  }

  //========================================================
  //  Class ChargeRegulatorMethodFeedback
  //========================================================
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ChargeRegulatorMethodFeedback::ChargeRegulatorMethodFeedback
    (ChargeRegulator *chargeRegulator, 
     ChargeRegulator::ChargeRegulatorParameters & p)
      : ChargeRegulatorMethod (chargeRegulator, p),
      controlNodes_(nodes_),
      influenceGroupBit_(p.groupBit)
  {
    nControlNodes_   = controlNodes_.size(); 
    sum_.resize(1); 
  }
  //--------------------------------------------------------
  //  Initialize
  //--------------------------------------------------------
  void ChargeRegulatorMethodFeedback::initialize(void)
  {
    ChargeRegulatorMethod::initialize();
    if (surfaceType_ != ChargeRegulator::CONDUCTOR) 
      throw ATC_Error("currently charge feedback can only mimic a conductor");
    set_influence();  
    set_influence_matrix(); 
    initialized_ = true;
  }
  //--------------------------------------------------------
  //  find measurement atoms and nodes 
  //--------------------------------------------------------
  void ChargeRegulatorMethodFeedback::set_influence(void)
  {

    // get nodes that overlap influence atoms & compact list of influence atoms
    boundary_ = 
      atc_->nodal_influence(influenceGroupBit_,influenceNodes_,influenceAtoms_);
    nInfluenceAtoms_ = influenceAtoms_.size(); // local
    nInfluenceNodes_ = influenceNodes_.size(); // global
    stringstream ss; ss << "control nodes: " << nControlNodes_ << " influence nodes: " << nInfluenceNodes_ << " local influence atoms: " << nInfluenceAtoms_ ;
    lammpsInterface_->print_msg(ss.str());
    if (nInfluenceNodes_ == 0) throw ATC_Error("no influence nodes");

    const Array<int> & map = (boundary_) ? atc_->ghost_to_atom_map() : atc_->internal_to_atom_map(); 
    for (set<int>::const_iterator itr = influenceAtoms_.begin(); itr != influenceAtoms_.end(); itr++) {
      influenceAtomsIds_.insert(map(*itr));
    }
  }
  //--------------------------------------------------------
  //  constuct a Green's submatrix 
  //--------------------------------------------------------
  void ChargeRegulatorMethodFeedback::set_influence_matrix(void)
  {
    // construct control-influence matrix bar{G}^-1: ds{p} = G{p,m}^-1 dphi{m}


//
    if (nInfluenceNodes_ < nControlNodes_) throw ATC_Error(" least square not implmented ");
    if (nInfluenceNodes_ > nControlNodes_) throw ATC_Error(" solve not possible ");
    DENS_MAT G(nInfluenceNodes_,nControlNodes_); 
    DENS_VEC G_I;
    set<int>::const_iterator itr,itr2,itr3;
    const Array<int> & nmap = atc_->fe_engine()->fe_mesh()->global_to_unique_map();
    int i = 0;
    for (itr = influenceNodes_.begin(); itr != influenceNodes_.end(); itr++) {
      poissonSolver_->greens_function(*itr, G_I);
      int j = 0;
      for (itr2 = controlNodes_.begin(); itr2 != controlNodes_.end(); itr2++) {
        int jnode = nmap(*itr2);
        G(i,j++) = G_I(jnode);  
      }
      i++;
    }
    invG_ = inv(G);

    // construct the prolong-restrict projector N N^T for influence nodes only

    InterscaleManager & interscaleManager(atc_->interscale_manager());
    const SPAR_MAT & N_Ia = (boundary_) ?
      (interscaleManager.per_atom_sparse_matrix("InterpolantGhost"))->quantity():
      (interscaleManager.per_atom_sparse_matrix("Interpolant"))->quantity();
    NT_.reset(nInfluenceAtoms_,nInfluenceNodes_);
    DENS_MAT NNT(nInfluenceNodes_,nInfluenceNodes_);
    int k = 0;
    for (itr3 = influenceAtoms_.begin(); itr3 != influenceAtoms_.end(); itr3++) {
      int katom = *itr3;
      int i = 0;
      for (itr = influenceNodes_.begin(); itr != influenceNodes_.end(); itr++) {
        int Inode = *itr;
        int j = 0;
        NT_(k,i) = N_Ia(katom,Inode);
        for (itr2 = influenceNodes_.begin(); itr2 != influenceNodes_.end(); itr2++) {
          int Jnode = *itr2;
          NNT(i,j++) += N_Ia(katom,Inode)*N_Ia(katom,Jnode);
        }
        i++;
      }
      k++;
    }
    // swap contributions across processors
    DENS_MAT localNNT = NNT;
    int count = NNT.nRows()*NNT.nCols(); 
    lammpsInterface_->allsum(localNNT.ptr(),NNT.ptr(),count);
    invNNT_ = inv(NNT);
  
    // total influence matrix
    if (nInfluenceAtoms_ > 0) { NTinvNNTinvG_ = NT_*invNNT_*invG_; }

  }
  
  //--------------------------------------------------------
  // change potential/charge pre-force calculation
  //--------------------------------------------------------
  void ChargeRegulatorMethodFeedback::apply_pre_force(double dt)
  {

    sum_ = 0; 
    if (nInfluenceAtoms_ == 0) return; // nothing to do
    apply_feedback_charges();
  }
  //--------------------------------------------------------
  // apply feedback charges to atoms                        
  //--------------------------------------------------------
  void ChargeRegulatorMethodFeedback::apply_feedback_charges()
  {
    double * q = lammpsInterface_->atom_charge();
    // calculate error in potential on the control nodes
    
    const DENS_MAT & phiField = (atc_->field(ELECTRIC_POTENTIAL)).quantity();
    DENS_MAT dphi(nControlNodes_,1);
    int i = 0;
    set<int>::const_iterator itr;
    for (itr = controlNodes_.begin(); itr != controlNodes_.end(); itr++) {
      dphi(i++,0) = targetPhi_ - phiField(*itr,0);
    }

    // construct the atomic charges consistent with the correction
    DENS_MAT dq = NTinvNNTinvG_*dphi;
    i = 0;
    for (itr = influenceAtomsIds_.begin(); itr != influenceAtomsIds_.end(); itr++) {
      sum_(0) += dq(i,0); 
      q[*itr] += dq(i++,0); 
    }
    
    (interscaleManager_->fundamental_atom_quantity(LammpsInterface::ATOM_CHARGE))->force_reset();
    (interscaleManager_->fundamental_atom_quantity(LammpsInterface::ATOM_CHARGE, GHOST))->force_reset();
  }

  //========================================================
  //  Class ChargeRegulatorMethodImageCharge
  //========================================================
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ChargeRegulatorMethodImageCharge::ChargeRegulatorMethodImageCharge
    (ChargeRegulator *chargeRegulator, 
     ChargeRegulator::ChargeRegulatorParameters & p)
      : ChargeRegulatorMethod (chargeRegulator, p),
      imageNodes_(nodes_)
  {
  }
  //--------------------------------------------------------
  //  Initialize
  //--------------------------------------------------------
  void ChargeRegulatorMethodImageCharge::initialize(void)
  {
    ChargeRegulatorMethod::initialize();
    if (surfaceType_ != ChargeRegulator::DIELECTRIC) throw ATC_Error("currently image charge can only mimic a dielectric");
    double eps1 = permittivity_;// dielectric
    double eps2 = lammpsInterface_->dielectric();// ambient
    permittivityRatio_ = (eps2-eps1)/(eps2+eps1);
#ifdef ATC_VERBOSE
    stringstream ss; ss << "permittivity ratio: " << permittivityRatio_;
    lammpsInterface_->print_msg_once(ss.str());
#endif
    set_greens_functions();

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
    initialized_ = true;
  }

  //--------------------------------------------------------
  // change potential/charge post-force calculation
  //--------------------------------------------------------
  void ChargeRegulatorMethodImageCharge::apply_post_force(double dt)
  {
    sum_ = 0;
    apply_local_forces();
    
    //correct_forces();
  }
 
  //--------------------------------------------------------
  //  apply local coulomb forces 
  //  -- due to image charges
  //--------------------------------------------------------
  void ChargeRegulatorMethodImageCharge::apply_local_forces()
  {

    int inum = lammpsInterface_->neighbor_list_inum();
    int * ilist = lammpsInterface_->neighbor_list_ilist();
    int * numneigh = lammpsInterface_->neighbor_list_numneigh();
    int ** firstneigh = lammpsInterface_->neighbor_list_firstneigh();

    const int *mask = lammpsInterface_->atom_mask();
///..............................................
    double ** x = lammpsInterface_->xatom(); 
    double ** f = lammpsInterface_->fatom(); 
    double *  q = lammpsInterface_->atom_charge();

    // loop over neighbor list
    for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii];
      double qi = q[i];
      if ((mask[i] & atomGroupBit_) && qi != 0.) {
        double* fi = f[i];
        DENS_VEC xi(x[i],nsd_); 
        // distance to surface
        double dn = reflect(xi);
        // all ions near the interface/wall
        // (a) self image
        if (dn < rC_) { // close enough to wall to have explicit image charges
          double factor_coul = 1; 
          double dx = 2.*dn; // distance to image charge
          double fn = factor_coul*qi*qi*permittivityRatio_/dx;
          fi[0] += fn*normal_[0];
          fi[1] += fn*normal_[1];
          fi[2] += fn*normal_[2];
          sum_ += fn*normal_;
        // (b) neighbor images
        int * jlist = firstneigh[i];
        int jnum = numneigh[i];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          // this changes j
          double factor_coul = lammpsInterface_->coulomb_factor(j);
          double qj = q[j];
          if (qj != 0.) { // all charged neighbors
            DENS_VEC xj(x[j],nsd_);
            dn = reflect(xj);
            DENS_VEC dx = xi-xj;
            double r2 = dx.norm_sq();
            // neighbor image j' inside cutoff from i 
            if (r2 < rCsq_) { 
              double fm = factor_coul*qi*qj*permittivityRatio_/r2;
              fi[0] += fm*dx(0);
              fi[1] += fm*dx(1);
              fi[2] += fm*dx(2);
              sum_ += fm*dx;
              }
            }
          }
        } // end i < rC if
      }
    }
    // update managed data
    (interscaleManager_->fundamental_atom_quantity(LammpsInterface::ATOM_FORCE))->force_reset();
  }

  //--------------------------------------------------------
  // correct charge densities
  //  - to reflect image charges 
  //--------------------------------------------------------
  void ChargeRegulatorMethodImageCharge::correct_charge_densities()
  {
  }


  //--------------------------------------------------------
  //  correct_forces
  //  - due to image charge density used in short-range solution
  //--------------------------------------------------------
  void ChargeRegulatorMethodImageCharge::correct_forces()
  {
  }


  //========================================================
  //  Class ChargeRegulatorMethodEffectiveCharge
  //========================================================
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  
  ChargeRegulatorMethodEffectiveCharge::ChargeRegulatorMethodEffectiveCharge( 
    ChargeRegulator *chargeRegulator, 
    ChargeRegulator::ChargeRegulatorParameters & p)
      : ChargeRegulatorMethod (chargeRegulator, p),
      chargeDensity_(p.value),
      useSlab_(false) 
  {
  }
  //--------------------------------------------------------
  //  add_charged_surface
  //--------------------------------------------------------
  void ChargeRegulatorMethodEffectiveCharge::initialize( ) 
  {
    ChargeRegulatorMethod::initialize();
    boundary_ = atc_->is_ghost_group(atomGroupBit_);
    // set face sources to all point at unit function for use in integration
    SURFACE_SOURCE faceSources;
    map<PAIR, Array<XT_Function*> > & fs(faceSources[ELECTRIC_POTENTIAL]);
    XT_Function * f = XT_Function_Mgr::instance()->constant_function(1.);
    set< PAIR >::const_iterator fsItr;
    for (fsItr = surface_.begin(); fsItr != surface_.end(); fsItr++) {
      Array < XT_Function * > & dof  = fs[*fsItr];
      dof.reset(1);
      dof(0) = f;
    }

    // computed integrals of nodal shape functions on face
    FIELDS nodalFaceWeights;
    Array<bool> fieldMask(NUM_FIELDS); fieldMask(ELECTRIC_POTENTIAL) = true;
    (atc_->fe_engine())->compute_fluxes(fieldMask,0.,faceSources,nodalFaceWeights);
    const DENS_MAT & w = (nodalFaceWeights[ELECTRIC_POTENTIAL].quantity());

    // Get coordinates of each node in face set
    for (set<int>::const_iterator n =nodes_.begin(); n != nodes_.end(); n++) {
      DENS_VEC x = atc_->fe_engine()->fe_mesh()->nodal_coordinates(*n);
      // compute effective charge at each node I
      // multiply charge density by integral of N_I over face
      double v = w(*n,0)*chargeDensity_;
      pair<DENS_VEC,double> p(x,v);
      nodeXFMap_[*n] = p;
    }

    // set up data structure holding charged faceset information
    FIELDS sources;
    double k = lammpsInterface_->coulomb_constant();
    string fname = "radial_power";
    double xtArgs[8];
    xtArgs[0] = 0; xtArgs[1] = 0; xtArgs[2] = 0;
    xtArgs[3] = 1; xtArgs[4] = 1; xtArgs[5] = 1;
    xtArgs[6] = k*chargeDensity_;
    xtArgs[7] = -1.;
    const DENS_MAT & s(sources[ELECTRIC_POTENTIAL].quantity());
    NODE_TO_XF_MAP::iterator XFitr;
    for (XFitr = nodeXFMap_.begin(); XFitr != nodeXFMap_.end(); XFitr++) {
      // evaluate voltage at each node I
      // set up X_T function for integration: k*chargeDensity_/||x_I - x_s||
      // integral is approximated in two parts:
      // 1) near part with all faces within r < rcrit evaluated as 2 * pi * rcrit * k sigma A/A0, A is area of this region and A0 = pi * rcrit^2, so 2 k sigma A / rcrit
      // 2) far part evaluated using Gaussian quadrature on faceset
      DENS_VEC x((XFitr->second).first);
      xtArgs[0] = x(0); xtArgs[1] = x(1); xtArgs[2] = x(2);
      f = XT_Function_Mgr::instance()->function(fname,8,xtArgs);
      for (fsItr = surface_.begin(); fsItr != surface_.end(); fsItr++) {
        fs[*fsItr] = f;
      }

      // perform integration to get quantities at nodes on facesets
      // V_J' = int_S N_J k*sigma/|x_I - x_s| dS
      (atc_->fe_engine())->compute_fluxes(fieldMask,0.,faceSources,sources);

      // sum over all nodes in faceset to get total potential:
      // V_I = sum_J VJ'
      int node = XFitr->first;
      nodalChargePotential_[node] = s(node,0);
      double totalPotential = 0.;
      for (set<int>::const_iterator n =nodes_.begin(); n != nodes_.end(); n++) {
        totalPotential += s(*n,0); }

      // assign an XT function per each node and
      // then call the prescribed data manager and fix each node individually.
      f = XT_Function_Mgr::instance()->constant_function(totalPotential);
      (atc_->prescribed_data_manager())->fix_field(node,ELECTRIC_POTENTIAL,0,f);
    }
    initialized_ = true;
  }

  //--------------------------------------------------------
  //  add effective forces post LAMMPS force call
  //--------------------------------------------------------
  void ChargeRegulatorMethodEffectiveCharge::apply_post_force(double dt) 
  {
    apply_local_forces();
  }

  //--------------------------------------------------------
  //  apply_charged_surfaces
  //--------------------------------------------------------
  void ChargeRegulatorMethodEffectiveCharge::apply_local_forces()
  {
    double * q = lammpsInterface_->atom_charge();
    _atomElectricalForce_.resize(nlocal(),nsd_);

    double penalty = poissonSolver_->penalty_coefficient();
    if (penalty <= 0.0) throw ATC_Error("ExtrinsicModelElectrostatic::apply_charged_surfaces expecting non zero penalty");

    double dx[3];
    const DENS_MAT & xa((interscaleManager_->per_atom_quantity("AtomicCoarseGrainingPositions"))->quantity());

// WORKSPACE - most are static
    SparseVector<double> dv(nNodes_); 
    vector<SparseVector<double> > derivativeVectors;
    derivativeVectors.reserve(nsd_);
    const SPAR_MAT_VEC & shapeFunctionDerivatives((interscaleManager_->vector_sparse_matrix("InterpolateGradient"))->quantity());

    DenseVector<INDEX> nodeIndices;
    DENS_VEC nodeValues;

    NODE_TO_XF_MAP::const_iterator inode;
    for (inode = nodeXFMap_.begin(); inode != nodeXFMap_.end(); inode++) {
      
      int node = inode->first;
      DENS_VEC xI = (inode->second).first;
      double qI = (inode->second).second;
      double phiI = nodalChargePotential_[node];
      for (int i = 0; i < nlocal(); i++) {
        int atom = (atc_->internal_to_atom_map())(i);
        double qa = q[atom];
        if (qa != 0) { 
          double dxSq = 0.;
          for (int j = 0; j < nsd_; j++) {
            dx[j] = xa(i,j) - xI(j);
            dxSq += dx[j]*dx[j];
          }
          if (dxSq < rCsq_) { 
            // first apply pairwise coulombic interaction
            if (!useSlab_) { 
              double coulForce = qqrd2e_*qI*qa/(dxSq*sqrtf(dxSq));
              for (int j = 0; j < nsd_; j++) {
                _atomElectricalForce_(i,j) += dx[j]*coulForce; }
            }
            
            // second correct for FE potential induced by BCs
            // determine shape function derivatives at atomic location
            // and construct sparse vectors to store derivative data
            
            
            for (int j = 0; j < nsd_; j++) {
              shapeFunctionDerivatives[j]->row(i,nodeValues,nodeIndices);
              derivativeVectors.push_back(dv);
              for (int k = 0; k < nodeIndices.size(); k++) {
                derivativeVectors[j](nodeIndices(k)) = nodeValues(k); }
            }
              
            // compute greens function from charge quadrature
            
            SparseVector<double> shortFePotential(nNodes_); 
            shortFePotential.add_scaled(greensFunctions_[node],penalty*phiI);
              
            // compute electric field induced by charge
            DENS_VEC efield(nsd_);
            for (int j = 0; j < nsd_; j++) {
              efield(j) = -.1*dot(derivativeVectors[j],shortFePotential); }
              
            // apply correction in atomic forces
            double c = qV2e_*qa;
            for (int j = 0; j < nsd_; j++) {
              if ((!useSlab_) || (j==nsd_)) { 
                _atomElectricalForce_(i,j) -= c*efield(j);
              }
            }
          }
        }
      }
    }
    
  }

}; // end namespace
