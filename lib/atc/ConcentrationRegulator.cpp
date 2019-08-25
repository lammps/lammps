#include "ConcentrationRegulator.h"
#include "LammpsInterface.h"
#include "ATC_Coupling.h"
#include "ATC_Error.h"

using ATC_Utility::to_string;
using ATC_Utility::rnd;
using std::map;
using std::string;
using std::pair;
using std::min;
using std::max;

namespace ATC {

const double kMinScale_ = 10000.;

  //========================================================
  //  Class ConcentrationRegulator
  //========================================================

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ConcentrationRegulator::ConcentrationRegulator(ATC_Coupling * atc) :
    AtomicRegulator(atc)
  {
    // do nothing
  }
  //--------------------------------------------------------
  // Destructor
  //--------------------------------------------------------
  ConcentrationRegulator::~ConcentrationRegulator()
  {
    if (regulators_.size()) {
      map<string,ConcentrationRegulatorMethod *>::iterator it;
      for (it = regulators_.begin(); it != regulators_.end(); ++it) {
        delete it->second;
      }
      regulators_.clear();
    }
    if (parameters_.size()) parameters_.clear();
  }
  
  //--------------------------------------------------------
  //  modify:
  //    parses and adjusts charge regulator state based on
  //    user input, in the style of LAMMPS user input
  //--------------------------------------------------------
  bool ConcentrationRegulator::modify(int narg, char **arg)
  {
    bool foundMatch = false;
    return foundMatch;
  }

  //--------------------------------------------------------
  //  construct_methods:
  //--------------------------------------------------------
  void ConcentrationRegulator::construct_methods()
  {
    AtomicRegulator::construct_methods();

    if (atc_->reset_methods()) {
      // eliminate existing methods
      delete_method();
      // consruct new ones
      map<string, ConcentrationRegulatorParameters>::iterator itr;
      for (itr = parameters_.begin();
           itr != parameters_.end(); itr++) { 
        string tag = itr->first;
        if (regulators_.find(tag) != regulators_.end()) delete regulators_[tag];
        ConcentrationRegulatorParameters & p = itr->second;
        switch (p.method) {
        case NONE: {
          regulators_[tag] = new ConcentrationRegulatorMethod(this);
          break;
        }
        case TRANSITION: {
          p.type = atc_->tag_to_type(tag);
          p.groupbit = LammpsInterface::instance()->type_to_groupbit(p.type);
          p.transitionType = atc_->tag_to_type(p.transitionTag);
          regulators_[tag] = new ConcentrationRegulatorMethodTransition(this,p);
          break;
        }
        default: 
          throw ATC_Error("ConcentrationRegulator::initialize unknown concentration regulator type");
        }
      }
    }
  }
  //--------------------------------------------------------
  //  initialize:
  //--------------------------------------------------------
  void ConcentrationRegulator::initialize()
  {
    
    map<string, ConcentrationRegulatorMethod *>::iterator itr;
    for (itr = regulators_.begin();
         itr != regulators_.end(); itr++) { itr->second->initialize(); }

    atc_->set_boundary_integration_type(boundaryIntegrationType_); 
    AtomicRegulator::reset_nlocal();
    AtomicRegulator::delete_unused_data();
    needReset_ = false;
  }

  //--------------------------------------------------------
  //  pre_exchange
  //--------------------------------------------------------
  void ConcentrationRegulator::pre_exchange()
  {
    map<string, ConcentrationRegulatorMethod *>::iterator itr;
    for (itr = regulators_.begin();
         itr != regulators_.end(); itr++) { itr->second->pre_exchange();}
  }

  //--------------------------------------------------------
  //  pre_force
  //--------------------------------------------------------
  void ConcentrationRegulator::pre_force()
  {
    map<string, ConcentrationRegulatorMethod *>::iterator itr;
    for (itr = regulators_.begin();
         itr != regulators_.end(); itr++) { itr->second->pre_force();}
  }

  //--------------------------------------------------------
  //  finish
  //--------------------------------------------------------
  void ConcentrationRegulator::finish()
  {
    map<string, ConcentrationRegulatorMethod *>::iterator itr;
    for (itr = regulators_.begin();
         itr != regulators_.end(); itr++) { itr->second->finish();}
  }

  //--------------------------------------------------------
  //  output
  //--------------------------------------------------------
  void ConcentrationRegulator::output(OUTPUT_LIST & outputData) const
  {
    map<string, ConcentrationRegulatorMethod *>::const_iterator itr;
    for (itr = regulators_.begin();
         itr != regulators_.end(); itr++) { itr->second->output(outputData);}
  }

  //--------------------------------------------------------
  //  compute vector
  //--------------------------------------------------------
  double ConcentrationRegulator::compute_vector(int n) const
  {
    int s = regulators_.size();
    if (s == 0) return 0;
    int m = n / s;
    n     = n % s;
    int c = 0;

    map<string, ConcentrationRegulatorMethod *>::const_iterator itr;
    for (itr = regulators_.begin();
         itr != regulators_.end(); itr++) { 
      if (c++ == n) { return itr->second->compute_vector(m); }
    }
    return 0.;

  }
  //--------------------------------------------------------
  //  size vector
  //--------------------------------------------------------
  int ConcentrationRegulator::size_vector(int i) const
  {
    int n = (regulators_.size())*5; 
    if (n==0) n = 20; 
    return n; 
  }

  //========================================================
  //  Class ConcentrationRegulatorMethodTransition
  //========================================================
 
  //--------------------------------------------------------
  //  Constructor
  //         Grab references to ATC and ConcentrationRegulator 
  //--------------------------------------------------------
  ConcentrationRegulatorMethodTransition::ConcentrationRegulatorMethodTransition
    (ConcentrationRegulator *concReg, 
     ConcentrationRegulator::ConcentrationRegulatorParameters & p)
      : ConcentrationRegulatorMethod(concReg), 
        concentrationRegulator_(concReg),
        interscaleManager_(NULL),
        lammpsInterface_(LammpsInterface::instance()),
        list_(NULL),
        targetConcentration_(p.value), 
        targetCount_(0), 
        elemset_(p.elemset),
        p_(NULL),
        randomNumberGenerator_(NULL), 
        q0_(0),
        controlType_(p.type), 
        controlIndex_(0),
        transitionType_(p.transitionType),
        transitionInterval_(p.transitionInterval),
        transitionCounter_(0),
        nInTransition_(0),
        transitionFactor_(0),
        controlMask_(p.groupbit),
        frequency_(p.frequency),
        maxEnergy_(p.maxEnergy),
        maxExchanges_(p.maxExchanges),
        maxAttempts_(p.maxAttempts),
        nexchanges_(0),
        initialized_(false),
        _rngUniformCounter_(0),
        _rngNormalCounter_(0)
  {
    controlIndex_ = atc_->type_index(controlType_);
    LammpsInterface * li = LammpsInterface::instance();
    q0_ = li->type_to_charge(controlType_);
    double kB = (li->boltz())/(li->mvv2e()); // E/T*m*v^2/E = m v^2/T
    double m = li->atom_mass(controlType_);
    sigma_ = sqrt(kB/m); // v / sqrt(T)
    randomNumberGenerator_ = li->random_number_generator();
  }

  //--------------------------------------------------------
  //  Initialize
  //--------------------------------------------------------
  void ConcentrationRegulatorMethodTransition::initialize(void)
  {
#ifdef ATC_VERBOSE
    lammpsInterface_->print_msg_once(
      "\ncontrol type: "+to_string(controlType_)+
      "\ntransistion type:"+to_string(transitionType_)+
      "\ncontrol mask:"+to_string(controlMask_)+
      "\nfrequency:"+to_string(frequency_)+
      "\nmax exchanges:"+to_string(maxExchanges_)+
      "\nmax attempts:"+to_string(maxAttempts_)+
      "\nmax energy:"+to_string(maxEnergy_)
    );
#endif
    interscaleManager_ = &(atc_->interscale_manager());

    PerAtomQuantity<int> * a2el = atc_->atom_to_element_map();
    list_ = new AtomInElementSet(atc_,a2el,elemset_,controlType_);
    
    nNodes_ = atc_->num_nodes();
    DENS_MAT conc(nNodes_,1); conc = targetConcentration_;
    DENS_VEC integral = atc_->fe_engine()->integrate(conc,elemset_);
    targetCount_ = rnd(integral(0)) ;

    volumes_.resize(elemset_.size());
    ESET::const_iterator itr;
    int i = 0;
    DENS_MAT c(nNodes_,1); c = 1; 
    V_ = 0.;
    for (itr = elemset_.begin(); itr != elemset_.end(); itr++, i++) {
      ESET e; e.insert(*itr);
      DENS_VEC v = atc_->fe_engine()->integrate(c,e);
      volumes_(i) = v(0);
      V_ += v(0);
    }
    volumes_ *= 1./V_;
    for (int i = 1; i < volumes_.size(); i++) {
      volumes_(i) += volumes_(i-1);
    } 

    // record orginal energetic properties
    int ntypes = lammpsInterface_->ntypes();
    epsilon0_.reset(ntypes);
    p_ = lammpsInterface_->potential();
    lammpsInterface_->epsilons(controlType_,p_,epsilon0_.ptr());
    
#ifdef ATC_VERBOSE
    string msg = "type "+to_string(controlType_)+" target count " + to_string(targetCount_);
    msg += " volume " + to_string(V_);
    msg += " current count " + to_string(count());
    ATC::LammpsInterface::instance()->print_msg_once(msg);
    msg = "WARNING: ensure neighboring happens at least every "+to_string(frequency_);
    ATC::LammpsInterface::instance()->print_msg_once(msg);
#endif
  }
  double ConcentrationRegulatorMethodTransition::uniform() const {
    _rngUniformCounter_++;
    return lammpsInterface_->random_uniform(randomNumberGenerator_); 
  }
  double ConcentrationRegulatorMethodTransition::normal() const {
    _rngNormalCounter_++;
    return lammpsInterface_->random_normal(randomNumberGenerator_); 
  }
  //--------------------------------------------------------
  //  pre exchange
  //--------------------------------------------------------
  void ConcentrationRegulatorMethodTransition::pre_exchange(void)
  {
    // return if should not be called on this timestep
    if ( ! lammpsInterface_->now(frequency_)) return;
    nexchanges_ = excess();
    int  n = abs(nexchanges_); 
    bool success = false;
    if      (nexchanges_ > 0) { success = delete_atoms(n); }
    else if (nexchanges_ < 0) { success = insert_atoms(n); }
    else return;
    if (!success) throw ATC_Error("insertions/deletions did not succeed");

    if (nexchanges_ !=0)  {
      nInTransition_ = -nexchanges_;
      lammpsInterface_->reset_ghosts(-nexchanges_);
      atc_->reset_atoms();
    }
    transitionCounter_=0;
    transition();
  } 
  //--------------------------------------------------------
  //  pre force
  //--------------------------------------------------------
  void ConcentrationRegulatorMethodTransition::pre_force(void)
  {
    transition();
  }
  //--------------------------------------------------------
  //  accept
  //--------------------------------------------------------
  bool ConcentrationRegulatorMethodTransition::accept(double energy, double T) const 
  { 
#ifdef ATC_VERBOSE2
    if (energy < maxEnergy_) lammpsInterface_->print_msg(" energy "+to_string(energy)+" "+to_string(rngCounter_));
#endif
    return (energy < maxEnergy_);
  } 
  //--------------------------------------------------------
  //  energy
  //--------------------------------------------------------
  double ConcentrationRegulatorMethodTransition::energy(int id) const 
  {
    double e = lammpsInterface_->shortrange_energy(id,maxEnergy_);
#ifdef ATC_VERBOSE
{    
    int * tag = lammpsInterface_->atom_tag();
    lammpsInterface_->print_msg(to_string(controlType_)+" deletion energy "+to_string(e)+" id "+to_string(tag[id])+" "+to_string(_rngUniformCounter_)+":"+to_string(_rngNormalCounter_));
}
#endif
    return e;
  }
  double ConcentrationRegulatorMethodTransition::energy(double * x) const 
  {
    double e = lammpsInterface_->shortrange_energy(x,controlType_,maxEnergy_);
#ifdef ATC_VERBOSE
{
    lammpsInterface_->print_msg(to_string(controlType_)+" insertion energy "+to_string(e)+" x "+to_string(x[0])+","+to_string(x[1])+","+to_string(x[2])+" "+to_string(_rngUniformCounter_)+":"+to_string(_rngNormalCounter_));
}
#endif
    return e;
  }
  //--------------------------------------------------------
  //  excess
  //--------------------------------------------------------
  int ConcentrationRegulatorMethodTransition::excess(void) const
  {
     int nexcess = count()-targetCount_;
     nexcess = max(min(nexcess,maxExchanges_),-maxExchanges_);
     return nexcess;
  }
  //--------------------------------------------------------
  //  count
  //--------------------------------------------------------
  int ConcentrationRegulatorMethodTransition::count(void) const
  {
    // integrate concentration over region
    const DENS_MAT & c = (atc_->field(SPECIES_CONCENTRATION)).quantity();
    DENS_VEC integral = atc_->fe_engine()->integrate(c,elemset_);
    return rnd(integral(controlIndex_)) ;
  }
  //--------------------------------------------------------
  //  delete atoms
  //--------------------------------------------------------
  bool ConcentrationRegulatorMethodTransition::delete_atoms(int n)
  {
    ID_PAIR idPair;
    
    deletionIds_.clear(); 
    int deletions = 0;
    int attempts = 0;
    while(deletions < n && attempts < maxAttempts_){
      if(accept(deletion_id(idPair))) {
        deletionIds_.push_back(idPair);
        deletions += 1;
      }
      deletions = lammpsInterface_->int_allmax(deletions);
      attempts++;
    }
    ID_LIST::iterator itr;
    for (itr = deletionIds_.begin(); itr != deletionIds_.end(); itr++) {
      lammpsInterface_->delete_atom(itr->second);
    }
#ifdef ATC_VERBOSE
      string c = to_string(controlType_);
      lammpsInterface_->all_print(attempts, c+"-attempts  ");
      lammpsInterface_->all_print(deletions,c+" deletions ");
      lammpsInterface_->all_print(_rngUniformCounter_,c+" RNG-uniform  ");
      lammpsInterface_->all_print(_rngNormalCounter_,c+" RNG-normal  ");
//    lammpsInterface_->all_print(uniform()," RANDOM ");
#endif
    return (n == deletions); // success
  }
  //--------------------------------------------------------
  //  pick id
  //--------------------------------------------------------
  double ConcentrationRegulatorMethodTransition::deletion_id(ID_PAIR & id) const
  {
    if (atc_->parallel_consistency()) return deletion_id_consistent(id); 
    else                              return deletion_id_free(id);
  }
  double ConcentrationRegulatorMethodTransition::deletion_id_consistent(ID_PAIR & id) const
  {
    id.first  = -1;
    id.second = -1;
    int ntotal = lammpsInterface_->natoms();
    double r = uniform();
    r *= ntotal;
    const ID_LIST & list = list_->quantity();
    ID_LIST::const_iterator itr;
    int i=0, idx = -1;
    double min = ntotal;
    int * tag = lammpsInterface_->atom_tag();
    for (itr = list.begin(); itr != list.end(); itr++) {
      int atag = tag[itr->second]; 
      double d = abs(atag-r);
      if (d < min) {
        min = d;
        idx = i;
      } 
      i++;
    }
    int imin = kMinScale_*min;
    if(imin == lammpsInterface_->int_allmin(imin)) {
      if (idx < 0) throw ATC_Error("deletion_id failed to find a suitable atom");
      id = list_->item(idx);
      // avoid repeats
      ID_LIST & l = list_->set_quantity();
      l.erase(l.begin()+idx);
      return energy(id.second);
    }
    else {
      return maxEnergy_;
    }
  }
  double ConcentrationRegulatorMethodTransition::deletion_id_free(ID_PAIR & id) const
  {
    id.first  = -1;
    id.second = -1;
    int n = list_->size();
    double nrank = lammpsInterface_->int_scansum(n);
    int   ntotal = lammpsInterface_->int_allsum(n);
    if (ntotal == 0) throw ATC_Error("control type "+to_string(controlType_)+" is depleted");
    double r = uniform();
    r *= ntotal;
    if ( (r >= nrank-n) && (r < nrank)) { // pick processor

      r = uniform(); 
      int idx = rnd(r*(n-1));
      id = list_->item(idx);
      // avoid repeats
      ID_LIST & l = list_->set_quantity();
      l.erase(l.begin()+idx);
      return energy(id.second);
    }
    else { 
      return maxEnergy_;
    }
  }
  //--------------------------------------------------------
  //  insert atoms
  //--------------------------------------------------------
  bool ConcentrationRegulatorMethodTransition::insert_atoms(int n)
  {

    insertionIds_.clear();
    DENS_VEC x(3); x = 0; 
    DENS_VEC v(3); v = 0;
    const DENS_MAN & T = atc_->field(TEMPERATURE);
    int additions = 0;
    int attempts = 0;
    while(additions < n && attempts < maxAttempts_){
      if(accept(insertion_location(x))) {
        DENS_VEC Tv = atc_->fe_engine()->interpolate_field(x,T);
Tv(0) = 300.; 
        pick_velocity(v,Tv(0)); // 3 normal
        int nlocal = lammpsInterface_->insert_atom(transitionType_,controlMask_,x.ptr(),v.ptr()); // no charge
        insertionIds_.push_back(pair<int,int>(-1,nlocal)); // atc id unknown
        additions += 1;
#ifdef ATC_VERBOSE2
        lammpsInterface_->print_msg(">>> insert x:"+to_string(x(0))+" "+to_string(x(1))+" "+to_string(x(2))+" v:"+to_string(v(0))+" "+to_string(v(1))+" "+to_string(v(2))+" "+to_string(rngCounter_));
#endif
      }
      attempts++;
      //lammpsInterface_->barrier();
      additions = lammpsInterface_->int_allmax(additions);
#ifdef ATC_VERBOSE
{
      string c = to_string(controlType_);
      lammpsInterface_->all_print(_rngUniformCounter_,c+" rng-uniform  ");
      lammpsInterface_->all_print(_rngNormalCounter_,c+" rng-normal  ");
//    lammpsInterface_->all_print(uniform()," random ");
}
#endif
      if (atc_->parallel_consistency()) { sync_random_number_generators(); }
#ifdef ATC_VERBOSE2
        lammpsInterface_->print_msg("attempts: "+to_string(attempts)+" additions "+to_string(additions)+" : "+to_string(rngCounter_));
#endif
#ifdef ATC_VERBOSE
{
      string c = to_string(controlType_);
      lammpsInterface_->all_print(attempts, c+"+attempts  ");
      lammpsInterface_->all_print(additions,c+" additions ");
      lammpsInterface_->all_print(_rngUniformCounter_,c+" RNG-uniform  ");
      lammpsInterface_->all_print(_rngNormalCounter_,c+" RNG-normal  ");
//    lammpsInterface_->all_print(uniform()," RANDOM ");
}
#endif
    }
    return (n == additions); // success
  }
  //--------------------------------------------------------
  //  sync random number generators
  //--------------------------------------------------------
  void ConcentrationRegulatorMethodTransition::sync_random_number_generators() const
  {
    // normal
    int n = lammpsInterface_->int_allmax(_rngNormalCounter_);
    int dn = n - _rngNormalCounter_;
    lammpsInterface_->advance_random_normal(randomNumberGenerator_,dn);
    _rngNormalCounter_ = n;
    // uniform
    int u = lammpsInterface_->int_allmax(_rngUniformCounter_);
    int du = u - _rngUniformCounter_;
    lammpsInterface_->advance_random_uniform(randomNumberGenerator_,du);
    _rngUniformCounter_ = u;
  }
  //--------------------------------------------------------
  //  pick location
  //--------------------------------------------------------
  double ConcentrationRegulatorMethodTransition::insertion_location(DENS_VEC & x) const
  {
     // pick random element  
     int elem = pick_element(); // 1 uniform
     // pick random local coordinate
     DENS_VEC xi(3);
     pick_coordinates(elem,xi,x); // 3 uniform
//   if (! lammpsInterface_->in_box(x.ptr())) { throw ATC_Error("new atom is not in box");}
     if (lammpsInterface_->in_my_processor_box(x.ptr())) {
#ifdef ATC_VERBOSE2
       lammpsInterface_->print_msg(">>> insertion_location e:" +to_string(elem)+" xi:"+to_string(xi(0))+" "+to_string(xi(1))+" "+to_string(xi(2))+" x:"+to_string(x(0))+" "+to_string(x(1))+" "+to_string(x(2))+ " energy "+to_string(energy(x.ptr()))+" "+true_false(accept(energy(x.ptr())))+" "+to_string(rngUniformCounter_));
#endif
       return energy(x.ptr());
     }
     else { 
       return maxEnergy_; 
     }
  }
  //--------------------------------------------------------
  //  pick element
  //--------------------------------------------------------
  int ConcentrationRegulatorMethodTransition::pick_element() const
  {
    double r = uniform();
    ESET::const_iterator itr = elemset_.begin(); // global?
    for (int i = 0; i < volumes_.size() ; ++i) {
      if (r < volumes_(i)) return *itr;
      itr++;
    }
    return *itr;
  }
  //--------------------------------------------------------
  //  pick coordinates
  //--------------------------------------------------------
  void ConcentrationRegulatorMethodTransition::pick_coordinates(const int elem,
                                                     DENS_VEC & xi, 
                                                     DENS_VEC & x) const
  {
    xi.reset(3);
    xi(0) = 2.*uniform()-1.;
    xi(1) = 2.*uniform()-1.;
    xi(2) = 2.*uniform()-1.;
    atc_->fe_engine()->fe_mesh()->position(elem,xi,x);

  }
  //--------------------------------------------------------
  //  pick velocity
  //--------------------------------------------------------
  void ConcentrationRegulatorMethodTransition::pick_velocity(DENS_VEC & v,double T) const
  {
    double s = sigma_*sqrt(T);
    v(0) = s*normal();
    v(1) = s*normal();
    v(2) = s*normal();
//v = 0;
  }
  //--------------------------------------------------------
  //  transition
  //--------------------------------------------------------
  void ConcentrationRegulatorMethodTransition::transition()
  {
    transitionCounter_++;
    //if (insertionIds_.size() == 0) return; // 
    if      (transitionCounter_> transitionInterval_) { 
      nInTransition_ = 0; 
      return;
    }
    else if (transitionCounter_==transitionInterval_) {
      nInTransition_ -= lammpsInterface_->change_type(transitionType_,controlType_);
    }
    else {
      transitionFactor_ = insertion_factor(transitionCounter_);
      if (nInTransition_ < 0) transitionFactor_ = 1-transitionFactor_;
      double q = 0; 
      lammpsInterface_->set_charge(transitionType_,q);
      DENS_VEC eps = epsilon0_;
      
      lammpsInterface_->set_epsilons(transitionType_,p_,eps.ptr());
      lammpsInterface_->pair_reinit(); // epsilon
    }
  }
  //--------------------------------------------------------
  //  diagnostics
  //--------------------------------------------------------
  double ConcentrationRegulatorMethodTransition::compute_vector(int n) const
  {
    if      (n==0) return count() - targetCount_; 
    else if (n==1) return count()/V_;
    else if (n==2) return (1.-transitionFactor_)*nInTransition_;
    else if (n==3) return _rngUniformCounter_;
    else if (n==4) return _rngNormalCounter_;
    else if (n==5) return lammpsInterface_->random_state(randomNumberGenerator_);
    else return 0;
  }
}; // end namespace
