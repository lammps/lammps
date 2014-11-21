#ifndef LAMMPS_INTERFACE_H
#define LAMMPS_INTERFACE_H

#include <iostream>
#include <stdlib.h>
#include <map>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include "mpi.h"
#include "lammps.h"
#include "modify.h"
#include "memory.h"
#include "random_park.h"
typedef LAMMPS_NS::RanPark* RNG_POINTER;
#include "lmptype.h"
#include "compute.h"
typedef const LAMMPS_NS::Compute* COMPUTE_POINTER;
#include "update.h"
#include "min.h"
#include "ATC_Error.h"
#include "ATC_TypeDefs.h"
#include "MatrixDef.h"
#include "MPI_Wrappers.h"

typedef LAMMPS_NS::Pair* POTENTIAL;

// Forward class declarations for LAMMPS_NS namespace
namespace LAMMPS_NS {
  class LAMMPS;
  class NeighList;
  class Compute;
  class ComputePEAtom;
  class ComputeStressAtom;
  class ComputeCentroAtom;
  class ComputeCNAAtom;
  class ComputeCoordAtom;
  class ComputeKEAtom;
  class Pair;
  class PairEAM;
  class Fix;
  class RanPark;
}

namespace ATC {

static const std::string atomPeNameBase_ = "atcPE";
static const double big_ = 1.e20;

/**
 *  @class LammpsInterface
 *  @brief Singleton class that handles all interfacing with the lammps code
 */

class LammpsInterface {

 public:

  // Enumeration of fundamental per-atom quantities, i.e. those defined directly by Lammps
  enum FundamentalAtomQuantity {
    ATOM_MASS = 0,
    ATOM_CHARGE,
    ATOM_POSITION,
    ATOM_VELOCITY,
    ATOM_FORCE,
    NUM_FUNDAMENTAL_ATOM_QUANTITIES
  };

  // Enumeration for lattice type. this MUST match the enum in src/lattice.cpp
  enum LatticeType {
    NONE=0,
    SC,
    BCC,
    FCC,
    HCP,
    DIAMOND,
    SQ,
    SQ2,
    HEX,
    CUSTOM
  };

  // Enumeration for units type. this is internal to ATC
  enum UnitsType {
    UNKNOWN=0,
    ATC,
    LJ,
    REAL,
    METAL,
    SI
  };

  // Provides a struct for easily passing/recovering data about SparseMats
  struct SparseMatInfo { 
    INDEX rows;
    INDEX cols;
    INDEX rowsCRS;
    INDEX size;
  };

  /** Static instance of this class */
  static LammpsInterface * instance();

  /** Destroy */
  static void Destroy();

  /** Set lammps pointer */
  void set_lammps(LAMMPS_NS::LAMMPS * lammps) 
  { 
    lammps_ = lammps; 
    MPI_Comm_rank(lammps_->world, & commRank_);
    MPI_Comm_size(lammps_->world, & commSize_);
  }

  /** \name Methods that interface with lammps base class */
  /*@{*/
// begin MPI --------------------------------------------------------------------
  MPI_Comm world() const;

  void broadcast(double * buf, int count = 1) const
  {
    MPI_Wrappers::broadcast(lammps_->world, buf, count);
  }

  void int_broadcast(int * buf, int count = 1) const
  {
    MPI_Wrappers::int_broadcast(lammps_->world, buf, count);
  }
  // send_buf is frequently a void* so MPI_IN_PLACE can be passed in
  void allsum(void * send_buf, double * rec_buf, int count = 1) const
  {
    MPI_Wrappers::allsum(lammps_->world, send_buf, rec_buf, count);
  }

  void int_allsum(void * send_buf, int * rec_buf, int count = 1) const
  {
    MPI_Wrappers::int_allsum(lammps_->world, send_buf, rec_buf, count);
  }
  int int_allsum(int & i) const
  {
    int j = 0;
    MPI_Wrappers::int_allsum(lammps_->world, &i, &j, 1);
    return j;
  }
  int int_scansum(int & i) const
  {
    int j = 0;
    MPI_Wrappers::int_scansum(lammps_->world, &i, &j, 1);
    return j;
  }
  int int_allmax(int & i) const
  {
    int j = 0;
    MPI_Wrappers::int_allmax(lammps_->world, &i, &j, 1);
    return j;
  }
  int int_allmin(int & i) const
  {
    int j = 0;
    MPI_Wrappers::int_allmin(lammps_->world, &i, &j, 1);
    return j;
  }
  double allmin(double & i) const
  {
    double j = 0;
    MPI_Wrappers::allmin(lammps_->world, &i, &j, 1);
    return j;
  }
  void sparse_allsum(SparseMatrix<double> &toShare) const
#ifdef ISOLATE_FE
  {
    MPI_Wrappers::sparse_allsum(lammps_->world, toShare); 
  }
#else
;
#endif
  void allmax(double * send_buf, double * rec_buf, int count = 1)
  {
    MPI_Wrappers::allmax(lammps_->world, send_buf, rec_buf, count);
  }
  void int_allmax(int * send_buf, int * rec_buf, int count = 1) const
  {
    MPI_Wrappers::int_allmax(lammps_->world, send_buf, rec_buf, count);
  }
  void allmin(double * send_buf, double * rec_buf, int count = 1) const
  {
    MPI_Wrappers::allmin(lammps_->world, send_buf, rec_buf, count);
  }
  void int_allmin(int * send_buf, int * rec_buf, int count = 1) const
  {
    MPI_Wrappers::int_allmin(lammps_->world, send_buf, rec_buf, count);
  }
  int rank_min(double * send_buf, double * rec_buf, int count = 1) const
  {
    return MPI_Wrappers::rank_min(lammps_->world, send_buf, rec_buf, count);
  }
  void int_recv(int * recv_buf, int max_size, int iproc) const
  {
    MPI_Wrappers::int_recv(lammps_->world, recv_buf, max_size, iproc);
  }
  void recv(double * recv_buf, int max_size, int iproc) const
  {
    MPI_Wrappers::recv(lammps_->world, recv_buf, max_size, iproc);
  }
  void int_send(int * send_buf, int send_size) const
  {
    MPI_Wrappers::int_send(lammps_->world, send_buf, send_size);
  }
  void send(double * send_buf, int send_size) const
  {
    MPI_Wrappers::send(lammps_->world, send_buf, send_size);
  }
  void allgatherv(double * send_buf, int send_count,
               double * rec_buf, int * rec_counts, int * displacements) const
  {
    MPI_Wrappers::allgatherv(lammps_->world, send_buf, send_count, rec_buf,
                             rec_counts, displacements);
  }
  void int_scatter(int * send_buf, int * rec_buf, int count = 1)
  {
    MPI_Wrappers::int_scatter(lammps_->world, send_buf, rec_buf, count);
  }
  void logical_or(void * send_buf, int * rec_buf, int count = 1) const
  {
    MPI_Wrappers::logical_or(lammps_->world, send_buf, rec_buf, count);
  }
  void barrier(void) const
  {
    MPI_Wrappers::barrier(lammps_->world);
  }
  void stop(std::string msg="") const
  {
    MPI_Wrappers::stop(lammps_->world, msg);
  }
// end MPI --------------------------------------------------------------------

  void print_debug(std::string msg="") const
  {
    std::cout << "rank " << comm_rank() << " " << msg << "\n" << std::flush; 
    barrier();
  }

  int comm_rank(void) const { return commRank_;}
  int comm_size(void) const { return commSize_;}
  bool rank_zero(void) const { return (commRank_==0);}
  bool serial(void) const { 
    int size = 1; MPI_Comm_size(lammps_->world,&size);
    return (size==1);
  }

  void print_msg(std::string msg) const
  {
    int me;
    MPI_Comm_rank(lammps_->world,&me);
    std::stringstream full_msg;
    if (serial()) {
      full_msg << " ATC: " << msg << "\n";
    } 
    else {
      full_msg << " ATC: P" << me << ", " << msg << "\n";
    }
    std::string mesg = full_msg.str();
    
    if (lammps_->screen)  fprintf(lammps_->screen, "%s",mesg.c_str());
    if (lammps_->logfile) fprintf(lammps_->logfile,"%s",mesg.c_str());
  }

  void print_msg_once(std::string msg,bool prefix=true, bool endline=true) const
  {
    int me;
    MPI_Comm_rank(lammps_->world,&me);
    if (me==0) {
      std::stringstream full_msg;
      if (prefix) full_msg << " ATC: ";
      full_msg << msg;
      if (endline) full_msg << "\n";
      std::string mesg = full_msg.str();
      if (lammps_->screen)  fprintf(lammps_->screen, "%s",mesg.c_str());
      if (lammps_->logfile) fprintf(lammps_->logfile,"%s",mesg.c_str());
    }
  }

  void all_print(double data, std::string tag ="") const
  {
    int me;
    MPI_Comm_rank(lammps_->world,&me);
    std::stringstream full_msg;
    if (serial()) {
      full_msg << " ATC: " << tag << data << "\n";
    } 
    else {
      int commSize = comm_size();
      double recv[commSize];
      MPI_Wrappers::gather(lammps_->world,data,recv);
      if (rank_zero()) {
        full_msg << " ATC:" << tag;
        for (int i = 0; i < commSize; i++) {
          full_msg << " P" << i << ": " << recv[i] ;
        }
        full_msg << "\n";
      }
    }
    if (rank_zero()) {
      std::string mesg = full_msg.str();
      if (lammps_->screen)  fprintf(lammps_->screen, "%s",mesg.c_str());
      if (lammps_->logfile) fprintf(lammps_->logfile,"%s",mesg.c_str());
    }
  }

  void stream_msg_once(std::string msg,bool prefix=true, bool endline=true) const
  {
    int me;
    MPI_Comm_rank(lammps_->world,&me);
    if (me==0) {
      if (prefix) std::cout << " ATC: ";
      std::cout << msg;
      if (endline) std::cout << "\n";
      std::cout << std::flush;
    }
  }


  void forward_comm_fix() const;
  void comm_borders() const;
  /*@}*/

  /** \name Methods that interface with Atom class */
  /*@{*/
  void set_fix_pointer(LAMMPS_NS::Fix * thisFix);
  std::string fix_id() const;
  bool atoms_sorted() const;
  LAMMPS_NS::bigint natoms() const;
  int nlocal() const;
  int nghost() const;
  int nmax() const;
  int ntypes() const;
  double ** xatom() const; 
  double ** vatom() const; 
  double ** fatom() const; 
  const int * atom_mask() const; 
  int * atom_mask();
  int * atom_type() const; 
  int * atom_tag() const; 
  int * atom_to_molecule() const;
  int * num_bond() const;
  int ** bond_atom() const;
  int * image() const;
  int bond_per_atom() const;
  int newton_bond() const;
  int local_to_global_map(int global) const;
  int type_to_charge(int t) const;
  //* Returns a pointer to the atom masses (NOT SAFE).
  double * atom_mass() const;
  //* Indexes an atomic mass by atom type (NO BOUNDS CHECK).
  double   atom_mass(int iType) const;
  double * atom_rmass() const;
  double * atom_charge() const;
  double *  atom_scalar(FundamentalAtomQuantity quantityType) const;
  double ** atom_vector(FundamentalAtomQuantity quantityType) const;
  int       atom_quantity_ndof(FundamentalAtomQuantity quantityType) const;
  double    atom_quantity_conversion(FundamentalAtomQuantity quantityType) const;
  void unwrap_coordinates(int iatom, double* xatom) const;
  /*@}*/

  /** \name Methods that interface with Domain class */
  /*@{*/
  int dimension() const;
  int nregion() const;
  void box_bounds(double & boxxlo, double & boxxhi,
                      double & boxylo, double & boxyhi,
                      double & boxzlo, double & boxzhi) const;
  bool in_box(double * x) const;
  bool in_my_processor_box(double * x) const;
  void sub_bounds(double & subxlo, double & subxhi,
                      double & subylo, double & subyhi,
                      double & subzlo, double & subzhi) const;
  int xperiodic() const;
  int yperiodic() const;
  int zperiodic() const;
  int nperiodic() const;
  void box_periodicity(int & xperiodic,
                           int & yperiodic,
                           int & zperiodic) const;
  void periodicity_correction(double * x) const;
  void set_reference_box(void) const; // const since called by perd_corr
  int region_id(const char * regionName) const;
  double domain_xprd() const;
  double domain_yprd() const;
  double domain_zprd() const;
  double domain_volume() const;
  double domain_xy() const;
  double domain_xz() const;
  double domain_yz() const;
  int domain_triclinic() const;
  bool region_bounds(const char * regionName,
                           double &xmin, double &xmax,
                           double &ymin, double & ymax,
                           double &zmin, double &zmax,
                           double &xscale,
                           double &yscale,
                           double &zscale) const;
  bool region_bounds(const char * regionName,
                           double &xmin, double &xmax,
                           double &ymin, double & ymax,
                           double &zmin, double &zmax) const {
    double xs,ys,zs;
    bool ifBlock = region_bounds(regionName,
      xmin,xmax,ymin,ymax,zmin,zmax,xs,ys,zs);
    xmin *= xs;
    xmax *= xs;
    ymin *= ys;
    ymax *= ys;
    zmin *= zs;
    zmax *= zs;
    return ifBlock;
  }
  /*@}*/
  void minimum_image(double & dx, double & dy, double & dz) const;
  void closest_image(const double * const xi, const double * const xj, double * const xjImage) const; 


  /** \name Methods that interface with Update class */
  UnitsType units_style() const;
  double convert_units(double value, UnitsType in, UnitsType out, int massExp, int lenExp, int timeExp, int engExp=0) const;
  //double minimize_energy() { return lammps_->update->minimize->ecurrent; }
  double minimize_energy() const { return lammps_->update->minimize->eprevious; }
  /*@}*/
  
  /** \name Methods that interface with Lattice class */
  /*@{*/
  double xlattice() const;
  double ylattice() const;
  double zlattice() const;
  LatticeType lattice_style() const;
  int n_basis() const;
  void basis_vectors(double **basis) const;
  double max_lattice_constant(void) const;
  double near_neighbor_cutoff(void) const;
  void unit_cell(double *a1, double *a2, double *a3) const;
  /** these functions are more than just simple pass throughs */
  int  num_atoms_per_cell(void) const;
  double  volume_per_atom(void) const;
  void lattice(Matrix<double> &N, Matrix<double> &B) const;
  /*@}*/

  /** \name Methods that interface with Force class */
  /*@{*/
  double boltz() const;
  double mvv2e() const;
  double ftm2v() const;
  double nktv2p() const;
  double qqr2e() const;
  double qe2f() const;
  double dielectric() const;
  double qqrd2e() const;
  double qv2e() const; // converts charge * voltage --> mass length^2 / time^2
  /*@}*/

  /** \name Methods that interface with pair class */
  /*@{*/
  // interface to "single"
  double pair_force(int i, int j, double rsq, double& fmag_over_rmag) const; // pair class
  double pair_force(int n, double rsq, double& fmag_over_rmag) const; // bond class
  double pair_force(std::map< std::pair< int,int >,int >::const_iterator itr, double rsq, double& fmag_over_rmag, int nbonds = 0) const; 
  double pair_force(std::pair< std::pair< int,int >,int > apair, double rsq, double& fmag_over_rmag, int nbonds = 0) const; 
  double pair_cutoff() const;
  void pair_reinit() const;
  int single_enable() const;
  LAMMPS_NS::PairEAM * pair_eam(void) const; 
  double bond_stiffness(int i, int j, double rsq) const;
  /*@}*/

  /** \name Methods for addition/deletion of atoms*/
  /*@{*/
  int delete_atom(int id) const;
  int insert_atom(int type, int mask, double* x, double* v, double q = 0) const;
  double shortrange_energy(double *x, int type, int id = -1, 
    double max = big_) const;
  int reset_ghosts(int dn) const;
  double shortrange_energy(int id, double max = big_) const;
  POTENTIAL potential(void) const;
  int type_to_groupbit(int itype) const; 
  int      change_type(int itype, int jtype) const;
  int      count_type(int itype) const;
  bool     epsilons(int type, POTENTIAL p, double * epsilons) const;
  bool set_epsilons(int type, POTENTIAL p, double * epsilons) const;
  int  set_charge(int type, double charge) const;
  /*@}*/

  /** \name interface to random number generator */
  /*@{*/
  RNG_POINTER random_number_generator() const;
  double random_uniform(RNG_POINTER p) const;
  double random_normal (RNG_POINTER p) const;
  int random_state (RNG_POINTER p) const;
  void set_random_state (RNG_POINTER p, int seed) const;
  void advance_random_generator (RNG_POINTER p, int n = 1) const;
  void advance_random_uniform (RNG_POINTER p, int n = 1) const;
  void advance_random_normal  (RNG_POINTER p, int n = 1) const;
  /*@}*/

  /** these functions are more than just simple pass throughs */
  /*@{*/
  /** Boltzmann's constant in M,L,T,t units */
  double kBoltzmann(void) const;
  /** Planck's constant (energy-time units)  */
  double hbar(void) const;
  /** Dulong-Petit heat capacity per volume in M,L,T,t units */
  double heat_capacity(void) const;
  /** mass per volume in reference configuraturation in M,L units */
  double mass_density(int* numPerType=NULL) const;
  /**  permittivity of free space, converts from LAMMPS potential units implied by the electric field units to LAMMPS charge units/LAMMPS length units (e.g., V to elemental charge/A) */
  double epsilon0(void) const;
  double coulomb_constant(void) const;
  double * special_coul() const;
  int newton_pair() const;
  double coulomb_factor(int & j) const {
    int n = nlocal() + nghost();
    double * sc =  special_coul();
    double factor_coul = 1.;
    if (j >= n) {
      factor_coul = sc[j/n];
      j %= n;
    }
    return factor_coul;
  } 
  /*@}*/

  /** \name Methods that interface with Group class */
  /*@{*/
  int ngroup() const;
  int group_bit(std::string name) const;
  int group_bit(int iGroup) const;
  int group_index(std::string name) const;
  int group_inverse_mask(int iGroup) const;
  char * group_name(int iGroup) const;
  void group_bounds(int iGroup, double * b) const;
  /*@}*/

  /** \name Methods that interface with Memory class */
  /*@{*/
  double * create_1d_double_array(int length, const char *name) const;
  double * grow_1d_double_array(double *array, int length, const char *name) const;
  void destroy_1d_double_array(double * d) const;
  double ** create_2d_double_array(int n1, int n2, const char *name) const;
  void destroy_2d_double_array(double **d) const;
  double **grow_2d_double_array(double **array, int n1, int n2, const char *name) const;
  int * create_1d_int_array(int length, const char *name) const;
  int * grow_1d_int_array(int *array, int length, const char *name) const;
  void destroy_1d_int_array(int * d) const;
  int ** create_2d_int_array(int n1, int n2, const char *name) const;
  void destroy_2d_int_array(int **i) const;
  int ** grow_2d_int_array(int **array, int n1, int n2, const char *name) const;
  template <typename T>
    T * grow_array(T *&array, int n, const char *name) const {return lammps_->memory->grow(array,n,name);};
  template <typename T>
    void destroy_array(T * array) {lammps_->memory->destroy(array);};
  template <typename T>
    T ** grow_array(T **&array, int n1, int n2, const char *name) const {return lammps_->memory->grow(array,n1,n2,name);};
  template <typename T>
    void destroy_array(T ** array) const {lammps_->memory->destroy(array);};
  /*@}*/

  /** \name Methods that interface with Update class */
  /*@{*/
  double dt() const;
  LAMMPS_NS::bigint ntimestep() const;
  int    nsteps() const;
  bool now(LAMMPS_NS::bigint f) { return (ntimestep() % f == 0); }
  /*@}*/

  /** \name Methods that interface with neighbor list */
  /*@{*/
  void neighbor_remap(int & j) const { j &= NEIGHMASK; }
  int sbmask(int j) const;
  void set_list(int id, LAMMPS_NS::NeighList *ptr) ;
  int   neighbor_list_inum() const;
  int * neighbor_list_numneigh() const;
  int * neighbor_list_ilist() const;
  int ** neighbor_list_firstneigh() const;
  int   neighbor_ago() const;
  int reneighbor_frequency() const;
  LAMMPS_NS::NeighList * neighbor_list(void) const { return list_;}
  /*@}*/

  /** \name Methods that interface with bond list */
  /*@{*/
  int   bond_list_length() const;
  int ** bond_list() const; // direct access
  int * bond_list(int n) const { return bond_list()[n];} 
  int   bond_list_i(int n) const    { return bond_list(n)[0];} 
  int   bond_list_j(int n) const    { return bond_list(n)[1];}
  int   bond_list_type(int n) const { return bond_list(n)[2];}
  /*@}*/

  /** \name Methods that interface with Region class */
  /*@{*/
  char * region_name(int iRegion) const;
  char * region_style(int iRegion) const;
  double region_xlo(int iRegion) const;
  double region_xhi(int iRegion) const;
  double region_ylo(int iRegion) const;
  double region_yhi(int iRegion) const;
  double region_zlo(int iRegion) const;
  double region_zhi(int iRegion) const;
  double region_xscale(int iRegion) const;
  double region_yscale(int iRegion) const;
  double region_zscale(int iRegion) const;
  int region_match(int iRegion, double x, double y, double z) const;
  /*@}*/

  /** \name Methods that interface with compute class */
  enum COMPUTE_INVOKED 
   {INVOKED_SCALAR=1,INVOKED_VECTOR=2,INVOKED_ARRAY=4,INVOKED_PERATOM=8};
  enum PER_ATOM_COMPUTE 
   {PE_ATOM,
    STRESS_ATOM,
    CENTRO_ATOM,
    CNA_ATOM,
    COORD_ATOM,
    KE_ATOM,
    NUM_PER_ATOM_COMPUTES};
  // computes owned by LAMMPS
  COMPUTE_POINTER compute_pointer(std::string tag) const;
  int      compute_ncols_peratom(COMPUTE_POINTER computePointer) const;
  double*  compute_vector_peratom(COMPUTE_POINTER computePointer) const;
  double** compute_array_peratom(COMPUTE_POINTER computePointer) const;

  void     computes_addstep(int step) const;
  void     compute_addstep(COMPUTE_POINTER computePointer, int step) const;
  int     compute_matchstep(COMPUTE_POINTER computePointer, int step) const;
  void     reset_invoked_flag(COMPUTE_POINTER computePointer) const;
  // computes owned by ATC
  int      create_compute_pe_peratom(void) const;
  double * compute_pe_peratom(void) const;
  std::string compute_pe_name(void) const {return atomPeNameBase_;};//  +fix_id();}; enables unique names, if desired
  void     computes_clearstep(void) const {lammps_->modify->clearstep_compute();};
  /*@}*/
 


  /** Return lammps pointer -- only as a last resort! */
  LAMMPS_NS::LAMMPS * lammps_pointer() const;

 protected:
  /** transfer a const compute pointer to a non-const computer pointer */
  LAMMPS_NS::Compute * const_to_active(const LAMMPS_NS::Compute* computePointer) const;

  LAMMPS_NS::LAMMPS * lammps_;

  LAMMPS_NS::Fix * fixPointer_;

  /** access to neighbor list */
  mutable LAMMPS_NS::NeighList *list_;

  /** constructor */
  LammpsInterface();

  /** comm rank */
  int commRank_;

  /** number of processes */
  int commSize_;

  /** compute pe/atom */
  mutable LAMMPS_NS::Compute * atomPE_;

  /** box info */
  mutable bool refBoxIsSet_;
  mutable double upper_[3],lower_[3],length_[3];

  /** registry of computer pointers */
  mutable std::set<LAMMPS_NS::Compute * > computePointers_;

  /** a random number generator from lammps */
  mutable LAMMPS_NS::RanPark * random_;
  mutable LAMMPS_NS::RanPark * globalrandom_;

 private:

  static LammpsInterface * myInstance_;

};

  class HeartBeat
  {
    public:
      HeartBeat(std::string name, int freq) :
        name_(name), freq_(freq), counter_(0) {};
      ~HeartBeat(){};
      void start() const
        { ATC::LammpsInterface::instance()->stream_msg_once(name_,true,false);}
      void next() const { if (counter_++ % freq_ == 0 )
          ATC::LammpsInterface::instance()->stream_msg_once(".",false,false);}
      void finish() const
        { ATC::LammpsInterface::instance()->stream_msg_once("done",false,true);}
    protected:
      std::string name_;
      int freq_;
      mutable int counter_;
    private:
      HeartBeat();
  };


} // end namespace ATC



#endif
