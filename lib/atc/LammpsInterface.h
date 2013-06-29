#ifndef LAMMPS_INTERFACE_H
#define LAMMPS_INTERFACE_H

#include <iostream>
#include <stdlib.h>
#include "mpi.h"
#include "../../src/lammps.h"
#include "../../src/lmptype.h"

#include "ATC_TypeDefs.h"

// Forward class declarations for LAMMPS_NS namespace
namespace LAMMPS_NS {
  class LAMMPS;
  class NeighList;
  class Compute;
}

namespace ATC {

/**
 *  @class LammpsInterface
 *  @brief Singleton class that handles all interfacing with the lammps code
 */

class LammpsInterface {

 public:

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
    LJ,
    REAL,
    METAL
  };

  /** Static instance of this class */
  static LammpsInterface * instance();

  /** Set lammps pointer */
  void set_lammps(LAMMPS_NS::LAMMPS * lammps) 
  { 
    lammps_ = lammps; 
    MPI_Comm_rank(lammps_->world, & commRank_);
  }

  /** \name Methods that interface with Lammps base class */
  /*@{*/
  MPI_Comm world();

  void allsum(double * send_buf, double * rec_buf, int count = 1)
  {
    MPI_Allreduce(send_buf, rec_buf, count, MPI_DOUBLE, MPI_SUM,
                  lammps_->world);
    MPI_Barrier(lammps_->world);
  }

  void int_allsum(int * send_buf, int * rec_buf, int count = 1)
  {
    MPI_Allreduce(send_buf, rec_buf, count, MPI_INT, MPI_SUM,
                  lammps_->world);
    MPI_Barrier(lammps_->world);
  }

  void allmax(double * send_buf, double * rec_buf, int count = 1)
  {
    MPI_Allreduce(send_buf, rec_buf, count, MPI_DOUBLE, MPI_MAX,
                  lammps_->world);
    MPI_Barrier(lammps_->world);
  }

  void int_allmax(int * send_buf, int * rec_buf, int count = 1)
  {
    MPI_Allreduce(send_buf, rec_buf, count, MPI_INT, MPI_MAX,
                  lammps_->world);
    MPI_Barrier(lammps_->world);
  }

  void logical_or(int * send_buf, int * rec_buf, int count = 1)
  {
    MPI_Allreduce(send_buf, rec_buf, count, MPI_INT, MPI_LOR,
                  lammps_->world);
    MPI_Barrier(lammps_->world);
  }

  int comm_rank(void) { return commRank_;}
  /*@}*/

  /** \name Methods that interface with Atom class */
  /*@{*/
  int nlocal();
  int nghost();
  int nmax();
  double natoms();
  double ** xatom(); 
  int ntypes();
  const double ** xatom() const;
  double ** vatom(); 
  double ** fatom(); 
  int * atom_mask(); 
  int * atom_type(); 
  int * atom_tag(); 
  double * atom_mass();
  double   atom_mass(int iType);
  double * atom_rmass();
  double * atom_charge();
  void unwrap_coordinates(int iatom, double* xatom);
  /*@}*/

  /** \name Methods that interface with Domain class */
  /*@{*/
  int dimension();
  int nregion();
  void get_box_bounds(double & boxxlo, double & boxxhi,
                      double & boxylo, double & boxyhi,
                      double & boxzlo, double &boxzhi);
  int xperiodic();
  int yperiodic();
  int zperiodic();
  int nperiodic();
  void get_box_periodicity(int & xperiodic,
                           int & yperiodic,
                           int & zperiodic);
  int get_region_id(const char * regionName);
  double domain_xprd();
  double domain_yprd();
  double domain_zprd();
  double domain_xy();
  double domain_xz();
  double domain_yz();
  int domain_triclinic();
  /*@}*/

  /** \name Methods that interface with Update class */
  UnitsType units_style();
  /*@}*/
  
  /** \name Methods that interface with Lattice class */
  /*@{*/
  double xlattice();
  double ylattice();
  double zlattice();
  LatticeType lattice_style();
  int get_n_basis();
  void get_basis(double **basis);
  void get_unit_cell(double *a1, double *a2, double *a3);
  /** these functions are more than just simple pass throughs */
  int  num_atoms_per_cell(void);
  double  volume_per_atom(void);
  void get_lattice(MATRIX &N, MATRIX &B);
  /*@}*/

  /** \name Methods that interface with Force class */
  /*@{*/
  double boltz();
  double mvv2e();
  double ftm2v();
  double nktv2p();
  double qqr2e();
  double qe2f();
  double dielectric();
  double qqrd2e();
  double pair_force(int i, int j, double rsq, double& fmag_over_rmag);
  double pair_cutoff();
  int single_enable();
  /** these functions are more than just simple pass throughs */
  /** Boltzmann's constant in M,L,T,t units */
  double kBoltzmann(void);
  /** Dulong-Petit heat capacity per volume in M,L,T,t units */
  double heat_capacity(void);
  /** mass per volume in reference configuraturation in M,L units */
  double mass_density(void);
  /** ratio of permittivity of free space over elemental charge in V, L units */
  double epsilon0(void) {return 0.00552635; } // [V A]^-1
  // NOTE this won't work for LJ/SI/CGS units where [L] != A
  /*@}*/

  /** \name Methods that interface with Group class */
  /*@{*/
  int ngroup();
  int group_bit(int iGroup);
  int find_group(const char * c);
  int group_inverse_mask(int iGroup);
  char * group_name(int iGroup);
  void group_bounds(int iGroup, double * b);
  /*@}*/

  /** \name Methods that interface with Memory class */
  /*@{*/
  double * create_1d_double_array(int nlo, int nhi, const char *name);
  void destroy_1d_double_array(double * d, int i);
  double ** create_2d_double_array(int n1, int n2, const char *name);
  void destroy_2d_double_array(double **d);
  double **grow_2d_double_array(double **array, int n1, int n2, const char *name);
  int ** create_2d_int_array(int n1, int n2, const char *name);
  void destroy_2d_int_array(int **i);
  int ** grow_2d_int_array(int **array, int n1, int n2, const char *name);
  /*@}*/

  /** \name Methods that interface with Update class */
  /*@{*/
  double dt();
  int    ntimestep();
  int    nsteps();
  /*@}*/

  /** \name Methods that interface with neighbor list */
  /*@{*/
  void init_list(int id, LAMMPS_NS::NeighList *ptr);
  int   neighbor_list_inum();
  int * neighbor_list_numneigh();
  int * neighbor_list_ilist();
  int ** neighbor_list_firstneigh();
  int   neighbor_ago();
  /*@}*/

  /** \name Methods that interface with Region class */
  /*@{*/
  char * region_name(int iRegion);
  char * region_style(int iRegion);
  double region_xlo(int iRegion);
  double region_xhi(int iRegion);
  double region_ylo(int iRegion);
  double region_yhi(int iRegion);
  double region_zlo(int iRegion);
  double region_zhi(int iRegion);
  double region_xscale(int iRegion);
  double region_yscale(int iRegion);
  double region_zscale(int iRegion);
  int region_match(int iRegion, double x, double y, double z);
  /*@}*/

  /** \name Methods that interface with compute class */
  enum COMPUTE_INVOKED 
   {DUMMY0,INVOKED_SCALAR,INVOKED_VECTOR,DUMMMY3,INVOKED_PERATOM};
  int find_compute(const char* tag);
  LAMMPS_NS::Compute*  get_compute(const char* tag);
  double*   compute_scalar_data(const char* tag);
  double**  compute_vector_data(const char* tag);
  int  compute_ncols(const char* tag);
  /*@}*/

  /** \name Methods that interface with compute pe/atom class */
  /*@{*/
  int  atomPE_create(void);
  void atomPE_init(void);
  void atomPE_addstep(LAMMPS_NS::bigint step);
  double * atomPE_compute(void);
  /*@}*/

  /** Return lammps pointer -- only as a last resort! */
  LAMMPS_NS::LAMMPS * get_lammps_ptr();

 protected:

  LAMMPS_NS::LAMMPS * lammps_;

  /** access to neighbor list */
  LAMMPS_NS::NeighList *list_;

  /** constructor */
  LammpsInterface();

  /** comm rank */
  int commRank_;

  /** compute pe/atom */
  LAMMPS_NS::Compute * atomPE_;

 private:

  static LammpsInterface * myInstance_;

};

} // end namespace ATC



#endif
