/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_H
#define LMP_PAIR_H

#include "pointers.h"  // IWYU pragma: export

namespace LAMMPS_NS {

class Pair : protected Pointers {
  friend class AngleSDK;
  friend class AngleSDKOMP;
  friend class BondQuartic;
  friend class BondQuarticOMP;
  friend class DihedralCharmm;
  friend class DihedralCharmmOMP;
  friend class FixGPU;
  friend class FixIntel;
  friend class FixOMP;
  friend class ThrOMP;
  friend class Info;

 public:
  static int instance_total;     // # of Pair classes ever instantiated

  double eng_vdwl,eng_coul;      // accumulated energies
  double virial[6];              // accumulated virial
  double *eatom,**vatom;         // accumulated per-atom energy/virial
  double **cvatom;               // accumulated per-atom centroid virial

  double cutforce;               // max cutoff for all atom pairs
  double **cutsq;                // cutoff sq for each atom pair
  int **setflag;                 // 0/1 = whether each i,j has been set

  int comm_forward;              // size of forward communication (0 if none)
  int comm_reverse;              // size of reverse communication (0 if none)
  int comm_reverse_off;          // size of reverse comm even if newton off

  int single_enable;             // 1 if single() routine exists
  int single_hessian_enable;     // 1 if single_hessian() routine exists
  int restartinfo;               // 1 if pair style writes restart info
  int respa_enable;              // 1 if inner/middle/outer rRESPA routines
  int one_coeff;                 // 1 if allows only one coeff * * call
  int manybody_flag;             // 1 if a manybody potential
  int no_virial_fdotr_compute;   // 1 if does not invoke virial_fdotr_compute()
  int writedata;                 // 1 if writes coeffs to data file
  int ghostneigh;                // 1 if pair style needs neighbors of ghosts
  double **cutghost;             // cutoff for each ghost pair

  int ewaldflag;                 // 1 if compatible with Ewald solver
  int pppmflag;                  // 1 if compatible with PPPM solver
  int msmflag;                   // 1 if compatible with MSM solver
  int dispersionflag;            // 1 if compatible with LJ/dispersion solver
  int tip4pflag;                 // 1 if compatible with TIP4P solver
  int dipoleflag;                // 1 if compatible with dipole solver
  int spinflag;                  // 1 if compatible with spin solver
  int reinitflag;                // 1 if compatible with fix adapt and alike

  int centroidstressflag;        // compatibility with centroid atomic stress
                                 // 1 if same as two-body atomic stress
                                 // 2 if implemented and different from two-body
                                 // 4 if not compatible/implemented

  int tail_flag;                 // pair_modify flag for LJ tail correction
  double etail,ptail;            // energy/pressure tail corrections
  double etail_ij,ptail_ij;

  int evflag;                    // energy,virial settings
  int eflag_either,eflag_global,eflag_atom;
  int vflag_either,vflag_global,vflag_atom,cvflag_atom;

  int ncoultablebits;            // size of Coulomb table, accessed by KSpace
  int ndisptablebits;            // size of dispersion table
  double tabinnersq;
  double tabinnerdispsq;
  double *rtable,*drtable,*ftable,*dftable,*ctable,*dctable;
  double *etable,*detable,*ptable,*dptable,*vtable,*dvtable;
  double *rdisptable, *drdisptable, *fdisptable, *dfdisptable;
  double *edisptable, *dedisptable;
  int ncoulshiftbits,ncoulmask;
  int ndispshiftbits, ndispmask;

  int nextra;                    // # of extra quantities pair style calculates
  double *pvector;               // vector of extra pair quantities

  int single_extra;              // number of extra single values calculated
  double *svector;               // vector of extra single quantities

  class NeighList *list;         // standard neighbor list used by most pairs
  class NeighList *listhalf;     // half list used by some pairs
  class NeighList *listfull;     // full list used by some pairs

  int allocated;                 // 0/1 = whether arrays are allocated
                                 //       public so external driver can check
  int compute_flag;              // 0 if skip compute()

  enum{GEOMETRIC,ARITHMETIC,SIXTHPOWER};   // mixing options

  int beyond_contact, nondefault_history_transfer;   // for granular styles

  // KOKKOS host/device flag and data masks

  ExecutionSpace execution_space;
  unsigned int datamask_read,datamask_modify;

  Pair(class LAMMPS *);
  virtual ~Pair();

  // top-level Pair methods

  void init();
  virtual void reinit();
  virtual void setup() {}
  double mix_energy(double, double, double, double);
  double mix_distance(double, double);
  void write_file(int, char **);
  void init_bitmap(double, double, int, int &, int &, int &, int &);
  virtual void modify_params(int, char **);
  void compute_dummy(int, int);

  // need to be public, so can be called by pair_style reaxc

  void v_tally(int, double *, double *);
  void ev_tally(int, int, int, int, double, double, double,
                double, double, double);
  void ev_tally3(int, int, int, double, double,
                 double *, double *, double *, double *);
  void v_tally3(int, int, int, double *, double *, double *, double *);
  void v_tally4(int, int, int, int, double *, double *, double *,
                double *, double *, double *);
  void ev_tally_xyz(int, int, int, int, double, double,
                    double, double, double, double, double, double);

  // general child-class methods

  virtual void compute(int, int) = 0;
  virtual void compute_inner() {}
  virtual void compute_middle() {}
  virtual void compute_outer(int, int) {}

  virtual double single(int, int, int, int,
                        double, double, double,
                        double& fforce) {
    fforce = 0.0;
    return 0.0;
  }

  void hessian_twobody(double fforce, double dfac, double delr[3], double phiTensor[6]);

  virtual double single_hessian(int, int, int, int,
                        double, double[3], double, double,
                        double& fforce, double d2u[6]) {
    fforce = 0.0;
    for (int i=0; i<6; i++) d2u[i] = 0;
    return 0.0;
  }

  virtual void settings(int, char **) = 0;
  virtual void coeff(int, char **) = 0;

  virtual void init_style();
  virtual void init_list(int, class NeighList *);
  virtual double init_one(int, int) {return 0.0;}

  virtual void init_tables(double, double *);
  virtual void init_tables_disp(double);
  virtual void free_tables();
  virtual void free_disp_tables();

  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual void write_restart_settings(FILE *) {}
  virtual void read_restart_settings(FILE *) {}
  virtual void write_data(FILE *) {}
  virtual void write_data_all(FILE *) {}

  virtual int pack_forward_comm(int, int *, double *, int, int *) {return 0;}
  virtual void unpack_forward_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_reverse_comm(int, int *, double *) {}
  virtual double memory_usage();

  void set_copymode(int value) {copymode = value;}

  // specific child-class methods for certain Pair styles

  virtual void *extract(const char *, int &) {return NULL;}
  virtual void swap_eam(double *, double **) {}
  virtual void reset_dt() {}
  virtual void min_xf_pointers(int, double **, double **) {}
  virtual void min_xf_get(int) {}
  virtual void min_x_set(int) {}
  virtual void transfer_history(double *, double*) {}

  // management of callbacks to be run from ev_tally()

 protected:
  int num_tally_compute;
  class Compute **list_tally_compute;
 public:
  virtual void add_tally_callback(class Compute *);
  virtual void del_tally_callback(class Compute *);

 protected:
  int instance_me;        // which Pair class instantiation I am

  int special_lj[4];           // copied from force->special_lj for Kokkos

  int suffix_flag;             // suffix compatibility flag

                                       // pair_modify settings
  int offset_flag,mix_flag;            // flags for offset and mixing
  double tabinner;                     // inner cutoff for Coulomb table
  double tabinner_disp;                 // inner cutoff for dispersion table


 public:
  // custom data type for accessing Coulomb tables

  typedef union {int i; float f;} union_int_float_t;

  // Accessor for the user-intel package to determine virial calc for hybrid

  inline int fdotr_is_set() const { return vflag_fdotr; }

 protected:
  int vflag_fdotr;
  int maxeatom,maxvatom,maxcvatom;

  int copymode;   // if set, do not deallocate during destruction
                  // required when classes are used as functors by Kokkos

  void ev_init(int eflag, int vflag, int alloc = 1) {
    if (eflag||vflag) ev_setup(eflag, vflag, alloc);
    else ev_unset();
  }
  virtual void ev_setup(int, int, int alloc = 1);
  void ev_unset();
  void ev_tally_full(int, double, double, double, double, double, double);
  void ev_tally_xyz_full(int, double, double,
                         double, double, double, double, double, double);
  void ev_tally4(int, int, int, int, double,
                 double *, double *, double *, double *, double *, double *);
  void ev_tally_tip4p(int, int *, double *, double, double);
  void v_tally2(int, int, double, double *);
  void v_tally_tensor(int, int, int, int,
                      double, double, double, double, double, double);
  void virial_fdotr_compute();

  // union data struct for packing 32-bit and 64-bit ints into double bufs
  // see atom_vec.h for documentation

  union ubuf {
    double d;
    int64_t i;
    ubuf(double arg) : d(arg) {}
    ubuf(int64_t arg) : i(arg) {}
    ubuf(int arg) : i(arg) {}
  };

  inline int sbmask(int j) const {
    return j >> SBBITS & 3;
  }
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Too many total bits for bitmapped lookup table

Table size specified via pair_modify command is too large.  Note that
a value of N generates a 2^N size table.

E: Cannot have both pair_modify shift and tail set to yes

These 2 options are contradictory.

E: Cannot use pair tail corrections with 2d simulations

The correction factors are only currently defined for 3d systems.

W: Using pair tail corrections with non-periodic system

This is probably a bogus thing to do, since tail corrections are
computed by integrating the density of a periodic system out to
infinity.

W: Using pair tail corrections with pair_modify compute no

The tail corrections will thus not be computed.

W: Using pair potential shift with pair_modify compute no

The shift effects will thus not be computed.

W: Using a manybody potential with bonds/angles/dihedrals and special_bond exclusions

This is likely not what you want to do.  The exclusion settings will
eliminate neighbors in the neighbor list, which the manybody potential
needs to calculated its terms correctly.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Fix adapt interface to this pair style not supported

New coding for the pair style would need to be done.

E: Pair style requires a KSpace style

No kspace style is defined.

E: BUG: restartinfo=1 but no restart support in pair style

The pair style has a bug, where it does not support reading
and writing information to a restart file, but does not set
the member variable restartinfo to 0 as required in that case.

E: Cannot yet use compute tally with Kokkos

This feature is not yet supported.

E: Pair style does not support pair_write

The pair style does not have a single() function, so it can
not be invoked by pair write.

E: Invalid atom types in pair_write command

Atom types must range from 1 to Ntypes inclusive.

E: Invalid style in pair_write command

Self-explanatory.  Check the input script.

E: Invalid cutoffs in pair_write command

Inner cutoff must be larger than 0.0 and less than outer cutoff.

E: Cannot open pair_write file

The specified output file for pair energies and forces cannot be
opened.  Check that the path and name are correct.

E: Bitmapped lookup tables require int/float be same size

Cannot use pair tables on this machine, because of word sizes.  Use
the pair_modify command with table 0 instead.

W: Table inner cutoff >= outer cutoff

You specified an inner cutoff for a Coulombic table that is longer
than the global cutoff.  Probably not what you wanted.

E: Too many exponent bits for lookup table

Table size specified via pair_modify command does not work with your
machine's floating point representation.

E: Too many mantissa bits for lookup table

Table size specified via pair_modify command does not work with your
machine's floating point representation.

E: Too few bits for lookup table

Table size specified via pair_modify command does not work with your
machine's floating point representation.

*/
