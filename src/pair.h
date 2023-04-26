/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_H
#define LMP_PAIR_H

#include "pointers.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class Pair : protected Pointers {
  friend class AngleSPICA;
  friend class AngleSPICAOMP;
  friend class BondQuartic;
  friend class BondQuarticOMP;
  friend class DihedralCharmm;
  friend class DihedralCharmmOMP;
  friend class FixGPU;
  friend class FixIntel;
  friend class FixOMP;
  friend class FixQEq;
  friend class PairHybrid;
  friend class PairHybridScaled;
  friend class ThrOMP;
  friend class Info;
  friend class Neighbor;

 public:
  static int instance_total;    // # of Pair classes ever instantiated

  double eng_vdwl, eng_coul;    // accumulated energies
  double virial[6];             // accumulated virial: xx,yy,zz,xy,xz,yz
  double *eatom, **vatom;       // accumulated per-atom energy/virial
  double **cvatom;              // accumulated per-atom centroid virial

  double cutforce;    // max cutoff for all atom pairs
  double **cutsq;     // cutoff sq for each atom pair
  int **setflag;      // 0/1 = whether each i,j has been set

  int comm_forward;        // size of forward communication (0 if none)
  int comm_reverse;        // size of reverse communication (0 if none)
  int comm_reverse_off;    // size of reverse comm even if newton off

  int single_enable;              // 1 if single() routine exists
  int born_matrix_enable;         // 1 if born_matrix() routine exists
  int single_hessian_enable;      // 1 if single_hessian() routine exists
  int restartinfo;                // 1 if pair style writes restart info
  int respa_enable;               // 1 if inner/middle/outer rRESPA routines
  int one_coeff;                  // 1 if allows only one coeff * * call
  int manybody_flag;              // 1 if a manybody potential
  int unit_convert_flag;          // value != 0 indicates support for unit conversion.
  int no_virial_fdotr_compute;    // 1 if does not invoke virial_fdotr_compute()
  int writedata;                  // 1 if writes coeffs to data file
  int finitecutflag;              // 1 if cut depends on finite atom size
  int ghostneigh;                 // 1 if pair style needs neighbors of ghosts
  double **cutghost;              // cutoff for each ghost pair

  int ewaldflag;         // 1 if compatible with Ewald solver
  int pppmflag;          // 1 if compatible with PPPM solver
  int msmflag;           // 1 if compatible with MSM solver
  int dispersionflag;    // 1 if compatible with LJ/dispersion solver
  int tip4pflag;         // 1 if compatible with TIP4P solver
  int dipoleflag;        // 1 if compatible with dipole solver
  int spinflag;          // 1 if compatible with spin solver
  int reinitflag;        // 1 if compatible with fix adapt and alike

  int centroidstressflag;    // centroid stress compared to two-body stress
                             // CENTROID_SAME = same as two-body stress
                             // CENTROID_AVAIL = different and implemented
                             // CENTROID_NOTAVAIL = different, not yet implemented

  int tail_flag;          // pair_modify flag for LJ tail correction
  double etail, ptail;    // energy/pressure tail corrections
  double etail_ij, ptail_ij;
  int trim_flag;    // pair_modify flag for trimming neigh list

  int evflag;    // energy,virial settings
  int eflag_either, eflag_global, eflag_atom;
  int vflag_either, vflag_global, vflag_atom, cvflag_atom;

  int ncoultablebits;    // size of Coulomb table, accessed by KSpace
  int ndisptablebits;    // size of dispersion table
  double tabinnersq;
  double tabinnerdispsq;
  double *rtable, *drtable, *ftable, *dftable, *ctable, *dctable;
  double *etable, *detable, *ptable, *dptable, *vtable, *dvtable;
  double *rdisptable, *drdisptable, *fdisptable, *dfdisptable;
  double *edisptable, *dedisptable;
  int ncoulshiftbits, ncoulmask;
  int ndispshiftbits, ndispmask;

  int nextra;         // # of extra quantities pair style calculates
  double *pvector;    // vector of extra pair quantities

  int single_extra;    // number of extra single values calculated
  double *svector;     // vector of extra single quantities

  class NeighList *list;        // standard neighbor list used by most pairs
  class NeighList *listhalf;    // half list used by some pairs
  class NeighList *listfull;    // full list used by some pairs

  int allocated;       // 0/1 = whether arrays are allocated
                       //       public so external driver can check
  int compute_flag;    // 0 if skip compute()
  int mixed_flag;      // 1 if all itype != jtype coeffs are from mixing
  bool did_mix;        // set to true by mix_energy() to indicate that mixing was performed

  enum { GEOMETRIC, ARITHMETIC, SIXTHPOWER };    // mixing options

  int beyond_contact, nondefault_history_transfer;    // for granular styles

  // KOKKOS host/device flag and data masks

  ExecutionSpace execution_space;
  unsigned int datamask_read, datamask_modify;
  int kokkosable;             // 1 if Kokkos pair
  int reverse_comm_device;    // 1 if reverse comm on Device

  Pair(class LAMMPS *);
  ~Pair() override;

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

  void ev_tally(int, int, int, int, double, double, double, double, double, double);
  void ev_tally3(int, int, int, double, double, double *, double *, double *, double *);
  void v_tally2_newton(int, double *, double *);
  void v_tally3(int, int, int, double *, double *, double *, double *);
  void v_tally4(int, int, int, int, double *, double *, double *, double *, double *, double *);

  // general child-class methods

  virtual void compute(int, int) = 0;
  virtual void compute_inner() {}
  virtual void compute_middle() {}
  virtual void compute_outer(int, int) {}

  virtual double single(int, int, int, int, double, double, double, double &fforce)
  {
    fforce = 0.0;
    return 0.0;
  }

  void hessian_twobody(double fforce, double dfac, double delr[3], double phiTensor[6]);

  virtual double single_hessian(int, int, int, int, double, double[3], double, double,
                                double &fforce, double d2u[6])
  {
    fforce = 0.0;
    for (int i = 0; i < 6; i++) d2u[i] = 0;
    return 0.0;
  }

  virtual void born_matrix(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/, double /*rsq*/,
                           double /*factor_coul*/, double /*factor_lj*/, double &du, double &du2)
  {
    du = du2 = 0.0;
  }

  virtual void finish() {}
  virtual void settings(int, char **) = 0;
  virtual void coeff(int, char **) = 0;

  virtual void init_style();
  virtual void init_list(int, class NeighList *);
  virtual double init_one(int, int) { return 0.0; }

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

  virtual int pack_forward_comm(int, int *, double *, int, int *) { return 0; }
  virtual void unpack_forward_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) { return 0; }
  virtual void unpack_reverse_comm(int, int *, double *) {}

  virtual void reset_grid() {}

  virtual void pack_forward_grid(int, void *, int, int *) {}
  virtual void unpack_forward_grid(int, void *, int, int *) {}
  virtual void pack_reverse_grid(int, void *, int, int *) {}
  virtual void unpack_reverse_grid(int, void *, int, int *) {}

  virtual double memory_usage();

  void set_copymode(int value) { copymode = value; }

  // specific child-class methods for certain Pair styles

  virtual void *extract(const char *, int &) { return nullptr; }
  virtual void *extract_peratom(const char *, int &) { return nullptr; }
  virtual void swap_eam(double *, double **) {}
  virtual void reset_dt() {}
  virtual void min_xf_pointers(int, double **, double **) {}
  virtual void min_xf_get(int) {}
  virtual void min_x_set(int) {}
  virtual void transfer_history(double *, double *, int, int) {}
  virtual double atom2cut(int) { return 0.0; }
  virtual double radii2cut(double, double) { return 0.0; }

  // management of callbacks to be run from ev_tally()

 protected:
  int num_tally_compute, did_tally_flag;
  class Compute **list_tally_compute;

 public:
  virtual void add_tally_callback(class Compute *);
  virtual void del_tally_callback(class Compute *);
  bool did_tally_callback() const { return did_tally_flag != 0; }

 protected:
  int instance_me;      // which Pair class instantiation I am
  int special_lj[4];    // copied from force->special_lj for Kokkos
  int suffix_flag;      // suffix compatibility flag

  // pair_modify settings
  int offset_flag, mix_flag;    // flags for offset and mixing
  double tabinner;              // inner cutoff for Coulomb table
  double tabinner_disp;         // inner cutoff for dispersion table

 protected:
  // for mapping of elements to atom types and parameters
  // mostly used for manybody potentials
  int nelements;        // # of unique elements
  char **elements;      // names of unique elements
  int *elem1param;      // mapping from elements to parameters
  int **elem2param;     // mapping from element pairs to parameters
  int ***elem3param;    // mapping from element triplets to parameters
  int *map;             // mapping from atom types to elements
  int nparams;          // # of stored parameter sets
  int maxparam;         // max # of parameter sets
  void map_element2type(int, char **, bool update_setflag = true);

 public:
  // custom data type for accessing Coulomb tables

  typedef union {
    int i;
    float f;
  } union_int_float_t;

  // Accessor for the INTEL package to determine virial calc for hybrid

  inline int fdotr_is_set() const { return vflag_fdotr; }

 protected:
  int vflag_fdotr;
  int maxeatom, maxvatom, maxcvatom;

  int copymode;    // if set, do not deallocate during destruction
                   // required when classes are used as functors by Kokkos

  void ev_init(int eflag, int vflag, int alloc = 1)
  {
    if (eflag || vflag)
      ev_setup(eflag, vflag, alloc);
    else
      ev_unset();
  }
  virtual void ev_setup(int, int, int alloc = 1);
  void ev_unset();
  void ev_tally_full(int, double, double, double, double, double, double);
  void ev_tally_xyz_full(int, double, double, double, double, double, double, double, double);
  void ev_tally4(int, int, int, int, double, double *, double *, double *, double *, double *,
                 double *);
  void ev_tally_tip4p(int, int *, double *, double, double);
  void ev_tally_xyz(int, int, int, int, double, double, double, double, double, double, double,
                    double);
  void v_tally2(int, int, double, double *);
  void v_tally_tensor(int, int, int, int, double, double, double, double, double, double);
  void virial_fdotr_compute();

  inline int sbmask(int j) const { return j >> SBBITS & 3; }
};

}    // namespace LAMMPS_NS

#endif
