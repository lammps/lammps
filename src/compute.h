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

#ifndef LMP_COMPUTE_H
#define LMP_COMPUTE_H

#include "pointers.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class Compute : protected Pointers {
  friend class Neighbor;

 public:
  // clang-format off
  enum {
    INVOKED_NONE    = 0,
    INVOKED_SCALAR  = 1<<0,
    INVOKED_VECTOR  = 1<<1,
    INVOKED_ARRAY   = 1<<2,
    INVOKED_PERATOM = 1<<3,
    INVOKED_LOCAL   = 1<<4,
    INVOKED_PERGRID = 1<<5,
  };
  // clang-format on
  static int instance_total;    // # of Compute classes ever instantiated

  char *id, *style;
  int igroup, groupbit;

  double scalar;           // computed global scalar
  double *vector;          // computed global vector
  double **array;          // computed global array
  double *vector_atom;     // computed per-atom vector
  double **array_atom;     // computed per-atom array
  double *vector_local;    // computed local vector
  double **array_local;    // computed local array

  int scalar_flag;                 // 0/1 if compute_scalar() function exists
  int vector_flag;                 // 0/1 if compute_vector() function exists
  int array_flag;                  // 0/1 if compute_array() function exists
  int thermo_modify_colname;       // 1 if fix has custom column names for output
  int size_vector;                 // length of global vector
  int size_array_rows;             // rows in global array
  int size_array_cols;             // columns in global array
  int size_vector_variable;        // 1 if vec length is unknown in advance
  int size_array_rows_variable;    // 1 if array rows is unknown in advance

  int peratom_flag;         // 0/1 if compute_peratom() function exists
  int size_peratom_cols;    // 0 = vector, N = columns in peratom array

  int local_flag;         // 0/1 if compute_local() function exists
  int size_local_rows;    // rows in local vector or array
  int size_local_cols;    // 0 = vector, N = columns in local array

  int pergrid_flag;    // 0/1 if compute_pergrid() function exists

  int extscalar;    // 0/1 if global scalar is intensive/extensive
  int extvector;    // 0/1/-1 if global vector is all int/ext/extlist
  int *extlist;     // list of 0/1 int/ext for each vec component
  int extarray;     // 0/1 if global array is all intensive/extensive

  int tempflag;            // 1 if Compute can be used as temperature
                           // must have both compute_scalar, compute_vector
  int pressflag;           // 1 if Compute can be used as pressure (uses virial)
                           // must have both compute_scalar, compute_vector
  int pressatomflag;       // 1 if Compute calculates per-atom virial
                           // 2 if Compute calculates per-atom centroid virial
                           // 3 if Compute calculates both
  int peflag;              // 1 if Compute calculates PE (uses Force energies)
  int peatomflag;          // 1 if Compute calculates per-atom PE
  int create_attribute;    // 1 if compute stores attributes that need
                           // setting when a new atom is created

  int tempbias;    // 0/1 if Compute temp includes self/extra bias

  int timeflag;     // 1 if Compute stores list of timesteps it's called on
  int ntime;        // # of entries in time list
  int maxtime;      // max # of entries time list can hold
  bigint *tlist;    // list of timesteps the Compute is called on

  int initialized_flag;      // 1 if compute is initialized, 0 if not
  int invoked_flag;          // non-zero if invoked or accessed this step, 0 if not
  bigint invoked_scalar;     // last timestep on which compute_scalar() was invoked
  bigint invoked_vector;     // ditto for compute_vector()
  bigint invoked_array;      // ditto for compute_array()
  bigint invoked_peratom;    // ditto for compute_peratom()
  bigint invoked_local;      // ditto for compute_local()
  bigint invoked_pergrid;    // ditto for compute_grid()

  double dof;    // degrees-of-freedom for temperature

  int comm_forward;           // size of forward communication (0 if none)
  int comm_reverse;           // size of reverse communication (0 if none)
  int dynamic_group_allow;    // 1 if can be used with dynamic group, else 0

  // KOKKOS host/device flag and data masks

  ExecutionSpace execution_space;
  unsigned int datamask_read, datamask_modify;

  int copymode, kokkosable;

  Compute(class LAMMPS *, int, char **);
  ~Compute() override;
  void modify_params(int, char **);
  virtual void reset_extra_dof();

  void init_flags();
  virtual void init() = 0;
  virtual void init_list(int, class NeighList *) {}
  virtual void setup() {}
  virtual double compute_scalar() { return 0.0; }
  virtual void compute_vector() {}
  virtual void compute_array() {}
  virtual void compute_peratom() {}
  virtual void compute_local() {}
  virtual void compute_pergrid() {}
  virtual void set_arrays(int) {}
  virtual std::string get_thermo_colname(int) { return {};  }

  virtual int pack_forward_comm(int, int *, double *, int, int *) { return 0; }
  virtual void unpack_forward_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) { return 0; }
  virtual void unpack_reverse_comm(int, int *, double *) {}

  virtual void reset_grid() {};

  virtual int get_grid_by_name(const std::string &, int &) { return -1; };
  virtual void *get_grid_by_index(int) { return nullptr; };
  virtual int get_griddata_by_name(int, const std::string &, int &) { return -1; };
  virtual void *get_griddata_by_index(int) { return nullptr; };

  virtual void dof_remove_pre() {}
  virtual int dof_remove(int) { return 0; }
  virtual void remove_bias(int, double *) {}
  virtual void remove_bias_thr(int, double *, double *) {}
  virtual void remove_bias_all() {}
  virtual void remove_bias_all_kk() {}
  virtual void reapply_bias_all() {}
  virtual void restore_bias(int, double *) {}
  virtual void restore_bias_thr(int, double *, double *) {}
  virtual void restore_bias_all() {}

  virtual void reset_extra_compute_fix(const char *);

  virtual void lock_enable() {}
  virtual void lock_disable() {}
  virtual int lock_length() { return 0; }
  virtual void lock(class Fix *, bigint, bigint) {}
  virtual void unlock(class Fix *) {}

  virtual void refresh() {}

  void addstep(bigint);
  int matchstep(bigint);
  void clearstep();

  bool is_initialized() const { return initialized_flag == 1; }

  virtual double memory_usage() { return 0.0; }

  virtual void pair_setup_callback(int, int) {}
  virtual void pair_tally_callback(int, int, int, int, double, double, double, double, double,
                                   double)
  {
  }

 protected:
  int instance_me;    // which Compute class instantiation I am

  double natoms_temp;    // # of atoms used for temperature calculation
  double extra_dof;      // extra DOF for temperature computes
  double fix_dof;        // DOF due to fixes
  int dynamic;           // recount atoms for temperature computes
  int dynamic_user;      // user request for temp compute to be dynamic

  double vbias[3];      // stored velocity bias for one atom
  double **vbiasall;    // stored velocity bias for all atoms
  int maxbias;          // size of vbiasall array

  inline int sbmask(int j) const { return j >> SBBITS & 3; }

  // private methods

  void adjust_dof_fix();
};

}    // namespace LAMMPS_NS

#endif
