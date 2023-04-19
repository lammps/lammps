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

#ifndef LMP_FIX_H
#define LMP_FIX_H

#include "pointers.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class Fix : protected Pointers {
  friend class Neighbor;

 public:
  static int instance_total;    // # of Fix classes ever instantiated

  char *id, *style;
  int igroup, groupbit;

  int restart_global;      // 1 if Fix saves global state, 0 if not
  int restart_peratom;     // 1 if Fix saves peratom state, 0 if not
  int restart_file;        // 1 if Fix writes own restart file, 0 if not
  int force_reneighbor;    // 1 if Fix forces reneighboring, 0 if not

  int box_change;    // >0 if Fix changes box size, shape, or sub-domains, 0 if not
  // clang-format off
  enum {
    NO_BOX_CHANGE     = 0,
    BOX_CHANGE_ANY    = 1 << 0,
    BOX_CHANGE_DOMAIN = 1 << 1,
    BOX_CHANGE_X      = 1 << 2,
    BOX_CHANGE_Y      = 1 << 3,
    BOX_CHANGE_Z      = 1 << 4,
    BOX_CHANGE_YZ     = 1 << 5,
    BOX_CHANGE_XZ     = 1 << 6,
    BOX_CHANGE_XY     = 1 << 7,
    BOX_CHANGE_SIZE   = BOX_CHANGE_X  | BOX_CHANGE_Y  | BOX_CHANGE_Z,
    BOX_CHANGE_SHAPE  = BOX_CHANGE_YZ | BOX_CHANGE_XZ | BOX_CHANGE_XY
  };
  // clang-format on

  bigint next_reneighbor;      // next timestep to force a reneighboring
  int nevery;                  // how often to call an end_of_step fix
  int thermo_energy;           // 1 if fix_modify energy enabled, 0 if not
  int thermo_virial;           // 1 if fix_modify virial enabled, 0 if not
  int energy_global_flag;      // 1 if contributes to global eng
  int energy_peratom_flag;     // 1 if contributes to peratom eng
  int virial_global_flag;      // 1 if contributes to global virial
  int virial_peratom_flag;     // 1 if contributes to peratom virial
  int ecouple_flag;            // 1 if thermostat fix outputs cumulative
                               //      reservoir energy via compute_scalar()
  int time_integrate;          // 1 if performs time integration, 0 if no
  int rigid_flag;              // 1 if integrates rigid bodies, 0 if not
  int no_change_box;           // 1 if cannot swap ortho <-> triclinic
  int time_depend;             // 1 if requires continuous timestepping
  int create_attribute;        // 1 if fix stores attributes that need
                               //      setting when a new atom is created
  int restart_pbc;             // 1 if fix moves atoms (except integrate)
                               //      so write_restart must remap to PBC
  int wd_header;               // # of header values fix writes to data file
  int wd_section;              // # of sections fix writes to data file
  int dynamic_group_allow;     // 1 if can be used with dynamic group, else 0
  int dof_flag;                // 1 if has dof() method (not min_dof())
  int special_alter_flag;      // 1 if has special_alter() meth for spec lists
  int enforce2d_flag;          // 1 if has enforce2d method
  int respa_level_support;     // 1 if fix supports fix_modify respa
  int respa_level;             // which respa level to apply fix (1-Nrespa)
  int maxexchange;             // max # of per-atom values for Comm::exchange()
  int maxexchange_dynamic;     // 1 if fix sets maxexchange dynamically
  int pre_exchange_migrate;    // 1 if fix migrates atoms in pre_exchange()
  int stores_ids;              // 1 if fix stores atom IDs

  int scalar_flag;                 // 0/1 if compute_scalar() function exists
  int vector_flag;                 // 0/1 if compute_vector() function exists
  int array_flag;                  // 0/1 if compute_array() function exists
  int size_vector;                 // length of global vector
  int size_array_rows;             // rows in global array
  int size_array_cols;             // columns in global array
  int size_vector_variable;        // 1 if vec length is unknown in advance
  int size_array_rows_variable;    // 1 if array rows is unknown in advance
  int global_freq;                 // frequency s/v data is available at

  int peratom_flag;         // 0/1 if per-atom data is stored
  int size_peratom_cols;    // 0 = vector, N = columns in peratom array
  int peratom_freq;         // frequency per-atom data is available at

  int local_flag;         // 0/1 if local data is stored
  int size_local_rows;    // rows in local vector or array
  int size_local_cols;    // 0 = vector, N = columns in local array
  int local_freq;         // frequency local data is available at

  int pergrid_flag;       // 0/1 if per-grid data is stored
  int pergrid_freq;       // frequency per-grid data is available at

  int extscalar;    // 0/1 if global scalar is intensive/extensive
  int extvector;    // 0/1/-1 if global vector is all int/ext/extlist
  int *extlist;     // list of 0/1 int/ext for each vec component
  int extarray;     // 0/1 if global array is intensive/extensive

  double *vector_atom;     // computed per-atom vector
  double **array_atom;     // computed per-atom array
  double *vector_local;    // computed local vector
  double **array_local;    // computed local array

  int comm_forward;    // size of forward communication (0 if none)
  int comm_reverse;    // size of reverse communication (0 if none)
  int comm_border;     // size of border communication (0 if none)

  double virial[6];          // virial for this timestep
  double *eatom, **vatom;    // per-atom energy/virial for this timestep
  double **cvatom;           // per-atom centroid virial for this timestep

  int centroidstressflag;    // centroid stress compared to two-body stress
                             // CENTROID_SAME = same as two-body stress
                             // CENTROID_AVAIL = different and implemented
                             // CENTROID_NOTAVAIL = different, not yet implemented

  int restart_reset;    // 1 if restart just re-initialized fix

  // KOKKOS host/device flag and data masks

  int kokkosable;             // 1 if Kokkos fix
  int forward_comm_device;    // 1 if forward comm on Device
  int exchange_comm_device;   // 1 if exchange comm on Device
  ExecutionSpace execution_space;
  unsigned int datamask_read, datamask_modify;

  Fix(class LAMMPS *, int, char **);
  ~Fix() override;
  void modify_params(int, char **);

  virtual int setmask() = 0;

  virtual void post_constructor() {}
  virtual void init() {}
  virtual void init_list(int, class NeighList *) {}
  virtual void setup(int) {}
  virtual void setup_pre_exchange() {}
  virtual void setup_pre_neighbor() {}
  virtual void setup_post_neighbor() {}
  virtual void setup_pre_force(int) {}
  virtual void setup_pre_reverse(int, int) {}
  virtual void min_setup(int) {}
  virtual void initial_integrate(int) {}
  virtual void post_integrate() {}
  virtual void pre_exchange() {}
  virtual void pre_neighbor() {}
  virtual void post_neighbor() {}
  virtual void pre_force(int) {}
  virtual void pre_reverse(int, int) {}
  virtual void post_force(int) {}
  virtual void final_integrate() {}
  virtual void end_of_step() {}
  virtual void post_run() {}
  virtual void write_restart(FILE *) {}
  virtual void write_restart_file(const char *) {}
  virtual void restart(char *) {}

  virtual void grow_arrays(int) {}
  virtual void copy_arrays(int, int, int) {}
  virtual void set_arrays(int) {}
  virtual void update_arrays(int, int) {}
  virtual void set_molecule(int, tagint, int, double *, double *, double *);
  virtual void clear_bonus() {}

  virtual int pack_border(int, int *, double *) { return 0; }
  virtual int unpack_border(int, int, double *) { return 0; }
  virtual int pack_exchange(int, double *) { return 0; }
  virtual int unpack_exchange(int, double *) { return 0; }
  virtual int pack_restart(int, double *) { return 0; }
  virtual void unpack_restart(int, int) {}
  virtual int size_restart(int) { return 0; }
  virtual int maxsize_restart() { return 0; }

  virtual void setup_pre_force_respa(int, int) {}
  virtual void initial_integrate_respa(int, int, int) {}
  virtual void post_integrate_respa(int, int) {}
  virtual void pre_force_respa(int, int, int) {}
  virtual void post_force_respa(int, int, int) {}
  virtual void final_integrate_respa(int, int) {}

  virtual void min_pre_exchange() {}
  virtual void min_pre_neighbor() {}
  virtual void min_post_neighbor() {}
  virtual void min_pre_force(int) {}
  virtual void min_pre_reverse(int, int) {}
  virtual void min_post_force(int) {}

  virtual double min_energy(double *) { return 0.0; }
  virtual void min_store() {}
  virtual void min_clearstore() {}
  virtual void min_pushstore() {}
  virtual void min_popstore() {}
  virtual int min_reset_ref() { return 0; }
  virtual void min_step(double, double *) {}
  virtual double max_alpha(double *) { return 0.0; }
  virtual int min_dof() { return 0; }

  virtual int pack_forward_comm(int, int *, double *, int, int *) { return 0; }
  virtual void unpack_forward_comm(int, int, double *) {}
  virtual int pack_reverse_comm_size(int, int) { return 0; }
  virtual int pack_reverse_comm(int, int, double *) { return 0; }
  virtual void unpack_reverse_comm(int, int *, double *) {}

  virtual void reset_grid(){};

  virtual void pack_forward_grid(int, void *, int, int *){};
  virtual void unpack_forward_grid(int, void *, int, int *){};
  virtual void pack_reverse_grid(int, void *, int, int *){};
  virtual void unpack_reverse_grid(int, void *, int, int *){};
  virtual void pack_remap_grid(int, void *, int, int *){};
  virtual void unpack_remap_grid(int, void *, int, int *){};
  virtual int unpack_read_grid(int, char *) {return 0;};
  virtual void pack_write_grid(int, void *){};
  virtual void unpack_write_grid(int, void *, int *){};

  virtual int get_grid_by_name(const std::string &, int &) { return -1; };
  virtual void *get_grid_by_index(int) { return nullptr; };
  virtual int get_griddata_by_name(int, const std::string &, int &) { return -1; };
  virtual void *get_griddata_by_index(int) { return nullptr; };

  virtual double compute_scalar() { return 0.0; }
  virtual double compute_vector(int) { return 0.0; }
  virtual double compute_array(int, int) { return 0.0; }

  virtual int dof(int) { return 0; }
  virtual void deform(int) {}
  virtual void reset_target(double) {}
  virtual void reset_dt() {}
  virtual void enforce2d() {}

  virtual void read_data_header(char *) {}
  virtual void read_data_section(char *, int, char *, tagint) {}
  virtual bigint read_data_skip_lines(char *) { return 0; }

  virtual void write_data_header(FILE *, int) {}
  virtual void write_data_section_size(int, int &, int &) {}
  virtual void write_data_section_pack(int, double **) {}
  virtual void write_data_section_keyword(int, FILE *) {}
  virtual void write_data_section(int, FILE *, int, double **, int) {}

  virtual void zero_momentum() {}
  virtual void zero_rotation() {}

  virtual void rebuild_special() {}

  virtual int image(int *&, double **&) { return 0; }

  virtual int modify_param(int, char **) { return 0; }
  virtual void *extract(const char *, int &) { return nullptr; }

  virtual double memory_usage() { return 0.0; }

 protected:
  int instance_me;    // which Fix class instantiation I am

  int evflag;
  int eflag_either, eflag_global, eflag_atom;
  int vflag_either, vflag_global, vflag_atom, cvflag_atom;
  int maxeatom, maxvatom, maxcvatom;

  int copymode;    // if set, do not deallocate during destruction
                   // required when classes are used as functors by Kokkos

  int dynamic;    // recount atoms for temperature computes

  void ev_init(int eflag, int vflag)
  {
    if ((eflag && thermo_energy) || (vflag && thermo_virial))
      ev_setup(eflag, vflag);
    else
      evflag = eflag_either = eflag_global = eflag_atom = vflag_either = vflag_global = vflag_atom =
          cvflag_atom = 0;
  }
  void ev_setup(int, int);
  void ev_tally(int, int *, double, double, double *);

  void v_init(int vflag)
  {
    if (vflag && thermo_virial)
      v_setup(vflag);
    else
      evflag = vflag_either = vflag_global = vflag_atom = cvflag_atom = 0;
  }
  void v_setup(int);
  void v_tally(int, int *, double, double *);
  void v_tally(int, int *, double, double *, int, int, int[][2], double *, double[][3]);
  void v_tally(int, int *, double, double *, double[][3], double[][3], double[]);
  void v_tally(int, double *);
  void v_tally(int, int, double);
};

namespace FixConst {
  enum {
    INITIAL_INTEGRATE = 1 << 0,
    POST_INTEGRATE = 1 << 1,
    PRE_EXCHANGE = 1 << 2,
    PRE_NEIGHBOR = 1 << 3,
    POST_NEIGHBOR = 1 << 4,
    PRE_FORCE = 1 << 5,
    PRE_REVERSE = 1 << 6,
    POST_FORCE = 1 << 7,
    FINAL_INTEGRATE = 1 << 8,
    END_OF_STEP = 1 << 9,
    POST_RUN = 1 << 10,
    INITIAL_INTEGRATE_RESPA = 1 << 11,
    POST_INTEGRATE_RESPA = 1 << 12,
    PRE_FORCE_RESPA = 1 << 13,
    POST_FORCE_RESPA = 1 << 14,
    FINAL_INTEGRATE_RESPA = 1 << 15,
    MIN_PRE_EXCHANGE = 1 << 16,
    MIN_PRE_NEIGHBOR = 1 << 17,
    MIN_POST_NEIGHBOR = 1 << 18,
    MIN_PRE_FORCE = 1 << 19,
    MIN_PRE_REVERSE = 1 << 20,
    MIN_POST_FORCE = 1 << 21,
    MIN_ENERGY = 1 << 22
  };
}

}    // namespace LAMMPS_NS

#endif
