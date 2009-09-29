/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FIX_H
#define FIX_H

#include "pointers.h"

namespace LAMMPS_NS {

class Fix : protected Pointers {
 public:
  char *id,*style;
  int igroup,groupbit;

  int restart_global;            // 1 if Fix saves global state, 0 if not
  int restart_peratom;           // 1 if Fix saves peratom state, 0 if not
  int force_reneighbor;          // 1 if Fix forces reneighboring, 0 if not
  int box_change;                // 1 if Fix changes box size, 0 if not
  int next_reneighbor;           // next timestep to force a reneighboring
  int thermo_energy;             // 1 if fix_modify enabled ThEng, 0 if not
  int nevery;                    // how often to call an end_of_step fix
  int rigid_flag;                // 1 if Fix integrates rigid bodies, 0 if not
  int virial_flag;               // 1 if Fix contributes to virial, 0 if not
  int no_change_box;             // 1 if cannot swap ortho <-> triclinic
  int time_integrate;            // 1 if fix performs time integration, 0 if no
  int time_depend;               // 1 if fix is timestep dependent, 0 if not
  int create_attribute;          // 1 if fix stores attributes that need
                                 //   setting when a new atom is created
  int restart_pbc;               // 1 if fix moves atoms (except integrate)
                                 //   so that write_restart must remap to PBC

  int scalar_flag;               // 0/1 if compute_scalar() function exists
  int vector_flag;               // 0/1 if compute_vector() function exists
  int size_vector;               // N = size of global vector
  int scalar_vector_freq;        // frequency s/v data is available at
  int extscalar;                 // 0/1 if scalar is intensive/extensive
  int extvector;                 // 0/1/-1 if vector is all int/ext/extlist
  int *extlist;                  // list of 0/1 int/ext for each vec component

  int peratom_flag;              // 0/1 if per-atom data is stored
  int size_peratom;              // 0 = scalar_atom, N = size of vector_atom
  double *scalar_atom;           // computed per-atom scalar
  double **vector_atom;          // computed per-atom vector
  int peratom_freq;              // frequency per-atom data is available at

  int comm_forward;              // size of forward communication (0 if none)
  int comm_reverse;              // size of reverse communication (0 if none)

  double virial[6];              // accumlated virial
  double **vatom;                // accumulated per-atom virial

  int INITIAL_INTEGRATE,POST_INTEGRATE;    // mask settings
  int PRE_EXCHANGE,PRE_NEIGHBOR;
  int PRE_FORCE,POST_FORCE,FINAL_INTEGRATE,END_OF_STEP,THERMO_ENERGY;
  int INITIAL_INTEGRATE_RESPA,POST_INTEGRATE_RESPA;
  int PRE_FORCE_RESPA,POST_FORCE_RESPA,FINAL_INTEGRATE_RESPA;
  int MIN_PRE_EXCHANGE,MIN_POST_FORCE,MIN_ENERGY;

  Fix(class LAMMPS *, int, char **);
  virtual ~Fix();
  void modify_params(int, char **);

  virtual int setmask() = 0;

  virtual void init() {}
  virtual void init_list(int, class NeighList *) {}
  virtual void setup(int) {}
  virtual void min_setup(int) {}
  virtual void initial_integrate(int) {}
  virtual void post_integrate() {}
  virtual void pre_exchange() {}
  virtual void pre_neighbor() {}
  virtual void pre_force(int) {}
  virtual void post_force(int) {}
  virtual void final_integrate() {}
  virtual void end_of_step() {}
  virtual void write_restart(FILE *) {}
  virtual void restart(char *) {}

  virtual void grow_arrays(int) {}
  virtual void copy_arrays(int, int) {}
  virtual void set_arrays(int) {}
  virtual int pack_exchange(int, double *) {return 0;}
  virtual int unpack_exchange(int, double *) {return 0;}
  virtual int pack_restart(int, double *) {return 0;}
  virtual void unpack_restart(int, int) {}
  virtual int size_restart(int) {return 0;}
  virtual int maxsize_restart() {return 0;}

  virtual void initial_integrate_respa(int, int, int) {}
  virtual void post_integrate_respa(int, int) {}
  virtual void pre_force_respa(int, int, int) {}
  virtual void post_force_respa(int, int, int) {}
  virtual void final_integrate_respa(int) {}

  virtual void min_pre_exchange() {}
  virtual void min_post_force(int) {}

  virtual double min_energy(double *) {return 0.0;}
  virtual void min_store() {}
  virtual void min_clearstore() {}
  virtual void min_pushstore() {}
  virtual void min_popstore() {}
  virtual void min_step(double, double *) {}
  virtual double max_alpha(double *) {return 0.0;}
  virtual int min_dof() {return 0;}

  virtual int pack_comm(int, int *, double *, int, int *) {return 0;}
  virtual void unpack_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_reverse_comm(int, int *, double *) {}

  virtual double compute_scalar() {return 0.0;}
  virtual double compute_vector(int) {return 0.0;}

  virtual int dof(int) {return 0;}
  virtual void deform(int) {}
  virtual void reset_target(double) {}
  virtual void reset_dt() {}

  virtual int modify_param(int, char **) {return 0;}

  virtual void *extract(char *) {return NULL;}
  virtual double memory_usage() {return 0.0;}

 protected:
  int evflag;
  int vflag_global,vflag_atom;
  int maxvatom;

  void v_setup(int);
  void v_tally(int, int *, double, double *);
};

}

#endif
