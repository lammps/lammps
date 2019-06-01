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

#ifndef LMP_MODIFY_H
#define LMP_MODIFY_H

#include <cstdio>
#include "pointers.h"
#include <map>
#include <string>

namespace LAMMPS_NS {

class Modify : protected Pointers {
  friend class Info;
  friend class FixSRP;
  friend class Respa;
  friend class RespaOMP;

 public:
  int nfix,maxfix;
  int n_initial_integrate,n_post_integrate,n_pre_exchange;
  int n_pre_neighbor,n_post_neighbor;
  int n_pre_force,n_pre_reverse,n_post_force;
  int n_final_integrate,n_end_of_step,n_thermo_energy,n_thermo_energy_atom;
  int n_initial_integrate_respa,n_post_integrate_respa;
  int n_pre_force_respa,n_post_force_respa,n_final_integrate_respa;
  int n_min_pre_exchange,n_min_pre_neighbor,n_min_post_neighbor;
  int n_min_pre_force,n_min_pre_reverse,n_min_post_force,n_min_energy;

  int restart_pbc_any;       // 1 if any fix sets restart_pbc
  int nfix_restart_global;   // stored fix global info from restart file
  int nfix_restart_peratom;  // stored fix peratom info from restart file

  class Fix **fix;           // list of fixes
  int *fmask;                // bit mask for when each fix is applied

  int ncompute,maxcompute;   // list of computes
  class Compute **compute;

  Modify(class LAMMPS *);
  virtual ~Modify();
  virtual void init();
  virtual void setup(int);
  virtual void setup_pre_exchange();
  virtual void setup_pre_neighbor();
  virtual void setup_post_neighbor();
  virtual void setup_pre_force(int);
  virtual void setup_pre_reverse(int, int);
  virtual void initial_integrate(int);
  virtual void post_integrate();
  virtual void pre_exchange();
  virtual void pre_neighbor();
  virtual void post_neighbor();
  virtual void pre_force(int);
  virtual void pre_reverse(int,int);
  virtual void post_force(int);
  virtual void final_integrate();
  virtual void end_of_step();
  virtual double thermo_energy();
  virtual void thermo_energy_atom(int, double *);
  virtual void post_run();
  virtual void create_attribute(int);

  virtual void setup_pre_force_respa(int, int);
  virtual void initial_integrate_respa(int, int, int);
  virtual void post_integrate_respa(int, int);
  virtual void pre_force_respa(int, int, int);
  virtual void post_force_respa(int, int, int);
  virtual void final_integrate_respa(int, int);

  virtual void min_pre_exchange();
  virtual void min_pre_neighbor();
  virtual void min_post_neighbor();
  virtual void min_pre_force(int);
  virtual void min_pre_reverse(int,int);
  virtual void min_post_force(int);

  virtual double min_energy(double *);
  virtual void min_store();
  virtual void min_step(double, double *);
  virtual void min_clearstore();
  virtual void min_pushstore();
  virtual void min_popstore();
  virtual double max_alpha(double *);
  virtual int min_dof();
  virtual int min_reset_ref();

  void add_fix(int, char **, int trysuffix=1);
  void modify_fix(int, char **);
  void delete_fix(const char *);
  void delete_fix(int);
  int find_fix(const char *);
  int find_fix_by_style(const char *);
  int check_package(const char *);
  int check_rigid_group_overlap(int);
  int check_rigid_region_overlap(int, class Region *);
  int check_rigid_list_overlap(int *);

  void add_compute(int, char **, int trysuffix=1);
  void modify_compute(int, char **);
  void delete_compute(const char *);
  int find_compute(const char *);

  void clearstep_compute();
  void addstep_compute(bigint);
  void addstep_compute_all(bigint);

  void write_restart(FILE *);
  int read_restart(FILE *);
  void restart_deallocate(int);

  bigint memory_usage();

 protected:

  // lists of fixes to apply at different stages of timestep

  int *list_initial_integrate,*list_post_integrate;
  int *list_pre_exchange,*list_pre_neighbor,*list_post_neighbor;
  int *list_pre_force,*list_pre_reverse,*list_post_force;
  int *list_final_integrate,*list_end_of_step,*list_thermo_energy;
  int *list_thermo_energy_atom;
  int *list_initial_integrate_respa,*list_post_integrate_respa;
  int *list_pre_force_respa,*list_post_force_respa;
  int *list_final_integrate_respa;
  int *list_min_pre_exchange,*list_min_pre_neighbor,*list_min_post_neighbor;
  int *list_min_pre_force,*list_min_pre_reverse,*list_min_post_force;
  int *list_min_energy;

  int *end_of_step_every;

  int n_timeflag;            // list of computes that store time invocation
  int *list_timeflag;

  char **id_restart_global;           // stored fix global info
  char **style_restart_global;        // from read-in restart file
  char **state_restart_global;
  int *used_restart_global;

  char **id_restart_peratom;          // stored fix peratom info
  char **style_restart_peratom;       // from read-in restart file
  int *index_restart_peratom;
  int *used_restart_peratom;

  int index_permanent;        // fix/compute index returned to library call

  void list_init(int, int &, int *&);
  void list_init_end_of_step(int, int &, int *&);
  void list_init_thermo_energy(int, int &, int *&);
  void list_init_thermo_energy_atom(int &, int *&);
  void list_init_dofflag(int &, int *&);
  void list_init_compute();

 public:
  typedef Compute *(*ComputeCreator)(LAMMPS *, int, char **);
  typedef std::map<std::string,ComputeCreator> ComputeCreatorMap;
  ComputeCreatorMap *compute_map;

  typedef Fix *(*FixCreator)(LAMMPS *, int, char **);
  typedef std::map<std::string,FixCreator> FixCreatorMap;
  FixCreatorMap *fix_map;

 protected:
  template <typename T> static Compute *compute_creator(LAMMPS *, int, char **);
  template <typename T> static Fix *fix_creator(LAMMPS *, int, char **);
};

}

#endif

/* ERROR/WARNING messages:

E: Fix %s does not allow use of dynamic group

Dynamic groups have not yet been enabled for this fix.

E: Compute %s does not allow use of dynamic group

Dynamic groups have not yet been enabled for this compute.

W: One or more atoms are time integrated more than once

This is probably an error since you typically do not want to
advance the positions or velocities of an atom more than once
per timestep.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix command before simulation box is defined

The fix command cannot be used before a read_data, read_restart, or
create_box command.

E: Could not find fix group ID

A group ID used in the fix command does not exist.

E: Replacing a fix, but new style != old style

A fix ID can be used a 2nd time, but only if the style matches the
previous fix.  In this case it is assumed you with to reset a fix's
parameters.  This error may mean you are mistakenly re-using a fix ID
when you do not intend to.

W: Replacing a fix, but new group != old group

The ID and style of a fix match for a fix you are changing with a fix
command, but the new group you are specifying does not match the old
group.

E: Unrecognized fix style %s

The choice of fix style is unknown.

E: Could not find fix_modify ID

A fix ID used in the fix_modify command does not exist.

E: Could not find fix ID to delete

Self-explanatory.

E: Reuse of compute ID

A compute ID cannot be used twice.

E: Unrecognized compute style %s

The choice of compute style is unknown.

E: Could not find compute_modify ID

Self-explanatory.

E: Could not find compute ID to delete

Self-explanatory.

U: Unknown fix style

The choice of fix style is unknown.

U: Unknown compute style

The choice of compute style is unknown.

*/
