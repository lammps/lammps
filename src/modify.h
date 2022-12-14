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

#ifndef LMP_MODIFY_H
#define LMP_MODIFY_H

#include "pointers.h"

#include <map>

namespace LAMMPS_NS {

class Compute;
class Fix;

class Modify : protected Pointers {
  friend class Info;
  friend class FixSRP;
  friend class Respa;
  friend class RespaOMP;

 public:
  int n_initial_integrate, n_post_integrate, n_pre_exchange;
  int n_pre_neighbor, n_post_neighbor;
  int n_pre_force, n_pre_reverse, n_post_force_any;
  int n_final_integrate, n_end_of_step;
  int n_energy_couple, n_energy_global, n_energy_atom;
  int n_initial_integrate_respa, n_post_integrate_respa;
  int n_pre_force_respa, n_post_force_respa_any, n_final_integrate_respa;
  int n_min_pre_exchange, n_min_pre_neighbor, n_min_post_neighbor;
  int n_min_pre_force, n_min_pre_reverse, n_min_post_force, n_min_energy;

  int restart_pbc_any;         // 1 if any fix sets restart_pbc
  int nfix_restart_global;     // stored fix global info from restart file
  int nfix_restart_peratom;    // stored fix peratom info from restart file

  int nfix, maxfix;
  Fix **fix;     // list of fixes
  int *fmask;    // bit mask for when each fix is applied

  int ncompute, maxcompute;
  Compute **compute;    // list of computes

  Modify(class LAMMPS *);
  ~Modify() override;
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
  virtual void pre_reverse(int, int);
  virtual void post_force(int);
  virtual void final_integrate();
  virtual void end_of_step();
  virtual double energy_couple();
  virtual double energy_global();
  virtual void energy_atom(int, double *);
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
  virtual void min_pre_reverse(int, int);
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

  void reset_grid();

  Fix *add_fix(int, char **, int trysuffix = 1);
  Fix *add_fix(const std::string &, int trysuffix = 1);
  Fix *replace_fix(const char *, int, char **, int trysuffix = 1);
  Fix *replace_fix(const std::string &, const std::string &, int trysuffix = 1);
  void modify_fix(int, char **);
  void delete_fix(const std::string &);
  void delete_fix(int);

  // deprecated API
  int find_fix(const std::string &);
  // new API
  Fix *get_fix_by_id(const std::string &) const;
  Fix *get_fix_by_index(int idx) const { return ((idx >= 0) && (idx < nfix)) ? fix[idx] : nullptr; }
  const std::vector<Fix *> get_fix_by_style(const std::string &) const;
  const std::vector<Fix *> &get_fix_list();

  Compute *add_compute(int, char **, int trysuffix = 1);
  Compute *add_compute(const std::string &, int trysuffix = 1);
  void modify_compute(int, char **);
  void delete_compute(const std::string &);
  void delete_compute(int);

  // deprecated API
  int find_compute(const std::string &);
  // new API
  Compute *get_compute_by_id(const std::string &) const;
  Compute *get_compute_by_index(int idx) const
  {
    return ((idx >= 0) && (idx < ncompute)) ? compute[idx] : nullptr;
  }
  const std::vector<Compute *> get_compute_by_style(const std::string &) const;
  const std::vector<Compute *> &get_compute_list();

  void clearstep_compute();
  void addstep_compute(bigint);
  void addstep_compute_all(bigint);

  int check_package(const char *);
  int check_rigid_group_overlap(int);
  int check_rigid_region_overlap(int, class Region *);
  int check_rigid_list_overlap(int *);

  void write_restart(FILE *);
  int read_restart(FILE *);
  void restart_deallocate(int);

  double memory_usage();

 protected:
  // internal fix counts

  int n_post_force, n_post_force_group, n_post_force_respa;

  // lists of fixes to apply at different stages of timestep

  int *list_initial_integrate, *list_post_integrate;
  int *list_pre_exchange, *list_pre_neighbor, *list_post_neighbor;
  int *list_pre_force, *list_pre_reverse;
  int *list_post_force, *list_post_force_group;
  int *list_final_integrate, *list_end_of_step;
  int *list_energy_couple, *list_energy_global, *list_energy_atom;
  int *list_initial_integrate_respa, *list_post_integrate_respa;
  int *list_pre_force_respa, *list_post_force_respa;
  int *list_final_integrate_respa;
  int *list_min_pre_exchange, *list_min_pre_neighbor, *list_min_post_neighbor;
  int *list_min_pre_force, *list_min_pre_reverse, *list_min_post_force;
  int *list_min_energy;

  int *end_of_step_every;

  int n_timeflag;    // list of computes that store time invocation
  int *list_timeflag;

  char **id_restart_global;       // stored fix global info
  char **style_restart_global;    // from read-in restart file
  char **state_restart_global;
  int *used_restart_global;

  char **id_restart_peratom;       // stored fix peratom info
  char **style_restart_peratom;    // from read-in restart file
  int *index_restart_peratom;
  int *used_restart_peratom;

  int index_permanent;    // fix/compute index returned to library call

  // vectors to be used for new-API accessors as wrapper
  std::vector<Fix *> fix_list;
  std::vector<Compute *> compute_list;

  void list_init(int, int &, int *&);
  void list_init_end_of_step(int, int &, int *&);
  void list_init_energy_couple(int &, int *&);
  void list_init_energy_global(int &, int *&);
  void list_init_energy_atom(int &, int *&);
  void list_init_post_force_group(int &, int *&);
  void list_init_post_force_respa_group(int &, int *&);
  void list_init_dofflag(int &, int *&);
  void list_init_compute();

 public:
  typedef Compute *(*ComputeCreator)(LAMMPS *, int, char **);
  typedef std::map<std::string, ComputeCreator> ComputeCreatorMap;
  ComputeCreatorMap *compute_map;

  typedef Fix *(*FixCreator)(LAMMPS *, int, char **);
  typedef std::map<std::string, FixCreator> FixCreatorMap;
  FixCreatorMap *fix_map;

 protected:
  void create_factories();
};

}    // namespace LAMMPS_NS

#endif
