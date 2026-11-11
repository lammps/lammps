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

#include "creator_registry.h"

#include <map>

namespace LAMMPS_NS {

class Compute;
class Fix;

class Modify : protected Pointers {
 public:
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
  virtual void fused_integrate(int) {}
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
  void image_flip(int, int, int);

  Fix *add_fix(int, char **, int trysuffix = 1);
  Fix *add_fix(const std::string &, int trysuffix = 1);
  Fix *replace_fix(const std::string &, int, char **, int trysuffix = 1);
  Fix *replace_fix(const std::string &, const std::string &, int trysuffix = 1);
  void modify_fix(int, char **);
  void delete_fix(const std::string &);
  void delete_fix(Fix *);

  [[nodiscard]] Fix *get_fix_by_id(const std::string &) const;
  [[nodiscard]] Fix *get_fix_by_index(int idx) const
  {
    return ((idx >= 0) && (idx < (int) fixes.size())) ? fixes[idx] : nullptr;
  }
  [[nodiscard]] std::vector<Fix *> get_fix_by_style(const std::string &) const;
  // the list must not be changed (fix added or deleted) while iterating over it
  [[nodiscard]] const std::vector<Fix *> &get_fix_list() const { return fixes; }
  [[nodiscard]] int get_fix_mask(Fix *) const;
  void set_fix_mask(Fix *, int);
  void clear_fix_mask(Fix *, int);

  // true if any fix is invoked at any of the stages in mask (bits from FixConst),
  // as of the last init(). used by integrators to skip stages without fixes

  [[nodiscard]] bool has_fixes(int mask) const { return (stagemask & mask) != 0; }

  // true if fixes or computes were added, deleted, or modified since the last init()

  [[nodiscard]] bool changed_since_init() const { return nchanges != nchanges_at_init; }

  Compute *add_compute(int, char **, int trysuffix = 1);
  Compute *add_compute(const std::string &, int trysuffix = 1);
  void modify_compute(int, char **);
  void delete_compute(const std::string &);
  void delete_compute(Compute *);

  [[nodiscard]] Compute *get_compute_by_id(const std::string &) const;
  [[nodiscard]] Compute *get_compute_by_index(int idx) const
  {
    return ((idx >= 0) && (idx < (int) computes.size())) ? computes[idx] : nullptr;
  }
  [[nodiscard]] std::vector<Compute *> get_compute_by_style(const std::string &) const;
  // the list must not be changed (compute added or deleted) while iterating over it
  [[nodiscard]] const std::vector<Compute *> &get_compute_list() const { return computes; }

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

  // true if any fix requires that atoms are remapped into the periodic box
  // when writing a restart file, as of the last init()

  [[nodiscard]] bool has_restart_pbc() const { return restart_pbc_any != 0; }

  // number of fixes with per-atom data read from a restart file that is
  // waiting for the fixes to be re-defined (reset by the next init())

  [[nodiscard]] int get_nfix_restart_peratom() const { return nfix_restart_peratom; }

  double memory_usage();

 protected:
  int restart_pbc_any;         // 1 if any fix sets restart_pbc
  int nfix_restart_global;     // stored fix global info from restart file
  int nfix_restart_peratom;    // stored fix peratom info from restart file

  // lists of fixes and computes in the order they were defined.
  // fixmask[i] is the mask of stages at which fixes[i] is invoked, as set by
  // Fix::setmask() and changed by set_fix_mask()/clear_fix_mask(), plus the
  // internal bits below.  A fix is replaced in place to keep its position.

  std::vector<Fix *> fixes;
  std::vector<int> fixmask;
  std::vector<Compute *> computes;

  // copy of fixmask taken by init() with the bits for the energy contributions
  // added: fixes are invoked according to this snapshot, so a fix added or a
  // mask changed after init() takes effect at the next init() only, as with
  // the lists of fix indices used before. stagemask is the union of callmask

  std::vector<int> callmask;
  int stagemask;

  // count of changes to the lists (add, delete, replace, fix_modify, compute_modify)
  // and its value at the last init(), for changed_since_init()

  int nchanges, nchanges_at_init;

  // invoke func(fix) for all fixes whose mask at init() has bit set, in list
  // order. the loop re-reads the list size and keeps no state in the class, so
  // func() may add or delete fixes and may itself call another stage (nested
  // calls from MC fixes) without disturbing the loop

  template <typename Func> void for_each_fix(int bit, Func &&func)
  {
    for (std::size_t i = 0; i < fixes.size(); ++i)
      if (callmask[i] & bit) func(fixes[i]);
  }
  void erase_fix_slot(std::size_t);

  char **id_restart_global;       // stored fix global info
  char **style_restart_global;    // from read-in restart file
  char **state_restart_global;
  int *used_restart_global;

  char **id_restart_peratom;       // stored fix peratom info
  char **style_restart_peratom;    // from read-in restart file
  int *index_restart_peratom;
  int *used_restart_peratom;

 public:
  using ComputeCreator = Compute *(*) (LAMMPS *, int, char **);
  using FixCreator = Fix *(*) (LAMMPS *, int, char **);

  // global registries of fix and compute style factory functions
  static CreatorRegistry<FixCreator> &fix_styles();
  static CreatorRegistry<ComputeCreator> &compute_styles();
};

}    // namespace LAMMPS_NS

#endif
