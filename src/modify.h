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

/** \class LAMMPS_NS::Modify
 * \brief Manager for fixes and computes in LAMMPS
 *
 * The Modify class manages all Fix and Compute instances in a LAMMPS
 * simulation. Fixes modify the simulation state (e.g., applying forces,
 * constraints, or thermostats), while Computes calculate properties from
 * the simulation (e.g., temperature, pressure, or per-atom quantities).
 *
 * Key responsibilities include:
 *
 * - Creating, storing, and deleting Fix and Compute instances
 * - Calling fixes at appropriate points during the timestep (initial_integrate,
 *   post_force, final_integrate, etc.)
 * - Managing fix execution order via bitmask flags
 * - Handling restart file read/write for fixes
 * - Managing energy contributions from fixes for thermodynamics
 * - Coordinating minimization callbacks for fixes
 *
 * The class maintains lists of fixes to be called at each stage of the
 * timestep, organized by callback type. It also provides style maps for
 * creating new fix and compute instances dynamically.
 *
 * \sa LAMMPS_NS::Fix, LAMMPS_NS::Compute, LAMMPS_NS::Pointers */
class Modify : protected Pointers {
  friend class Info;
  friend class FixSRP;
  friend class Respa;
  friend class RespaOMP;

 public:
  int n_initial_integrate;    /**< Number of fixes called at initial_integrate */
  int n_post_integrate;       /**< Number of fixes called at post_integrate */
  int n_pre_exchange;         /**< Number of fixes called at pre_exchange */
  int n_pre_neighbor;         /**< Number of fixes called at pre_neighbor */
  int n_post_neighbor;        /**< Number of fixes called at post_neighbor */
  int n_pre_force;            /**< Number of fixes called at pre_force */
  int n_pre_reverse;          /**< Number of fixes called at pre_reverse */
  int n_post_force_any;       /**< Number of fixes called at any post_force stage */
  int n_final_integrate;      /**< Number of fixes called at final_integrate */
  int n_end_of_step;          /**< Number of fixes called at end_of_step */
  int n_energy_couple;        /**< Number of fixes contributing coupled energy */
  int n_energy_global;        /**< Number of fixes contributing global energy */
  int n_energy_atom;          /**< Number of fixes contributing per-atom energy */
  int n_initial_integrate_respa;    /**< Number of fixes for rRESPA initial_integrate */
  int n_post_integrate_respa;       /**< Number of fixes for rRESPA post_integrate */
  int n_pre_force_respa;            /**< Number of fixes for rRESPA pre_force */
  int n_post_force_respa_any;       /**< Number of fixes for any rRESPA post_force */
  int n_final_integrate_respa;      /**< Number of fixes for rRESPA final_integrate */
  int n_min_pre_exchange;     /**< Number of fixes for minimization pre_exchange */
  int n_min_pre_neighbor;     /**< Number of fixes for minimization pre_neighbor */
  int n_min_post_neighbor;    /**< Number of fixes for minimization post_neighbor */
  int n_min_pre_force;        /**< Number of fixes for minimization pre_force */
  int n_min_pre_reverse;      /**< Number of fixes for minimization pre_reverse */
  int n_min_post_force;       /**< Number of fixes for minimization post_force */
  int n_min_energy;           /**< Number of fixes for minimization energy */

  int restart_pbc_any;        /**< 1 if any fix sets restart_pbc flag */
  int nfix_restart_global;    /**< Number of stored fix global info from restart */
  int nfix_restart_peratom;   /**< Number of stored fix per-atom info from restart */

  int nfix;       /**< Current number of defined fixes */
  int maxfix;     /**< Allocated size of fix array */
  Fix **fix;      /**< Array of pointers to Fix instances */
  int *fmask;     /**< Bitmask for when each fix is applied */

  int ncompute;       /**< Current number of defined computes */
  int maxcompute;     /**< Allocated size of compute array */
  Compute **compute;  /**< Array of pointers to Compute instances */

  /** Constructor for Modify class
   * \param lmp Pointer to the main LAMMPS instance */
  Modify(class LAMMPS *);

  /** Destructor - cleans up all fixes, computes, and style maps */
  ~Modify() override;

  /** Initialize all fixes and computes before a run */
  virtual void init();

  /** Setup all fixes at start of a run
   * \param vflag Virial flag */
  virtual void setup(int);

  /** Setup fixes for pre_exchange callback */
  virtual void setup_pre_exchange();

  /** Setup fixes for pre_neighbor callback */
  virtual void setup_pre_neighbor();

  /** Setup fixes for post_neighbor callback */
  virtual void setup_post_neighbor();

  /** Setup fixes for pre_force callback
   * \param vflag Virial flag */
  virtual void setup_pre_force(int);

  /** Setup fixes for pre_reverse callback
   * \param eflag Energy flag
   * \param vflag Virial flag */
  virtual void setup_pre_reverse(int, int);

  /** Call initial_integrate on all applicable fixes
   * \param vflag Virial flag */
  virtual void initial_integrate(int);

  /** Call post_integrate on all applicable fixes */
  virtual void post_integrate();

  /** Call pre_exchange on all applicable fixes */
  virtual void pre_exchange();

  /** Call pre_neighbor on all applicable fixes */
  virtual void pre_neighbor();

  /** Call post_neighbor on all applicable fixes */
  virtual void post_neighbor();

  /** Call pre_force on all applicable fixes
   * \param vflag Virial flag */
  virtual void pre_force(int);

  /** Call pre_reverse on all applicable fixes
   * \param eflag Energy flag
   * \param vflag Virial flag */
  virtual void pre_reverse(int, int);

  /** Call post_force on all applicable fixes
   * \param vflag Virial flag */
  virtual void post_force(int);

  /** Call final_integrate on all applicable fixes */
  virtual void final_integrate();

  /** Fused integration for performance
   * \param vflag Virial flag */
  virtual void fused_integrate(int) {}

  /** Call end_of_step on all applicable fixes */
  virtual void end_of_step();

  /** Compute coupled energy contribution from fixes
   * \return Coupled energy from all contributing fixes */
  virtual double energy_couple();

  /** Compute global energy contribution from fixes
   * \return Global energy from all contributing fixes */
  virtual double energy_global();

  /** Compute per-atom energy contribution from fixes
   * \param flag Energy flag
   * \param energy Array to store per-atom energies */
  virtual void energy_atom(int, double *);

  /** Call post_run on all applicable fixes */
  virtual void post_run();

  /** Create per-atom attributes for fixes
   * \param n Number of atoms */
  virtual void create_attribute(int);

  /** Setup for rRESPA pre_force
   * \param vflag Virial flag
   * \param ilevel rRESPA level */
  virtual void setup_pre_force_respa(int, int);

  /** Call initial_integrate for rRESPA
   * \param vflag Virial flag
   * \param ilevel rRESPA level
   * \param iloop Loop counter */
  virtual void initial_integrate_respa(int, int, int);

  /** Call post_integrate for rRESPA
   * \param ilevel rRESPA level
   * \param iloop Loop counter */
  virtual void post_integrate_respa(int, int);

  /** Call pre_force for rRESPA
   * \param vflag Virial flag
   * \param ilevel rRESPA level
   * \param iloop Loop counter */
  virtual void pre_force_respa(int, int, int);

  /** Call post_force for rRESPA
   * \param vflag Virial flag
   * \param ilevel rRESPA level
   * \param iloop Loop counter */
  virtual void post_force_respa(int, int, int);

  /** Call final_integrate for rRESPA
   * \param ilevel rRESPA level
   * \param iloop Loop counter */
  virtual void final_integrate_respa(int, int);

  /** Call pre_exchange for minimization */
  virtual void min_pre_exchange();

  /** Call pre_neighbor for minimization */
  virtual void min_pre_neighbor();

  /** Call post_neighbor for minimization */
  virtual void min_post_neighbor();

  /** Call pre_force for minimization
   * \param vflag Virial flag */
  virtual void min_pre_force(int);

  /** Call pre_reverse for minimization
   * \param eflag Energy flag
   * \param vflag Virial flag */
  virtual void min_pre_reverse(int, int);

  /** Call post_force for minimization
   * \param vflag Virial flag */
  virtual void min_post_force(int);

  /** Compute energy for minimization
   * \param fextra Extra forces array
   * \return Total energy */
  virtual double min_energy(double *);

  /** Store current state for line search */
  virtual void min_store();

  /** Take a step along search direction
   * \param alpha Step size
   * \param hextra Extra forces array */
  virtual void min_step(double, double *);

  /** Clear stored state */
  virtual void min_clearstore();

  /** Push current state onto stack */
  virtual void min_pushstore();

  /** Pop state from stack */
  virtual void min_popstore();

  /** Calculate maximum allowed step size
   * \param hextra Extra forces array
   * \return Maximum step size */
  virtual double max_alpha(double *);

  /** Get degrees of freedom for minimization
   * \return Number of degrees of freedom */
  virtual int min_dof();

  /** Reset reference state for minimization
   * \return Status flag */
  virtual int min_reset_ref();

  /** Reset all grid-based data structures */
  void reset_grid();

  /** Add a new fix from command-line arguments
   * \param narg Number of arguments
   * \param arg Argument array
   * \param trysuffix Flag to try accelerated suffix versions
   * \return Pointer to the newly created Fix */
  Fix *add_fix(int, char **, int trysuffix = 1);

  /** Add a new fix from a single string
   * \param fixcmd Command string containing fix definition
   * \param trysuffix Flag to try accelerated suffix versions
   * \return Pointer to the newly created Fix */
  Fix *add_fix(const std::string &, int trysuffix = 1);

  /** Replace an existing fix with a new one
   * \param id ID of the fix to replace
   * \param narg Number of arguments
   * \param arg Argument array
   * \param trysuffix Flag to try accelerated suffix versions
   * \return Pointer to the replacement Fix */
  Fix *replace_fix(const std::string &, int, char **, int trysuffix = 1);

  /** Replace an existing fix with a new one from a string
   * \param id ID of the fix to replace
   * \param fixcmd Command string for new fix
   * \param trysuffix Flag to try accelerated suffix versions
   * \return Pointer to the replacement Fix */
  Fix *replace_fix(const std::string &, const std::string &, int trysuffix = 1);

  /** Modify settings of an existing fix
   * \param narg Number of arguments
   * \param arg Argument array */
  void modify_fix(int, char **);

  /** Delete a fix by ID
   * \param id ID of the fix to delete */
  void delete_fix(const std::string &);

  /** Delete a fix by index
   * \param idx Index of the fix in the array */
  void delete_fix(int);

  // deprecated API

  /** Find a fix by ID (deprecated)
   * \param id ID of the fix to find
   * \return Index of fix in array, or -1 if not found
   * \deprecated Use get_fix_by_id() instead */
  int find_fix(const std::string &);

  // new API

  /** Get a fix by its ID
   * \param id ID of the fix to find
   * \return Pointer to Fix, or nullptr if not found */
  Fix *get_fix_by_id(const std::string &) const;

  /** Get a fix by its index in the array
   * \param idx Index of the fix
   * \return Pointer to Fix, or nullptr if index out of range */
  Fix *get_fix_by_index(int idx) const { return ((idx >= 0) && (idx < nfix)) ? fix[idx] : nullptr; }

  /** Get all fixes matching a style pattern
   * \param style Style name to match
   * \return Vector of pointers to matching Fix instances */
  const std::vector<Fix *> get_fix_by_style(const std::string &) const;

  /** Get the list of all fixes
   * \return Reference to vector containing all Fix pointers */
  const std::vector<Fix *> &get_fix_list();

  /** Get the execution mask for a fix
   * \param ifix Pointer to the Fix
   * \return Bitmask for the fix, or 0 if not found */
  int get_fix_mask(Fix *ifix) const
  {
    for (int i = 0; i < nfix; ++i) {
      if (fix[i] == ifix) return fmask[i];
    }
    return 0;
  }

  /** Add a new compute from command-line arguments
   * \param narg Number of arguments
   * \param arg Argument array
   * \param trysuffix Flag to try accelerated suffix versions
   * \return Pointer to the newly created Compute */
  Compute *add_compute(int, char **, int trysuffix = 1);

  /** Add a new compute from a single string
   * \param computecmd Command string containing compute definition
   * \param trysuffix Flag to try accelerated suffix versions
   * \return Pointer to the newly created Compute */
  Compute *add_compute(const std::string &, int trysuffix = 1);

  /** Modify settings of an existing compute
   * \param narg Number of arguments
   * \param arg Argument array */
  void modify_compute(int, char **);

  /** Delete a compute by ID
   * \param id ID of the compute to delete */
  void delete_compute(const std::string &);

  /** Delete a compute by index
   * \param idx Index of the compute in the array */
  void delete_compute(int);

  // deprecated API

  /** Find a compute by ID (deprecated)
   * \param id ID of the compute to find
   * \return Index of compute in array, or -1 if not found
   * \deprecated Use get_compute_by_id() instead */
  int find_compute(const std::string &);

  // new API

  /** Get a compute by its ID
   * \param id ID of the compute to find
   * \return Pointer to Compute, or nullptr if not found */
  Compute *get_compute_by_id(const std::string &) const;

  /** Get a compute by its index in the array
   * \param idx Index of the compute
   * \return Pointer to Compute, or nullptr if index out of range */
  Compute *get_compute_by_index(int idx) const
  {
    return ((idx >= 0) && (idx < ncompute)) ? compute[idx] : nullptr;
  }

  /** Get all computes matching a style pattern
   * \param style Style name to match
   * \return Vector of pointers to matching Compute instances */
  const std::vector<Compute *> get_compute_by_style(const std::string &) const;

  /** Get the list of all computes
   * \return Reference to vector containing all Compute pointers */
  const std::vector<Compute *> &get_compute_list();

  /** Clear the timestep for all computes */
  void clearstep_compute();

  /** Add a timestep to compute invocation list
   * \param step Timestep number */
  void addstep_compute(bigint);

  /** Add a timestep to all compute invocation lists
   * \param step Timestep number */
  void addstep_compute_all(bigint);

  /** Check if a package is being used
   * \param pkg Package name
   * \return 1 if package is in use, 0 otherwise */
  int check_package(const char *);

  /** Check for overlapping rigid body groups
   * \param groupbit Group bitmask
   * \return 1 if overlap detected, 0 otherwise */
  int check_rigid_group_overlap(int);

  /** Check for rigid body overlap with a region
   * \param groupbit Group bitmask
   * \param region Pointer to Region
   * \return 1 if overlap detected, 0 otherwise */
  int check_rigid_region_overlap(int, class Region *);

  /** Check for rigid body overlap with a list of atoms
   * \param list List of atom indices
   * \return 1 if overlap detected, 0 otherwise */
  int check_rigid_list_overlap(int *);

  /** Write fix info to restart file
   * \param fp File pointer for restart file */
  void write_restart(FILE *);

  /** Read fix info from restart file
   * \param fp File pointer for restart file
   * \return Number of fixes read */
  int read_restart(FILE *);

  /** Deallocate restart storage
   * \param flag Flag indicating what to deallocate */
  void restart_deallocate(int);

  /** Calculate memory usage
   * \return Memory usage in bytes */
  double memory_usage();

 protected:
  // internal fix counts

  int n_post_force;           /**< Number of fixes for post_force */
  int n_post_force_group;     /**< Number of fixes for post_force by group */
  int n_post_force_respa;     /**< Number of fixes for rRESPA post_force */

  // lists of fixes to apply at different stages of timestep

  int *list_initial_integrate;    /**< Fix indices for initial_integrate */
  int *list_post_integrate;       /**< Fix indices for post_integrate */
  int *list_pre_exchange;         /**< Fix indices for pre_exchange */
  int *list_pre_neighbor;         /**< Fix indices for pre_neighbor */
  int *list_post_neighbor;        /**< Fix indices for post_neighbor */
  int *list_pre_force;            /**< Fix indices for pre_force */
  int *list_pre_reverse;          /**< Fix indices for pre_reverse */
  int *list_post_force;           /**< Fix indices for post_force */
  int *list_post_force_group;     /**< Fix indices for post_force by group */
  int *list_final_integrate;      /**< Fix indices for final_integrate */
  int *list_end_of_step;          /**< Fix indices for end_of_step */
  int *list_energy_couple;        /**< Fix indices for coupled energy */
  int *list_energy_global;        /**< Fix indices for global energy */
  int *list_energy_atom;          /**< Fix indices for per-atom energy */
  int *list_initial_integrate_respa;    /**< Fix indices for rRESPA initial_integrate */
  int *list_post_integrate_respa;       /**< Fix indices for rRESPA post_integrate */
  int *list_pre_force_respa;            /**< Fix indices for rRESPA pre_force */
  int *list_post_force_respa;           /**< Fix indices for rRESPA post_force */
  int *list_final_integrate_respa;      /**< Fix indices for rRESPA final_integrate */
  int *list_min_pre_exchange;     /**< Fix indices for minimization pre_exchange */
  int *list_min_pre_neighbor;     /**< Fix indices for minimization pre_neighbor */
  int *list_min_post_neighbor;    /**< Fix indices for minimization post_neighbor */
  int *list_min_pre_force;        /**< Fix indices for minimization pre_force */
  int *list_min_pre_reverse;      /**< Fix indices for minimization pre_reverse */
  int *list_min_post_force;       /**< Fix indices for minimization post_force */
  int *list_min_energy;           /**< Fix indices for minimization energy */

  int *end_of_step_every;    /**< Frequency of end_of_step calls for each fix */

  int n_timeflag;       /**< Number of computes with time invocation tracking */
  int *list_timeflag;   /**< List of compute indices with time tracking */

  char **id_restart_global;       /**< IDs of fixes with global restart info */
  char **style_restart_global;    /**< Styles of fixes with global restart info */
  char **state_restart_global;    /**< State strings for global restart */
  int *used_restart_global;       /**< Flags for used global restart info */

  char **id_restart_peratom;       /**< IDs of fixes with per-atom restart info */
  char **style_restart_peratom;    /**< Styles of fixes with per-atom restart info */
  int *index_restart_peratom;      /**< Indices for per-atom restart info */
  int *used_restart_peratom;       /**< Flags for used per-atom restart info */

  int index_permanent;    /**< Fix/compute index for library access */

  // vectors to be used for new-API accessors as wrapper
  std::vector<Fix *> fix_list;          /**< Vector wrapper for fix array */
  std::vector<Compute *> compute_list;  /**< Vector wrapper for compute array */

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
  using ComputeCreator = Compute *(*) (LAMMPS *, int, char **);     /**< Function pointer type for compute creation */
  using ComputeCreatorMap = std::map<std::string, ComputeCreator>;  /**< Map type for compute style factory */
  ComputeCreatorMap *compute_map;    /**< Factory map for creating compute styles */

  using FixCreator = Fix *(*) (LAMMPS *, int, char **);    /**< Function pointer type for fix creation */
  using FixCreatorMap = std::map<std::string, FixCreator>; /**< Map type for fix style factory */
  FixCreatorMap *fix_map;    /**< Factory map for creating fix styles */

 protected:
  void create_factories();
};

}    // namespace LAMMPS_NS

#endif
