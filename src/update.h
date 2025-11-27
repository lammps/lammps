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

#ifndef LMP_UPDATE_H
#define LMP_UPDATE_H

#include "pointers.h"

#include <map>

namespace LAMMPS_NS {

/** \class LAMMPS_NS::Update
 * \brief Time integration and minimization management for LAMMPS
 *
 * The Update class manages time integration (molecular dynamics) and energy
 * minimization in LAMMPS. It holds the simulation timestep and tracks the
 * current step number. Key responsibilities include:
 *
 * - Managing the timestep size (dt) and current timestep number (ntimestep)
 * - Creating and managing time integrators (e.g., Verlet, rRESPA)
 * - Creating and managing energy minimizers (e.g., CG, FIRE, SD)
 * - Tracking simulation time and run boundaries
 * - Managing energy and virial tallying flags
 * - Handling unit style settings that affect the timestep
 *
 * The class maintains factories (style maps) for creating integrator and
 * minimizer styles dynamically based on user input.
 *
 * \sa LAMMPS_NS::Integrate, LAMMPS_NS::Min, LAMMPS_NS::Pointers */
class Update : protected Pointers {
 public:
  /** Constructor for Update class
   * \param lmp Pointer to the main LAMMPS instance */
  Update(class LAMMPS *);

  /** Destructor - cleans up integrator, minimizer, and style maps */
  ~Update() override;

  /** Initialize the update before a run */
  void init();

  /** Set unit style and adjust timestep accordingly
   * \param style Name of the unit style (e.g., "real", "metal", "lj") */
  void set_units(const char *);

  /** Create a new time integrator
   * \param narg Number of arguments
   * \param arg Argument array (style name and options)
   * \param trysuffix Flag to try accelerated suffix versions */
  void create_integrate(int, char **, int);

  /** Create a new energy minimizer
   * \param narg Number of arguments
   * \param arg Argument array (style name and options)
   * \param trysuffix Flag to try accelerated suffix versions */
  void create_minimize(int, char **, int);

  /** Reset the timestep counter
   * \param narg Number of arguments
   * \param arg Argument array containing new timestep value */
  void reset_timestep(int, char **);

  /** Reset the timestep counter to a specific value
   * \param newstep New timestep value
   * \param update_all Flag to update all timestep-dependent values */
  void reset_timestep(bigint, bool);

  /** Update the simulation time based on current timestep */
  void update_time();

  /** Calculate memory usage of the Update class
   * \return Memory usage in bytes */
  double memory_usage();

  double dt;              /**< Timestep size for dynamics */
  double etol;            /**< Energy tolerance for minimization convergence */
  double ftol;            /**< Force tolerance for minimization convergence */
  bigint ntimestep;       /**< Current timestep number (dynamics or minimization iteration) */
  int nsteps;             /**< Number of steps/iterations to run */
  int whichflag;          /**< Run mode: 0 = unset, 1 = dynamics, 2 = minimization */
  double atime;           /**< Simulation time at atimestep */
  bigint atimestep;       /**< Timestep at which atime was last updated */
  bigint firststep;       /**< First step of current run */
  bigint laststep;        /**< Last step of current run */
  bigint beginstep;       /**< First step of multiple concatenated runs */
  bigint endstep;         /**< Last step of multiple concatenated runs */
  int first_update;       /**< 0 before first update, 1 after (used for setup) */
  int max_eval;           /**< Maximum force evaluations for minimizer */
  int restrict_output;    /**< 1 if output should not write dump/restart files */
  int setupflag;          /**< 1 when setup() is computing forces, 0 otherwise */
  int multireplica;       /**< 1 if minimization across replicas, 0 otherwise */
  int dt_default;         /**< 1 if dt is at default value for unit style, 0 if user-set */

  bigint eflag_global;    /**< Timestep when global energy was last tallied */
  bigint eflag_atom;      /**< Timestep when per-atom energy was last tallied */
  bigint vflag_global;    /**< Timestep when global virial was last tallied */
  bigint vflag_atom;      /**< Timestep when per-atom virial was last tallied */

  char *unit_style;       /**< Current unit style name */

  class Integrate *integrate;    /**< Pointer to current time integrator instance */
  char *integrate_style;         /**< Current integrator style name */

  class Min *minimize;           /**< Pointer to current minimizer instance */
  char *minimize_style;          /**< Current minimizer style name */

  using IntegrateCreator = Integrate *(*) (LAMMPS *, int, char **);    /**< Function pointer type for integrator creation */
  using MinimizeCreator = Min *(*) (LAMMPS *);                         /**< Function pointer type for minimizer creation */

  using IntegrateCreatorMap = std::map<std::string, IntegrateCreator>;    /**< Map type for integrator style factory */
  using MinimizeCreatorMap = std::map<std::string, MinimizeCreator>;      /**< Map type for minimizer style factory */

  IntegrateCreatorMap *integrate_map;    /**< Factory map for creating integrator styles */
  MinimizeCreatorMap *minimize_map;      /**< Factory map for creating minimizer styles */

 private:
  void new_integrate(char *, int, char **, int, int &);
  void new_minimize(char *, int, char **, int, int &);
};

}    // namespace LAMMPS_NS

#endif
