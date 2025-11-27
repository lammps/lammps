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

#ifndef LMP_OUTPUT_H
#define LMP_OUTPUT_H

#include "pointers.h"

#include <map>

namespace LAMMPS_NS {

class Dump;
struct json_metadata;

/** \class LAMMPS_NS::Output
 * \brief Output management for LAMMPS simulations
 *
 * The Output class manages all forms of output from LAMMPS simulations,
 * including thermodynamic data, dump files, and restart files. It coordinates
 * the timing and execution of output operations during runs.
 *
 * Key responsibilities include:
 *
 * - Managing thermodynamic output (thermo) and its frequency
 * - Managing multiple dump files with different output frequencies
 * - Managing restart file output with single or double file modes
 * - Coordinating output timing to determine next output timestep
 * - Supporting variable-based output frequencies
 * - Providing a factory for creating dump styles
 *
 * The class tracks the next timestep for each type of output and provides
 * methods to force immediate output of dumps or restart files.
 *
 * \sa LAMMPS_NS::Dump, LAMMPS_NS::Thermo, LAMMPS_NS::WriteRestart */
class Output : protected Pointers {
 public:
  bigint next;    /**< Next timestep for any kind of output */

  bigint next_thermo;      /**< Next timestep for thermodynamic output */
  int thermo_every;        /**< Output frequency for thermo (0 = first/last only) */
  bigint last_thermo;      /**< Last timestep thermo was output */
  char *var_thermo;        /**< Variable name for thermo frequency (null if fixed) */
  int ivar_thermo;         /**< Variable index for thermo frequency */
  class Thermo *thermo;    /**< Pointer to Thermo instance */

  int ndump;                    /**< Number of defined dumps */
  int max_dump;                 /**< Allocated size of dump arrays */
  bigint next_dump_any;         /**< Next timestep for any dump */
  bigint next_time_dump_any;    /**< Next timestep for any time-based dump with computes */
  int any_time_dumps;           /**< 1 if any time-based dump defined, 0 otherwise */
  int *mode_dump;               /**< 0 = every N timesteps, 1 = every Delta simulation time */
  int *every_dump;              /**< Dump frequency in timesteps (0 = variable) */
  double *every_time_dump;      /**< Dump frequency in simulation time (0.0 = variable) */
  bigint *next_dump;            /**< Next timestep for each dump */
  double *next_time_dump;       /**< Next simulation time for each dump (mode = 1) */
  bigint *last_dump;            /**< Last timestep each dump was output */
  char **var_dump;              /**< Variable names for dump frequencies */
  int *ivar_dump;               /**< Variable indices for dump frequencies */
  Dump **dump;                  /**< Array of Dump instance pointers */

  int restart_flag;               /**< 1 if any restart files written, 0 otherwise */
  int restart_flag_single;        /**< 1 if single restart files written */
  int restart_flag_double;        /**< 1 if double restart files written */
  bigint next_restart;            /**< Next timestep for any restart output */
  bigint next_restart_single;     /**< Next timestep for single restart file */
  bigint next_restart_double;     /**< Next timestep for double restart file */
  int restart_every_single;       /**< Single restart frequency (0 = variable) */
  int restart_every_double;       /**< Double restart frequency (0 = variable) */
  bigint last_restart;            /**< Last timestep any restart was output */
  int restart_toggle;             /**< 0 = use restart2a, 1 = use restart2b */
  char *var_restart_single;       /**< Variable name for single restart frequency */
  char *var_restart_double;       /**< Variable name for double restart frequency */
  int ivar_restart_single;        /**< Variable index for single restart frequency */
  int ivar_restart_double;        /**< Variable index for double restart frequency */
  char *restart1;                 /**< Filename for single restart file */
  char *restart2a, *restart2b;    /**< Filenames for double restart files */
  class WriteRestart *restart;    /**< Pointer to WriteRestart instance */

  using DumpCreator = Dump *(*) (LAMMPS *, int, char **);     /**< Function pointer type for dump creation */
  using DumpCreatorMap = std::map<std::string, DumpCreator>;  /**< Map type for dump style factory */
  DumpCreatorMap *dump_map;    /**< Factory map for dump styles */

  /** Create MPI datatype for particle data
   * \return MPI_Datatype for particle structure */
  MPI_Datatype createParticleStructType();

  /** Constructor for Output class
   * \param lmp Pointer to the main LAMMPS instance */
  Output(class LAMMPS *);

  /** Destructor - cleans up thermo, dumps, restart, and style map */
  ~Output() override;

  /** Initialize output before a run */
  void init();

  /** Setup output at start of run
   * \param memflag 1 to print memory usage, 0 to skip (default 1) */
  void setup(int memflag = 1);

  /** Perform all output for current timestep
   * \param ntimestep Current timestep */
  void write(bigint);

  /** Force output of all dump files
   * \param ntimestep Current timestep */
  void write_dump(bigint);

  /** Force output of restart file
   * \param ntimestep Current timestep */
  void write_restart(bigint);

  /** Write molecule data in JSON format
   * \param fp File pointer
   * \param nmolecules Number of molecules
   * \param natoms Number of atoms
   * \param types Atom types array
   * \param metadata JSON metadata structure */
  void write_molecule_json(FILE *, int, int, int *, json_metadata *);

  /** Reset output after timestep is reset
   * \param ntimestep New timestep value */
  void reset_timestep(bigint);

  /** Reset output after timestep size (dt) changes */
  void reset_dt();

  /** Add a new dump
   * \param narg Number of arguments
   * \param arg Argument array
   * \return Pointer to new Dump instance */
  Dump *add_dump(int, char **);

  /** Modify settings of an existing dump
   * \param narg Number of arguments
   * \param arg Argument array */
  void modify_dump(int, char **);

  /** Delete a dump by ID
   * \param id Dump ID to delete */
  void delete_dump(const std::string &);

  /** Get a dump by its ID
   * \param id Dump ID
   * \return Pointer to Dump, or nullptr if not found */
  Dump *get_dump_by_id(const std::string &) const;

  /** Get a dump by its index
   * \param idx Index in dump array
   * \return Pointer to Dump, or nullptr if out of range */
  Dump *get_dump_by_index(int idx) const
  {
    return ((idx >= 0) && (idx < ndump)) ? dump[idx] : nullptr;
  }

  /** Get list of all dumps
   * \return Reference to vector of Dump pointers */
  const std::vector<Dump *> &get_dump_list();

  /** Check if any time-based dump outputs this timestep
   * \param ntimestep Current timestep
   * \return 1 if output needed, 0 otherwise */
  int check_time_dumps(bigint);

  /** Set thermodynamic output frequency
   * \param narg Number of arguments
   * \param arg Argument array */
  void set_thermo(int, char **);

  /** Create a new thermo style
   * \param narg Number of arguments
   * \param arg Argument array */
  void create_thermo(int, char **);

  /** Create restart file output
   * \param narg Number of arguments
   * \param arg Argument array */
  void create_restart(int, char **);

  /** Print memory usage information */
  void memory_usage();

 private:
  std::vector<Dump *> dump_list;    /**< Vector wrapper for dump array */
  void calculate_next_dump(int, int, bigint);
};

}    // namespace LAMMPS_NS

#endif
