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
  bigint next;    // next timestep for any kind of output

  bigint next_thermo;      // next timestep for thermo output
  int thermo_every;        // output freq for thermo, 0 if first/last only
  bigint last_thermo;      // last timestep thermo was output
  char *var_thermo;        // variable name for thermo freq, null pointer if every
  int ivar_thermo;         // variable index for thermo frequency
  class Thermo *thermo;    // Thermodynamic computations

  int ndump;                    // # of Dumps defined
  int max_dump;                 // max size of Dump list
  bigint next_dump_any;         // next timestep for any dump
  bigint next_time_dump_any;    // next timestep for any time dump with computes
  int any_time_dumps;           // 1 if any time dump defined
  int *mode_dump;               // 0/1 if write every N timesteps or Delta in sim time
  int *every_dump;              // dump every N timesteps, 0 if variable
  double *every_time_dump;      // dump every Delta of sim time, 0.0 if variable
  bigint *next_dump;            // next timestep to perform dump
  double *next_time_dump;       // next simulation time to perform dump (mode = 1)
  bigint *last_dump;            // last timestep each snapshot was output
  char **var_dump;              // variable name for next dump (steps or sim time)
  int *ivar_dump;               // variable index of var_dump name
  Dump **dump;                  // list of defined Dumps

  int restart_flag;               // 1 if any restart files are written
  int restart_flag_single;        // 1 if single restart files are written
  int restart_flag_double;        // 1 if double restart files are written
  bigint next_restart;            // next timestep to write any restart file
  bigint next_restart_single;     // next timestep to write a single restart file
  bigint next_restart_double;     // next timestep to write a double restart file
  int restart_every_single;       // single restart file write freq, 0 if var
  int restart_every_double;       // double restart file write freq, 0 if var
  bigint last_restart;            // last timestep any restart file was output
  int restart_toggle;             // 0 if use restart2a as prefix, 1 if restart2b
  char *var_restart_single;       // variable name for single restart freq
  char *var_restart_double;       // variable name for double restart freq
  int ivar_restart_single;        // index of var_restart_single
  int ivar_restart_double;        // index of var_restart_double
  char *restart1;                 // name single restart file
  char *restart2a, *restart2b;    // names of double restart files
  class WriteRestart *restart;    // class for writing restart files

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
  [[nodiscard]] Dump *get_dump_by_id(const std::string &) const;

  /** Get a dump by its index
   * \param idx Index in dump array
   * \return Pointer to Dump, or nullptr if out of range */
  [[nodiscard]] Dump *get_dump_by_index(int idx) const
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
  std::vector<Dump *> dump_list;
  void calculate_next_dump(int, int, bigint);
};

}    // namespace LAMMPS_NS

#endif
