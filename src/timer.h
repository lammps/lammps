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

#ifndef LMP_TIMER_H
#define LMP_TIMER_H

#include "pointers.h"

namespace LAMMPS_NS {

/** \class LAMMPS_NS::Timer
 * \brief CPU and wall-clock timing for LAMMPS simulations
 *
 * The Timer class provides performance timing for various stages of a LAMMPS
 * simulation. It tracks both CPU time and wall-clock time for different
 * components of the simulation loop.
 *
 * Key features include:
 *
 * - Timing of pair, bond, kspace, neighbor, communication operations
 * - Multiple levels of detail (off, loop, normal, full)
 * - Optional MPI synchronization before timing
 * - Timeout support for limiting simulation wall time
 * - Separate timing for replica-exchange and NEB operations
 *
 * Timer levels control the amount of timing overhead:
 * - OFF: No timing
 * - LOOP: Only total loop time
 * - NORMAL: Major simulation components
 * - FULL: All detailed timing including sync points
 *
 * \sa LAMMPS_NS::Pointers */
class Timer : protected Pointers {
 public:
  /** Timer type identifiers for different simulation stages */
  enum ttype {
    RESET = -2,      /**< Reset all timers */
    START = -1,      /**< Start timing */
    TOTAL = 0,       /**< Total simulation time */
    PAIR,            /**< Pair force computation */
    BOND,            /**< Bond/angle/dihedral/improper force computation */
    KSPACE,          /**< Long-range (kspace) force computation */
    NEIGH,           /**< Neighbor list construction */
    COMM,            /**< MPI communication */
    MODIFY,          /**< Fix and compute operations */
    OUTPUT,          /**< Output (thermo, dump, restart) */
    SYNC,            /**< MPI synchronization */
    ALL,             /**< Total for all-inclusive timing */
    DEPHASE,         /**< Dephasing time for hyper dynamics */
    DYNAMICS,        /**< MD dynamics time for hyper dynamics */
    QUENCH,          /**< Quench time for hyper dynamics */
    NEB,             /**< NEB (nudged elastic band) time */
    REPCOMM,         /**< Replica communication time */
    REPOUT,          /**< Replica output time */
    NUM_TIMER        /**< Number of timer types */
  };

  /** Timer detail levels */
  enum tlevel {
    OFF = 0,    /**< No timing */
    LOOP,       /**< Only loop timing */
    NORMAL,     /**< Normal detail level */
    FULL,       /**< Full detail with sync */
    NUMLVL      /**< Number of levels */
  };

  /** Constructor for Timer class
   * \param lmp Pointer to the main LAMMPS instance */
  Timer(class LAMMPS *);

  /** Initialize timers before a run */
  void init();

  // inline function to reduce overhead if we want no detailed timings

  /** Record a timestamp for timing
   * \param which Timer type (default START to begin timing) */
  void stamp(enum ttype which = START)
  {
    if (_level > LOOP) _stamp(which);
  }

  /** Start barrier timing (with MPI sync) */
  void barrier_start();

  /** Stop barrier timing */
  void barrier_stop();

  // accessor methods for supported level of detail

  /** Check if loop timing is enabled
   * \return true if level >= LOOP */
  bool has_loop() const { return (_level >= LOOP); }

  /** Check if normal timing is enabled
   * \return true if level >= NORMAL */
  bool has_normal() const { return (_level >= NORMAL); }

  /** Check if full timing is enabled
   * \return true if level >= FULL */
  bool has_full() const { return (_level >= FULL); }

  /** Check if synchronization is enabled
   * \return true if sync is enabled */
  bool has_sync() const { return (_sync != OFF); }

  /** Check if timeout is enabled
   * \return true if timeout is set */
  bool has_timeout() const { return (_timeout >= 0.0); }

  /** Check if timeout has expired
   * \return true if timeout == 0.0 (expired) */
  bool is_timeout() const { return (_timeout == 0.0); }

  /** Get elapsed wall-clock time for a timer
   * \param which Timer type
   * \return Elapsed wall-clock time in seconds */
  double elapsed(enum ttype);

  /** Get elapsed CPU time for a timer
   * \param which Timer type
   * \return Elapsed CPU time in seconds */
  double cpu(enum ttype);

  /** Get stored CPU time for a timer
   * \param which Timer type
   * \return CPU time in seconds */
  double get_cpu(enum ttype which) const { return cpu_array[which]; };

  /** Get stored wall-clock time for a timer
   * \param which Timer type
   * \return Wall-clock time in seconds */
  double get_wall(enum ttype which) const { return wall_array[which]; };

  /** Set wall-clock time for a timer
   * \param which Timer type
   * \param value Time value in seconds */
  void set_wall(enum ttype, double);

  /** Initialize timeout timer from command-line or environment */
  void init_timeout();

  /** Force an immediate timeout */
  void force_timeout() { _timeout = 0.0; }

  /** Restore original timeout after forced timeout */
  void reset_timeout() { _timeout = _s_timeout; }

  /** Get remaining time before timeout
   * \return Seconds remaining (0.0 if inactive, negative if expired) */
  double get_timeout_remain();

  /** Print timeout message to file
   * \param fp File pointer for output */
  void print_timeout(FILE *);

  // check for timeout. inline wrapper around internal
  // function to reduce overhead in case there is no check.

  /** Check if timeout has occurred
   * \param step Current timestep
   * \return true if timeout has occurred */
  bool check_timeout(int step)
  {
    if (_timeout == 0.0) return true;
    if (_nextcheck != step)
      return false;
    else
      return _check_timeout();
  }

  /** Modify timer parameters
   * \param narg Number of arguments
   * \param arg Argument array */
  void modify_params(int, char **);

 private:
  double cpu_array[NUM_TIMER];     /**< CPU time array for each timer type */
  double wall_array[NUM_TIMER];    /**< Wall-clock time array for each timer type */
  double previous_cpu;             /**< Previous CPU time for delta calculation */
  double previous_wall;            /**< Previous wall time for delta calculation */
  double timeout_start;            /**< Wall time when timeout started */
  double _timeout;                 /**< Max allowed wall time (negative = infinite) */
  double _s_timeout;               /**< Saved timeout for restoration after force */
  int _level;                      /**< Timer detail level (OFF, LOOP, NORMAL, FULL) */
  int _sync;                       /**< Sync flag (non-zero = synchronize before timing) */
  int _checkfreq;                  /**< Frequency of timeout checking */
  int _nextcheck;                  /**< Next step to check timeout */

  // update one specific timer array
  void _stamp(enum ttype);

  // check for timeout
  bool _check_timeout();
};

}    // namespace LAMMPS_NS

#endif
