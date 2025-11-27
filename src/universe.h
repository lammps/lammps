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

#ifndef LMP_UNIVERSE_H
#define LMP_UNIVERSE_H

#include "pointers.h"

namespace LAMMPS_NS {

/** \class LAMMPS_NS::Universe
 * \brief Class for managing multi-partition simulations in LAMMPS
 *
 * The Universe class manages the parallel environment for LAMMPS simulations,
 * particularly when running multi-partition jobs (e.g., replica exchange,
 * parallel tempering, or nudged elastic band calculations). It handles:
 *
 * - The overall MPI communicator spanning all processors (universe)
 * - Division of processors into separate partitions (worlds)
 * - Communication between partitions
 * - Screen and log file output for the universe
 *
 * When LAMMPS is run with the `-partition` command-line flag, the Universe
 * class divides the total MPI communicator into separate worlds, each running
 * an independent simulation. These worlds can communicate with each other
 * through the universe communicator for coordinated multi-replica methods.
 *
 * For single-partition runs (the common case), the universe communicator
 * is equivalent to the world communicator, and nworlds = 1.
 *
 * \sa LAMMPS_NS::Pointers, LAMMPS_NS::Comm */
class Universe : protected Pointers {
 public:
  MPI_Comm uworld;    /**< MPI communicator for entire universe (all processors) */
  int me;             /**< Rank of this processor in the universe communicator */
  int nprocs;         /**< Total number of processors in the universe */

  FILE *uscreen;      /**< File pointer for universe screen output (or NULL) */
  FILE *ulogfile;     /**< File pointer for universe log file (or NULL) */

  int existflag;          /**< 1 if universe exists due to -partition flag, 0 otherwise */
  int nworlds;            /**< Number of separate worlds (partitions) in the universe */
  int iworld;             /**< Which world (partition) this processor belongs to (0 to nworlds-1) */
  int *procs_per_world;   /**< Array of processor counts for each world [nworlds] */
  int *root_proc;         /**< Array of root processor ranks for each world [nworlds] */

  MPI_Comm uorig;         /**< Original MPI communicator passed to LAMMPS instance */
  int *uni2orig;          /**< Mapping from universe rank to original communicator rank */

  /** Constructor for Universe class
   * \param lmp Pointer to the main LAMMPS instance
   * \param communicator MPI communicator for the universe */
  Universe(class LAMMPS *, MPI_Comm);

  /** Destructor - cleans up allocated memory and closes universe files */
  ~Universe() override;

  /** Reorder processor assignment among partitions
   * \param style Style of reordering (e.g., "nth" or "custom")
   * \param arg Argument for the reordering style */
  void reorder(char *, char *);

  /** Add a new world (partition) to the universe
   * \param arg String specifying number of processors for the new world */
  void add_world(char *);

  /** Check if all partitions have consistent settings
   * \return 1 if consistent, 0 if inconsistent */
  int consistent();
};

}    // namespace LAMMPS_NS

#endif
