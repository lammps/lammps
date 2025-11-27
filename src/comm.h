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

#ifndef LMP_COMM_H
#define LMP_COMM_H

#include "pointers.h"    // IWYU pragma: export

namespace LAMMPS_NS {

/** \class LAMMPS_NS::Comm
 * \brief Base class for inter-processor communication in LAMMPS
 *
 * The Comm class manages all inter-processor communication for parallel
 * LAMMPS simulations. It handles the decomposition of the simulation box
 * across MPI processes and the exchange of atom data between processors.
 *
 * Key responsibilities include:
 *
 * - Domain decomposition across MPI processes (brick or tiled layouts)
 * - Forward communication of atom coordinates to ghost atoms
 * - Reverse communication of forces from ghost atoms to owning processors
 * - Exchange of atoms that migrate between processor domains
 * - Setting up communication patterns for ghost atom acquisition
 * - Communication for Pair, Bond, Fix, Compute, and Dump styles
 *
 * This is an abstract base class with two concrete implementations:
 * - CommBrick: Traditional 6-way stencil communication for regular grids
 * - CommTiled: Communication for irregular domain decompositions (from RCB)
 *
 * \sa LAMMPS_NS::CommBrick, LAMMPS_NS::CommTiled, LAMMPS_NS::Pointers */
class Comm : protected Pointers {
 public:
  /** Communication style: BRICK = regular grid, TILED = irregular domains */
  enum { BRICK, TILED };
  int style;    /**< Communication style: BRICK or TILED */

  /** Domain layout types */
  enum { LAYOUT_UNIFORM, LAYOUT_NONUNIFORM, LAYOUT_TILED };
  int layout;    /**< Domain layout: UNIFORM, NONUNIFORM, or TILED */

  /** Communication mode for neighbor cutoffs */
  enum { SINGLE, MULTI };
  int mode;    /**< Cutoff mode: SINGLE or MULTI (per-collection) */

  /** Standard buffer size for fixed-size per-atom data */
  enum { BUFEXTRA = 1024 };

  int me;                     /**< MPI rank of this processor in world communicator */
  int nprocs;                 /**< Total number of MPI processes */
  int ghost_velocity;         /**< 1 if ghost atoms have velocity, 0 if not */
  double cutghost[3];         /**< Ghost cutoff distances in each dimension */
  double cutghostuser;        /**< User-specified ghost cutoff (mode == SINGLE) */
  double *cutusermulti;       /**< Per-collection user ghost cutoffs (mode == MULTI) */
  int ncollections;           /**< Number of atom collections known to comm */
  int ncollections_cutoff;    /**< Number of collections with stored cutoffs */
  int recv_from_partition;    /**< Partition to receive layout from (-1 if none) */
  int send_to_partition;      /**< Partition to send layout to (-1 if none) */
  int other_partition_style;  /**< Style of layout dependency on other partition */

  int nthreads;    /**< Number of OpenMP threads per MPI process */

  // public settings specific to layout = UNIFORM, NONUNIFORM

  int procgrid[3];          /**< Processor count in each dimension of 3D grid */
  int user_procgrid[3];     /**< User-requested processor counts per dimension */
  int myloc[3];             /**< This processor's location in 3D grid (0 to N-1) */
  int procneigh[3][2];      /**< 6 neighboring processor ranks (left/right in x,y,z) */
  double *xsplit;           /**< Fractional (0-1) x-boundaries of subdomains */
  double *ysplit;           /**< Fractional (0-1) y-boundaries of subdomains */
  double *zsplit;           /**< Fractional (0-1) z-boundaries of subdomains */
  int ***grid2proc;         /**< Mapping from 3D grid location to processor rank */

  // public settings specific to layout = TILED

  int rcbnew;              /**< 1 if just reset by rebalance, 0 otherwise */
  double mysplit[3][2];    /**< Fractional (0-1) bounds of this processor's subdomain */
  double rcbcutfrac;       /**< Fractional RCB cut position */
  int rcbcutdim;           /**< Dimension of RCB cut */

  // methods

  /** Constructor for Comm class
   * \param lmp Pointer to the main LAMMPS instance */
  Comm(class LAMMPS *);

  /** Destructor - cleans up communication buffers and data structures */
  ~Comm() override;

  // NOTE: copy_arrays is called from a constructor and must not be made virtual

  /** Copy arrays from another Comm instance
   * \param oldcomm Comm instance to copy from */
  void copy_arrays(class Comm *);

  /** Initialize communication before a run */
  virtual void init();

  /** Modify communication parameters
   * \param narg Number of arguments
   * \param arg Argument array */
  void modify_params(int, char **);

  /** Set 3D processor grid attributes
   * \param narg Number of arguments
   * \param arg Argument array */
  void set_processors(int, char **);

  /** Setup the 3D processor grid
   * \param outflag 1 to output grid info, 0 for silent (default 1) */
  virtual void set_proc_grid(int outflag = 1);

  /** Determine the communication cutoff distance
   * \return Communication cutoff in distance units */
  double get_comm_cutoff();

  /** Setup 3D communication pattern (pure virtual)
   * Must be implemented by derived classes */
  virtual void setup() = 0;

  /** Forward communication of atom coordinates (pure virtual)
   * \param dummy Unused parameter for interface compatibility */
  virtual void forward_comm(int dummy = 0) = 0;

  /** Reverse communication of forces (pure virtual) */
  virtual void reverse_comm() = 0;

  /** Exchange atoms that have moved to other processors (pure virtual) */
  virtual void exchange() = 0;

  /** Setup list of atoms to communicate as ghosts (pure virtual) */
  virtual void borders() = 0;

  // forward/reverse comm from a Pair, Bond, Fix, Compute, Dump

  /** Forward communication for Pair styles
   * \param pair Pointer to Pair instance
   * \param size Number of values per atom (0 = default from pair) */
  virtual void forward_comm(class Pair *, int size = 0) = 0;

  /** Reverse communication for Pair styles
   * \param pair Pointer to Pair instance
   * \param size Number of values per atom (0 = default from pair) */
  virtual void reverse_comm(class Pair *, int size = 0) = 0;

  /** Forward communication for Bond styles
   * \param bond Pointer to Bond instance
   * \param size Number of values per atom */
  virtual void forward_comm(class Bond *, int size = 0) = 0;

  /** Reverse communication for Bond styles
   * \param bond Pointer to Bond instance
   * \param size Number of values per atom */
  virtual void reverse_comm(class Bond *, int size = 0) = 0;

  /** Forward communication for Fix styles
   * \param fix Pointer to Fix instance
   * \param size Number of values per atom */
  virtual void forward_comm(class Fix *, int size = 0) = 0;

  /** Reverse communication for Fix styles
   * \param fix Pointer to Fix instance
   * \param size Number of values per atom */
  virtual void reverse_comm(class Fix *, int size = 0) = 0;

  /** Reverse communication for Fix with variable-size data
   * \param fix Pointer to Fix instance */
  virtual void reverse_comm_variable(class Fix *) = 0;

  /** Forward communication for Compute styles
   * \param compute Pointer to Compute instance
   * \param size Number of values per atom */
  virtual void forward_comm(class Compute *, int size = 0) = 0;

  /** Reverse communication for Compute styles
   * \param compute Pointer to Compute instance
   * \param size Number of values per atom */
  virtual void reverse_comm(class Compute *, int size = 0) = 0;

  /** Forward communication for Dump styles
   * \param dump Pointer to Dump instance
   * \param size Number of values per atom */
  virtual void forward_comm(class Dump *, int size = 0) = 0;

  /** Reverse communication for Dump styles
   * \param dump Pointer to Dump instance
   * \param size Number of values per atom */
  virtual void reverse_comm(class Dump *, int size = 0) = 0;

  // forward comm of an array

  /** Forward communication of an array
   * \param ncol Number of columns in array
   * \param array 2D array to communicate */
  virtual void forward_comm_array(int, double **) = 0;

  // map a point to a processor, based on current decomposition

  /** Setup for coordinate-to-processor mapping */
  virtual void coord2proc_setup() {}

  /** Map a point to the owning processor
   * \param x Coordinates of the point
   * \param igx Output: grid x-index
   * \param igy Output: grid y-index
   * \param igz Output: grid z-index
   * \return Rank of owning processor */
  virtual int coord2proc(double *, int &, int &, int &);

  // memory usage

  /** Calculate memory usage
   * \return Memory usage in bytes */
  virtual double memory_usage() = 0;

  // non-virtual functions common to all Comm styles

  /** Ring communication pattern
   * \param n Number of data items
   * \param nper Size of each data item
   * \param inbuf Input buffer
   * \param messtag Message tag
   * \param callback Function to call for processing
   * \param ptr User pointer for callback
   * \param outbuf Output buffer
   * \param self 1 to include self in ring (default 1) */
  void ring(int, int, void *, int, void (*)(int, char *, void *), void *, void *, int self = 1);

  /** Rendezvous communication pattern
   * \param mode Communication mode (1 or 2)
   * \param n Number of data items
   * \param inbuf Input buffer
   * \param insize Size of input data per item
   * \param inorder 1 if input is ordered, 0 otherwise
   * \param procs Array of destination processors
   * \param callback Function to call at rendezvous
   * \param outorder 1 if output should be ordered, 0 otherwise
   * \param outbuf Output buffer
   * \param outsize Size of output data per item
   * \param ptr User pointer for callback
   * \param statflag Statistics flag (default 0)
   * \return Number of output items */
  int rendezvous(int, int, char *, int, int, int *,
                 int (*)(int, char *, int &, int *&, char *&, void *), int, char *&, int, void *,
                 int statflag = 0);

  // extract data useful to other classes

  /** Extract internal data
   * \param name Name of data to extract
   * \param len Output: length of data
   * \return Pointer to data, or nullptr if not found */
  virtual void *extract(const char *, int &) { return nullptr; }

 protected:
  int bordergroup;    /**< Only communicate atoms in this group for borders() */

  int triclinic;              /**< 0 if orthogonal domain, 1 if triclinic */
  int map_style;              /**< Atom ID mapping style (non-0 if mapping enabled) */
  int comm_x_only;            /**< 1 if only exchange x in forward comm */
  int comm_f_only;            /**< 1 if only exchange f in reverse comm */

  int size_forward;    /**< Number of per-atom datums in forward comm */
  int size_reverse;    /**< Number of per-atom datums in reverse comm */
  int size_border;     /**< Number of per-atom datums in border comm */

  int maxforward;              /**< Max datums in forward comm buffer */
  int maxreverse;              /**< Max datums in reverse comm buffer */
  int maxexchange;             /**< Max size of one exchanged atom */
  int maxexchange_atom;        /**< Contribution to maxexchange from AtomVec */
  int maxexchange_fix;         /**< Static contribution to maxexchange from Fixes */
  int maxexchange_fix_dynamic; /**< 1 if a fix has dynamic exchange contribution */
  int bufextra;                /**< Augment send buffer size for exchange */
  int bufextra_max;            /**< Maximum bufextra value */

  int gridflag;        /**< Option for creating 3D processor grid */
  int mapflag;         /**< Option for mapping processors to grid */
  char xyz[4];         /**< xyz ordering for processor to grid mapping */
  char *customfile;    /**< File with custom processor map */
  char *outfile;       /**< File for processor grid/map output */
  int numa_nodes;      /**< Number of NUMA domains per socket */

  int otherflag;            /**< 1 if this partition depends on another */
  int other_style;          /**< Style of dependency on other partition */
  int other_procgrid[3];    /**< Processor layout of other partition */
  int other_coregrid[3];    /**< Core layout of other partition */
  int ncores;               /**< Number of cores per node */
  int coregrid[3];          /**< 3D grid of cores within a node */
  int user_coregrid[3];     /**< User-requested cores per dimension */
  int multi_reduce;         /**< 1 if multi cutoff is intra-collection cutoff */

  void init_exchange();
  int rendezvous_irregular(int, char *, int, int, int *,
                           int (*)(int, char *, int &, int *&, char *&, void *), int, char *&, int,
                           void *, int);
  int rendezvous_all2all(int, char *, int, int, int *,
                         int (*)(int, char *, int &, int *&, char *&, void *), int, char *&, int,
                         void *, int);
  void rendezvous_stats(int, int, int, int, int, int, bigint);

 public:
  enum { MULTIPLE };    /**< Flag for multiple send/recv operations */
};

}    // namespace LAMMPS_NS

#endif
