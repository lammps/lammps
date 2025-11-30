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

#ifndef LMP_COMM_BRICK_H
#define LMP_COMM_BRICK_H

#include "comm.h"

namespace LAMMPS_NS {

/** \class LAMMPS_NS::CommBrick
 * \brief Communication class for regular brick-style domain decomposition
 *
 * CommBrick implements inter-processor communication for LAMMPS simulations
 * using a traditional 6-way stencil communication pattern. This is the default
 * communication style and is optimal for uniform or non-uniform (but still
 * logically regular) domain decompositions.
 *
 * Key features include:
 *
 * - 6-way stencil communication in 3D (2 neighbors per dimension)
 * - Support for both uniform and non-uniform processor grids
 * - Efficient swap-based communication for ghost atom updates
 * - Multi-collection support for different cutoff distances
 *
 * This class is used when the simulation domain is decomposed into a regular
 * 3D grid of processor subdomains. Each processor exchanges atoms and data
 * only with its 6 nearest neighbors (or fewer in lower dimensions).
 *
 * \sa LAMMPS_NS::Comm, LAMMPS_NS::CommTiled */
class CommBrick : public Comm {
 public:
  /** Constructor for CommBrick class
   * \param lmp Pointer to the main LAMMPS instance */
  CommBrick(class LAMMPS *);

  /** Copy constructor from another Comm instance
   * \param lmp Pointer to the main LAMMPS instance
   * \param oldcomm Pointer to existing Comm instance to copy from */
  CommBrick(class LAMMPS *, class Comm *);

  /** Destructor - cleans up communication buffers and swap arrays */
  ~CommBrick() override;

  /** Initialize communication before a run */
  void init() override;

  /** Setup 3D communication pattern based on processor grid */
  void setup() override;

  /** Forward communication of atom coordinates to ghost atoms
   * \param dummy Unused parameter for interface compatibility */
  void forward_comm(int dummy = 0) override;

  /** Reverse communication of forces from ghost atoms to owners */
  void reverse_comm() override;

  /** Exchange atoms that have migrated to other processor domains */
  void exchange() override;

  /** Setup list of atoms to communicate as ghost atoms */
  void borders() override;

  /** Forward communication for Pair styles
   * \param pair Pointer to Pair instance
   * \param size Number of values per atom (0 = use default) */
  void forward_comm(class Pair *, int size = 0) override;

  /** Reverse communication for Pair styles
   * \param pair Pointer to Pair instance
   * \param size Number of values per atom (0 = use default) */
  void reverse_comm(class Pair *, int size = 0) override;

  /** Forward communication for Bond styles
   * \param bond Pointer to Bond instance
   * \param size Number of values per atom */
  void forward_comm(class Bond *, int size = 0) override;

  /** Reverse communication for Bond styles
   * \param bond Pointer to Bond instance
   * \param size Number of values per atom */
  void reverse_comm(class Bond *, int size = 0) override;

  /** Forward communication for Fix styles
   * \param fix Pointer to Fix instance
   * \param size Number of values per atom */
  void forward_comm(class Fix *, int size = 0) override;

  /** Reverse communication for Fix styles
   * \param fix Pointer to Fix instance
   * \param size Number of values per atom */
  void reverse_comm(class Fix *, int size = 0) override;

  /** Reverse communication with variable-size data for Fix styles
   * \param fix Pointer to Fix instance */
  void reverse_comm_variable(class Fix *) override;

  /** Forward communication for Compute styles
   * \param compute Pointer to Compute instance
   * \param size Number of values per atom */
  void forward_comm(class Compute *, int size = 0) override;

  /** Reverse communication for Compute styles
   * \param compute Pointer to Compute instance
   * \param size Number of values per atom */
  void reverse_comm(class Compute *, int size = 0) override;

  /** Forward communication for Dump styles
   * \param dump Pointer to Dump instance
   * \param size Number of values per atom */
  void forward_comm(class Dump *, int size = 0) override;

  /** Reverse communication for Dump styles
   * \param dump Pointer to Dump instance
   * \param size Number of values per atom */
  void reverse_comm(class Dump *, int size = 0) override;

  /** Forward communication of a 2D array
   * \param n Number of columns in the array
   * \param array Pointer to 2D array to communicate */
  void forward_comm_array(int, double **) override;

  /** Extract internal data
   * \param name Name of data to extract
   * \param dim Output: dimensionality of data
   * \return Pointer to requested data, or NULL if name not recognized */
  void *extract(const char *, int &) override;

  /** Calculate memory usage
   * \return Memory usage in bytes */
  double memory_usage() override;

 protected:
  int nswap;              /**< Number of swaps to perform = sum of maxneed */
  int recvneed[3][2];     /**< Number of procs away to receive atoms from per dim */
  int sendneed[3][2];     /**< Number of procs away to send atoms to per dim */
  int maxneed[3];         /**< Max procs away any proc needs, per dimension */
  int maxswap;            /**< Max number of swaps memory is allocated for */
  int *sendnum;           /**< Number of atoms to send in each swap */
  int *recvnum;           /**< Number of atoms to receive in each swap */
  int *sendproc;          /**< Processor to send to at each swap */
  int *recvproc;          /**< Processor to receive from at each swap */
  int *size_forward_recv; /**< Number of values to receive in each forward comm */
  int *size_reverse_send; /**< Number of values to send in each reverse comm */
  int *size_reverse_recv; /**< Number of values to receive in each reverse comm */
  double *slablo;         /**< Lower bound of slab to send at each swap */
  double *slabhi;         /**< Upper bound of slab to send at each swap */
  double **multilo;       /**< Lower bounds of slabs for multi-collection swap */
  double **multihi;       /**< Upper bounds of slabs for multi-collection swap */
  double **cutghostmulti; /**< Ghost cutoff on a per-collection basis */
  int *pbc_flag;          /**< General flag for sending atoms through PBC */
  int **pbc;              /**< Dimension flags for PBC adjustments */

  int *firstrecv;        /**< Where to put first received atom in each swap */
  int **sendlist;        /**< List of atoms to send in each swap */
  int *localsendlist;    /**< Indexed list of local sendlist atoms */
  int *maxsendlist;      /**< Max size of send list for each swap */

  double *buf_send;      /**< Send buffer for all communication */
  double *buf_recv;      /**< Receive buffer for all communication */
  int maxsend;           /**< Current size of send buffer */
  int maxrecv;           /**< Current size of receive buffer */
  int smax;              /**< Max size in atoms of single borders send */
  int rmax;              /**< Max size in atoms of single borders recv */

  // NOTE: init_buffers is called from a constructor and must not be made virtual
  void init_buffers();

  int updown(int, int, int, double, int, double *);
  // compare cutoff to procs
  virtual void grow_send(int, int);       /**< Reallocate send buffer */
  virtual void grow_recv(int);            /**< Free/allocate receive buffer */
  virtual void grow_list(int, int);       /**< Reallocate one sendlist */
  virtual void grow_swap(int);            /**< Grow swap and multi arrays */
  virtual void allocate_swap(int);        /**< Allocate swap arrays */
  virtual void allocate_multi(int);       /**< Allocate multi arrays */
  virtual void free_swap();               /**< Free swap arrays */
  virtual void free_multi();              /**< Free multi arrays */
};

}    // namespace LAMMPS_NS

#endif
