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

#ifndef LMP_COMM_TILED_H
#define LMP_COMM_TILED_H

#include "comm.h"

namespace LAMMPS_NS {

/** \class LAMMPS_NS::CommTiled
 * \brief Communication class for irregular tiled domain decomposition
 *
 * CommTiled implements inter-processor communication for LAMMPS simulations
 * using irregular (tiled) domain decompositions. This communication style
 * is used when the simulation domain has been decomposed using recursive
 * coordinate bisection (RCB) load balancing, which can result in non-rectangular
 * processor subdomains.
 *
 * Key features include:
 *
 * - Support for arbitrary processor subdomain shapes from RCB decomposition
 * - Communication with potentially many neighboring processors
 * - Efficient handling of overlapping ghost regions
 * - Multi-collection support for different cutoff distances
 *
 * Unlike CommBrick which uses a fixed 6-way stencil, CommTiled must dynamically
 * determine which processors own ghost atoms based on the irregular subdomain
 * boundaries. This allows for better load balancing but with some additional
 * communication overhead.
 *
 * \sa LAMMPS_NS::Comm, LAMMPS_NS::CommBrick */
class CommTiled : public Comm {
 public:
  /** Constructor for CommTiled class
   * \param lmp Pointer to the main LAMMPS instance */
  CommTiled(class LAMMPS *);

  /** Copy constructor from another Comm instance
   * \param lmp Pointer to the main LAMMPS instance
   * \param oldcomm Pointer to existing Comm instance to copy from */
  CommTiled(class LAMMPS *, class Comm *);

  /** Destructor - cleans up communication buffers and swap arrays */
  ~CommTiled() override;

  /** Initialize communication before a run */
  void init() override;

  /** Setup communication pattern for tiled decomposition */
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

  /** Setup for coordinate-to-processor mapping */
  void coord2proc_setup() override;

  /** Map a point to the owning processor for tiled decomposition
   * \param x Coordinates of the point
   * \param igx Output: grid x-index (not used for tiled)
   * \param igy Output: grid y-index (not used for tiled)
   * \param igz Output: grid z-index (not used for tiled)
   * \return Rank of owning processor */
  int coord2proc(double *, int &, int &, int &) override;

  /** Calculate memory usage
   * \return Memory usage in bytes */
  double memory_usage() override;

 protected:
  int nswap;      /**< Number of swaps to perform = 2*dim */
  int maxswap;    /**< Largest nswap can be = 6 */

  // forward/reverse comm info, proc lists include self

  int *nsendproc;          /**< Number of procs to send to per swap */
  int *nrecvproc;          /**< Number of procs to recv from per swap */
  int *sendother;          /**< 1 if send to other proc per swap */
  int *recvother;          /**< 1 if recv from other proc per swap */
  int *sendself;           /**< 1 if send to self per swap */
  int *nprocmax;           /**< Current max number of send procs per swap */
  int **sendproc;          /**< Procs to send to per swap */
  int **recvproc;          /**< Procs to recv from per swap */
  int **sendnum;           /**< Number of atoms to send per swap/proc */
  int **recvnum;           /**< Number of atoms to recv per swap/proc */
  int **size_forward_recv; /**< Number of values to recv in each forward swap/proc */
  int **firstrecv;         /**< Where to put first recv atom per swap/proc */
  int **size_reverse_send; /**< Number of values to send in each reverse swap/proc */
  int **size_reverse_recv; /**< Number of values to recv in each reverse swap/proc */
  int **forward_recv_offset;  /**< Forward comm offsets in buf_recv per swap/proc */
  int **reverse_recv_offset;  /**< Reverse comm offsets in buf_recv per swap/proc */
  int ***sendlist;         /**< List of atoms to send per swap/proc */
  int **maxsendlist;       /**< Max size of send list per swap/proc */
  int **pbc_flag;          /**< General flag for sending atoms through PBC */
  int ***pbc;              /**< Dimension flags for PBC adjustments */

  double ***sendbox;       /**< Bounding box of atoms to send per swap/proc */

  double **cutghostmulti;     /**< Ghost cutoff on a per-collection basis */
  double ****sendbox_multi;   /**< Bounding box of atoms to send per swap/proc for multi comm */

  // exchange comm info, proc lists do not include self

  int *nexchproc;       /**< Number of procs to exchange with in each dim */
  int *nexchprocmax;    /**< Current max number of exchange procs for each dim */
  int **exchproc;       /**< Procs to exchange with per dim */
  int **exchnum;        /**< Number of values received per dim/proc */

  double *buf_send;     /**< Send buffer for all communication */
  double *buf_recv;     /**< Receive buffer for all communication */
  int maxsend;          /**< Current size of send buffer */
  int maxrecv;          /**< Current size of receive buffer */
  int smaxone;          /**< Max size in atoms of single borders send */
  int rmaxone;          /**< Max size in atoms of single borders recv */
  int smaxall;          /**< Max size in atoms of any borders send to all procs in one swap */
  int rmaxall;          /**< Max size in atoms of any borders recv from all procs in one swap */

  int maxrequest;       /**< Max size of MPI Request vector */
  MPI_Request *requests;  /**< Array of MPI request handles */

  /** Structure to hold RCB decomposition info for one processor */
  struct RCBinfo {
    double mysplit[3][2];    /**< Fractional RCB bounding box for one proc */
    double cutfrac;          /**< Fractional position of cut this proc owns */
    int dim;                 /**< Dimension = 0/1/2 of cut */
  };

  RCBinfo *rcbinfo;    /**< List of RCB info for all procs */

  int noverlap;      /**< Number of overlapping procs */
  int maxoverlap;    /**< Current max length of overlap list */
  int *overlap;      /**< List of overlapping procs */

  double *prd;       /**< Local pointer to Domain periodic dimensions */
  double *boxlo;     /**< Local pointer to Domain box lower bounds */
  double *boxhi;     /**< Local pointer to Domain box upper bounds */
  double *sublo;     /**< Local pointer to subdomain lower bounds */
  double *subhi;     /**< Local pointer to subdomain upper bounds */
  int dimension;     /**< Simulation dimensionality (2 or 3) */

  void init_pointers();
  void init_buffers();
  int init_buffers_flag;  /**< Flag indicating if buffers have been initialized */

  // box drop and other functions

  using BoxDropPtr = void (CommTiled::*)(int, double *, double *, int &);
  BoxDropPtr box_drop;
  void box_drop_brick(int, double *, double *, int &);
  void box_drop_tiled(int, double *, double *, int &);
  void box_drop_tiled_recurse(double *, double *, int, int, int &);

  using BoxOtherPtr = void (CommTiled::*)(int, int, int, double *, double *);
  BoxOtherPtr box_other;
  void box_other_brick(int, int, int, double *, double *);
  void box_other_tiled(int, int, int, double *, double *);

  using BoxTouchPtr = int (CommTiled::*)(int, int, int);
  BoxTouchPtr box_touch;
  int box_touch_brick(int, int, int);
  int box_touch_tiled(int, int, int);

  using PointDropPtr = int (CommTiled::*)(int, double *);
  PointDropPtr point_drop;
  int point_drop_brick(int, double *);
  int point_drop_tiled(int, double *);
  int point_drop_tiled_recurse(double *, int, int);
  int closer_subbox_edge(int, double *);

  virtual void grow_send(int, int);              /**< Reallocate send buffer */
  virtual void grow_recv(int, int flag = 0);     /**< Free/allocate receive buffer */
  virtual void grow_list(int, int, int);         /**< Reallocate sendlist for one swap/proc */
  void allocate_swap(int);                       /**< Allocate swap arrays */
  virtual void grow_swap_send(int, int, int);    /**< Grow swap arrays for send and recv */
  void grow_swap_send_multi(int, int);           /**< Grow multi swap arrays for send and recv */
  void grow_swap_recv(int, int);                 /**< Grow swap recv arrays */
  void deallocate_swap(int);                     /**< Deallocate swap arrays */
};

}    // namespace LAMMPS_NS

#endif
