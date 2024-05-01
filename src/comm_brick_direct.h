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

#ifndef LMP_COMM_BRICK_DIRECT_H
#define LMP_COMM_BRICK_DIRECT_H

#include "comm_brick.h"

namespace LAMMPS_NS {

class CommBrickDirect : public CommBrick {
 public:
  CommBrickDirect(class LAMMPS *);
  CommBrickDirect(class LAMMPS *, class Comm *);
  ~CommBrickDirect() override;

  void init() override;                         // init error checks
  void setup() override;                        // setup direct comm data structs
  void forward_comm(int dummy = 0) override;    // forward comm of atom coords
  void reverse_comm() override;                 // reverse comm of forces
  void borders() override;                      // setup list of atoms to comm

  void forward_comm(class Pair *) override;                 // forward comm from a Pair
  void reverse_comm(class Pair *) override;                 // reverse comm from a Pair
  void forward_comm(class Bond *) override;                 // forward comm from a Bond
  void reverse_comm(class Bond *) override;                 // reverse comm from a Bond
  void forward_comm(class Fix *, int size = 0) override;    // forward comm from a Fix
  void reverse_comm(class Fix *, int size = 0) override;    // reverse comm from a Fix
  void reverse_comm_variable(class Fix *) override;         // variable size reverse comm from a Fix
  void forward_comm(class Compute *) override;              // forward from a Compute
  void reverse_comm(class Compute *) override;              // reverse from a Compute
  void forward_comm(class Dump *) override;                 // forward comm from a Dump
  void reverse_comm(class Dump *) override;                 // reverse comm from a Dump
  
  void forward_comm_array(int, double **) override;         // forward comm of array

 protected:
  // per-swap data
  // swap = exchange of data between me and another proc in stencil, including self

  int ndirect;                            // # of direct swaps with nearby procs, including self
  int maxdirect;                          // max size which all swap-length data is allocated for
  int nself_direct;                       // # of swaps with self, non-empty or empty

  int **swaporder;                        // ordering (ijk indices) of swaps within 3d stencil

  int *send_indices_direct;               // indices of non-empty swap sends to other procs
  int *recv_indices_direct;               // indices of non-empty swap recvs from other procs
  int *self_indices_direct;               // indices of non-empty swaps with self

  int *proc_direct;                       // proc to send/recv to/from for each swap, can be me
  int *pbc_flag_direct;                   // overall flag for sending atoms thru PBC
  int **pbc_direct;                       // 6 dimension flags for PBC adjusts, including triclinc
  int *sendtag, *recvtag;                 // MPI tags for send/recv in each swap

  int *sendnum_direct;                    // # of atoms to send in each swap
  int *recvnum_direct;                    // # of atoms to recv in each swap

  int *size_forward_recv_direct;          // max # of values to recv in each forward comm
  int *size_reverse_send_direct;          // max # of values to send in each reverse comm
  int *size_reverse_recv_direct;          // max # of values to recv in each reverse comm
  int *size_border_recv_direct;           // max # of values to recv in each border comm

  int *swap2list;                         // index to list of atoms each swap uses
  int **sendlist_direct;                  // ptrs to sendatoms_list for each swap
  int *firstrecv_direct;                  // index of first received ghost atom in each swap

  int *recv_offset_forward_direct;  // offsets into buf_recv_direct for forward comm receives
  int *recv_offset_reverse_direct;  // offsets into buf_recv_direct for reverse comm receives
  int *recv_offset_border_direct;   // offsets into buf_recv_direct for border comm receives
  int *recv_offset_forward_atoms;   // offsets in atom counts for forward comm receives
  int *recv_offset_reverse_atoms;   // offsets in atom counts for reverse comm receives

  // per-list data
  // list = indices of atom to send in a swap
  // only 27 (3d) or 9 (2d) possible lists
  // each may be used in multiple swaps or not used (or defined)

  int maxlist;                  // max possible lists
  int *active_list;             // 1 if each list is generated and used in a swap
  int **check_list;             // clist[I][J} = 1 if list I requires bounds check in dim J
  double ***bounds_list;        // blist[I][J][K] = lo/hi bounds K=0/1 in dim J for list I
  int *sendnum_list;            // # of atom indices in each list
  int **sendatoms_list;         // list of owned atom indices
  int *maxsendatoms_list;       // max size of each allocated list

  double cutxlo, cutxhi;    // cutoffs for sending owned atoms to procs on 6 faces of stencil
  double cutylo, cutyhi;
  double cutzlo, cutzhi;

  // communication buffers for MPI sends and receives as well as self data copies

  int smax_direct,rmax_direct;    // send/recv buf sizes in atom counts
  int ssum_direct,rsum_direct;    // max = max for one swap, sum = sum over all swaps
  
  double *buf_send_direct;  // send buffer used for every swap (large enough for any)
  double *buf_recv_direct;  // recv buffer used for all swaps (large enough for all)

  int maxsend_direct;       // size of buf_send_direct
  int maxrecv_direct;       // size of buf_recv_direct

  MPI_Request *requests;    // list of requests, length = ndirect

  // private methods
  // init_pointers and init_buffers_direct are called from a constructor
  //   so must not be made virtual

  void init_pointers();
  void init_buffers_direct();

  void order_swaps(int, int, int, int, int, int);
  void allocate_direct();
  void allocate_lists();
  void deallocate_direct();
  void deallocate_lists(int);

  void check_buffer_sizes();
  void grow_send_direct(int, int);
  void grow_recv_direct(int);
  void grow_list_direct(int, int);
};

}    // namespace LAMMPS_NS

#endif
