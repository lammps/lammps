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

  void setup() override;                        // setup 3d comm pattern
  void forward_comm(int dummy = 0) override;    // forward comm of atom coords
  void reverse_comm() override;                 // reverse comm of forces
  void borders() override;                      // setup list of atoms to comm

 protected:
  
  struct DirectSwap {
    int proc;
    int allflag;
    int xcheck,ycheck,zcheck;
    double xlo,xhi;
    double ylo,yhi;
    double zlo,zhi;
    int pbc_flag;
    int pbc[6];
  };

  DirectSwap *dswap;
  int ndirect;                            // # of DirectSwaps with nearby procs, including self
  int maxdirect;                          // max # of DirectSwaps dswap is allocated for
  
  int nself_direct;                       // # of non-empty swaps with self
  
  int *send_indices_direct;               // indices of non-empty swap sends to other procs
  int *recv_indices_direct;               // indices of non-empty swap recvs with other procs
  int *self_indices_direct;               // indices of non-empty swaps with self

  int *sendnum_direct, *recvnum_direct;   // # of atoms to send/recv in each swap
  int *proc_direct;                       // proc to send/recv to/from at each swap
  int *size_forward_recv_direct;          // # of values to recv in each forward comm
  int *size_reverse_send_direct;          // # of values to send in each reverse comm
  int *size_reverse_recv_direct;          // # of values to recv in each reverse comm
  int *pbc_flag_direct;                   // general flag for sending atoms thru PBC
  int **pbc_direct;                       // dimension flags for PBC adjustments

  int *firstrecv_direct;    // index of where to put 1st ghost atom in each swap
  int **sendlist_direct;    // list of owned atoms to send in each swap
  int *maxsendlist_direct;  // max size of each sendlist_direct list

  double *buf_send_direct;  // send buffer used for every swap (large enough for any)
  double **buf_recv_direct; // lsit of recv buffers for all swaps (large enough for each)

  MPI_Request *requests;    // list of requests, length = ndirect

  void grow_list_direct(int, int);
};

}    // namespace LAMMPS_NS

#endif
