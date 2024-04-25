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
    int proc;                   // proc ID to perform direct swap with, can be self
    int allflag;                // 1 if sending all my owned atoms
    int xcheck,ycheck,zcheck;   // which coord dims to check to send subset
    int pbc_flag;               // 1 if any PBC coord shifts are needed when sender packs
    int pbc[6];                 // coord shift for each dim plus triclinic shifts
    int sendtag,recvtag;        // MPI send/recv tags for matching swap on sender/receiver
    double xlo,xhi;             // lo/hi bounds in each coord dim when sending subset
    double ylo,yhi;
    double zlo,zhi;
  };

  DirectSwap *dswap;
  int ndirect;                            // # of DirectSwaps with nearby procs, including self
  int maxdirect;                          // max # of DirectSwaps dswap is allocated for
  
  int nself_direct;                       // # of non-empty swaps with self
  
  int *send_indices_direct;               // indices of non-empty swap sends to other procs
  int *recv_indices_direct;               // indices of non-empty swap recvs with other procs
  int *self_indices_direct;               // indices of non-empty swaps with self

  int *proc_direct;                       // proc to send/recv to/from for each swap
  int *pbc_flag_direct;                   // general flag for sending atoms thru PBC
  int **pbc_direct;                       // dimension flags for PBC adjustments
  int *sendtag, *recvtag;                 // MPI tags for send and recv in each swap
  
  int *sendnum_direct;                    // # of atoms to send in each swap
  int *recvnum_direct;                    // # of atoms to recv in each swap
  int *size_forward_recv_direct;          // # of values to recv in each forward comm
  int *size_reverse_send_direct;          // # of values to send in each reverse comm
  int *size_reverse_recv_direct;          // # of values to recv in each reverse comm

  int *firstrecv_direct;    // index of where to put 1st ghost atom in each swap
  int **sendlist_direct;    // list of owned atoms to send in each swap
  int *maxsendlist_direct;  // max size of each sendlist_direct list

  double *buf_send_direct;  // send buffer used for every swap (large enough for any)
  double **buf_recv_direct; // list of recv buffers for all swaps (large enough for each)

  MPI_Request *requests;    // list of requests, length = ndirect

  void grow_list_direct(int, int);
};

}    // namespace LAMMPS_NS

#endif
