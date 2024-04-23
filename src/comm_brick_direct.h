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
  int ndirect;                            // # of swaps to perform with direct neighbors
  int nrecv_direct;
  int nself_direct;
  int maxdirect;                          // max # of swaps memory is allocated for
  int *sendnum_direct, *recvnum_direct;   // # of atoms to send/recv in each swap
  int *sendproc_direct, *recvproc_direct; // proc to send/recv to/from at each swap
  int *size_forward_recv_direct;          // # of values to recv in each forward comm
  int *size_reverse_send_direct;          // # to send in each reverse comm
  int *size_reverse_recv_direct;          // # to recv in each reverse comm
  int *pbc_flag_direct;                   // general flag for sending atoms thru PBC
  int **pbc_direct;                       // dimension flags for PBC adjustments

  int *self_direct;
  int *firstrecv_direct;    // where to put 1st recv atom in each swap
  int **sendlist_direct;    // list of atoms to send in each swap

  double *buf_send_direct;  // send buffer for all comm
  double **buf_recv_direct; // recv buffer for all comm

  MPI_Request *requests;

};

}    // namespace LAMMPS_NS

#endif
