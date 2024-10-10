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

class CommTiled : public Comm {
 public:
  CommTiled(class LAMMPS *);
  CommTiled(class LAMMPS *, class Comm *);

  ~CommTiled() override;

  void init() override;
  void setup() override;                        // setup comm pattern
  void forward_comm(int dummy = 0) override;    // forward comm of atom coords
  void reverse_comm() override;                 // reverse comm of forces
  void exchange() override;                     // move atoms to new procs
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

  void forward_comm_array(int, double **) override;            // forward comm of array

  void coord2proc_setup() override;
  int coord2proc(double *, int &, int &, int &) override;

  double memory_usage() override;

 protected:
  int nswap;      // # of swaps to perform = 2*dim
  int maxswap;    // largest nswap can be = 6

  // forward/reverse comm info, proc lists include self

  int *nsendproc, *nrecvproc;    // # of procs to send/recv to/from per swap
  int *sendother, *recvother;    // 1 if send/recv to/from other proc per swap
  int *sendself;                 // 1 if send to self per swap
  int *nprocmax;                 // current max # of send procs per swap
  int **sendproc, **recvproc;    // procs to send/recv to/from per swap
  int **sendnum, **recvnum;      // # of atoms to send/recv per swap/proc
  int **size_forward_recv;       // # of values to recv in each forward swap/proc
  int **firstrecv;               // where to put 1st recv atom per swap/proc
  int **size_reverse_send;       // # of values to send in each reverse swap/proc
  int **size_reverse_recv;       // # of values to recv in each reverse swap/proc
  int **forward_recv_offset;     // forward comm offsets in buf_recv per swap/proc
  int **reverse_recv_offset;     // reverse comm offsets in buf_recv per swap/proc
  int ***sendlist;               // list of atoms to send per swap/proc
  int **maxsendlist;             // max size of send list per swap/proc
  int **pbc_flag;                // general flag for sending atoms thru PBC
  int ***pbc;                    // dimension flags for PBC adjustments

  double ***sendbox;    // bounding box of atoms to send per swap/proc

  double **cutghostmulti;         // cutghost on a per-collection basis
  double **cutghostmultiold;      // cutghost on a per-type basis
  double ****sendbox_multi;       // bounding box of atoms to send
                                  //   per swap/proc for multi comm
  double ****sendbox_multiold;    // bounding box of atoms to send
                                  //   per swap/proc for multi/old comm

  // exchange comm info, proc lists do not include self

  int *nexchproc;       // # of procs to send/recv to/from in each dim
  int *nexchprocmax;    // current max # of exch procs for each dim
  int **exchproc;       // procs to exchange with per dim
  int **exchnum;        // # of values received per dim/proc

  double *buf_send;        // send buffer for all comm
  double *buf_recv;        // recv buffer for all comm
  int maxsend, maxrecv;    // current size of send/recv buffer
  int smaxone, rmaxone;    // max size in atoms of single borders send/recv
  int smaxall, rmaxall;    // max size in atoms of any borders send/recv
                           //   for comm to all procs in one swap

  int maxrequest;    // max size of Request vector
  MPI_Request *requests;

  struct RCBinfo {
    double mysplit[3][2];    // fractional RCB bounding box for one proc
    double cutfrac;          // fractional position of cut this proc owns
    int dim;                 // dimension = 0/1/2 of cut
  };

  RCBinfo *rcbinfo;    // list of RCB info for all procs

  int noverlap;      // # of overlapping procs
  int maxoverlap;    // current max length of overlap
  int *overlap;      // list of overlapping procs

  double *prd;    // local ptrs to Domain attributes
  double *boxlo, *boxhi;
  double *sublo, *subhi;
  int dimension;

  // NOTE: init_pointers and init_buffers are called from a constructor
  //  and must not be made virtual

  void init_pointers();
  void init_buffers();
  int init_buffers_flag;

  // box drop and other functions

  typedef void (CommTiled::*BoxDropPtr)(int, double *, double *, int &);
  BoxDropPtr box_drop;
  void box_drop_brick(int, double *, double *, int &);
  void box_drop_tiled(int, double *, double *, int &);
  void box_drop_tiled_recurse(double *, double *, int, int, int &);

  typedef void (CommTiled::*BoxOtherPtr)(int, int, int, double *, double *);
  BoxOtherPtr box_other;
  void box_other_brick(int, int, int, double *, double *);
  void box_other_tiled(int, int, int, double *, double *);

  typedef int (CommTiled::*BoxTouchPtr)(int, int, int);
  BoxTouchPtr box_touch;
  int box_touch_brick(int, int, int);
  int box_touch_tiled(int, int, int);

  typedef int (CommTiled::*PointDropPtr)(int, double *);
  PointDropPtr point_drop;
  int point_drop_brick(int, double *);
  int point_drop_tiled(int, double *);
  int point_drop_tiled_recurse(double *, int, int);
  int closer_subbox_edge(int, double *);

  virtual void grow_send(int, int);               // reallocate send buffer
  virtual void grow_recv(int, int flag = 0);      // free/allocate recv buffer
  virtual void grow_list(int, int, int);          // reallocate sendlist for one swap/proc
  void allocate_swap(int);                // allocate swap arrays
  virtual void grow_swap_send(int, int, int);     // grow swap arrays for send and recv
  void grow_swap_send_multi(int, int);    // grow multi swap arrays for send and recv
  void grow_swap_recv(int, int);
  void deallocate_swap(int);    // deallocate swap arrays
};

}    // namespace LAMMPS_NS

#endif
