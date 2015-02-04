/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_COMM_KOKKOS_H
#define LMP_COMM_KOKKOS_H

#include "comm_brick.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class CommKokkos : public CommBrick {
 public:
  class AtomKokkos *atomKK;

  bool exchange_comm_classic;
  bool forward_comm_classic;
  bool exchange_comm_on_host;
  bool forward_comm_on_host;

  CommKokkos(class LAMMPS *);
  ~CommKokkos();
  void init();

  void forward_comm(int dummy = 0);    // forward comm of atom coords
  void reverse_comm();              // reverse comm of atom coords
  void exchange();                     // move atoms to new procs
  void borders();                      // setup list of atoms to comm

  void forward_comm_pair(class Pair *);    // forward comm from a Pair
  void reverse_comm_pair(class Pair *);    // reverse comm from a Pair
  void forward_comm_fix(class Fix *, int size=0);      // forward comm from a Fix
  void reverse_comm_fix(class Fix *, int size=0);      // reverse comm from a Fix
  void forward_comm_compute(class Compute *);  // forward from a Compute
  void reverse_comm_compute(class Compute *);  // reverse from a Compute
  void forward_comm_dump(class Dump *);    // forward comm from a Dump
  void reverse_comm_dump(class Dump *);    // reverse comm from a Dump

  template<class DeviceType> void forward_comm_device(int dummy);
  template<class DeviceType> void forward_comm_pair_device(Pair *pair);
  template<class DeviceType> void exchange_device();
  template<class DeviceType> void borders_device();

 protected:
  DAT::tdual_int_2d k_sendlist;
  DAT::tdual_xfloat_2d k_buf_send,k_buf_recv;
  DAT::tdual_int_1d k_exchange_sendlist,k_exchange_copylist,k_sendflag;
  DAT::tdual_int_1d k_count;
  //double *buf_send;                 // send buffer for all comm
  //double *buf_recv;                 // recv buffer for all comm

  void grow_send(int, int);
  void grow_recv(int);
  void grow_send_kokkos(int, int, ExecutionSpace space = Host);
  void grow_recv_kokkos(int, ExecutionSpace space = Host);
  void grow_list(int, int);
  void grow_swap(int);
};

}

#endif

/* ERROR/WARNING messages:

*/
