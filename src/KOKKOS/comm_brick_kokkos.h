// clang-format off
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

#ifndef LMP_COMM_BRICK_KOKKOS_H
#define LMP_COMM_BRICK_KOKKOS_H

#include "comm_brick.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class CommBrickKokkos : public CommBrick {
 public:


  bool exchange_comm_classic;
  bool forward_comm_classic;
  bool forward_pair_comm_classic;
  bool reverse_pair_comm_classic;
  bool forward_fix_comm_classic;
  bool reverse_comm_classic;
  bool exchange_comm_on_host;
  bool forward_comm_on_host;
  bool reverse_comm_on_host;

  CommBrickKokkos(class LAMMPS *);
  ~CommBrickKokkos() override;
  void init() override;

  using CommBrick::forward_comm;
  using CommBrick::reverse_comm;
  void forward_comm(int dummy = 0) override;    // forward comm of atom coords
  void reverse_comm() override;                 // reverse comm of atom coords
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

  template<class DeviceType> void forward_comm_device();
  template<class DeviceType> void reverse_comm_device();
  template<class DeviceType> void forward_comm_device(Pair *pair);
  template<class DeviceType> void reverse_comm_device(Pair *pair);
  template<class DeviceType> void forward_comm_device(Fix *fix, int size=0);
  template<class DeviceType> void exchange_device();
  template<class DeviceType> void borders_device();

 protected:
  DAT::tdual_int_2d k_sendlist;
  DAT::tdual_int_scalar k_total_send;
  DAT::tdual_xfloat_2d k_buf_send,k_buf_recv;
  DAT::tdual_int_1d k_exchange_sendlist,k_exchange_copylist,k_indices;
  DAT::tdual_int_scalar k_count;

  DAT::tdual_int_2d k_swap;
  DAT::tdual_int_2d k_swap2;
  DAT::tdual_int_2d k_pbc;
  DAT::tdual_int_1d k_pbc_flag;
  DAT::tdual_int_1d k_g2l;
  DAT::tdual_int_1d k_firstrecv;
  DAT::tdual_int_1d k_sendnum_scan;
  int totalsend;

  int max_buf_pair,max_buf_fix;
  DAT::tdual_xfloat_1d k_buf_send_pair, k_buf_send_fix;
  DAT::tdual_xfloat_1d k_buf_recv_pair, k_buf_recv_fix;
  void grow_buf_pair(int);
  void grow_buf_fix(int);

  void grow_send(int, int) override;
  void grow_recv(int) override;
  void grow_send_kokkos(int, int, ExecutionSpace space = Host);
  void grow_recv_kokkos(int, ExecutionSpace space = Host);
  void grow_list(int, int) override;
  void grow_swap(int) override;
  void copy_swap_info();
};

}    // namespace LAMMPS_NS

#endif
