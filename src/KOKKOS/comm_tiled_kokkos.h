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

#ifndef LMP_COMM_TILED_KOKKOS_H
#define LMP_COMM_TILED_KOKKOS_H

#include "comm_tiled.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class CommTiledKokkos : public CommTiled {
 public:
  CommTiledKokkos(class LAMMPS *);
  CommTiledKokkos(class LAMMPS *, class Comm *);

  ~CommTiledKokkos() override;

  bool exchange_comm_classic;
  bool forward_comm_classic;
  bool forward_pair_comm_classic;
  bool reverse_pair_comm_classic;
  bool forward_fix_comm_classic;
  bool reverse_comm_classic;
  bool exchange_comm_on_host;
  bool forward_comm_on_host;
  bool reverse_comm_on_host;

  using CommTiled::forward_comm;
  using CommTiled::reverse_comm;

  void init() override;
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

  void forward_comm_array(int, double **) override;          // forward comm of array

  template<class DeviceType> void forward_comm_device();
  template<class DeviceType> void reverse_comm_device();

 protected:

  DAT::tdual_int_3d k_sendlist;
  //DAT::tdual_int_scalar k_total_send;
  DAT::tdual_xfloat_2d k_buf_send,k_buf_recv;
  //DAT::tdual_int_scalar k_count;

  void grow_send(int, int) override;
  void grow_recv(int, int flag = 0) override;
  void grow_send_kokkos(int, int, ExecutionSpace space = Host);
  void grow_recv_kokkos(int, int, ExecutionSpace space = Host);
  void grow_list(int, int, int) override;
  void grow_swap_send(int, int, int) override;     // grow swap arrays for send and recv
};

}    // namespace LAMMPS_NS

#endif
