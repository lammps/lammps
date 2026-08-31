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

#ifndef LMP_COMM_BRICK_DIRECT_KOKKOS_H
#define LMP_COMM_BRICK_DIRECT_KOKKOS_H

#include "comm_brick_direct.h"
#include "comm_brick_kokkos.h"

namespace LAMMPS_NS {

class CommBrickDirectKokkos : public CommBrickDirect {
 public:
  CommBrickDirectKokkos(class LAMMPS *);
  CommBrickDirectKokkos(class LAMMPS *, class Comm *);
  ~CommBrickDirectKokkos() override;
  void setup() override;                        // setup direct comm data structs

  using CommBrick::forward_comm;
  using CommBrick::reverse_comm;
  void forward_comm(int dummy = 0) override;    // forward comm of atom coords
  void reverse_comm() override;                 // reverse comm of forces
  void exchange() override;                     // move atoms to new procs
  void borders() override;                      // setup list of atoms to comm

  template<class DeviceType> void forward_comm_device();
  template<class DeviceType> void reverse_comm_device();
  template<class DeviceType> void build_lists_device();
  template<class DeviceType> void borders_comm_device();

 private:
  DAT::tdual_double_2d_lr k_buf_send_direct,k_buf_recv_direct;
  DAT::tdual_double_2d_lr k_buf_send_reverse,k_buf_recv_reverse;
  DAT::tdual_double_2d_lr k_buf_send_border,k_buf_recv_border;
  DAT::tdual_int_scalar k_total_send;
  int border_device_flag;
  DAT::tdual_int_2d_lr k_sendatoms_list;
  DAT::tdual_int_1d k_swap2list;
  DAT::tdual_int_2d k_pbc_direct;
  DAT::tdual_int_1d k_pbc_flag_direct;
  DAT::tdual_int_1d k_firstrecv_direct;
  DAT::tdual_int_1d k_sendnum_scan_direct;
  DAT::tdual_int_1d k_self_flag;
  int totalsend;

  void grow_send_direct(int, int) override;
  void grow_recv_direct(int) override;

  // the send lists are the device's copy; sendatoms_list aliases its host side

  void build_lists() override;
  void borders_comm() override;

  void allocate_lists() override;
  void deallocate_lists(int) override;
  void grow_list_direct(int, int) override;
  void lists_to_kokkos();
};

}    // namespace LAMMPS_NS

#endif
