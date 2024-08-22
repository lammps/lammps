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
  void forward_comm(int dummy = 0) override;    // forward comm of atom coords
  void borders() override;                      // setup list of atoms to comm

  template<class DeviceType> void forward_comm_device();

 private:
  DAT::tdual_xfloat_2d k_buf_send_direct,k_buf_recv_direct;
  DAT::tdual_int_2d k_sendatoms_list;
  DAT::tdual_int_1d k_swap2list;
  DAT::tdual_int_2d k_pbc_direct;
  DAT::tdual_int_1d k_pbc_flag_direct;
  DAT::tdual_int_1d k_firstrecv_direct;
  DAT::tdual_int_1d k_sendnum_scan_direct;
  DAT::tdual_int_1d k_self_flags;
  int totalsend;

  void grow_send_direct(int, int) override;
  void grow_recv_direct(int) override;
};

}    // namespace LAMMPS_NS

#endif
