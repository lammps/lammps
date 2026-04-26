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

#ifdef FIX_CLASS
// clang-format off
FixStyle(rigid/small/kk,FixRigidSmallKokkos<LMPDeviceType>);
FixStyle(rigid/small/kk/device,FixRigidSmallKokkos<LMPDeviceType>);
FixStyle(rigid/small/kk/host,FixRigidSmallKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_RIGID_SMALL_KOKKOS_H
#define LMP_FIX_RIGID_SMALL_KOKKOS_H

#include "fix_rigid_small.h"
#include "fix_rigid_base_kokkos.h"


namespace LAMMPS_NS {

template<class DeviceType>
class FixRigidSmallKokkos : public FixRigidSmall, public FixRigidBaseKokkos<DeviceType, FixRigidSmall> {
 public:
  FixRigidSmallKokkos(class LAMMPS *, int, char **);
  ~FixRigidSmallKokkos() override;

  void post_constructor() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void post_force(int) override;
  void final_integrate() override;
  void setup_pre_neighbor() override;
  void pre_neighbor() override;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  void set_molecule(int, tagint, int, double *, double *, double *) override;

  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  bigint dof(int) override;
  void zero_momentum() override;
  void zero_rotation() override;
  double compute_scalar() override;

 protected:
  class AtomKokkos *atomKK;
  ExecutionSpace execution_space;

  void grow_body() override;

};

}    // namespace LAMMPS_NS

#endif // !LMP_FIX_RIGID_SMALL_KOKKOS_H
#endif
