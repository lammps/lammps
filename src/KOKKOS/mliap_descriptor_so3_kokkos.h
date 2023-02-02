/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Matt Bettencourt (NVIDIA)
 ------------------------------------------------------------------------- */

#ifndef LMP_MLIAP_DESCRIPTOR_SO3_KOKKOS_H
#define LMP_MLIAP_DESCRIPTOR_SO3_KOKKOS_H

#include "mliap_descriptor_kokkos.h"
#include "mliap_descriptor_so3.h"
#include "mliap_so3_kokkos.h"

namespace LAMMPS_NS {
template <class DeviceType>
class MLIAPDescriptorSO3Kokkos :
    public MLIAPDescriptorSO3,
    public MLIAPDescriptorKokkos<DeviceType> {
 public:
  MLIAPDescriptorSO3Kokkos(LAMMPS *, char *);

  void compute_descriptors(class MLIAPData *) override;
  void compute_forces(class MLIAPData *) override;
  void compute_force_gradients(class MLIAPData *) override;
  void compute_descriptor_gradients(class MLIAPData *) override;
  void init() override;
  double memory_usage() override;

 protected:
  template <typename ViewType>
  KOKKOS_FUNCTION static void v_tally(int vflag_either, int vflag_global, int vflag_atom, int i,
                                      int j, int ij, double *fij, ViewType rij,
                                      Kokkos::View<double[6], DeviceType> virial, ViewType vatom);
  class MLIAP_SO3Kokkos<DeviceType> *so3ptr_kokkos;

  // inherited from non-Kokkos class
};
}    // namespace LAMMPS_NS
#endif
