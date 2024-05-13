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

#ifndef LMP_MLIAP_UNIFIED_KOKKOS_H
#define LMP_MLIAP_UNIFIED_KOKKOS_H

#include "mliap_unified.h"
#include "mliap_descriptor_kokkos.h"
#include "mliap_model_kokkos.h"
#include "mliap_data_kokkos.h"

#include <Python.h>

namespace LAMMPS_NS {
template <class DeviceType>
class MLIAPDummyDescriptorKokkos : public MLIAPDummyDescriptor, public MLIAPDescriptorKokkos<DeviceType>{
 public:
  MLIAPDummyDescriptorKokkos(LAMMPS *);
  ~MLIAPDummyDescriptorKokkos() override;
  void compute_descriptors(class MLIAPData *) override;
  void compute_forces(class MLIAPData *) override;
  void compute_force_gradients(class MLIAPData *) override;
  void compute_descriptor_gradients(class MLIAPData *) override;
  void init() override;
  void set_elements(char **, int);
};
template <class DeviceType>
class MLIAPDummyModelKokkos : public MLIAPDummyModel, public MLIAPModelKokkos<DeviceType> {
 public:
  MLIAPDummyModelKokkos(LAMMPS *, char * = nullptr);
  ~MLIAPDummyModelKokkos() override;
  int get_nparams() override;
  int get_gamma_nnz(class MLIAPData *) override;
  void compute_gradients(class MLIAPData *) override;
  void compute_gradgrads(class MLIAPData *) override;
  void compute_force_gradients(class MLIAPData *) override;
  double memory_usage() override;

 protected:
  void read_coeffs(char *) override;
};

template <class DeviceType>
struct MLIAPBuildUnifiedKokkos_t {
  MLIAPDataKokkos<DeviceType> *data;
  MLIAPDummyDescriptorKokkos<DeviceType> *descriptor;
  MLIAPDummyModelKokkos<DeviceType> *model;
};
template <class DeviceType>
MLIAPBuildUnifiedKokkos_t<DeviceType> build_unified(char *, MLIAPDataKokkos<DeviceType> *, LAMMPS *, char * = NULL);
void update_pair_energy(MLIAPDataKokkosDevice *, double *);
void update_pair_forces(MLIAPDataKokkosDevice *, double *);
void update_atom_energy(MLIAPDataKokkosDevice *, double *);

}    // namespace LAMMPS_NS

#endif
