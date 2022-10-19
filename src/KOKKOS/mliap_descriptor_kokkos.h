/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing author: Matt Bettencourt (NVIDIA)
------------------------------------------------------------------------- */

#ifndef SRC_KOKKOS_MLIAP_DESCRIPTOR_KOKKOS_H_
#define SRC_KOKKOS_MLIAP_DESCRIPTOR_KOKKOS_H_



#include "pointers.h"
#include "kokkos_type.h"
#include "memory_kokkos.h"
#include "mliap_descriptor.h"

namespace LAMMPS_NS {

template<class DeviceType>
class MLIAPDescriptorKokkos : virtual protected Pointers {
public:
  MLIAPDescriptorKokkos(LAMMPS *lmp, MLIAPDescriptor *descriptor_in) : Pointers(lmp), descriptor(descriptor_in) {
    memoryKK->destroy_kokkos(k_wjelem);
  }

  void init_data() {
    int num_elems = descriptor->nelements;
    memoryKK->destroy_kokkos(k_wjelem);
    memoryKK->create_kokkos(k_wjelem,num_elems,"MLIAPDescriptorKokkos::k_wjelem");
    for (int i=0;i<num_elems; ++i)
      k_wjelem.h_view(i) = descriptor->wjelem[i];
    k_wjelem.modify<LMPHostType>();
    k_wjelem.sync<LMPDeviceType>();
  }
  virtual ~MLIAPDescriptorKokkos() {
//    memoryKK->destroy_kokkos(k_coeffelem);
//    model->coeffelem=nullptr;
    memoryKK->destroy_kokkos(k_wjelem);
  }


  MLIAPDescriptor *descriptor;
  DAT::tdual_float_1d k_wjelem;

};

}// namespace




#endif /* SRC_KOKKOS_MLIAP_DESCRIPTOR_KOKKOS_H_ */
