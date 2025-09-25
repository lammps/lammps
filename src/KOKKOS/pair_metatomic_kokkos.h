/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#ifdef PAIR_CLASS
// clang-format off
PairStyle(metatomic/kk, PairMetatomicKokkos<LMPDeviceType>);

// additional versions of pair_style metatomic/kk, to be used with pair_style
// hybrid when combining multiple metatomic potentials together
PairStyle(metatomic_1/kk, PairMetatomicKokkos<LMPDeviceType>);
PairStyle(metatomic_2/kk, PairMetatomicKokkos<LMPDeviceType>);
PairStyle(metatomic_3/kk, PairMetatomicKokkos<LMPDeviceType>);
PairStyle(metatomic_4/kk, PairMetatomicKokkos<LMPDeviceType>);
PairStyle(metatomic_5/kk, PairMetatomicKokkos<LMPDeviceType>);
PairStyle(metatomic_6/kk, PairMetatomicKokkos<LMPDeviceType>);
PairStyle(metatomic_7/kk, PairMetatomicKokkos<LMPDeviceType>);
PairStyle(metatomic_8/kk, PairMetatomicKokkos<LMPDeviceType>);
PairStyle(metatomic_9/kk, PairMetatomicKokkos<LMPDeviceType>);
// clang-format on
#else

#ifndef LMP_PAIR_METATOMIC_KOKKOS_H
#define LMP_PAIR_METATOMIC_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_metatomic.h"

namespace LAMMPS_NS {

template<class DeviceType>
class MetatomicSystemAdaptorKokkos;

template<class DeviceType>
struct PairMetatomicDataKokkos;

template<class DeviceType>
class PairMetatomicKokkos : public PairMetatomic {
public:
    PairMetatomicKokkos(class LAMMPS *);
    ~PairMetatomicKokkos();

    void init_style() override;

    void pre_compute() override;
    void store_forces(const at::Tensor& forces_tensor) override;

private:
    void pick_device(c10::Device& device, const char* requested) override;

    Kokkos::View<int32_t*, Kokkos::LayoutRight, DeviceType> type_mapping_kk;
};

}    // namespace LAMMPS_NS

#endif
#endif
