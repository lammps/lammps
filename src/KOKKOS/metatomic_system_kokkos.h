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

#ifndef LMP_METATENSOR_SYSTEM_KOKKOS_H
#define LMP_METATENSOR_SYSTEM_KOKKOS_H

#include "metatomic_system.h"

#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

/* ---------------------------------------------------------------------- */

// See https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html for a
// list of execution spaces.
template<class DeviceType>
struct KokkosDeviceToTorch {};

#if defined(KOKKOS_ENABLE_SERIAL)
template<> struct KokkosDeviceToTorch<Kokkos::Serial> {
    static torch::Device convert() {
        return torch::Device(torch::kCPU);
    }
};
#endif

#if defined(KOKKOS_ENABLE_CUDA)
template<> struct KokkosDeviceToTorch<Kokkos::Cuda> {
    static torch::Device convert() {
        return torch::Device(torch::kCUDA, Kokkos::device_id());
    }
};
#endif

#if defined(KOKKOS_ENABLE_HIP)
template<> struct KokkosDeviceToTorch<Kokkos::HIP> {
    static torch::Device convert() {
        return torch::Device(torch::kHIP, Kokkos::device_id());
    }
};
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
template<> struct KokkosDeviceToTorch<Kokkos::OpenMP> {
    static torch::Device convert() {
        return torch::Device(torch::kCPU);
    }
};
#endif

#if defined(KOKKOS_ENABLE_THREADS)
template<> struct KokkosDeviceToTorch<Kokkos::Threads> {
    static torch::Device convert() {
        return torch::Device(torch::kCPU);
    }
};
#endif

// Kokkos::SYCL, Kokkos::OpenMPTarget, Kokkos::HPX don't have a matching device
// in torch.

/* ---------------------------------------------------------------------- */

// data for metatomic neighbors lists
struct MetatomicNeighborsDataKokkos {
    // single neighbors sample containing [i, j, S_a, S_b, S_c]
    using sample_t = std::array<int32_t, 5>;

    // cutoff for this NL in LAMMPS units
    double cutoff;
    // options of the NL as requested by the model
    metatomic_torch::NeighborListOptions options;

    // Below are cached allocations for the LAMMPS -> metatomic NL translation
    // TODO: report memory usage for these?

    Kokkos::View<int32_t**, Kokkos::LayoutRight, LMPDeviceType> samples;
    // pairs distances vectors
    Kokkos::View<double**, Kokkos::LayoutRight, LMPDeviceType> distances_f64;
    Kokkos::View<float**, Kokkos::LayoutRight, LMPDeviceType> distances_f32;
};

template<class DeviceType>
class MetatomicSystemAdaptorKokkos : public MetatomicSystemAdaptor {
public:
    MetatomicSystemAdaptorKokkos(LAMMPS *lmp, MetatomicSystemOptions options);
    ~MetatomicSystemAdaptorKokkos() override {}

    void add_nl_request(
        double cutoff, metatomic_torch::NeighborListOptions request
    ) override;

    // Create a metatensor system matching the LAMMPS-Kokkos system data
    metatomic_torch::System system_from_lmp(
        NeighList* list,
        bool do_virial,
        torch::ScalarType dtype,
        torch::Device device
    ) override;

    void setup_neighbors_kk(metatomic_torch::System& system, NeighListKokkos<DeviceType>* list);

private:
    /// Torch device corresponding to the kokkos `DeviceType`
    torch::Device device_;

    // allocations caches for all the NL requested by the model
    std::vector<MetatomicNeighborsDataKokkos> nl_requests_kk_;
};

}    // namespace LAMMPS_NS

#endif
