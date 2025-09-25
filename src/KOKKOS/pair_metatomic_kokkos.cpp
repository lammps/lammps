/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Guillaume Fraux <guillaume.fraux@epfl.ch>
                         Filippo Bigi <filippo.bigi@epfl.ch>
------------------------------------------------------------------------- */
#include "pair_metatomic_kokkos.h"

#include "error.h"
#include "neigh_request.h"
#include "atom_masks.h"

#include "atom_kokkos.h"

#include "metatomic_system_kokkos.h"
#include "metatomic_types.h"

#include <algorithm>
#include <cctype>

using namespace LAMMPS_NS;

// LAMMPS uses `LAMMPS_NS::tagint` and `int` for tags and neighbor lists, respectively.
// For the moment, we require both to be int32_t for this interface
static_assert(std::is_same_v<LAMMPS_NS::tagint, int32_t>, "Error: LAMMPS_NS::tagint must be int32_t to compile metatensor/kk");
static_assert(std::is_same_v<int, int32_t>, "Error: int must be int32_t to compile metatensor/kk");

template<typename T, class DeviceType>
using UnmanagedView = Kokkos::View<T, Kokkos::LayoutRight, DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template<class DeviceType>
PairMetatomicKokkos<DeviceType>::PairMetatomicKokkos(LAMMPS* lmp): PairMetatomic(lmp) {
    respa_enable = 0;

    kokkosable = 1;
    execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

    datamask_read = X_MASK | F_MASK | TYPE_MASK | TAG_MASK | ENERGY_MASK | VIRIAL_MASK;
    datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

template<class DeviceType>
PairMetatomicKokkos<DeviceType>::~PairMetatomicKokkos() {}

template<class DeviceType>
void PairMetatomicKokkos<DeviceType>::init_style() {
    PairMetatomic::init_style();

    auto request = neighbor->find_request(this);
    request->set_kokkos_host(
        std::is_same_v<DeviceType, LMPHostType> &&
        !std::is_same_v<DeviceType, LMPDeviceType>
    );
    request->set_kokkos_device(std::is_same_v<DeviceType, LMPDeviceType>);

    // copy type mapping from host to device, to be able to give a device pointer
    // to MetatomicSystemAdaptorKokkos
    auto type_mapping_kk_host = UnmanagedView<int32_t*, LMPHostType>(this->type_mapping, atom->ntypes + 1);
    this->type_mapping_kk = Kokkos::View<int32_t*, Kokkos::LayoutRight, DeviceType>("type_mapping_kk", atom->ntypes + 1);
    Kokkos::deep_copy(this->type_mapping_kk, type_mapping_kk_host);

    auto options = MetatomicSystemOptions{
        this->type_mapping_kk.data(),
        mta_data->max_cutoff,
        mta_data->check_consistency,
        !(mta_data->non_conservative),
    };

    // override the system adaptor with the kokkos version
    this->system_adaptor = std::make_unique<MetatomicSystemAdaptorKokkos<DeviceType>>(lmp, options);

    // request NL with the new adaptor
    auto requested_nl = mta_data->model->run_method("requested_neighbor_lists");
    for (const auto& ivalue: requested_nl.toList()) {
        auto options = ivalue.get().toCustomClass<metatomic_torch::NeighborListOptionsHolder>();
        auto cutoff = options->engine_cutoff(mta_data->evaluation_options->length_unit());
        assert(cutoff <= mta_data->max_cutoff);

        this->system_adaptor->add_nl_request(cutoff, options);
    }
}

template<class DeviceType>
void PairMetatomicKokkos<DeviceType>::pick_device(torch::Device& device, const char* requested) {
    device = KokkosDeviceToTorch<DeviceType>::convert();

    if (requested != nullptr) {
        auto requested_str = std::string(requested);
        std::transform(requested_str.begin(), requested_str.end(), requested_str.begin(), ::tolower);
        if (c10::DeviceTypeName(device.type(), /*lower_case=*/true) != requested_str) {
            error->all(FLERR,
                "requested device '{}' does not match the device being used by kokkos '{}', "
                "use the non-kokkos version of this pair style to use a different "
                "device for the model and LAMMPS",
                requested, device.str()
            );
        }
    }
}

template<class DeviceType>
void PairMetatomicKokkos<DeviceType>::pre_compute() {
    /// Declare what we need to read from the atomKK object and what we will modify
    this->atomKK->sync(ExecutionSpaceFromDevice<DeviceType>::space, datamask_read);
    this->atomKK->modified(ExecutionSpaceFromDevice<DeviceType>::space, datamask_modify);
}

template<class DeviceType>
void PairMetatomicKokkos<DeviceType>::store_forces(const at::Tensor& forces_tensor) {
    assert(forces_tensor.scalar_type() == torch::kFloat64);
    auto forces = forces_tensor.contiguous();

    auto forces_lammps_kk = this->atomKK->k_f.template view<DeviceType>();
    auto forces_metatensor_kk = UnmanagedView<double**, DeviceType>(
        forces.template data_ptr<double>(),
        forces.size(0), 3
    );

    int num_forces_to_update;
    if (mta_data->non_conservative) {
        num_forces_to_update = atomKK->nlocal;
    } else {
        num_forces_to_update = atomKK->nlocal + atomKK->nghost;
    }

    double scale = this->scale;  // the GPU can't access the `this` pointer
    Kokkos::parallel_for(
        num_forces_to_update,
        KOKKOS_LAMBDA(size_t i) {
            forces_lammps_kk(i, 0) += scale * forces_metatensor_kk(i, 0);
            forces_lammps_kk(i, 1) += scale * forces_metatensor_kk(i, 1);
            forces_lammps_kk(i, 2) += scale * forces_metatensor_kk(i, 2);
        }
    );
}


namespace LAMMPS_NS {
template class PairMetatomicKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairMetatomicKokkos<LMPHostType>;
#endif
}
