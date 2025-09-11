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
#include "metatomic_timer.h"

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
void PairMetatomicKokkos<DeviceType>::pick_device(torch::Device* device, const char* requested) {
    *device = KokkosDeviceToTorch<DeviceType>::convert();

    if (requested != nullptr) {
        auto requested_str = std::string(requested);
        std::transform(requested_str.begin(), requested_str.end(), requested_str.begin(), ::tolower);
        if (c10::DeviceTypeName(device->type(), /*lower_case=*/true) != requested_str) {
            error->all(FLERR,
                "requested device '{}' does not match the device being used by kokkos '{}', "
                "use the non-kokkos version of this pair style to use a different "
                "device for the model and LAMMPS",
                requested, device->str()
            );
        }
    }
}

template<class DeviceType>
void PairMetatomicKokkos<DeviceType>::compute(int eflag, int vflag) {
    if (std::getenv("LAMMPS_METATENSOR_PROFILE") != nullptr) {
        MetatomicTimer::enable(true);
    } else {
        MetatomicTimer::enable(false);
    }

    auto _ = MetatomicTimer("PairMetatomicKokkos::compute");

    /// Declare what we need to read from the atomKK object and what we will modify
    this->atomKK->sync(ExecutionSpaceFromDevice<DeviceType>::space, datamask_read);
    this->atomKK->modified(ExecutionSpaceFromDevice<DeviceType>::space, datamask_modify);

    if (eflag || vflag) {
        ev_setup(eflag, vflag);
    } else {
        evflag = vflag_fdotr = eflag_global = eflag_atom = 0;
    }

    if (eflag_atom) {
        mta_data->evaluation_options->outputs.at("energy")->per_atom = true;
    } else {
        mta_data->evaluation_options->outputs.at("energy")->per_atom = false;
    }

    auto dtype = torch::kFloat64;
    if (mta_data->capabilities->dtype() == "float64") {
        dtype = torch::kFloat64;
    } else if (mta_data->capabilities->dtype() == "float32") {
        dtype = torch::kFloat32;
    } else {
        error->all(FLERR, "the model requested an unsupported dtype '{}'", mta_data->capabilities->dtype());
    }

    // transform from LAMMPS to metatensor System
    auto system = this->system_adaptor->system_from_lmp(
        mta_list,
        static_cast<bool>(vflag_global),
        mta_data->remap_pairs,
        dtype,
        mta_data->device
    );

    // only run the calculation for atoms actually in the current domain
    mta_data->selected_atoms_values.resize_({atom->nlocal, 2});
    mta_data->selected_atoms_values.index_put_({torch::indexing::Slice(), 0}, 0);
    auto tensor_options = mta_data->selected_atoms_values.options();
    mta_data->selected_atoms_values.index_put_(
        {torch::indexing::Slice(), 1},
        torch::arange(atom->nlocal, tensor_options)
    );

    auto selected_atoms = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        std::vector<std::string>{"system", "atom"},
        mta_data->selected_atoms_values,
        metatensor::assume_unique{}
    );
    mta_data->evaluation_options->set_selected_atoms(selected_atoms);

    torch::IValue result_ivalue;
    try {
        auto _ = MetatomicTimer("running Model::forward");
        result_ivalue = mta_data->model->forward({
            std::vector<metatomic_torch::System>{system},
            mta_data->evaluation_options,
            mta_data->check_consistency
        });
    } catch (const std::exception& e) {
        error->all(FLERR, "error evaluating the torch model: {}", e.what());
    }

    auto result = result_ivalue.toGenericDict();
    auto energy = result.at("energy").toCustomClass<metatensor_torch::TensorMapHolder>();
    auto energy_block = metatensor_torch::TensorMapHolder::block_by_id(energy, 0);
    auto energy_tensor = energy_block->values();

    torch::Tensor forces_tensor;
    torch::Tensor virial_tensor;

    if (mta_data->non_conservative) {
        auto forces = result.at("non_conservative_forces").toCustomClass<metatensor_torch::TensorMapHolder>();;
        auto forces_block = metatensor_torch::TensorMapHolder::block_by_id(forces, 0);
        forces_tensor = forces_block->values().squeeze(-1);
        forces_tensor = forces_tensor.to(torch::kFloat64);
        auto stress = result.at("non_conservative_stress").toCustomClass<metatensor_torch::TensorMapHolder>();;
        auto stress_block = metatensor_torch::TensorMapHolder::block_by_id(stress, 0);
        auto stress_tensor = stress_block->values().squeeze(0).squeeze(-1);
        virial_tensor = - stress_tensor * torch::abs(torch::det(system->cell()));
        virial_tensor = virial_tensor.to(torch::kCPU).to(torch::kFloat64);
    } else {
        // compute forces/virial on device with backward propagation
        // reset gradients to zero before calling backward
        this->system_adaptor->positions.mutable_grad() = torch::Tensor();
        this->system_adaptor->strain.mutable_grad() = torch::Tensor();

        auto _ = MetatomicTimer("running Model::backward");
        energy_tensor.backward(-torch::ones_like(energy_tensor));

        forces_tensor = this->system_adaptor->positions.grad();
        virial_tensor = this->system_adaptor->strain.grad().to(torch::kCPU);
    }

    {
        auto _ = MetatomicTimer("storing model output in LAMMPS data structures");

        // move results to cpu for storing
        auto energy_detached = energy_tensor.detach().to(torch::kCPU).to(torch::kFloat64);
        auto energy_samples = energy_block->samples();

        // store the energy returned by the model
        torch::Tensor global_energy;
        if (eflag_atom) {
            assert(energy_samples->size() == 2);
            assert(energy_samples->names()[0] == "system");
            assert(energy_samples->names()[1] == "atom");

            auto samples_values = energy_samples->values().to(torch::kCPU);
            auto samples = samples_values.accessor<int32_t, 2>();

            int64_t n_atoms = atom->nlocal + atom->nghost;
            assert(samples_values.sizes() == mta_data->selected_atoms_values.sizes());

            auto energies = energy_detached.accessor<double, 2>();
            for (int64_t i=0; i<energy_samples->count(); i++) {
                assert(samples[i][0] == 0);
                // handle potentially out of order samples in
                // the per-atom energy tensor
                auto atom_i = samples[i][1];
                assert(atom_i < n_atoms);
                eatom[atom_i] += this->scale * energies[i][0];
            }

            global_energy = energy_detached.sum(0);
            assert(energy_detached.sizes() == std::vector<int64_t>({1}));
        } else {
            assert(energy_samples->size() == 1);
            assert(energy_samples->names()[0] == "system");

            assert(energy_detached.sizes() == std::vector<int64_t>({1, 1}));
            global_energy = energy_detached.reshape({1});
        }

        if (eflag_global) {
            eng_vdwl += this->scale * global_energy.item<double>();
        }

        // store forces/virial
        assert(forces_tensor.scalar_type() == torch::kFloat64);
        forces_tensor = forces_tensor.contiguous();

        auto forces_lammps_kk = this->atomKK->k_f.template view<DeviceType>();
        auto forces_metatensor_kk = UnmanagedView<double**, DeviceType>(
            forces_tensor.template data_ptr<double>(),
            forces_tensor.size(0), 3
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

        assert(!vflag_fdotr);

        if (vflag_global) {
            assert(virial_tensor.is_cpu() && virial_tensor.scalar_type() == torch::kFloat64);
            auto predicted_virial = virial_tensor.template accessor<double, 2>();

            virial[0] += this->scale * predicted_virial[0][0];
            virial[1] += this->scale * predicted_virial[1][1];
            virial[2] += this->scale * predicted_virial[2][2];

            virial[3] += this->scale * 0.5 * (predicted_virial[1][0] + predicted_virial[0][1]);
            virial[4] += this->scale * 0.5 * (predicted_virial[2][0] + predicted_virial[0][2]);
            virial[5] += this->scale * 0.5 * (predicted_virial[2][1] + predicted_virial[1][2]);
        }

        if (vflag_atom) {
            error->all(FLERR, "per atom virial is not implemented");
        }
    }
}

namespace LAMMPS_NS {
template class PairMetatomicKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairMetatomicKokkos<LMPHostType>;
#endif
}
