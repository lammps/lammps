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
#include "metatomic_system_kokkos.h"

#include "metatomic_timer.h"

#include "domain.h"
#include "error.h"

#include "atom_kokkos.h"

#include <torch/cuda.h>

using namespace LAMMPS_NS;

template<typename T, class DeviceType>
using UnmanagedView = Kokkos::View<T, Kokkos::LayoutRight, DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template<class DeviceType>
MetatomicSystemAdaptorKokkos<DeviceType>::MetatomicSystemAdaptorKokkos(LAMMPS *lmp, MetatomicSystemOptions options):
    MetatomicSystemAdaptor(lmp, options),
    device_(KokkosDeviceToTorch<DeviceType>::convert())
{
    // MetatomicSystemAdaptor allocate on CPU, move to the right device
    this->atomic_types_ = this->atomic_types_.to(this->device_);

    auto tensor_options = torch::TensorOptions()
        .dtype(torch::kFloat64)
        .device(this->device_)
        .requires_grad(options_.requires_grad);

    this->strain = torch::eye(3, tensor_options);
}

#include "comm.h"

template<class DeviceType>
void MetatomicSystemAdaptorKokkos<DeviceType>::setup_neighbors_remap_kk(metatomic_torch::System& system, NeighListKokkos<DeviceType>* list) {
    auto _ = MetatomicTimer("converting kokkos neighbors with ghosts remapping");
    auto dtype = system->positions().scalar_type();

    auto total_n_atoms = atomKK->nlocal + atomKK->nghost;

    {
        auto _ = MetatomicTimer("identifying ghosts and real atoms");
        /*-------------- this will be done on CPU for now ------------------------*/
        // The hashmap in the following code is not easy to implement in either Kokkos or torch
        // The cost of this section seems to be very low anyway

        // Collect the local atom id of all local & ghosts atoms, mapping ghosts
        // atoms which are periodic images of local atoms back to the local atoms.
        //
        // Metatomic expects pairs corresponding to periodic atoms to be between
        // the main atoms, but using the actual distance vector between the atom and
        // the ghost.
        original_atom_id_.clear();
        original_atom_id_.reserve(total_n_atoms);

        // identify all local atom by their LAMMPS atom tag.
        local_atoms_tags_.clear();
        for (int i=0; i<atom->nlocal; i++) {
            original_atom_id_.emplace_back(i);
            local_atoms_tags_.emplace(atom->tag[i], i);
        }

        // now loop over ghosts & map them back to the main cell if needed
        ghost_atoms_tags_.clear();
        for (int i=atom->nlocal; i<total_n_atoms; i++) {
            auto tag = atom->tag[i];
            auto it = local_atoms_tags_.find(tag);
            if (it != local_atoms_tags_.end()) {
                // this is the periodic image of an atom already owned by this domain
                original_atom_id_.emplace_back(it->second);
            } else {
                // this can either be a periodic image of an atom owned by another
                // domain, or directly an atom from another domain. Since we can not
                // really distinguish between these, we take the first atom as the
                // "main" one and remap all atoms with the same tag to the first one
                auto it = ghost_atoms_tags_.find(tag);
                if (it != ghost_atoms_tags_.end()) {
                    // we already found this atom elsewhere in the system
                    original_atom_id_.emplace_back(it->second);
                } else {
                    // this is the first time we are seeing this atom
                    original_atom_id_.emplace_back(i);
                    ghost_atoms_tags_.emplace(tag, i);
                }
            }
        }
    }
    /*----------- end of "this will be done on CPU for now" --------------*/

    auto original_id = torch::from_blob(
        original_atom_id_.data(),
        {total_n_atoms},
        torch::TensorOptions().dtype(torch::kInt32).device(torch::kCPU)
    ).to(this->device_);

    auto neighbors_kk = list->d_neighbors;
    auto max_number_of_neighbors = list->maxneighs;

    auto neighbors = torch::zeros(
        {total_n_atoms, max_number_of_neighbors},
        torch::TensorOptions().dtype(torch::kInt32).device(this->device_)
    );
    // mask neighbors_kk with NEIGHMASK. Torch doesn't have this functionality, we do it in Kokkos
    auto neighbors_kk_masked = UnmanagedView<int32_t**, DeviceType>(
        neighbors.template data_ptr<int32_t>(),
        total_n_atoms,
        max_number_of_neighbors
    );
    Kokkos::parallel_for(
        Kokkos::MDRangePolicy({0, 0}, {total_n_atoms, max_number_of_neighbors}),
        KOKKOS_LAMBDA(size_t i, size_t j) {
            neighbors_kk_masked(i, j) = neighbors_kk(i, j) & NEIGHMASK;
        }
    );

    // Convert NL-related data to torch tensors
    auto numneigh = torch::from_blob(
        list->d_numneigh.data(),
        {total_n_atoms},
        torch::TensorOptions().dtype(torch::kInt32).device(this->device_)
    );
    auto ilist = torch::from_blob(
        list->d_ilist.data(),
        {total_n_atoms},
        torch::TensorOptions().dtype(torch::kInt32).device(this->device_)
    );

    auto x = system->positions().detach();
    auto cell_inverse = system->cell().detach().inverse();

    // convert from LAMMPS NL format to metatomic NL format
    auto expanded_arange = torch::arange(
        max_number_of_neighbors,
        torch::TensorOptions().dtype(torch::kInt32).device(this->device_)
    ).unsqueeze(0).expand({total_n_atoms, -1});
    auto neighbor_2d_mask = expanded_arange < numneigh.unsqueeze(1);

    auto expanded_arange_other_dim = torch::arange(
        total_n_atoms,
        torch::TensorOptions().dtype(torch::kInt32).device(this->device_)
    ).unsqueeze(1).expand({-1, max_number_of_neighbors});
    auto index_for_ilist = expanded_arange_other_dim.masked_select(neighbor_2d_mask);

    auto centers_id = ilist.index_select(0, index_for_ilist);
    auto neighbors_id = neighbors.masked_select(neighbor_2d_mask);

    // change centers and neighbors to the original atom ids
    auto centers_original_id = original_id.index_select(0, centers_id);
    auto neighbors_original_id = original_id.index_select(0, neighbors_id);

    // The following code is a direct translation of the code in the non-Kokkos
    // version (MetatomicSystemAdaptor::setup_neighbors_remap), but rewritten
    // in torch to use the GPU
    for (auto& cache: caches_) {
        // current values of various tensors, these change depending on full/half setting
        torch::Tensor centers_id_cur;
        torch::Tensor neighbors_id_cur;
        torch::Tensor centers_original_id_cur;
        torch::Tensor neighbors_original_id_cur;

        // filtered tensors, i.e. only containing pairs actually below the cutoff
        torch::Tensor centers_original_id_filt_cur;
        torch::Tensor neighbors_original_id_filt_cur;
        torch::Tensor distances_filt_cur;
        torch::Tensor cell_shifts_cur;

        // other tensors that need to live across multiple timed sections
        torch::Tensor samples_indices;
        torch::Tensor samples_values;
        {
            auto _ = MetatomicTimer("filtering LAMMPS neighbor list");
            // half list mask, if necessary
            auto full_list = cache.options->full_list();

            if (full_list) {
                centers_id_cur = centers_id;
                neighbors_id_cur = neighbors_id;
                centers_original_id_cur = centers_original_id;
                neighbors_original_id_cur = neighbors_original_id;
            } else {
                auto half_list_mask = centers_original_id <= neighbors_original_id;
                centers_id_cur = centers_id.masked_select(half_list_mask);
                neighbors_id_cur = neighbors_id.masked_select(half_list_mask);
                centers_original_id_cur = centers_original_id.masked_select(half_list_mask);
                neighbors_original_id_cur = neighbors_original_id.masked_select(half_list_mask);
            }

            // distance mask
            auto distances = x.index_select(0, neighbors_id_cur) - x.index_select(0, centers_id_cur);
            auto cutoff_mask = torch::sum(distances.pow(2), 1) < cache.cutoff*cache.cutoff;

            // index everything with the mask
            auto centers_original_id_filt = centers_original_id_cur.masked_select(cutoff_mask);
            auto neighbors_original_id_filt = neighbors_original_id_cur.masked_select(cutoff_mask);
            auto distances_filt = distances.index({cutoff_mask, torch::indexing::Slice()});

            // find filtered interatomic vectors using the original atoms
            auto original_distances_filtered =
                x.index_select(0, neighbors_original_id_filt)
                - x.index_select(0, centers_original_id_filt);

            // cell shifts
            auto pair_shifts = distances_filt - original_distances_filtered;
            auto cell_shifts = pair_shifts.matmul(cell_inverse);
            cell_shifts = torch::round(cell_shifts).to(torch::kInt32);

            if (full_list) {
                centers_original_id_filt_cur = centers_original_id_filt;
                neighbors_original_id_filt_cur = neighbors_original_id_filt;
                distances_filt_cur = distances_filt;
                cell_shifts_cur = cell_shifts;
            } else {
                auto half_list_cell_mask = centers_original_id_filt > neighbors_original_id_filt;
                auto pair_with_image_mask = centers_original_id_filt == neighbors_original_id_filt;
                auto negative_half_space_mask = torch::sum(cell_shifts, 1) < 0;
                // reproduce this mask (from MetatomicSystemAdaptor::setup_neighbors_remap) with torch:
                // if ((shift[0] + shift[1] + shift[2] == 0) && (shift[2] < 0 || (shift[2] == 0 && shift[1] < 0)))
                auto edge_mask = (
                    (torch::sum(cell_shifts, 1) == 0) & (
                        (cell_shifts.index({torch::indexing::Slice(), 2}) < 0) | (
                            cell_shifts.index({torch::indexing::Slice(), 2}) == 0 &
                            cell_shifts.index({torch::indexing::Slice(), 1}) < 0
                        )
                    )
                );
                auto final_mask = torch::logical_not(
                    half_list_cell_mask | (
                        pair_with_image_mask & (negative_half_space_mask | edge_mask)
                    )
                );
                centers_original_id_filt_cur = centers_original_id_filt.masked_select(final_mask);
                neighbors_original_id_filt_cur = neighbors_original_id_filt.masked_select(final_mask);
                distances_filt_cur = distances_filt.index({final_mask, torch::indexing::Slice()});
                cell_shifts_cur = cell_shifts.index({final_mask, torch::indexing::Slice()});
            }

            // make sure all the sample are unique
            samples_values = torch::concatenate({
                centers_original_id_filt_cur.unsqueeze(-1),
                neighbors_original_id_filt_cur.unsqueeze(-1),
                cell_shifts_cur
            }, /*dim=*/1);

            auto [samples_values_unique, samples_inverse, _counts] = torch::unique_dim(
                samples_values, /*dim=*/0, /*sorted=*/true, /*return_inverse=*/true, /*return_counts=*/false
            );
            samples_values = samples_values_unique;

            auto permutation = torch::arange(samples_inverse.size(0), samples_inverse.options());
            samples_inverse = samples_inverse.flip({0});
            permutation = permutation.flip({0});

            samples_indices = torch::empty(samples_values.size(0), samples_inverse.options());
            samples_indices.scatter_(0, samples_inverse, permutation);
        }

        // wrap into metatensor data structures
        torch::intrusive_ptr<metatensor_torch::LabelsHolder> samples;
        {
            auto n_pairs = samples_values.size(0);
            auto _ = MetatomicTimer("creating samples Labels (" +  std::to_string(n_pairs) + " pairs)");
            samples = torch::make_intrusive<metatensor_torch::LabelsHolder>(
                std::vector<std::string>{"first_atom", "second_atom", "cell_shift_a", "cell_shift_b", "cell_shift_c"},
                samples_values,
                metatensor::assume_unique{}
            );
        }

        torch::intrusive_ptr<metatensor_torch::TensorBlockHolder> neighbors;
        {
            auto _ = MetatomicTimer("creating neighbors TensorBlock");

            neighbors = torch::make_intrusive<metatensor_torch::TensorBlockHolder>(
                distances_filt_cur.index_select(0, samples_indices).unsqueeze(-1),
                samples,
                std::vector<metatensor_torch::Labels>{
                    metatensor_torch::LabelsHolder::create({"xyz"}, {{0}, {1}, {2}})->to(this->device_),
                },
                metatensor_torch::LabelsHolder::create({"distance"}, {{0}})->to(this->device_)
            );
        }

        metatomic_torch::register_autograd_neighbors(system, neighbors, options_.check_consistency);
        system->add_neighbor_list(cache.options, neighbors);
    }
}


template<class DeviceType>
metatomic_torch::System MetatomicSystemAdaptorKokkos<DeviceType>::system_from_lmp(
    NeighList* list,
    bool do_virial,
    bool remap_pairs,
    torch::ScalarType dtype,
    torch::Device device
) {
    auto _ = MetatomicTimer("creating System from LAMMPS-kokkos data");
    assert(device == this->device_);

    auto total_n_atoms = atomKK->nlocal + atomKK->nghost;

    atomic_types_.resize_({total_n_atoms});
    auto atomic_types_kk = UnmanagedView<int32_t*, DeviceType>(atomic_types_.data_ptr<int32_t>(), total_n_atoms);
    auto type_mapping_kk = UnmanagedView<int32_t*, DeviceType>(options_.types_mapping, atomKK->ntypes + 1);
    auto types_kk = atomKK->k_type.view<DeviceType>();
    Kokkos::parallel_for(total_n_atoms, KOKKOS_LAMBDA(int i) {
        atomic_types_kk(i) = type_mapping_kk(types_kk(i));
    });

    // atomKK->k_x contains "real" and then ghost atoms, in that order
    auto k_x = atomKK->k_x.view<DeviceType>();
    auto tensor_options = torch::TensorOptions().dtype(torch::kFloat64).device(this->device_);

    this->positions = torch::from_blob(
        k_x.data(), {total_n_atoms, 3},
        // requires_grad=true since we always need gradients w.r.t. positions
        tensor_options.requires_grad(options_.requires_grad)
    );

    auto cell = torch::zeros({3, 3}, tensor_options);
    cell[0][0] = domain->xprd;

    cell[1][0] = domain->xy;
    cell[1][1] = domain->yprd;

    cell[2][0] = domain->xz;
    cell[2][1] = domain->yz;
    cell[2][2] = domain->zprd;

    auto system_positions = this->positions.to(dtype);
    cell = cell.to(dtype);

    if (do_virial) {
        auto model_strain = this->strain.to(dtype);

        // pretend to scale positions/cell by the strain so that
        // it enters the computational graph.
        system_positions = system_positions.matmul(model_strain);
        cell = cell.matmul(model_strain);
    }

    // Periodic boundary conditions handling.
    //
    // While Metatomic models can support mixed PBC settings, we currently
    // assume that the system is fully periodic and we throw an error otherwise
    if (!domain->xperiodic || !domain->yperiodic || !domain->zperiodic) {
        error->one(FLERR, "metatomic/kk currently requires a fully periodic system");
    }
    auto pbc = torch::tensor(
        {domain->xperiodic, domain->yperiodic, domain->zperiodic},
        torch::TensorOptions().dtype(torch::kBool).device(this->device_)
    );

    auto system = torch::make_intrusive<metatomic_torch::SystemHolder>(
        atomic_types_,
        system_positions,
        cell,
        pbc
    );

    if (remap_pairs) {
        auto* kk_list = dynamic_cast<NeighListKokkos<DeviceType>*>(list);
        assert(kk_list != nullptr);
        this->setup_neighbors_remap_kk(system, kk_list);
    } else {
        error->one(FLERR, "the kokkos version of metatomic requires remap_pairs to be true");
    }

    return system;
}

namespace LAMMPS_NS {
template class MetatomicSystemAdaptorKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class MetatomicSystemAdaptorKokkos<LMPHostType>;
#endif
}
