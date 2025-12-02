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
#include "memory_kokkos.h"
#include "metatomic_timer.h"

#include "atom_masks.h"
#include "atom_kokkos.h"
#include "comm.h"
#include "domain.h"
#include "error.h"

#include <torch/cuda.h>

using namespace LAMMPS_NS;

template<typename T, class DeviceType>
using UnmanagedView = Kokkos::View<T, Kokkos::LayoutRight, DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

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

template<class DeviceType>
void MetatomicSystemAdaptorKokkos<DeviceType>::add_nl_request(double cutoff, metatomic_torch::NeighborListOptions request) {
    if (cutoff > options_.interaction_range) {
        error->one(FLERR,
            "Invalid metatomic model: one of the requested neighbor lists "
            "has a cutoff ({}) larger than the model interaction range ({})",
            cutoff, options_.interaction_range
        );
    } else if (cutoff < 0 || !std::isfinite(cutoff)) {
        error->one(FLERR,
            "model requested an invalid cutoff for neighbors list: {} "
            "(cutoff in model units is {})",
            cutoff, request->cutoff()
        );
    }

    nl_requests_kk_.push_back({
        cutoff,
        request,
        /*samples = */ {},
        /*distances_f64 = */ {},
        /*distances_f32 = */ {},
    });
}


template<class DeviceType>
void MetatomicSystemAdaptorKokkos<DeviceType>::setup_neighbors_kk(metatomic_torch::System& system, NeighListKokkos<DeviceType>* list) {
    auto _ = MetatomicTimer("converting kokkos neighbors list");

    auto dtype = system->positions().scalar_type();
    auto total_n_atoms = atomKK->nlocal + atomKK->nghost;
    auto max_number_of_neighbors = list->maxneighs;

    auto original_id = UnmanagedView<int*, LMPHostType>(
        original_atom_id_.data(),
        original_atom_id_.size()
    );
    auto d_original_id = Kokkos::View<int*, Kokkos::LayoutRight, LMPDeviceType>("", original_atom_id_.size());
    Kokkos::deep_copy(d_original_id, original_id);

    auto lmp_to_mta = UnmanagedView<int*, LMPHostType>(
        lmp_to_mta_.data(),
        lmp_to_mta_.size()
    );
    auto d_lmp_to_mta = Kokkos::View<int*, Kokkos::LayoutRight, LMPDeviceType>("", lmp_to_mta_.size());
    Kokkos::deep_copy(d_lmp_to_mta, lmp_to_mta);


    auto cell_inv = this->cell_inverse();
    auto d_cell_inv = Kokkos::View<double**, Kokkos::LayoutRight, LMPDeviceType>("cell_inv", 3, 3);
    Kokkos::deep_copy(
        d_cell_inv,
        UnmanagedView<double**, LMPHostType>(reinterpret_cast<double*>(cell_inv.data()), 3, 3)
    );

    auto x = atomKK->k_x.view<DeviceType>();
    auto cell_shifts = KOKKOS_LAMBDA(
        const Kokkos::View<double**, Kokkos::LayoutRight, LMPDeviceType>& cell_inv,
        const double pair_shift[3],
        int32_t shift_out[3]
    ) {
        shift_out[0] = static_cast<int32_t>(std::round(
            cell_inv(0, 0) * pair_shift[0] +
            cell_inv(0, 1) * pair_shift[1] +
            cell_inv(0, 2) * pair_shift[2]
        ));
        shift_out[1] = static_cast<int32_t>(std::round(
            cell_inv(1, 0) * pair_shift[0] +
            cell_inv(1, 1) * pair_shift[1] +
            cell_inv(1, 2) * pair_shift[2]
        ));
        shift_out[2] = static_cast<int32_t>(std::round(
            cell_inv(2, 0) * pair_shift[0] +
            cell_inv(2, 1) * pair_shift[1] +
            cell_inv(2, 2) * pair_shift[2]
        ));
    };

    for (auto& nl: nl_requests_kk_) {
        auto cutoff2 = nl.cutoff * nl.cutoff;
        auto full_list = nl.options->full_list();

        auto cutoff_ratio = nl.cutoff / options_.interaction_range;
        size_t max_n_pairs = total_n_atoms * list->maxneighs;
        // Allocate for a much smaller number of pairs to avoid wasting memory
        // (especially when the interaction range is much larger than the NL
        // cutoff)
        auto pairs_capacity = std::max(
            std::max(
                static_cast<size_t>(max_n_pairs * cutoff_ratio / 1000),
                100ul
            ),
            nl.samples.extent(0)
        );

        // actual number of pairs we found
        auto n_pairs = Kokkos::DualView<int64_t, LMPDeviceType>("n_pairs");
        auto d_n_pairs = n_pairs.view_device();
        auto h_n_pairs = n_pairs.view_host();

        auto processed_all_pairs = Kokkos::DualView<bool, LMPDeviceType>("processed_all_pairs");
        auto d_processed_all_pairs = processed_all_pairs.view_device();
        auto h_processed_all_pairs = processed_all_pairs.view_host();

        Kokkos::deep_copy(h_processed_all_pairs, false);
        processed_all_pairs.template modify<LMPHostType>();
        processed_all_pairs.template sync<LMPDeviceType>();

        while (!h_processed_all_pairs()) {
            {
                auto _ = MetatomicTimer("allocating caches for kokkos neighbors list (capacity for " + std::to_string(pairs_capacity) + " pairs)");
                if (nl.samples.extent(0) < pairs_capacity) {
                    MemoryKokkos::realloc_kokkos(nl.samples, "metatomic:samples", pairs_capacity, 5);
                }

                if (dtype == torch::kFloat64) {
                    if (nl.distances_f64.extent(0) < pairs_capacity) {
                        MemoryKokkos::realloc_kokkos(nl.distances_f64, "metatomic:distances_f64", pairs_capacity, 3);
                    }
                } else if (dtype == torch::kFloat32) {
                    if (nl.distances_f32.extent(0) < pairs_capacity) {
                        MemoryKokkos::realloc_kokkos(nl.distances_f32, "metatomic:distances_f32", pairs_capacity, 3);
                    }
                } else {
                    // should be unreachable
                    error->one(FLERR, "invalid dtype, this is a bug");
                }
            }

            auto _ = MetatomicTimer("filtering kokkos neighbors list");

            Kokkos::deep_copy(h_n_pairs, 0);
            n_pairs.template modify<LMPHostType>();
            n_pairs.template sync<LMPDeviceType>();

            Kokkos::deep_copy(h_processed_all_pairs, true);
            processed_all_pairs.template modify<LMPHostType>();
            processed_all_pairs.template sync<LMPDeviceType>();

            Kokkos::parallel_for(
                Kokkos::MDRangePolicy<DeviceType, Kokkos::Rank<2>>(
                    {0, 0},
                    {list->inum + list->gnum, max_number_of_neighbors}
                ),
                KOKKOS_LAMBDA(size_t ii, size_t jj) {
                    if (jj >= list->d_numneigh[ii]) {
                        return;
                    }

                    auto atom_i = list->d_ilist[ii];
                    auto original_atom_i = d_original_id[atom_i];
                    auto i_is_original = (atom_i == original_atom_i);

                    auto atom_j = list->d_neighbors(ii, jj) & NEIGHMASK;
                    auto original_atom_j = d_original_id[atom_j];
                    auto j_is_original = (atom_j == original_atom_j);

                    if (!full_list && original_atom_i > original_atom_j) {
                        // Remove extra pairs if the model requested half-lists
                        return;
                    }

                    if (!i_is_original && !j_is_original) {
                        // both atoms are periodic ghosts, skip the pair
                        return;
                    }

                    if (!i_is_original && j_is_original) {
                        // this pair will be accounted for when we will process
                        // atom_j as the central atom
                        return;
                    }

                    double distance[3] = {
                        x(atom_j, 0) - x(atom_i, 0),
                        x(atom_j, 1) - x(atom_i, 1),
                        x(atom_j, 2) - x(atom_i, 2),
                    };

                    auto distance2 = (
                        distance[0] * distance[0] +
                        distance[1] * distance[1] +
                        distance[2] * distance[2]
                    );
                    if (distance2 > cutoff2) {
                        // LAMMPS neighbors list contains some pairs after the
                        // cutoff, we filter them here
                        return;
                    }

                    // Compute the cell shift for the pair.
                    double shift_i[3] = {
                        x(atom_i, 0) - x(original_atom_i, 0),
                        x(atom_i, 1) - x(original_atom_i, 1),
                        x(atom_i, 2) - x(original_atom_i, 2),
                    };
                    double shift_j[3] = {
                        x(atom_j, 0) - x(original_atom_j, 0),
                        x(atom_j, 1) - x(original_atom_j, 1),
                        x(atom_j, 2) - x(original_atom_j, 2),
                    };
                    double pair_shift[3] = {
                        shift_j[0] - shift_i[0],
                        shift_j[1] - shift_i[1],
                        shift_j[2] - shift_i[2],
                    };

                    int32_t shift[3] = {0, 0, 0};
                    if (pair_shift[0] != 0 || pair_shift[1] != 0 || pair_shift[2] != 0) {
                        cell_shifts(d_cell_inv, pair_shift, shift);

                        if (!full_list && original_atom_i == original_atom_j) {
                            // If a half neighbors list has been requested, do
                            // not include the same pair between an atom and
                            // it's periodic image twice with opposite cell
                            // shifts (e.g. [1, -1, 1] and [-1, 1, -1]).
                            //
                            // Instead we pick pairs in the positive plan of
                            // shifts.
                            if (shift[0] + shift[1] + shift[2] < 0) {
                                // drop shifts on the negative half-space
                                return;
                            }

                            if ((shift[0] + shift[1] + shift[2] == 0)
                                && (shift[2] < 0 || (shift[2] == 0 && shift[1] < 0))) {
                                // drop shifts in the negative half plane or the
                                // negative shift[1] axis.
                                //
                                // See below for a graphical representation: we are
                                // keeping the shifts indicated with `O` and
                                // dropping the ones indicated with `X`
                                //
                                //  O O O │ O O O
                                //  O O O │ O O O
                                //  O O O │ O O O
                                // ─X─X─X─┼─O─O─O─
                                //  X X X │ X X X
                                //  X X X │ X X X
                                //  X X X │ X X X
                                return;
                            }
                        }
                    }

                    auto pair_i = Kokkos::atomic_fetch_add(&d_n_pairs(), 1);
                    if (pair_i >= nl.samples.extent(0)) {
                        // stop and re-allocate larger arrays
                        d_processed_all_pairs() = false;
                        return;
                    }

                    nl.samples(pair_i, 0) = d_lmp_to_mta[original_atom_i];
                    nl.samples(pair_i, 1) = d_lmp_to_mta[original_atom_j];
                    nl.samples(pair_i, 2) = shift[0];
                    nl.samples(pair_i, 3) = shift[1];
                    nl.samples(pair_i, 4) = shift[2];
                    if (dtype == torch::kFloat64) {
                        nl.distances_f64(pair_i, 0) = distance[0];
                        nl.distances_f64(pair_i, 1) = distance[1];
                        nl.distances_f64(pair_i, 2) = distance[2];
                    } else {
                        assert(dtype == torch::kFloat32);
                        nl.distances_f32(pair_i, 0) = static_cast<float>(distance[0]);
                        nl.distances_f32(pair_i, 1) = static_cast<float>(distance[1]);
                        nl.distances_f32(pair_i, 2) = static_cast<float>(distance[2]);
                    }
                }
            );

            processed_all_pairs.template modify<LMPDeviceType>();
            processed_all_pairs.template sync<LMPHostType>();
            if (!h_processed_all_pairs()) {
                pairs_capacity *= 2;
            }
        }

        n_pairs.template modify<LMPDeviceType>();
        n_pairs.template sync<LMPHostType>();

        auto samples_values = torch::from_blob(
            nl.samples.data(), {h_n_pairs(), 5}, torch::TensorOptions().dtype(torch::kInt32).device(this->device_)
        );

        torch::Tensor distances;
        if (dtype == torch::kFloat64) {
            distances = torch::from_blob(
                nl.distances_f64.data(), {h_n_pairs(), 3, 1}, torch::TensorOptions().dtype(dtype).device(this->device_)
            );
        } else if (dtype == torch::kFloat32) {
            distances = torch::from_blob(
                nl.distances_f32.data(), {h_n_pairs(), 3, 1}, torch::TensorOptions().dtype(dtype).device(this->device_)
            );
        } else {
            // should be unreachable
            error->one(FLERR, "invalid dtype, this is a bug");
        }

        // wrap into metatensor data structures
        torch::intrusive_ptr<metatensor_torch::LabelsHolder> samples;
        {
            auto _ = MetatomicTimer("creating samples Labels (" +  std::to_string(h_n_pairs()) + " pairs)");

            if (options_.check_consistency) {
                samples = torch::make_intrusive<metatensor_torch::LabelsHolder>(
                    std::vector<std::string>{"first_atom", "second_atom", "cell_shift_a", "cell_shift_b", "cell_shift_c"},
                    samples_values
                );
            } else {
                samples = torch::make_intrusive<metatensor_torch::LabelsHolder>(
                    std::vector<std::string>{"first_atom", "second_atom", "cell_shift_a", "cell_shift_b", "cell_shift_c"},
                    samples_values,
                    metatensor::assume_unique{}
                );
            }
        }

        torch::intrusive_ptr<metatensor_torch::TensorBlockHolder> neighbors;
        {
            auto _ = MetatomicTimer("creating neighbors TensorBlock");

            neighbors = torch::make_intrusive<metatensor_torch::TensorBlockHolder>(
                distances,
                samples,
                std::vector<metatensor_torch::Labels>{
                    metatensor_torch::LabelsHolder::create({"xyz"}, {{0}, {1}, {2}})->to(this->device_),
                },
                metatensor_torch::LabelsHolder::create({"distance"}, {{0}})->to(this->device_)
            );
        }

        {
            auto _ = MetatomicTimer("registering neighbors in System");
            metatomic_torch::register_autograd_neighbors(system, neighbors, options_.check_consistency);
            system->add_neighbor_list(nl.options, neighbors);
        }
    }
}


template<class DeviceType>
metatomic_torch::System MetatomicSystemAdaptorKokkos<DeviceType>::system_from_lmp(
    NeighList* list,
    bool do_virial,
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

    this->positions = torch::from_blob(k_x.data(), {total_n_atoms, 3}, tensor_options);

    auto cell = torch::zeros({3, 3}, tensor_options);
    cell[0][0] = domain->xprd;

    cell[1][0] = domain->xy;
    cell[1][1] = domain->yprd;

    cell[2][0] = domain->xz;
    cell[2][1] = domain->yz;
    cell[2][2] = domain->zprd;

    // Periodic boundary conditions handling.
    auto pbc = torch::tensor(
        {domain->xperiodic, domain->yperiodic, domain->zperiodic},
        torch::TensorOptions().dtype(torch::kBool).device(this->device_)
    );

    cell.index_put_(
        {torch::logical_not(pbc)},
        torch::tensor({0.0}, tensor_options)
    );

    // make sure to sync the updated tags to host
    atomKK->sync(ExecutionSpaceFromDevice<LMPHostType>::space, TAG_MASK);
    this->guess_periodic_ghosts();

    // Only keep the atoms which are not periodic images of other atoms
    auto mta_to_lmp_tensor = torch::from_blob(
        mta_to_lmp.data(),
        {static_cast<int64_t>(mta_to_lmp.size())},
        torch::TensorOptions().dtype(torch::kInt).device(torch::kCPU)
    ).to(this->device_);
    this->atomic_types_ = this->atomic_types_.index_select(0, mta_to_lmp_tensor);
    this->positions = this->positions.index_select(0, mta_to_lmp_tensor);

    this->positions.set_requires_grad(options_.requires_grad);

    auto system_positions = this->positions.to(dtype).to(device);
    cell = cell.to(dtype);

    if (do_virial) {
        auto model_strain = this->strain.to(dtype);

        // scale positions/cell by the strain so that it enters the
        // computational graph.
        system_positions = system_positions.matmul(model_strain);
        cell = cell.matmul(model_strain);
    }

    auto system = torch::make_intrusive<metatomic_torch::SystemHolder>(
        atomic_types_,
        system_positions,
        cell,
        pbc
    );

    auto* kk_list = dynamic_cast<NeighListKokkos<DeviceType>*>(list);
    assert(kk_list != nullptr);
    this->setup_neighbors_kk(system, kk_list);

    return system;
}

namespace LAMMPS_NS {
template class MetatomicSystemAdaptorKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class MetatomicSystemAdaptorKokkos<LMPHostType>;
#endif
}
