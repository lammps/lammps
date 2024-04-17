/* ----------------------------------------------------------------------
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
   Contributing authors: Guillaume Fraux <guillaume.fraux@epfl.ch>
------------------------------------------------------------------------- */
#include "metatensor_system.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "neighbor.h"

#include "neigh_list.h"
#include "neigh_request.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MetatensorSystemAdaptor::MetatensorSystemAdaptor(LAMMPS *lmp, Pair* requestor, MetatensorSystemOptions options):
    Pointers(lmp),
    list_(nullptr),
    options_(std::move(options)),
    caches_(),
    atomic_types_(torch::zeros({0}, torch::TensorOptions().dtype(torch::kInt32).device(torch::kCPU)))
{
    // We ask LAMMPS for a full neighbor lists because we need to know about
    // ALL pairs, even if options->full_list() is false. We will then filter
    // the pairs to only include each pair once where needed.
    auto request = neighbor->add_request(requestor, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);
    request->set_id(0);
    request->set_cutoff(options_.interaction_range);

    this->strain = torch::eye(3, torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU).requires_grad(true));
}

MetatensorSystemAdaptor::MetatensorSystemAdaptor(LAMMPS *lmp, Compute* requestor, MetatensorSystemOptions options):
    Pointers(lmp),
    list_(nullptr),
    options_(std::move(options)),
    caches_(),
    atomic_types_(torch::zeros({0}, torch::TensorOptions().dtype(torch::kInt32).device(torch::kCPU)))
{
    auto request = neighbor->add_request(requestor, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);
    request->set_id(0);
    request->set_cutoff(options_.interaction_range);

    this->strain = torch::eye(3, torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU).requires_grad(true));
}

MetatensorSystemAdaptor::~MetatensorSystemAdaptor() {

}

void MetatensorSystemAdaptor::init_list(int id, NeighList* ptr) {
    assert(id == 0);
    list_ = ptr;
}

void MetatensorSystemAdaptor::add_nl_request(double cutoff, metatensor_torch::NeighborsListOptions request) {
    if (cutoff > options_.interaction_range) {
        error->all(FLERR,
            "Invalid metatensor model: one of the requested neighbor lists "
            "has a cutoff ({}) larger than the model interaction range ({})",
            cutoff, options_.interaction_range
        );
    } else if (cutoff < 0 || !std::isfinite(cutoff)) {
        error->all(FLERR,
            "model requested an invalid cutoff for neighbors list: {} "
            "(cutoff in model units is {})",
            cutoff, request->cutoff()
        );
    }

    caches_.push_back({
        cutoff,
        request,
        /*known_samples = */ {},
        /*samples = */ {},
        /*distances_f64 = */ {},
        /*distances_f32 = */ {},
    });
}


static std::array<int32_t, 3> cell_shifts(
    const std::array<std::array<double, 3>, 3>& cell_inv,
    const std::array<double, 3>& pair_shift
) {
    auto shift_a = static_cast<int32_t>(std::round(
        cell_inv[0][0] * pair_shift[0] +
        cell_inv[0][1] * pair_shift[1] +
        cell_inv[0][2] * pair_shift[2]
    ));
    auto shift_b = static_cast<int32_t>(std::round(
        cell_inv[1][0] * pair_shift[0] +
        cell_inv[1][1] * pair_shift[1] +
        cell_inv[1][2] * pair_shift[2]
    ));
    auto shift_c = static_cast<int32_t>(std::round(
        cell_inv[2][0] * pair_shift[0] +
        cell_inv[2][1] * pair_shift[1] +
        cell_inv[2][2] * pair_shift[2]
    ));

    return {shift_a, shift_b, shift_c};
}


void MetatensorSystemAdaptor::setup_neighbors(metatensor_torch::System& system) {
    auto dtype = system->positions().scalar_type();
    auto device = system->positions().device();

    double** x = atom->x;
    auto total_n_atoms = atom->nlocal + atom->nghost;

    auto cell_inv_tensor = system->cell().inverse().t().to(torch::kCPU).to(torch::kFloat64);
    auto cell_inv_accessor = cell_inv_tensor.accessor<double, 2>();
    auto cell_inv = std::array<std::array<double, 3>, 3>{{
        {{cell_inv_accessor[0][0], cell_inv_accessor[0][1], cell_inv_accessor[0][2]}},
        {{cell_inv_accessor[1][0], cell_inv_accessor[1][1], cell_inv_accessor[1][2]}},
        {{cell_inv_accessor[2][0], cell_inv_accessor[2][1], cell_inv_accessor[2][2]}},
    }};

    // Collect the local atom id of all local & ghosts atoms, mapping ghosts
    // atoms which are periodic images of local atoms back to the local atoms.
    //
    // Metatensor expects pairs corresponding to periodic atoms to be between
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

    for (auto& cache: caches_) {
        auto cutoff2 = cache.cutoff * cache.cutoff;
        auto full_list = cache.options->full_list();

        // convert from LAMMPS neighbors list to metatensor format
        cache.known_samples.clear();
        cache.samples.clear();
        cache.distances_f32.clear();
        cache.distances_f64.clear();
        for (int ii=0; ii<(list_->inum + list_->gnum); ii++) {
            auto atom_i = list_->ilist[ii];
            auto original_atom_i = original_atom_id_[atom_i];

            auto neighbors = list_->firstneigh[ii];
            for (int jj=0; jj<list_->numneigh[ii]; jj++) {
                auto atom_j = neighbors[jj];
                auto original_atom_j = original_atom_id_[atom_j];

                if (!full_list && original_atom_i > original_atom_j) {
                    // Remove extra pairs if the model requested half-lists
                    continue;
                }

                auto distance = std::array<double, 3>{
                    x[atom_j][0] - x[atom_i][0],
                    x[atom_j][1] - x[atom_i][1],
                    x[atom_j][2] - x[atom_i][2],
                };

                auto distance2 = (
                    distance[0] * distance[0] +
                    distance[1] * distance[1] +
                    distance[2] * distance[2]
                );
                if (distance2 > cutoff2) {
                    // LAMMPS neighbors list contains some pairs after the
                    // cutoff, we filter them here
                    continue;
                }

                // Compute the cell shift for the pair.
                auto shift_i = std::array<double, 3>{
                    x[atom_i][0] - x[original_atom_i][0],
                    x[atom_i][1] - x[original_atom_i][1],
                    x[atom_i][2] - x[original_atom_i][2],
                };
                auto shift_j = std::array<double, 3>{
                    x[atom_j][0] - x[original_atom_j][0],
                    x[atom_j][1] - x[original_atom_j][1],
                    x[atom_j][2] - x[original_atom_j][2],
                };
                auto pair_shift = std::array<double, 3>{
                    shift_j[0] - shift_i[0],
                    shift_j[1] - shift_i[1],
                    shift_j[2] - shift_i[2],
                };

                auto shift = std::array<int32_t, 3>{0, 0, 0};
                if (pair_shift[0] != 0 || pair_shift[1] != 0 || pair_shift[2] != 0) {
                    shift = cell_shifts(cell_inv, pair_shift);

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
                            continue;
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
                            continue;
                        }
                    }
                }

                auto sample = std::array<int32_t, 5>{
                    original_atom_i,
                    original_atom_j,
                    shift[0],
                    shift[1],
                    shift[2],
                };

                // only add the pair if it is not already known. The same pair
                // can occur multiple time between two periodic ghosts shifted
                // around by the same amount, but we only want one of these pairs.
                if (cache.known_samples.insert(sample).second) {
                    cache.samples.push_back(sample);

                    if (dtype == torch::kFloat64) {
                        cache.distances_f64.push_back(distance);
                    } else if (dtype == torch::kFloat32) {
                        cache.distances_f32.push_back({
                            static_cast<float>(distance[0]),
                            static_cast<float>(distance[1]),
                            static_cast<float>(distance[2])
                        });
                    } else {
                        // should be unreachable
                        error->all(FLERR, "invalid dtype, this is a bug");
                    }
                }
            }
        }

        int64_t n_pairs = cache.samples.size();
        auto samples_values = torch::from_blob(
            reinterpret_cast<int32_t*>(cache.samples.data()),
            {n_pairs, 5},
            torch::TensorOptions().dtype(torch::kInt32).device(torch::kCPU)
        );

        auto samples = torch::make_intrusive<metatensor_torch::LabelsHolder>(
            std::vector<std::string>{"first_atom", "second_atom", "cell_shift_a", "cell_shift_b", "cell_shift_c"},
            samples_values
        );

        auto distances_vectors = torch::Tensor();
        if (dtype == torch::kFloat64) {
            distances_vectors = torch::from_blob(
                cache.distances_f64.data(),
                {n_pairs, 3, 1},
                torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU)
            );
        } else if (dtype == torch::kFloat32) {
            distances_vectors = torch::from_blob(
                cache.distances_f32.data(),
                {n_pairs, 3, 1},
                torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU)
            );
        } else {
            // should be unreachable
            error->all(FLERR, "invalid dtype, this is a bug");
        }

        auto neighbors = torch::make_intrusive<metatensor_torch::TensorBlockHolder>(
            distances_vectors.to(dtype).to(device),
            samples->to(device),
            std::vector<metatensor_torch::TorchLabels>{
                metatensor_torch::LabelsHolder::create({"xyz"}, {{0}, {1}, {2}})->to(device),
            },
            metatensor_torch::LabelsHolder::create({"distance"}, {{0}})->to(device)
        );

        metatensor_torch::register_autograd_neighbors(system, neighbors, options_.check_consistency);
        system->add_neighbors_list(cache.options, neighbors);
    }
}


metatensor_torch::System MetatensorSystemAdaptor::system_from_lmp(
    bool do_virial,
    torch::ScalarType dtype,
    torch::Device device
) {
    double** x = atom->x;
    auto total_n_atoms = atom->nlocal + atom->nghost;

    atomic_types_.resize_({total_n_atoms});
    for (int i=0; i<total_n_atoms; i++) {
        atomic_types_[i] = options_.types_mapping[atom->type[i]];
    }

    auto tensor_options = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);

    // atom->x contains "real" and then ghost atoms, in that order
    this->positions = torch::from_blob(
        *x, {total_n_atoms, 3},
        // requires_grad=true since we always need gradients w.r.t. positions
        tensor_options.requires_grad(true)
    );

    auto cell = torch::zeros({3, 3}, tensor_options);
    cell[0][0] = domain->xprd;

    cell[1][0] = domain->xy;
    cell[1][1] = domain->yprd;

    cell[2][0] = domain->xz;
    cell[2][1] = domain->yz;
    cell[2][2] = domain->zprd;

    auto system_positions = this->positions.to(dtype).to(device);
    cell = cell.to(dtype).to(device);

    if (do_virial) {
        auto model_strain = this->strain.to(dtype).to(device);

        // pretend to scale positions/cell by the strain so that
        // it enters the computational graph.
        system_positions = system_positions.matmul(model_strain);
        cell = cell.matmul(model_strain);
    }

    auto system = torch::make_intrusive<metatensor_torch::SystemHolder>(
        atomic_types_.to(device),
        system_positions,
        cell
    );

    this->setup_neighbors(system);
    return system;
}
