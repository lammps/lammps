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
------------------------------------------------------------------------- */
#include "metatomic_system.h"
#include "metatomic_timer.h"

#include "atom.h"
#include "domain.h"
#include "error.h"

#include "neigh_list.h"

#include <metatensor/torch.hpp>

using namespace LAMMPS_NS;

using vector_t = std::array<double, 3>;
using matrix_t = std::array<std::array<double, 3>, 3>;

static vector_t cross(vector_t a, vector_t b) {
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

static double dot(vector_t a, vector_t b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static vector_t normalize(vector_t a) {
    double norm = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    return {a[0] / norm, a[1] / norm, a[2] / norm};
}

static double determinant(matrix_t a) {
    return a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2])
         - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
         + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

matrix_t inverse(matrix_t a) {
    auto det = determinant(a);

    if (std::abs(det) < 1e-10) {
        throw std::runtime_error("this matrix is not invertible");
    }

    auto inverse = matrix_t();
    inverse[0][0] = (a[1][1] * a[2][2] - a[2][1] * a[1][2]) / det;
    inverse[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
    inverse[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
    inverse[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
    inverse[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
    inverse[1][2] = (a[1][0] * a[0][2] - a[0][0] * a[1][2]) / det;
    inverse[2][0] = (a[1][0] * a[2][1] - a[2][0] * a[1][1]) / det;
    inverse[2][1] = (a[2][0] * a[0][1] - a[0][0] * a[2][1]) / det;
    inverse[2][2] = (a[0][0] * a[1][1] - a[1][0] * a[0][1]) / det;
    return inverse;
}

/// Compute the inverse of the cell matrix of the system, accounting for
/// non-periodic directions by setting the corresponding rows to an unit vector
/// orthogonal to the periodic directions. This is used to compute the cell
/// shifts of neighbor pairs.
static std::array<std::array<double, 3>, 3> cell_inverse(Domain* domain) {
    auto periodic = std::array<bool, 3>{
        static_cast<bool>(domain->xperiodic),
        static_cast<bool>(domain->yperiodic),
        static_cast<bool>(domain->zperiodic),
    };

    auto cell = std::array<std::array<double, 3>, 3>{{0}};
    cell[0][0] = domain->xprd;
    cell[1][0] = domain->xy;
    cell[1][1] = domain->yprd;
    cell[2][0] = domain->xz;
    cell[2][1] = domain->yz;
    cell[2][2] = domain->zprd;

    // find number of periodic directions and their indices
    int n_periodic = 0;
    int periodic_idx_1 = -1;
    int periodic_idx_2 = -1;
    for (int i = 0; i < 3; ++i) {
        if (periodic[i]) {
            n_periodic += 1;
            if (periodic_idx_1 == -1) {
                periodic_idx_1 = i;
            } else if (periodic_idx_2 == -1) {
                periodic_idx_2 = i;
            }
        }
    }

    // adjust the box matrix to have a simple orthogonal dimension along
    // non-periodic directions
    if (n_periodic == 0) {
        return {
            std::array<double, 3>{1, 0, 0},
            std::array<double, 3>{0, 1, 0},
            std::array<double, 3>{0, 0, 1},
        };
    } else if (n_periodic == 1) {
        assert(periodic_idx_1 != -1);
        // Make the two non-periodic directions orthogonal to the periodic one
        auto a = cell[periodic_idx_1];
        auto b = std::array<double, 3>{0, 1, 0};
        if (std::abs(dot(normalize(a), b)) > 0.9) {
            b = std::array<double, 3>{0, 0, 1};
        }
        auto c = normalize(cross(a, b));
        b = normalize(cross(c, a));

        // Assign back to the cell picking the "non-periodic" indices without ifs
        cell[(periodic_idx_1 + 1) % 3] = b;
        cell[(periodic_idx_1 + 2) % 3] = c;
    } else if (n_periodic == 2) {
        assert(periodic_idx_1 != -1 && periodic_idx_2 != -1);
        // Make the one non-periodic direction orthogonal to the two periodic ones
        auto a = cell[periodic_idx_1];
        auto b = cell[periodic_idx_2];
        auto c = normalize(cross(a, b));

        // Assign back to the matrix picking the "non-periodic" index without ifs
        cell[(3 - periodic_idx_1 - periodic_idx_2)] = c;
    }

    return inverse(cell);
}

MetatomicSystemAdaptor::MetatomicSystemAdaptor(LAMMPS *lmp, MetatomicSystemOptions options):
    Pointers(lmp),
    options_(std::move(options)),
    nl_requests_(),
    atomic_types_(torch::zeros({0}, torch::TensorOptions().dtype(torch::kInt32).device(torch::kCPU)))
{
    auto tensor_options = torch::TensorOptions()
        .dtype(torch::kFloat64)
        .device(torch::kCPU)
        .requires_grad(options_.requires_grad);

    this->strain = torch::eye(3, tensor_options);
}

MetatomicSystemAdaptor::~MetatomicSystemAdaptor() {}

void MetatomicSystemAdaptor::add_nl_request(double cutoff, metatomic_torch::NeighborListOptions request) {
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

    nl_requests_.push_back({
        cutoff,
        request,
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

void MetatomicSystemAdaptor::guess_periodic_ghosts() {
    auto _ = MetatomicTimer("identifying periodic ghosts");
    auto total_n_atoms = atom->nlocal + atom->nghost;

    // Collect the local atom id of all local & ghosts atoms, mapping ghosts
    // atoms which are periodic images of local atoms back to the local atoms.
    //
    // metatomic expects pairs corresponding to periodic atoms to be between
    // the main atoms, but using the actual distance vector between the atom and
    // the ghost.
    original_atom_id_.clear();
    original_atom_id_.reserve(total_n_atoms);

    lmp_to_mta_.clear();
    lmp_to_mta_.reserve(total_n_atoms);

    mta_to_lmp.clear();
    mta_to_lmp.reserve(total_n_atoms);

    // identify all local atom by their LAMMPS atom tag.
    local_atoms_tags_.clear();
    for (int i=0; i<atom->nlocal; i++) {
        original_atom_id_.emplace_back(i);
        lmp_to_mta_.emplace_back(i);
        mta_to_lmp.emplace_back(i);
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
            lmp_to_mta_.emplace_back(-1);
        } else {
            // this can either be a periodic image of an atom owned by another
            // domain, or directly an atom from another domain. Since we can not
            // really distinguish between these, we take the first atom as the
            // "main" one and remap all atoms with the same tag to the first one
            auto it = ghost_atoms_tags_.find(tag);
            if (it != ghost_atoms_tags_.end()) {
                // we already found this atom elsewhere in the system
                original_atom_id_.emplace_back(it->second);
                lmp_to_mta_.emplace_back(-1);
            } else {
                // this is the first time we are seeing this atom
                original_atom_id_.emplace_back(i);
                ghost_atoms_tags_.emplace(tag, i);

                lmp_to_mta_.emplace_back(mta_to_lmp.size());
                mta_to_lmp.emplace_back(i);
            }
        }
    }
}


void MetatomicSystemAdaptor::setup_neighbors(metatomic_torch::System& system, NeighList *list) {
    auto _ = MetatomicTimer("converting neighbors list");
    auto dtype = system->positions().scalar_type();
    auto device = system->positions().device();

    double** x = atom->x;
    auto total_n_atoms = atom->nlocal + atom->nghost;
    auto cell_inv = cell_inverse(domain);

    for (auto& nl: nl_requests_) {
        {
            auto _ = MetatomicTimer("filtering LAMMPS neighbor list");

            auto cutoff2 = nl.cutoff * nl.cutoff;
            auto full_list = nl.options->full_list();

            // convert from LAMMPS neighbors list to metatomic format
            nl.samples.clear();
            nl.distances_f32.clear();
            nl.distances_f64.clear();
            for (int ii=0; ii<(list->inum + list->gnum); ii++) {
                auto atom_i = list->ilist[ii];
                assert(atom_i < total_n_atoms);
                auto original_atom_i = original_atom_id_[atom_i];
                auto i_is_original = (atom_i == original_atom_i);

                auto neighbors = list->firstneigh[ii];
                for (int jj=0; jj<list->numneigh[ii]; jj++) {
                    auto atom_j = neighbors[jj] & NEIGHMASK;
                    assert(atom_j < total_n_atoms);
                    auto original_atom_j = original_atom_id_[atom_j];
                    auto j_is_original = (atom_j == original_atom_j);

                    if (!full_list && original_atom_i > original_atom_j) {
                        // Remove extra pairs if the model requested half-lists
                        continue;
                    }

                    if (!i_is_original && !j_is_original) {
                        // both atoms are periodic ghosts, skip the pair
                        continue;
                    }

                    if (!i_is_original && j_is_original) {
                        // this pair will be accounted for when we will process
                        // atom_j as the central atom
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
                        lmp_to_mta_[original_atom_i],
                        lmp_to_mta_[original_atom_j],
                        shift[0],
                        shift[1],
                        shift[2],
                    };

                    nl.samples.push_back(sample);
                    if (dtype == torch::kFloat64) {
                        nl.distances_f64.push_back(distance);
                    } else if (dtype == torch::kFloat32) {
                        nl.distances_f32.push_back({
                            static_cast<float>(distance[0]),
                            static_cast<float>(distance[1]),
                            static_cast<float>(distance[2])
                        });
                    } else {
                        // should be unreachable
                        error->one(FLERR, "invalid dtype, this is a bug");
                    }
                }
            }
        }

        int64_t n_pairs = nl.samples.size();
        auto samples_values = torch::from_blob(
            reinterpret_cast<int32_t*>(nl.samples.data()),
            {n_pairs, 5},
            torch::TensorOptions().dtype(torch::kInt32).device(torch::kCPU)
        );

        torch::intrusive_ptr<metatensor_torch::LabelsHolder> samples;
        {
            auto _ = MetatomicTimer("creating samples Labels (" +  std::to_string(n_pairs) + " pairs)");
            if (options_.check_consistency) {
                // pairs should be unique, but I'm not 100% sure yet
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

        auto distances_vectors = torch::Tensor();
        if (dtype == torch::kFloat64) {
            distances_vectors = torch::from_blob(
                nl.distances_f64.data(),
                {n_pairs, 3, 1},
                torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU)
            );
        } else if (dtype == torch::kFloat32) {
            distances_vectors = torch::from_blob(
                nl.distances_f32.data(),
                {n_pairs, 3, 1},
                torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU)
            );
        } else {
            // should be unreachable
            error->one(FLERR, "invalid dtype, this is a bug");
        }

        {
            auto _ = MetatomicTimer("moving neighbor data to dtype/device");
            distances_vectors = distances_vectors.to(dtype).to(device);
            samples = samples->to(device);
        }

        torch::intrusive_ptr<metatensor_torch::TensorBlockHolder> neighbors;
        {
            auto _ = MetatomicTimer("creating neighbors TensorBlock");
            neighbors = torch::make_intrusive<metatensor_torch::TensorBlockHolder>(
                distances_vectors,
                samples,
                std::vector<metatensor_torch::Labels>{
                    metatensor_torch::LabelsHolder::create({"xyz"}, {{0}, {1}, {2}})->to(device),
                },
                metatensor_torch::LabelsHolder::create({"distance"}, {{0}})->to(device)
            );
        }

        metatomic_torch::register_autograd_neighbors(system, neighbors, options_.check_consistency);
        system->add_neighbor_list(nl.options, neighbors);
    }
}

metatomic_torch::System MetatomicSystemAdaptor::system_from_lmp(
    NeighList* list,
    bool do_virial,
    torch::ScalarType dtype,
    torch::Device device
) {
    auto _ = MetatomicTimer("creating System from LAMMPS data");

    double** x = atom->x;
    auto total_n_atoms = atom->nlocal + atom->nghost;

    atomic_types_.resize_({total_n_atoms});
    for (int i=0; i<total_n_atoms; i++) {
        atomic_types_[i] = options_.types_mapping[atom->type[i]];
    }

    auto tensor_options = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);

    // atom->x contains "real" and then ghost atoms, in that order
    this->positions = torch::from_blob(*x, {total_n_atoms, 3}, tensor_options);

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
        torch::TensorOptions().dtype(torch::kBool).device(torch::kCPU)
    );

    cell.index_put_({torch::logical_not(pbc)}, torch::tensor({0.0}, tensor_options));

    this->guess_periodic_ghosts();

    // Only keep the atoms which are not periodic images of other atoms
    auto mta_to_lmp_tensor = torch::from_blob(
        mta_to_lmp.data(),
        {static_cast<int64_t>(mta_to_lmp.size())},
        torch::TensorOptions().dtype(torch::kInt).device(torch::kCPU)
    );
    this->atomic_types_ = this->atomic_types_.index_select(0, mta_to_lmp_tensor);
    this->positions = this->positions.index_select(0, mta_to_lmp_tensor);

    this->positions.set_requires_grad(options_.requires_grad);

    auto system_positions = this->positions.to(dtype).to(device);
    cell = cell.to(dtype).to(device);

    if (do_virial) {
        auto model_strain = this->strain.to(dtype).to(device);

        // scale positions/cell by the strain so that it enters the
        // computational graph.
        system_positions = system_positions.matmul(model_strain);
        cell = cell.matmul(model_strain);
    }

    auto system = torch::make_intrusive<metatomic_torch::SystemHolder>(
        atomic_types_.to(device),
        system_positions,
        cell,
        pbc.to(device)
    );

    this->setup_neighbors(system, list);

    return system;
}
