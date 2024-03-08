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
   Contributing authors: TODO
------------------------------------------------------------------------- */
#include "pair_metatensor.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
#include "citeme.h"
#include "comm.h"

#include "neigh_list.h"
#include "neigh_request.h"

#include <torch/script.h>

#include <metatensor/torch.hpp>
#include <metatensor/torch/atomistic.hpp>

using namespace LAMMPS_NS;

static metatensor_torch::System system_from_lmp(Atom* atom);

/* ---------------------------------------------------------------------- */

PairMetatensor::PairMetatensor(LAMMPS *lmp) : Pair(lmp) {
    auto options = torch::TensorOptions().dtype(torch::kInt32);
    this->selected_atoms_values = torch::zeros({0, 2}, options);
    this->atomic_types = torch::zeros({0}, options);

    // Initialize evaluation_options
    this->evaluation_options = torch::make_intrusive<metatensor_torch::ModelEvaluationOptionsHolder>();
    auto energy_unit = std::string();
    if (strcmp(update->unit_style, "real") == 0) {
        this->evaluation_options->set_length_unit("angstrom");
        energy_unit = "kcal/mol";
    } else if (strcmp(update->unit_style, "metal") == 0) {
        this->evaluation_options->set_length_unit("angstrom");
        energy_unit = "eV";
    } else if (strcmp(update->unit_style, "si") == 0) {
        this->evaluation_options->set_length_unit("meter");
        energy_unit = "joule";
    } else if (strcmp(update->unit_style, "cgs") == 0) {
        this->evaluation_options->set_length_unit("centimeter");
        energy_unit = "erg";
    } else if (strcmp(update->unit_style, "electron") == 0) {
        this->evaluation_options->set_length_unit("Bohr");
        energy_unit = "Hartree";
    } else {
        error->all(FLERR, "unsupported units '{}' for pair metatensor ", update->unit_style);
    }

    auto output = torch::make_intrusive<metatensor_torch::ModelOutputHolder>();
    output->explicit_gradients = {};
    output->set_quantity("energy");
    output->set_unit(std::move(energy_unit));
    output->per_atom = false; // TODO, make this configurable

    this->evaluation_options->outputs.insert("energy", output);

    // we might not be running a pure pair potential,
    // so we can not compute virial as fdotr
    this->no_virial_fdotr_compute = 1;
}

PairMetatensor::~PairMetatensor() {
    delete this->lammps_to_metatensor_types;

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
    }
}


// called when finding `pair_style metatensor` in the input
void PairMetatensor::settings(int argc, char ** argv) {
    if (argc != 1) {
        error->all(FLERR, "expected 1 argument to pair_style metatensor, got {}", argc);
    }
    this->load_torch_model(argv[0]);

    if (!allocated) {
        allocate();
    }
}


void PairMetatensor::allocate() {
    allocated = 1;

    // setflags stores whether the coeff for a given pair of atom types are known
    setflag = memory->create(
        setflag,
        atom->ntypes + 1,
        atom->ntypes + 1,
        "pair:setflag"
    );

    for (int i = 1; i <= atom->ntypes; i++) {
        for (int j = i; j <= atom->ntypes; j++) {
            setflag[i][j] = 0;
        }
    }

    // cutsq stores the squared cutoff for each pair
    cutsq = memory->create(
        cutsq,
        atom->ntypes + 1,
        atom->ntypes + 1,
        "pair:cutsq"
    );

    // lammps_types_to_species stores the mapping from lammps atom types to
    // the metatensor model species
    lammps_to_metatensor_types = memory->create(
        lammps_to_metatensor_types,
        atom->ntypes + 1,
        "PairMetatensor:lammps_types_to_species"
    );

    for (int i = 1; i <= atom->ntypes; i++) {
        lammps_to_metatensor_types[i] = -1;
    }
}


// called on pair_coeff
void PairMetatensor::coeff(int argc, char ** argv) {
    if (argc < 3 || strcmp(argv[0], argv[1]) != 0) {
        error->all(FLERR, "pair_coeff can only be called with the same atom type twice (i.e. pair_coeff 1 1   6)");
    }

    int lammps_type = utils::inumeric(FLERR, argv[0], true, lmp);
    int species = utils::inumeric(FLERR, argv[2], true, lmp);

    lammps_to_metatensor_types[lammps_type] = species;

    // mark all pairs coeffs between the new atom type and already set ones as known
    for (int i = 1; i <= atom->ntypes; i++) {
        if (lammps_to_metatensor_types[i] >= 0) {
            setflag[i][lammps_type] = 1;
            setflag[lammps_type][i] = 1;
        }
    }
}


// called when the run starts
void PairMetatensor::init_style() {
    // Require newton pair on since we need to communicate forces accumulated on
    // ghost atoms to neighboring domains. These forces contributions come from
    // gradient of a local descriptor w.r.t. domain ghosts (periodic images
    // ghosts are handled separately).
    if (force->newton_pair != 1) {
        error->all(FLERR, "Pair style metatensor requires newton pair on");
    }

    // Translate from the metatensor neighbors lists requests to LAMMPS
    // neighbors lists requests.
    int n_requests = 0;
    auto requested_nl = this->torch_model->run_method("requested_neighbors_lists");
    for (const auto& ivalue: requested_nl.toList()) {
        auto options = ivalue.get().toCustomClass<metatensor_torch::NeighborsListOptionsHolder>();
        this->neigh_options.emplace_back(options);

        // We ask LAMMPS for a full neighbor lists because we need to know about
        // ALL pairs, even if options->full_list() is false. We will then filter
        // the pairs to only include each pair once where needed.
        auto request = neighbor->add_request(this, NeighConst::REQ_FULL);
        request->set_id(n_requests);
        request->set_cutoff(options->engine_cutoff(this->evaluation_options->length_unit()));

        n_requests += 1;
    }

    this->neigh_lists.resize(n_requests);
    this->neigh_data_cache.resize(n_requests);
}


void PairMetatensor::init_list(int id, NeighList *ptr) {
    this->neigh_lists[id] = ptr;
}


double PairMetatensor::init_one(int /*i*/, int /*j*/) {
    auto range = this->capabilities->engine_interaction_range(this->evaluation_options->length_unit());

    if (range < 0) {
        error->all(FLERR, "interaction_range is negative for this model");
    } else if (!std::isfinite(range)) {
        error->all(FLERR, "interaction_range is infinite for this model, this is not yet supported");
    } else if (range < 1) {
        return 1.0;
    } else {
        return range;
    }
}


void PairMetatensor::compute(int eflag, int vflag) {
    if (eflag || vflag) {
        ev_setup(eflag, vflag);
    } else {
        evflag = vflag_fdotr = eflag_global = eflag_atom = 0;
    }

    // transform from LAMMPS to metatensor System
    auto system = this->system_from_lmp();

    // only run the calculation for atoms actually in the current domain
    selected_atoms_values.resize_({atom->nlocal, 2});
    for (int i=0; i<atom->nlocal; i++) {
        selected_atoms_values[i][0] = 0;
        selected_atoms_values[i][1] = i;
    }
    auto selected_atoms = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        std::vector<std::string>{"system", "atom"}, selected_atoms_values
    );
    evaluation_options->set_selected_atoms(selected_atoms);

    torch::IValue result_ivalue;
    try {
        result_ivalue = this->torch_model->forward({
            std::vector<metatensor_torch::System>{system},
            evaluation_options,
            this->check_consistency
        });
    } catch (const std::exception& e) {
        error->all(FLERR, "error evaluating the torch model: {}", e.what());
    }

    auto result = result_ivalue.toGenericDict();
    auto energy = result.at("energy").toCustomClass<metatensor_torch::TensorMapHolder>();
    auto energy_tensor = metatensor_torch::TensorMapHolder::block_by_id(energy, 0)->values();

    // store the global energy returned by the model
    if (eflag_global) {
        eng_vdwl += energy_tensor.item<double>();
    }

    if (eflag_atom) {
        error->all(FLERR, "per atom energy is not implemented yet");
    }

    if (system->positions().size(0) != 0) {
        // if we have 0 atoms, we lost all of them! We'll let LAMMPS check for
        // this and report it later.

        energy_tensor.backward(-torch::ones_like(energy_tensor));

        auto forces_tensor = system->positions().grad().reshape({-1, 3});
        auto forces = forces_tensor.accessor<double, 2>();
        for (int i=0; i<atom->nlocal + atom->nghost; i++) {
            atom->f[i][0] += forces[i][0];
            atom->f[i][1] += forces[i][1];
            atom->f[i][2] += forces[i][2];
        }
    }

    assert(!vflag_fdotr);

    if (vflag_global) {
        error->warning(FLERR, "virial is not implemented yet");
    }

    if (vflag_atom) {
        error->all(FLERR, "per atom virial is not implemented");
    }
}


void PairMetatensor::load_torch_model(const char* path) {
    if (this->torch_model != nullptr) {
        error->all(FLERR, "torch model is already loaded");
    }

    // allow the user to request other devices
    auto device = torch::kCPU;

    try {
        metatensor_torch::check_atomistic_model(path);
        this->torch_model = std::make_unique<torch::jit::script::Module>(torch::jit::load(path, device));
    } catch (const c10::Error& e) {
        error->all(FLERR, "failed to load metatensor model at '{}': {}", path, e.what());
    }

    auto capabilities_ivalue = this->torch_model->run_method("capabilities");

    this->capabilities = capabilities_ivalue.toCustomClass<metatensor_torch::ModelCapabilitiesHolder>();

    bool found_valid_device = false;
    for (const auto& device: this->capabilities->supported_devices) {
        if (device == "cpu") {
            found_valid_device = true;
            break;
        }
    }

    if (!found_valid_device) {
        error->all(FLERR,
            "failed to find a valid device for the model at '{}': "
            "only 'cpu' is currently supported by LAMMPS, and the model supports {}",
            path, at::str(this->capabilities->supported_devices)
        );
    }

    // move all data to the right device
    this->selected_atoms_values = selected_atoms_values.to(device);
    this->atomic_types = atomic_types.to(device);

    if (lmp->comm->me == 0) {
        auto metadata_ivalue = this->torch_model->run_method("metadata");
        auto metadata = metadata_ivalue.toCustomClass<metatensor_torch::ModelMetadataHolder>();
        auto to_print = metadata->print();

        if (screen) {
            fprintf(screen, "\n%s\n", to_print.c_str());
        }
        if (logfile) {
            fprintf(logfile,"\n%s\n", to_print.c_str());
        }

        // add the model references to LAMMPS citation handling mechanism
        for (const auto& it: metadata->references) {
            for (const auto& ref: it.value()) {
                lmp->citeme->add(ref + "\n");
            }
        }
    }
}


metatensor_torch::System PairMetatensor::system_from_lmp() {
    double** x = atom->x;

    auto n_atoms = atom->nlocal + atom->nghost;

    this->atomic_types.resize_({n_atoms});
    for (int i=0; i<n_atoms; i++) {
        this->atomic_types[i] = this->lammps_to_metatensor_types[atom->type[i]];
    }

    // atom->x contains "real" and then ghost atoms, in that order
    auto positions = torch::from_blob(
        *x, {n_atoms, 3},
        // requires_grad=true since we always need gradients w.r.t. positions
        torch::TensorOptions().dtype(torch::kDouble).requires_grad(true)
    );

    auto cell = torch::zeros(
        {3, 3},
        torch::TensorOptions().dtype(torch::kDouble)
    );

    cell[0][0] = domain->xprd;

    cell[1][0] = domain->xy;
    cell[1][1] = domain->yprd;

    cell[2][0] = domain->xz;
    cell[2][1] = domain->yz;
    cell[2][2] = domain->zprd;

    auto cell_inv = cell.inverse().t();

    auto system = torch::make_intrusive<metatensor_torch::SystemHolder>(atomic_types, positions, cell);

    // Collect the local atom id of all local & ghosts atoms, mapping ghosts
    // atoms which are periodic images of local atoms back to the local atoms.
    //
    // Metatensor expects pairs corresponding to periodic atoms to be between
    // the main atoms, but using the actual distance vector between the atom and
    // the ghost.
    auto original_atom_id = std::vector<int>();
    original_atom_id.reserve(n_atoms);

    // identify all local atom by their LAMMPS atom tag.
    auto local_atoms_tags = std::unordered_map<tagint, int>();
    for (int i=0; i<atom->nlocal; i++) {
        original_atom_id.emplace_back(i);
        local_atoms_tags.emplace(atom->tag[i], i);
    }

    // now loop over ghosts & map them back to the main cell if needed
    for (int i=atom->nlocal; i<n_atoms; i++) {
        auto it = local_atoms_tags.find(atom->tag[i]);
        if (it != local_atoms_tags.end()) {
            original_atom_id.emplace_back(it->second);
        } else {
            original_atom_id.emplace_back(i);
        }
    }

    for (size_t i=0; i<this->neigh_options.size(); i++) {
        auto options = this->neigh_options[i];
        auto* list = this->neigh_lists[i];
        auto& cache = this->neigh_data_cache[i];

        if (cache.cutoff < 0) {
            cache.cutoff = options->engine_cutoff(this->evaluation_options->length_unit());
            if (cache.cutoff < 0 || !std::isfinite(cache.cutoff)) {
                error->all(FLERR,
                    "model requested an invalid cutoff for neighbors list: {} "
                    "(cutoff in model units is {})",
                    cache.cutoff, options->cutoff()
                );
            }
        }

        auto cutoff2 = cache.cutoff * cache.cutoff;

        // convert from LAMMPS neighbors list to metatensor format
        cache.samples.clear();
        cache.distances.clear();
        for (int ii=0; ii<list->inum; ii++) {
            auto atom_i = list->ilist[ii];
            auto original_atom_i = original_atom_id[atom_i];

            auto neighbors = list->firstneigh[ii];
            for (int jj=0; jj<list->numneigh[ii]; jj++) {
                auto atom_j = neighbors[jj];
                auto original_atom_j = original_atom_id[atom_j];

                if (!options->full_list() && original_atom_i > original_atom_j) {
                    // Remove extra pairs if the model requested half-lists
                    continue;
                }

                auto xx = x[atom_j][0] - x[atom_i][0];
                auto yy = x[atom_j][1] - x[atom_i][1];
                auto zz = x[atom_j][2] - x[atom_i][2];

                auto distance2 = xx * xx + yy * yy + zz * zz;
                if (distance2 > cutoff2) {
                    // LAMMPS neighbors list contains some pairs after the
                    // cutoff, we filter them here
                    continue;
                }

                // Compute the cell shift for the pair. This is non-zero only if
                // the pair is between an atom and a periodic image. With
                // LAMMPS, the second atoms in the pair is the one outside the
                // main cell, so we only need to check atom_j versus
                // original_atom_j
                double periodic_shift[3];
                periodic_shift[0] = x[atom_j][0] - x[original_atom_j][0];
                periodic_shift[1] = x[atom_j][1] - x[original_atom_j][1];
                periodic_shift[2] = x[atom_j][2] - x[original_atom_j][2];

                int32_t shift_a = 0;
                int32_t shift_b = 0;
                int32_t shift_c = 0;
                if (periodic_shift[0] != 0 || periodic_shift[1] != 0 || periodic_shift[2] != 0) {
                    auto periodic_shift_tensor = torch::tensor(
                        {periodic_shift[0], periodic_shift[1], periodic_shift[2]},
                        torch::TensorOptions().dtype(torch::kFloat64)
                    );
                    auto cell_shift_tensor = cell_inv.matmul(periodic_shift_tensor);
                    auto cell_shift = cell_shift_tensor.accessor<double, 1>();

                    shift_a = static_cast<int32_t>(std::round(cell_shift[0]));
                    shift_b = static_cast<int32_t>(std::round(cell_shift[1]));
                    shift_c = static_cast<int32_t>(std::round(cell_shift[2]));

                    if (!options->full_list() && original_atom_i == original_atom_j) {
                        // If a half neighbors list has been requested, do
                        // not include the same pair between an atom and
                        // it's periodic image twice with opposite cell
                        // shifts (e.g. [1, -1, 1] and [-1, 1, -1]).
                        //
                        // Instead we pick pairs in the positive plan of
                        // shifts.
                        if (shift_a + shift_b + shift_c < 0) {
                            // drop shifts on the negative half-space
                            continue;
                        }

                        if ((shift_a + shift_b + shift_c == 0)
                            && (shift_c < 0 || (shift_c == 0 && shift_b < 0))) {
                            // drop shifts in the negative half plane or the
                            // negative shift[1] axis. See below for a
                            // graphical representation: we are keeping the
                            // shifts indicated with `O` and dropping the
                            // ones indicated with `X`
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

                cache.distances.emplace_back(xx);
                cache.distances.emplace_back(yy);
                cache.distances.emplace_back(zz);

                cache.samples.emplace_back(original_atom_i);
                cache.samples.emplace_back(original_atom_j);
                cache.samples.emplace_back(shift_a);
                cache.samples.emplace_back(shift_b);
                cache.samples.emplace_back(shift_c);
            }
        }

        assert(cache.distances.size() % 3 == 0);
        int64_t n_pairs = cache.distances.size() / 3;
        assert(cache.samples.size() / 5 == n_pairs);
        auto samples_values = torch::from_blob(
            cache.samples.data(),
            {n_pairs, 5},
            torch::TensorOptions().dtype(torch::kInt32).device(torch::kCPU)
        );

        auto samples = torch::make_intrusive<metatensor_torch::LabelsHolder>(
            std::vector<std::string>{"first_atom", "second_atom", "cell_shift_a", "cell_shift_b", "cell_shift_c"},
            samples_values
        );

        auto distances_vectors = torch::from_blob(
            cache.distances.data(),
            {n_pairs, 3, 1},
            torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU)
        );

        auto neighbors = torch::make_intrusive<metatensor_torch::TensorBlockHolder>(
            distances_vectors,
            samples,
            std::vector<metatensor_torch::TorchLabels>{
                metatensor_torch::LabelsHolder::create({"xyz"}, {{0}, {1}, {2}}),
            },
            metatensor_torch::LabelsHolder::create({"distance"}, {{0}})
        );

        metatensor_torch::register_autograd_neighbors(system, neighbors, this->check_consistency);
        system->add_neighbors_list(options, neighbors);
    }

    return system;
}
