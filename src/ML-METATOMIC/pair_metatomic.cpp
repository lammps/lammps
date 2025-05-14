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
#include "pair_metatomic.h"
#include "metatomic_types.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
#include "comm.h"

#include "neigh_list.h"
#include "neigh_request.h"

#include <torch/version.h>
#include <torch/script.h>
#include <torch/cuda.h>

#if TORCH_VERSION_MAJOR >= 2
    #include <torch/mps.h>
#endif

#include <memory>

#include <metatensor/torch.hpp>
#include <metatomic/torch.hpp>

#include "metatomic_system.h"
#include "metatomic_timer.h"

using namespace LAMMPS_NS;

PairMetatomic::PairMetatomic(LAMMPS *lmp):
    Pair(lmp),
    type_mapping(nullptr),
    system_adaptor(nullptr)
{
    std::string energy_unit;
    std::string length_unit;
    if (strcmp(update->unit_style, "real") == 0) {
        length_unit = "angstrom";
        energy_unit = "kcal/mol";
    } else if (strcmp(update->unit_style, "metal") == 0) {
        length_unit = "angstrom";
        energy_unit = "eV";
    } else if (strcmp(update->unit_style, "si") == 0) {
        length_unit = "meter";
        energy_unit = "joule";
    } else if (strcmp(update->unit_style, "electron") == 0) {
        length_unit = "Bohr";
        energy_unit = "Hartree";
    } else {
        error->all(FLERR, "unsupported units '{}' for pair metatomic ", update->unit_style);
    }

    // we might not be running a pure pair potential,
    // so we can not compute virial as fdotr
    this->no_virial_fdotr_compute = 1;

    this->mta_data = new PairMetatomicData(std::move(length_unit), std::move(energy_unit));

    // settings for metatomic pair style
    this->single_enable = 0;
    this->restartinfo = 0;
    this->one_coeff = 1;
    this->manybody_flag = 1;
}

PairMetatomic::~PairMetatomic() {
    delete this->mta_data;

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(type_mapping);
    }
}

// called when finding `pair_style metatensor` in the input
void PairMetatomic::settings(int argc, char ** argv) {
    if (argc == 0) {
        error->all(FLERR, "expected at least 1 argument to pair_style metatomic, got {}", argc);
    }

    const char* model_path = argv[0];
    const char* extensions_directory = nullptr;
    const char* requested_device = nullptr;
    for (int i=1; i<argc; i++) {
        if (strcmp(argv[i], "check_consistency") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected <on/off> after 'check_consistency' in pair_style metatomic, got nothing");
            } else if (strcmp(argv[i + 1], "on") == 0) {
                mta_data->check_consistency = true;
            } else if (strcmp(argv[i + 1], "off") == 0) {
                mta_data->check_consistency = false;
            } else {
                error->all(FLERR, "expected <on/off> after 'check_consistency' in pair_style metatomic, got '{}'", argv[i + 1]);
            }

            i += 1;
        } else if (strcmp(argv[i], "remap_pairs") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected <on/off> after 'remap_pairs' in pair_style metatomic, got nothing");
            } else if (strcmp(argv[i + 1], "on") == 0) {
                mta_data->remap_pairs = true;
            } else if (strcmp(argv[i + 1], "off") == 0) {
                mta_data->remap_pairs = false;
            } else {
                error->all(FLERR, "expected <on/off> after 'remap_pairs' in pair_style metatomic, got '{}'", argv[i + 1]);
            }

            i += 1;
        } else if (strcmp(argv[i], "extensions") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected <path> after 'extensions' in pair_style metatomic, got nothing");
            }
            extensions_directory = argv[i + 1];
            i += 1;
        } else if (strcmp(argv[i], "device") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected string after 'device' in pair_style metatomic, got nothing");
            }
            requested_device = argv[i + 1];
            i += 1;
        } else {
            error->all(FLERR, "unexpected argument to pair_style metatomic: '{}'", argv[i]);
        }
    }

    // load the model and get it's capabilities (including supported devices)
    mta_data->load_model(this->lmp, model_path, extensions_directory);

    // Select the device to use based on the model's preference, the user choice
    // and what's available.
    this->pick_device(&mta_data->device, requested_device);

    // move all data to the correct device
    mta_data->model->to(mta_data->device);
    mta_data->selected_atoms_values = mta_data->selected_atoms_values.to(mta_data->device);

    auto message = "Running simulation on " + mta_data->device.str() + " device with " + mta_data->capabilities->dtype() + " data";
    if (screen) {
        fprintf(screen, "%s\n", message.c_str());
    }
    if (logfile) {
        fprintf(logfile,"%s\n", message.c_str());
    }

    if (!allocated) {
        allocate();
    }
}

std::vector<torch::DeviceType> PairMetatomic::available_devices() {
    auto devices = std::vector<torch::DeviceType>();
    for (const auto& supported: this->mta_data->capabilities->supported_devices) {
        if (supported == "cpu") {
            devices.push_back(torch::kCPU);
        } else if (supported == "cuda" && torch::cuda::is_available()) {
            devices.push_back(torch::kCUDA);
        } else if (supported == "mps") {
            #if TORCH_VERSION_MAJOR >= 2
            if (torch::mps::is_available()) {
                devices.push_back(torch::kMPS);
            }
            #endif
        } else {
            error->warning(FLERR,
                "the model declared support for unknown device '{}', it will be ignored", supported
            );
        }
    }

    if (devices.empty()) {
        error->all(FLERR,
            "failed to find a valid device for this model: "
            "the model supports {}, none of these where available",
            torch::str(this->mta_data->capabilities->supported_devices)
        );
    }

    return devices;
}

void PairMetatomic::pick_device(torch::Device* device, const char* requested) {
    auto available_devices = this->available_devices();

    auto picked_device_type = torch::kCPU;
    if (requested == nullptr) {
        // no user request, pick the device the model prefers
        picked_device_type = available_devices[0];
    } else {
        bool found_requested_device = false;
        for (const auto& device_type: available_devices) {
            if (device_type == torch::kCPU && strcmp(requested, "cpu") == 0) {
                picked_device_type = device_type;
                found_requested_device = true;
                break;
            } else if (device_type == torch::kCUDA && strcmp(requested, "cuda") == 0) {
                picked_device_type = device_type;
                found_requested_device = true;
                break;
            } else if (device_type == torch::kMPS && strcmp(requested, "mps") == 0) {
                picked_device_type = device_type;
                found_requested_device = true;
                break;
            }
        }

        if (!found_requested_device) {
            error->all(FLERR,
                "failed to find requested device ({}): it is either "
                "not supported by this model or not available on this machine",
                requested
            );
        }
    }

    if (picked_device_type == torch::kCUDA) {
        // distribute GPUs between multiple MPI processes on the same node

        // (1) get a MPI communicator for all processes on the current node
        MPI_Comm local;
        MPI_Comm_split_type(world, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &local);
        // (2) get the rank of this MPI process on the current node
        int local_rank;
        MPI_Comm_rank(local, &local_rank);

        int size;
        MPI_Comm_size(local, &size);
        if (size < torch::cuda::device_count()) {
            if (comm->me == 0) {
                error->warning(FLERR,
                    "found {} CUDA-capable GPUs, but only {} MPI processes on the current node; the remaining GPUs will not be used",
                    torch::cuda::device_count(), size
                );
            }
        }

        // (3) split GPUs between node-local processes using round-robin allocation
        int gpu_to_use = local_rank % torch::cuda::device_count();
        *device = torch::Device(picked_device_type, gpu_to_use);
    } else {
        *device = torch::Device(picked_device_type);
    }
}


void PairMetatomic::allocate() {
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
    // the metatomic model species
    type_mapping = memory->create(
        type_mapping,
        atom->ntypes + 1,
        "PairMetatomic:type_mapping"
    );

    for (int i = 1; i <= atom->ntypes; i++) {
        type_mapping[i] = -1;
    }
}

double PairMetatomic::init_one(int, int) {
    return mta_data->max_cutoff;
}


// called on pair_coeff
void PairMetatomic::coeff(int argc, char ** argv) {
    if (argc < 3 || strcmp(argv[0], "*") != 0 || strcmp(argv[1], "*") != 0) {
        error->all(FLERR, "invalid pair_coeff, expected `pair_coeff * * <list of types>`");
    }

    if (atom->ntypes != argc - 2) {
        error->all(FLERR,
            "invalid pair_coeff, expected `pair_coeff * * <list of types>` with {} types",
            atom->ntypes
        );
    }

    for (int lammps_type=1; lammps_type<argc - 1; lammps_type++) {
        int type = utils::inumeric(FLERR, argv[lammps_type + 1], true, lmp);
        type_mapping[lammps_type] = type;
    }

    // mark all pairs coeffs as known
    for (int i = 1; i <= atom->ntypes; i++) {
        for (int j = 1; j <= atom->ntypes; j++) {
            setflag[i][j] = 1;
            setflag[j][i] = 1;
        }
    }
}


// called when the run starts
void PairMetatomic::init_style() {
    // Require newton pair on since we need to communicate forces accumulated on
    // ghost atoms to neighboring domains. These forces contributions come from
    // gradient of a local descriptor w.r.t. domain ghosts (periodic images
    // ghosts are handled separately).
    if (force->newton_pair != 1) {
        error->all(FLERR, "Pair style metatomic requires newton pair on");
    }

    // get the model's interaction range
    auto range = mta_data->capabilities->engine_interaction_range(mta_data->evaluation_options->length_unit());
    if (range < 0) {
        error->all(FLERR, "interaction_range is negative for this model");
    } else if (!std::isfinite(range)) {
        if (comm->nprocs > 1) {
            error->all(FLERR,
                "interaction_range is infinite for this model, "
                "using multiple MPI domains is not supported"
            );
        }

        // determine the maximal cutoff in the NL
        auto requested_nl = mta_data->model->run_method("requested_neighbor_lists");
        for (const auto& ivalue: requested_nl.toList()) {
            auto options = ivalue.get().toCustomClass<metatomic_torch::NeighborListOptionsHolder>();
            auto cutoff = options->engine_cutoff(mta_data->evaluation_options->length_unit());

            mta_data->max_cutoff = std::max(mta_data->max_cutoff, cutoff);
        }
    } else {
        mta_data->max_cutoff = range;
    }

    if (!std::isfinite(mta_data->max_cutoff)) {
        error->all(FLERR,
            "the largest cutoff of this model is infinite, "
            "we can't compute the corresponding neighbor list"
        );
    }

    // create system adaptor
    auto options = MetatomicSystemOptions{
        this->type_mapping,
        mta_data->max_cutoff,
        mta_data->check_consistency,
    };
    this->system_adaptor = std::make_unique<MetatomicSystemAdaptor>(lmp, options);

    // We ask LAMMPS for a full neighbor lists because we need to know about
    // ALL pairs, even if options->full_list() is false. We will then filter
    // the pairs to only include each pair once where needed.
    auto request = neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);
    request->set_cutoff(mta_data->max_cutoff);

    // Translate from the metatomic neighbor lists requests to LAMMPS neighbor
    // lists requests.
    auto requested_nl = mta_data->model->run_method("requested_neighbor_lists");
    for (const auto& ivalue: requested_nl.toList()) {
        auto options = ivalue.get().toCustomClass<metatomic_torch::NeighborListOptionsHolder>();
        auto cutoff = options->engine_cutoff(mta_data->evaluation_options->length_unit());
        assert(cutoff <= mts_data->max_cutoff);

        this->system_adaptor->add_nl_request(cutoff, options);
    }
}


void PairMetatomic::init_list(int id, NeighList *ptr) {
    this->mta_list = ptr;
}


void PairMetatomic::compute(int eflag, int vflag) {
    if (std::getenv("LAMMPS_METATOMIC_PROFILE") != nullptr) {
        MetatomicTimer::enable(true);
    } else {
        MetatomicTimer::enable(false);
    }

    auto _ = MetatomicTimer("PairMetatomic::compute");

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

    // transform from LAMMPS to metatomic System
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
    auto options = mta_data->selected_atoms_values.options();
    mta_data->selected_atoms_values.index_put_(
        {torch::indexing::Slice(), 1},
        torch::arange(atom->nlocal, options)
    );

    auto selected_atoms = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        std::vector<std::string>{"system", "atom"}, mta_data->selected_atoms_values
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

    // compute forces/virial on device with backward propagation
    {
        // reset gradients to zero before calling backward
        this->system_adaptor->positions.mutable_grad() = torch::Tensor();
        this->system_adaptor->strain.mutable_grad() = torch::Tensor();

        auto _ = MetatomicTimer("running Model::backward");
        energy_tensor.backward(-torch::ones_like(energy_tensor));
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
            assert(samples_values.sizes() == mts_data->selected_atoms_values.sizes());

            auto energies = energy_detached.accessor<double, 2>();
            for (int64_t i=0; i<energy_samples->count(); i++) {
                assert(samples[i][0] == 0);
                // handle potentially out of order samples in
                // the per-atom energy tensor
                auto atom_i = samples[i][1];
                assert(atom_i < n_atoms);
                eatom[atom_i] += energies[i][0];
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
            eng_vdwl += global_energy.item<double>();
        }

        // store forces/virial
        auto forces_tensor = this->system_adaptor->positions.grad();
        assert(forces_tensor.is_cpu() && forces_tensor.scalar_type() == torch::kFloat64);

        auto forces = forces_tensor.accessor<double, 2>();
        for (int i=0; i<atom->nlocal + atom->nghost; i++) {
            atom->f[i][0] += forces[i][0];
            atom->f[i][1] += forces[i][1];
            atom->f[i][2] += forces[i][2];
        }

        assert(!vflag_fdotr);

        if (vflag_global) {
            auto virial_tensor = this->system_adaptor->strain.grad();
            assert(virial_tensor.is_cpu() && virial_tensor.scalar_type() == torch::kFloat64);
            auto predicted_virial = virial_tensor.accessor<double, 2>();

            virial[0] += predicted_virial[0][0];
            virial[1] += predicted_virial[1][1];
            virial[2] += predicted_virial[2][2];

            virial[3] += 0.5 * (predicted_virial[1][0] + predicted_virial[0][1]);
            virial[4] += 0.5 * (predicted_virial[2][0] + predicted_virial[0][2]);
            virial[5] += 0.5 * (predicted_virial[2][1] + predicted_virial[1][2]);
        }

        if (vflag_atom) {
            error->all(FLERR, "per atom virial is not implemented");
        }
    }
}
