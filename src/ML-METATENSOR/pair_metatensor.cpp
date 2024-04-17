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
#include "pair_metatensor.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
#include "citeme.h"
#include "comm.h"

#include "neigh_list.h"

#include <torch/version.h>
#include <torch/script.h>
#include <torch/cuda.h>

#if TORCH_VERSION_MAJOR >= 2
    #include <torch/mps.h>
#endif

#include <metatensor/torch.hpp>
#include <metatensor/torch/atomistic.hpp>

#include "metatensor_system.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMetatensor::PairMetatensor(LAMMPS *lmp):
    Pair(lmp),
    system_adaptor(nullptr),
    device(torch::kCPU)
{
    auto options = torch::TensorOptions().dtype(torch::kInt32);
    this->selected_atoms_values = torch::zeros({0, 2}, options);

    // default to true for now, this will be changed to false later
    this->check_consistency = true;

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
    output->per_atom = false;

    this->evaluation_options->outputs.insert("energy", output);

    // we might not be running a pure pair potential,
    // so we can not compute virial as fdotr
    this->no_virial_fdotr_compute = 1;
}

PairMetatensor::~PairMetatensor() {
    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(type_mapping);
    }
}

// called when finding `pair_style metatensor` in the input
void PairMetatensor::settings(int argc, char ** argv) {
    if (argc == 0) {
        error->all(FLERR, "expected at least 1 argument to pair_style metatensor, got {}", argc);
    }

    const char* model_path = argv[0];
    const char* extensions_directory = nullptr;
    const char* requested_device = nullptr;
    for (int i=1; i<argc; i++) {
        if (strcmp(argv[i], "check_consistency") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected <on/off> after 'check_consistency' in pair_style metatensor, got nothing");
            } else if (strcmp(argv[i + 1], "on") == 0) {
                this->check_consistency = true;
            } else if (strcmp(argv[i + 1], "off") == 0) {
                this->check_consistency = false;
            } else {
                error->all(FLERR, "expected <on/off> after 'check_consistency' in pair_style metatensor, got '{}'", argv[i + 1]);
            }

            i += 1;
        } else if (strcmp(argv[i], "extensions") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected <path> after 'extensions' in pair_style metatensor, got nothing");
            }
            extensions_directory = argv[i + 1];
            i += 1;
        } else if (strcmp(argv[i], "device") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected string after 'device' in pair_style metatensor, got nothing");
            }
            requested_device = argv[i + 1];
            i += 1;
        } else {
            error->all(FLERR, "unexpected argument to pair_style metatensor: '{}'", argv[i]);
        }
    }

    this->load_torch_model(model_path, extensions_directory);

    // Select the device to use based on the model's preference, the user choice
    // and what's available.
    auto available_devices = std::vector<torch::Device>();
    for (const auto& device: this->capabilities->supported_devices) {
        if (device == "cpu") {
            available_devices.push_back(torch::kCPU);
        } else if (device == "cuda") {
            if (torch::cuda::is_available()) {
                available_devices.push_back(torch::Device("cuda"));
            }
        } else if (device == "mps") {
            #if TORCH_VERSION_MAJOR >= 2
            if (torch::mps::is_available()) {
                available_devices.push_back(torch::Device("mps"));
            }
            #endif
        } else {
            error->warning(FLERR,
                "the model declared support for unknown device '{}', it will be ignored", device
            );
        }
    }

    if (available_devices.empty()) {
        error->all(FLERR,
            "failed to find a valid device for the model at '{}': "
            "the model supports {}, none of these where available",
            model_path, torch::str(this->capabilities->supported_devices)
        );
    }

    if (requested_device == nullptr) {
        // no user request, pick the device the model prefers
        this->device = available_devices[0];
    } else {
        bool found_requested_device = false;
        for (const auto& device: available_devices) {
            if (device.is_cpu() && strcmp(requested_device, "cpu") == 0) {
                this->device = device;
                found_requested_device = true;
                break;
            } else if (device.is_cuda() && strcmp(requested_device, "cuda") == 0) {
                this->device = device;
                found_requested_device = true;
                break;
            } else if (device.is_mps() && strcmp(requested_device, "mps") == 0) {
                this->device = device;
                found_requested_device = true;
                break;
            }
        }

        if (!found_requested_device) {
            error->all(FLERR,
                "failed to find requested device ({}): it is either "
                "not supported by this model or not available on this machine",
                requested_device
            );
        }
    }

    this->torch_model->to(device);

    auto message = "Running simulation on " + this->device.str() + " device with " + this->capabilities->dtype() + " data";
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
    type_mapping = memory->create(
        type_mapping,
        atom->ntypes + 1,
        "PairMetatensor:type_mapping"
    );

    for (int i = 1; i <= atom->ntypes; i++) {
        type_mapping[i] = -1;
    }
}

double PairMetatensor::init_one(int, int) {
    return this->interaction_range;
}


// called on pair_coeff
void PairMetatensor::coeff(int argc, char ** argv) {
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
void PairMetatensor::init_style() {
    // Require newton pair on since we need to communicate forces accumulated on
    // ghost atoms to neighboring domains. These forces contributions come from
    // gradient of a local descriptor w.r.t. domain ghosts (periodic images
    // ghosts are handled separately).
    if (force->newton_pair != 1) {
        error->all(FLERR, "Pair style metatensor requires newton pair on");
    }

    // get the model's interaction range
    auto range = this->capabilities->engine_interaction_range(this->evaluation_options->length_unit());
    if (range < 0) {
        error->all(FLERR, "interaction_range is negative for this model");
    } else if (!std::isfinite(range)) {
        error->all(FLERR, "interaction_range is infinite for this model, this is not yet supported");
    } else {
        this->interaction_range = range;
    }

    // create system adaptor
    auto options = MetatensorSystemOptions{
        type_mapping,
        this->interaction_range,
        check_consistency,
    };
    this->system_adaptor = std::make_unique<MetatensorSystemAdaptor>(lmp, this, options);

    // Translate from the metatensor neighbors lists requests to LAMMPS
    // neighbors lists requests.
    auto requested_nl = this->torch_model->run_method("requested_neighbors_lists");
    for (const auto& ivalue: requested_nl.toList()) {
        auto options = ivalue.get().toCustomClass<metatensor_torch::NeighborsListOptionsHolder>();
        auto cutoff = options->engine_cutoff(this->evaluation_options->length_unit());

        this->system_adaptor->add_nl_request(cutoff, options);
    }
}


void PairMetatensor::init_list(int id, NeighList *ptr) {
    system_adaptor->init_list(id, ptr);
}


void PairMetatensor::compute(int eflag, int vflag) {
    if (eflag || vflag) {
        ev_setup(eflag, vflag);
    } else {
        evflag = vflag_fdotr = eflag_global = eflag_atom = 0;
    }

    if (eflag_atom) {
        this->evaluation_options->outputs.at("energy")->per_atom = true;
    } else {
        this->evaluation_options->outputs.at("energy")->per_atom = false;
    }

    auto dtype = torch::kFloat64;
    if (this->capabilities->dtype() == "float64") {
        dtype = torch::kFloat64;
    } else if (this->capabilities->dtype() == "float32") {
        dtype = torch::kFloat32;
    } else {
        error->all(FLERR, "the model requested an unsupported dtype '{}'", this->capabilities->dtype());
    }

    // transform from LAMMPS to metatensor System
    auto system = system_adaptor->system_from_lmp(static_cast<bool>(vflag_global), dtype, this->device);

    // only run the calculation for atoms actually in the current domain
    selected_atoms_values.resize_({atom->nlocal, 2});
    for (int i=0; i<atom->nlocal; i++) {
        selected_atoms_values[i][0] = 0;
        selected_atoms_values[i][1] = i;
    }
    auto selected_atoms = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        std::vector<std::string>{"system", "atom"}, selected_atoms_values
    );
    evaluation_options->set_selected_atoms(selected_atoms->to(this->device));

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
    auto energy_detached = energy_tensor.detach().to(torch::kCPU).to(torch::kFloat64);

    // store the energy returned by the model
    torch::Tensor global_energy;
    if (eflag_atom) {
        auto energies = energy_detached.accessor<double, 2>();
        for (int i=0; i<atom->nlocal + atom->nghost; i++) {
            // TODO: handle out of order samples
            eatom[i] += energies[i][0];
        }

        global_energy = energy_detached.sum(0);
        assert(energy_detached.sizes() == std::vector<int64_t>({1}));
    } else {
        assert(energy_detached.sizes() == std::vector<int64_t>({1, 1}));
        global_energy = energy_detached.reshape({1});
    }

    if (eflag_global) {
        eng_vdwl += global_energy.item<double>();
    }

    // compute forces/virial with backward propagation
    energy_tensor.backward(-torch::ones_like(energy_tensor));
    auto forces_tensor = system->positions().grad().to(torch::kCPU).to(torch::kFloat64);
    auto forces = forces_tensor.accessor<double, 2>();
    for (int i=0; i<atom->nlocal + atom->nghost; i++) {
        atom->f[i][0] += forces[i][0];
        atom->f[i][1] += forces[i][1];
        atom->f[i][2] += forces[i][2];
    }

    assert(!vflag_fdotr);

    if (vflag_global) {
        auto virial_tensor = system_adaptor->strain.grad().to(torch::kCPU).to(torch::kFloat64);
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


void PairMetatensor::load_torch_model(
    const char* path,
    const char* extensions_directory
) {
    // TODO: seach for the model & extensions inside `$LAMMPS_POTENTIALS`?

    if (this->torch_model != nullptr) {
        error->all(FLERR, "torch model is already loaded");
    }

    torch::optional<std::string> extensions = torch::nullopt;
    if (extensions_directory != nullptr) {
        extensions = std::string(extensions_directory);
    }

    try {
        this->torch_model = std::make_unique<torch::jit::Module>(
            metatensor_torch::load_atomistic_model(path, extensions)
        );
    } catch (const c10::Error& e) {
        error->all(FLERR, "failed to load metatensor model at '{}': {}", path, e.what());
    }

    auto capabilities_ivalue = this->torch_model->run_method("capabilities");
    this->capabilities = capabilities_ivalue.toCustomClass<metatensor_torch::ModelCapabilitiesHolder>();

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
