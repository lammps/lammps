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
#include "domain.h"

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

static double compute_volume(Domain* domain) {
    // from Thermo::compute_vol
    if (domain->dimension == 3) {
        return domain->xprd * domain->yprd * domain->zprd;
    } else {
        return domain->xprd * domain->yprd;
    }
}

PairMetatomic::PairMetatomic(LAMMPS *lmp):
    Pair(lmp),
    type_mapping(nullptr),
    system_adaptor(nullptr),
    scale(1.0)
{
    if (strcmp(update->unit_style, "real") == 0) {
        this->length_unit = "angstrom";
        this->energy_unit = "kcal/mol";
    } else if (strcmp(update->unit_style, "metal") == 0) {
        this->length_unit = "angstrom";
        this->energy_unit = "eV";
    } else if (strcmp(update->unit_style, "si") == 0) {
        this->length_unit = "meter";
        this->energy_unit = "joule";
    } else if (strcmp(update->unit_style, "electron") == 0) {
        this->length_unit = "Bohr";
        this->energy_unit = "Hartree";
    } else {
        error->all(FLERR, "unsupported units '{}' for pair metatomic ", update->unit_style);
    }

    // we might not be running a pure pair potential,
    // so we can not compute virial as fdotr
    this->no_virial_fdotr_compute = 1;

    this->mta_data = new PairMetatomicData(this->length_unit);
    // use a default uncertainty threshold of 100 meV/atom
    this->mta_data->uncertainty_threshold = 0.1 * metatomic_torch::unit_conversion_factor("energy", "eV", energy_unit);

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

// called when finding `pair_style metatomic` in the input
void PairMetatomic::settings(int argc, char ** argv) {
    if (argc == 0) {
        error->all(FLERR, "expected at least 1 argument to pair_style metatomic, got {}", argc);
    }

    const char* model_path = argv[0];
    const char* extensions_directory = nullptr;
    const char* requested_device = nullptr;
    bool do_uncertainty = true;
    const char* variant = nullptr;
    const char* variant_energy = nullptr;
    const char* variant_energy_uq = nullptr;
    const char* variant_nc_forces = nullptr;
    const char* variant_nc_stress = nullptr;
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
        } else if (strcmp(argv[i], "non_conservative") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected <on/off> after 'non_conservative' in pair_style metatomic, got nothing");
            } else if (strcmp(argv[i + 1], "on") == 0) {
                mta_data->non_conservative = true;
            } else if (strcmp(argv[i + 1], "off") == 0) {
                mta_data->non_conservative = false;
            } else {
                error->all(FLERR, "expected <on/off> after 'non_conservative' in pair_style metatomic, got '{}'", argv[i + 1]);
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
        } else if (strcmp(argv[i], "scale") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected a number after 'scale' in pair_style metatomic, got nothing");
            }
            this->scale = utils::numeric(FLERR, argv[i + 1], false, lmp);
            i += 1;
        } else if (strcmp(argv[i], "uncertainty_threshold") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected a number or off after 'uncertainty_threshold' in pair_style metatomic, got nothing");
            } else if (strcmp(argv[i + 1], "off") == 0) {
                do_uncertainty = false;
            } else {
                mta_data->uncertainty_threshold = utils::numeric(FLERR, argv[i + 1], false, lmp);
            }

            if (mta_data->uncertainty_threshold <= 0) {
                error->all(FLERR, "'uncertainty_threshold' in pair_style metatomic must be positive");
            }
            i += 1;
        } else if (strcmp(argv[i], "variant") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected a name after 'variant' in pair_style metatomic, got nothing");
            }
            variant = argv[i + 1];
            i += 1;
        } else if (strcmp(argv[i], "variant/energy") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected a name or 'off' after 'variant/energy' in pair_style metatomic, got nothing");
            }
            variant_energy = argv[i + 1];
            i += 1;
        } else if (strcmp(argv[i], "variant/energy_uncertainty") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected a name or 'off' after 'variant/energy_uncertainty' in pair_style metatomic, got nothing");
            }
            variant_energy_uq = argv[i + 1];
            i += 1;
        } else if (strcmp(argv[i], "variant/non_conservative_forces") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected a name or 'off' after 'variant/non_conservative_forces' in pair_style metatomic, got nothing");
            }
            variant_nc_forces = argv[i + 1];
            i += 1;
        } else if (strcmp(argv[i], "variant/non_conservative_stress") == 0) {
            if (i == argc - 1) {
                error->all(FLERR, "expected a name or 'off' after 'variant/non_conservative_stress' in pair_style metatomic, got nothing");
            }
            variant_nc_stress = argv[i + 1];
            i += 1;
        } else {
            error->all(FLERR, "unexpected argument to pair_style metatomic: '{}'", argv[i]);
        }
    }

    // Load the model and get it's capabilities (including supported devices)
    mta_data->load_model(this->lmp, model_path, extensions_directory);

    // Set and resolve the variants to use
    mta_data->energy_key = "energy";
    mta_data->energy_uq_key = "energy_uncertainty";
    mta_data->nc_forces_key = "non_conservative_forces";
    mta_data->nc_stress_key = "non_conservative_stress";

    // Apply global variant (applies to all)
    if (variant != nullptr) {
        mta_data->energy_key += "/" + std::string(variant);
        mta_data->energy_uq_key += "/" + std::string(variant);
        mta_data->nc_forces_key += "/" + std::string(variant);
        mta_data->nc_stress_key += "/" + std::string(variant);
    }

    // Apply variant/energy
    if (variant_energy != nullptr) {
        if (strcmp(variant_energy, "off") == 0) {
            mta_data->energy_key = "energy";
        } else {
            mta_data->energy_key = "energy/" + std::string(variant_energy);
        }
    }

    // Apply variant/energy_uncertainty
    if (variant_energy_uq != nullptr) {
        if (strcmp(variant_energy_uq, "off") == 0) {
            mta_data->energy_uq_key = "energy_uncertainty";
        } else {
            mta_data->energy_uq_key = "energy_uncertainty/" + std::string(variant_energy_uq);
        }
    }

    // Handle non-conservative variants
    bool has_nc_forces = variant_nc_forces != nullptr;
    bool has_nc_stress = variant_nc_stress != nullptr;

    if (has_nc_forces && has_nc_stress) {
        bool forces_none = strcmp(variant_nc_forces, "off") == 0;
        bool stress_none = strcmp(variant_nc_stress, "off") == 0;
        if (forces_none != stress_none) {
            error->all(FLERR,
                "if both 'variant/non_conservative_stress' and "
                "'variant/non_conservative_forces' are set, they must either "
                "both be 'off' or both not be 'off'");
        }
    } else if (has_nc_forces && !has_nc_stress) {
        if (strcmp(variant_nc_forces, "off") != 0) {
            error->all(FLERR,
                "'variant/non_conservative_forces' is set but "
                "'variant/non_conservative_stress' is not; "
                "both must be set together or both be 'off'");
        }
    } else if (!has_nc_forces && has_nc_stress) {
        if (strcmp(variant_nc_stress, "off") != 0) {
            error->all(FLERR,
                "'variant/non_conservative_stress' is set but "
                "'variant/non_conservative_forces' is not; "
                "both must be set together or both be 'off'");
        }
    }

    if (has_nc_forces) {
        if (strcmp(variant_nc_forces, "off") == 0) {
            mta_data->nc_forces_key = "non_conservative_forces";
        } else {
            mta_data->nc_forces_key = "non_conservative_forces/" + std::string(variant_nc_forces);
        }
    }

    if (has_nc_stress) {
        if (strcmp(variant_nc_stress, "off") == 0) {
            mta_data->nc_stress_key = "non_conservative_stress";
        } else {
            mta_data->nc_stress_key = "non_conservative_stress/" + std::string(variant_nc_stress);
        }
    }

    // Check that the model has the required outputs
    const auto& outputs = mta_data->capabilities->outputs();
    auto energy_output = outputs.find(mta_data->energy_key);
    // LAMMPS assume that an energy will be available
    if (energy_output == outputs.end()) {
        lmp->error->all(FLERR,
            "the model at '{}' does not have an 'energy' output, "
            "we can not use it with pair_style metatomic.",
            model_path
        );
    }

    mta_data->is_energy_output_per_atom = energy_output->value()->per_atom;
    mta_data->energy_output = torch::make_intrusive<metatomic_torch::ModelOutputHolder>();
    mta_data->energy_output->set_quantity("energy");
    mta_data->energy_output->set_unit(this->energy_unit);

    auto uncertainty_output = outputs.find(mta_data->energy_uq_key);
    if (uncertainty_output != outputs.end()) {
        if (do_uncertainty && uncertainty_output->value()->per_atom) {
            // TODO: maybe if there is a global uncertainty output we should use
            // that as a fallback?

            mta_data->uncertainty_output = torch::make_intrusive<metatomic_torch::ModelOutputHolder>();
            mta_data->uncertainty_output->set_quantity("energy");
            mta_data->uncertainty_output->set_unit(this->energy_unit);
            mta_data->uncertainty_output->per_atom = true;

            if (comm->me == 0) {
                auto message = "Found '{}' output, we will check for atoms with high uncertainty on the energy predictions";
                if (screen) {
                    fprintf(screen, "%s\n", fmt::format(message, mta_data->energy_uq_key).c_str());
                }
                if (logfile) {
                    fprintf(logfile,"%s\n", fmt::format(message, mta_data->energy_uq_key).c_str());
                }
            }
        }
    }

    if (mta_data->non_conservative) {
        auto nc_forces = outputs.find(mta_data->nc_forces_key);
        if (nc_forces == outputs.end()) {
            error->all(FLERR,
                "the model at '{}' does not have a '{}' output, "
                "we can not enable non_conservative simulations",
                model_path, mta_data->nc_forces_key
            );
        }

        if (!nc_forces->value()->per_atom) {
            error->all(FLERR,
                "the '{}' output of the model at '{}' "
                "can not produce per-atom output, we can not enable non_conservative simulations",
                mta_data->nc_forces_key, model_path
            );
        }

        mta_data->nc_forces_output = torch::make_intrusive<metatomic_torch::ModelOutputHolder>();
        mta_data->nc_forces_output->set_quantity("force");
        mta_data->nc_forces_output->set_unit(this->energy_unit + "/" + this->length_unit);
        mta_data->nc_forces_output->per_atom = true;

        auto nc_stress = outputs.find(mta_data->nc_stress_key);
        if (nc_stress != outputs.end()) {
            mta_data->nc_stress_output = torch::make_intrusive<metatomic_torch::ModelOutputHolder>();
            mta_data->nc_stress_output->set_quantity("pressure");
            mta_data->nc_stress_output->set_unit(this->energy_unit + "/" + this->length_unit + "^3");
            mta_data->nc_stress_output->per_atom = false;
        } else {
            mta_data->nc_stress_output = nullptr;
        }
    }

    // Select the device to use based on the model's preference, the user choice
    // and what's available.
    this->pick_device(mta_data->device, requested_device);

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

void PairMetatomic::pick_device(torch::Device& device, const char* requested) {

    torch::optional<std::string> requested_string;
    std::string device_string;

    if (requested != nullptr) {
        requested_string = std::string(requested);
    } else {
        requested_string = torch::nullopt;
    }

    try {
        device_string = metatomic_torch::pick_device(
            this->mta_data->capabilities->supported_devices,
            requested_string
        );
    } catch (const c10::Error& e) {
        error->all(FLERR, "pair_style metatomic: {}", e.what());
    }

    if (device_string == "cuda") {
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
        device = torch::Device("cuda:" + std::to_string(gpu_to_use));
    } else {
        device = torch::Device(device_string);
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
        !(mta_data->non_conservative),
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
        assert(cutoff <= mta_data->max_cutoff);

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

    ev_init(eflag, vflag);

    mta_data->evaluation_options->outputs.clear();
    // we need an energy output if the energy was explicitly requested (through
    // `eflag_either`), or when running in standard/conservative mode, because
    // we'll get the forces as the gradient of the energy through autodiff.
    if (eflag_either || !mta_data->non_conservative) {
        if (eflag_atom) {
            if (!mta_data->is_energy_output_per_atom) {
                error->all(FLERR,
                    "the model at '{}' does not support per-atom 'energy' output",
                    mta_data->model_path
                );
            }
            mta_data->energy_output->per_atom = true;
        } else {
            assert(eflag_global);
            mta_data->energy_output->per_atom = false;
        }

        mta_data->evaluation_options->outputs.insert(mta_data->energy_key, mta_data->energy_output);

        if (mta_data->uncertainty_output != nullptr) {
            mta_data->evaluation_options->outputs.insert(mta_data->energy_uq_key, mta_data->uncertainty_output);
        }
    }

    if (mta_data->non_conservative) {
        mta_data->evaluation_options->outputs.insert(mta_data->nc_forces_key, mta_data->nc_forces_output);
        if (vflag_global) {
            if (mta_data->nc_stress_output == nullptr) {
                error->all(FLERR,
                    "the model at '{}' does not have a '{}' output, "
                    "we can not run non_conservative simulations that require computing the stress/virial",
                    mta_data->model_path, mta_data->nc_stress_key
                );
            }
            mta_data->evaluation_options->outputs.insert(mta_data->nc_stress_key, mta_data->nc_stress_output);
        }
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
        std::vector<std::string>{"system", "atom"},
        mta_data->selected_atoms_values,
        metatensor::assume_unique{}
    );
    mta_data->evaluation_options->set_selected_atoms(selected_atoms);

    torch::IValue results_ivalue;
    try {
        auto _ = MetatomicTimer("running Model::forward");
        results_ivalue = mta_data->model->forward({
            std::vector<metatomic_torch::System>{system},
            mta_data->evaluation_options,
            mta_data->check_consistency
        });
    } catch (const std::exception& e) {
        error->all(FLERR, "error evaluating the torch model: {}", e.what());
    }

    auto results = results_ivalue.toGenericDict();

    // check the max uncertainty
    if (mta_data->uncertainty_output != nullptr) {
        auto uncertainty = results.at(mta_data->energy_uq_key).toCustomClass<metatensor_torch::TensorMapHolder>();
        auto uncertainty_block = metatensor_torch::TensorMapHolder::block_by_id(uncertainty, 0);
        assert(uncertainty_block->values().sizes().size() == 2);
        assert(uncertainty_block->values().size(1) == 1);

        auto atoms_above_thresholds = uncertainty_block->values().reshape({-1}) > mta_data->uncertainty_threshold;
        if (torch::any(atoms_above_thresholds).to(torch::kCPU).item<bool>()) {
            auto atoms = uncertainty_block->samples()->column("atom").index({atoms_above_thresholds});
            std::ostringstream atoms_message;
            atoms_message << "atoms at index [";

            // only print the first 10 atoms above the threshold to avoid
            // flooding the output
            for (size_t i=0; i<std::min(static_cast<int64_t>(10), atoms.size(0)); i++) {
                if (i > 0) {
                    atoms_message << ", ";
                }
                atoms_message << atoms[i].item<int32_t>();
            }
            atoms_message << "]";

            if (atoms.size(0) > 10) {
                atoms_message << " and " << (atoms.size(0) - 10) << " more";
            }

            error->warning(FLERR,
                "The uncertainty on atomic energies for {} are larger than "
                "the threshold of {}. Be careful when analyzing the results, "
                "and consider retraining the model to better describe these "
                "configurations.",
                atoms_message.str(), mta_data->uncertainty_threshold
            );
        }
    }

    torch::Tensor energy_tensor;
    metatensor_torch::Labels energy_samples;

    // get the energy if we need to compute the energy, or if we are using it to
    // get the forces/virial with autograd
    if (eflag_either || !mta_data->non_conservative) {
        auto energy = results.at(mta_data->energy_key).toCustomClass<metatensor_torch::TensorMapHolder>();
        auto energy_block = metatensor_torch::TensorMapHolder::block_by_id(energy, 0);
        energy_tensor = energy_block->values();
        energy_samples = energy_block->samples();
    }

    torch::Tensor forces_tensor;
    torch::Tensor virial_tensor;

    if (mta_data->non_conservative) {
        auto forces = results.at(mta_data->nc_forces_key).toCustomClass<metatensor_torch::TensorMapHolder>();
        auto forces_block = metatensor_torch::TensorMapHolder::block_by_id(forces, 0);
        forces_tensor = forces_block->values().squeeze(-1);
        forces_tensor = forces_tensor.to(torch::kCPU).to(torch::kFloat64);

        if (vflag_global) {
            auto stress = results.at(mta_data->nc_stress_key).toCustomClass<metatensor_torch::TensorMapHolder>();
            auto stress_block = metatensor_torch::TensorMapHolder::block_by_id(stress, 0);
            auto stress_tensor = stress_block->values().squeeze(0).squeeze(-1);
            virial_tensor = - stress_tensor * compute_volume(domain);
            virial_tensor = virial_tensor.to(torch::kCPU).to(torch::kFloat64);
        }
    } else {
        // compute forces/virial on device with backward propagation
        // reset gradients to zero before calling backward
        this->system_adaptor->positions.mutable_grad() = torch::Tensor();
        this->system_adaptor->strain.mutable_grad() = torch::Tensor();

        auto _ = MetatomicTimer("running Model::backward");
        energy_tensor.backward(-torch::ones_like(energy_tensor));

        forces_tensor = this->system_adaptor->positions.grad();
        virial_tensor = this->system_adaptor->strain.grad();
    }

    {
        auto _ = MetatomicTimer("storing model output in LAMMPS data structures");

        // store the energy if requested
        if (eflag_either) {
            // move results to cpu for storing
            auto energy_detached = energy_tensor.detach().to(torch::kCPU).to(torch::kFloat64);

            // store the energy returned by the model
            if (eflag_atom) {
                assert(mta_data->energy_output->per_atom);
                assert(energy_samples->size() == 2);
                assert(energy_samples->names()[0] == "system");
                assert(energy_samples->names()[1] == "atom");

                auto samples_values = energy_samples->values().to(torch::kCPU);
                auto samples = samples_values.accessor<int32_t, 2>();

                assert(samples_values.sizes() == mta_data->selected_atoms_values.sizes());

                auto energies = energy_detached.accessor<double, 2>();
                for (int64_t i=0; i<energy_samples->count(); i++) {
                    assert(samples[i][0] == 0);
                    // handle potentially out of order samples in
                    // the per-atom energy tensor
                    auto atom_i = samples[i][1];
                    assert(atom_i < atom->nlocal + atom->nghost);
                    eatom[atom_i] += this->scale * energies[i][0];
                }
            }

            if (eflag_global) {
                torch::Tensor global_energy;
                if (mta_data->energy_output->per_atom) {
                    global_energy = energy_detached.sum(0);
                    assert(energy_detached.sizes() == std::vector<int64_t>({1}));
                } else {
                    assert(energy_samples->size() == 1);
                    assert(energy_samples->names()[0] == "system");

                    assert(energy_detached.sizes() == std::vector<int64_t>({1, 1}));
                    global_energy = energy_detached.reshape({1});
                }

                eng_vdwl += this->scale * global_energy.item<double>();
            }
        }

        // store forces/virial
        this->store_forces(forces_tensor);

        assert(!vflag_fdotr);

        if (vflag_global) {
            auto virial_cpu = virial_tensor.to(torch::kCPU);
            assert(virial_cpu.is_cpu() && virial_cpu.scalar_type() == torch::kFloat64);

            auto predicted_virial = virial_cpu.template accessor<double, 2>();

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

void PairMetatomic::pre_compute() {}

void PairMetatomic::store_forces(const at::Tensor& forces_tensor) {
    assert(forces_tensor.is_cpu() && forces_tensor.scalar_type() == torch::kFloat64);

    int num_forces_to_update;
    if (mta_data->non_conservative) {
        num_forces_to_update = atom->nlocal;
    } else {
        num_forces_to_update = atom->nlocal + atom->nghost;
    }

    auto forces = forces_tensor.accessor<double, 2>();
    for (int i=0; i<num_forces_to_update; i++) {
        atom->f[i][0] += this->scale * forces[i][0];
        atom->f[i][1] += this->scale * forces[i][1];
        atom->f[i][2] += this->scale * forces[i][2];
    }
}
