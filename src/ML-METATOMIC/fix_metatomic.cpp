// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Fix metatomic: ML-driven position and momentum prediction

   This fix implements machine learning-driven molecular dynamics where a
   trained model predicts atomic positions and momenta at each timestep.

   The integration scheme:
   1. initial_integrate: ML model predicts new positions and momenta
   2. post_force: Snapshot forces (includes e.g. any added stochastic forces)
   3. final_integrate: Apply force corrections to velocities
------------------------------------------------------------------------- */
#include "metatomic_types.h"
#include "metatomic_system.h"

#include "fix_metatomic.h"

#include "atom.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "group.h"
#include "force.h"
#include "update.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"

#include <vector>
#include <algorithm>

#include <metatomic/torch.hpp>
#include <metatensor/torch.hpp>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMetatomic::FixMetatomic(LAMMPS *lmp, int narg, char **arg): Fix(lmp, narg, arg) {
    time_integrate = 1;  // this tells LAMMPS that this fix advances simulation time
    dynamic_group_allow = 0;  // we don't allow dynamic groups for now

    // Check for multiple MPI processes - not currently supported
    if (comm->nprocs > 1) {
        error->all(FLERR, "fix metatomic does not support multiple MPI processes yet");
    }

    // Determine unit system for the ML model
    // Currently only 'metal' units are fully supported for momenta
    std::string energy_unit;
    std::string length_unit;
    if (strcmp(update->unit_style, "metal") == 0) {
        length_unit = "angstrom";
        this->momentum_conversion_factor = 10.1805057179 / 1000.0;
    } else if (strcmp(update->unit_style, "real") == 0) {
        length_unit = "angstrom";
        this->momentum_conversion_factor = 10.1805057179;
    } else if (strcmp(update->unit_style, "si") == 0) {
        length_unit = "m";
        this->momentum_conversion_factor = 10.1805057179 / 1.6605390666e-22;
    } else {
        error->all(FLERR, "unsupported units '{}' for fix metatomic", update->unit_style);
    }

    // For now, only metal units are fully tested and supported
    if (strcmp(update->unit_style, "metal") != 0) {
        error->all(FLERR, "fix metatomic currently only supports 'metal' units");
    }

    if (narg < 4) {
        error->all(FLERR,
            "Illegal fix metatomic command: expected at least 4 arguments "
            "(fix ID group-ID metatomic model_path ...); got %d", narg
        );
    }

    bool types_are_set = false;
    this->model_path = arg[3];
    this->requested_device = std::nullopt;
    this->extensions_directory = std::nullopt;
    std::vector<int> parsed_types;

    int iarg = 4;
    while (iarg < narg) {
        if (strcmp(arg[iarg], "types") == 0) {
            types_are_set = true;
            // Require exactly atom->ntypes integer values after the "types" keyword.
            iarg++;
            if (iarg + atom->ntypes > narg) {
                error->all(FLERR,
                    "Illegal fix metatomic command: expected %d type values "
                    "after 'types'", atom->ntypes
                );
            }
            for (int ti = 0; ti < atom->ntypes; ++ti) {
                int type = -1;
                const char *argstr = arg[iarg + ti];
                try {
                    type = std::stoi(argstr);
                } catch (const std::invalid_argument &) {
                    error->all(FLERR,
                        "Illegal fix metatomic command: expected integer for type %d, "
                        "got '%s'", ti + 1, argstr
                    );
                } catch (const std::out_of_range &) {
                    error->all(FLERR,
                        "Illegal fix metatomic command: type value out of range "
                        "for argument '%s'", argstr
                    );
                }
                if (type <= 0) {
                    error->all(FLERR, "Illegal fix metatomic command: type %d should be > 0", type);
                }
                parsed_types.push_back(type);
            }
            iarg += atom->ntypes;
        } else if (strcmp(arg[iarg], "device") == 0) {
            if (iarg + 1 >= narg) {
                error->all(FLERR,
                    "Illegal fix metatomic command: 'device' expects an argument "
                    "specifying the device (e.g. cpu, cuda, mps)"
                );
            }
            this->requested_device = std::string(arg[iarg + 1]);
            iarg += 2;
        } else if (strcmp(arg[iarg], "extensions_directory") == 0) {
            if (iarg + 1 >= narg) {
                error->all(FLERR,
                    "Illegal fix metatomic command: 'extensions_directory' expects "
                    "an argument specifying the directory path"
                );
            }
            this->extensions_directory = std::string(arg[iarg + 1]);
            iarg += 2;
        } else {
            error->all(FLERR,
                "Illegal fix metatomic command: unrecognized option '%s' (expected "
                "'types', 'device', or `extensions_directory`)", arg[iarg]
            );
        }
    }

    if (!types_are_set) {
        error->all(FLERR, "Illegal fix metatomic command: no types specified");
    }

    // Allocate and fill the type-mapping (1-based indexing)
    type_mapping = memory->create(type_mapping, atom->ntypes + 1, "fix_metatomic:type_mapping");
    for (int i = 1; i <= atom->ntypes; i++) {
        type_mapping[i] = parsed_types[i - 1];
    }

    this->mta_data = new FixMetatomicData(std::move(length_unit));

    // FlashMD needs position change delta-q and momenta p
    auto positions = torch::make_intrusive<metatomic_torch::ModelOutputHolder>(
        /*quantity =*/ "length",
        /*unit =*/ length_unit,
        /*per_atom =*/ true,
        /*explicit_gradients =*/ std::vector<std::string>{},
        /*description =*/ ""
    );
    this->mta_data->evaluation_options->outputs.insert("positions", positions);

    auto momenta = torch::make_intrusive<metatomic_torch::ModelOutputHolder>(
        /*quantity =*/ "momentum",
        /*unit =*/ "(eV*u)^(1/2)",
        /*per_atom =*/ true,
        /*explicit_gradients =*/ std::vector<std::string>{},
        /*description =*/ ""
    );
    this->mta_data->evaluation_options->outputs.insert("momenta", momenta);
}

FixMetatomic::~FixMetatomic() {
    memory->destroy(type_mapping);
}

/* ---------------------------------------------------------------------- */

int FixMetatomic::setmask() {
    return INITIAL_INTEGRATE | POST_FORCE | FINAL_INTEGRATE;
}

/* ---------------------------------------------------------------------- */

void FixMetatomic::init() {
    int fix_metatomic_index = -1;
    const auto& fixes = modify->get_fix_list();
    auto it = std::find(fixes.begin(), fixes.end(), this);
    if (it != fixes.end()) {
        fix_metatomic_index = int(it - fixes.begin());
    }
    if (fix_metatomic_index != 0) {
        error->all(FLERR, "fix metatomic should be defined as the first fix (before any other fix)");
    }

    if (comm->nprocs > 1) {
        error->all(FLERR,"fix metatomic currently does not support multiple processes");
    }

    if (!type_mapping) {
        error->all(FLERR, "fix metatomic internal error: type_mapping not initialized");
    }

    mta_data->load_model(
        this->lmp,
        this->model_path.c_str(),
        this->extensions_directory ? this->extensions_directory->c_str() : nullptr
    );

    double model_timestep = mta_data->model->attr("module").toModule().attr("timestep").toTensor().item<double>();
    model_timestep = model_timestep * 1e-3;  // fs to ps (metal units)
    if (std::abs(update->dt - model_timestep) > 1e-5 * model_timestep) {
        error->all(FLERR,
            "fix metatomic timestep (dt = {}) does not match the model's expected timestep ({}). "
            "Please set the timestep to match the model.",
            update->dt, model_timestep);
    }

    // Select the device to use based on the model's preference, the user choice
    // and what's available.
    this->pick_device(
        mta_data->device,
        this->requested_device ? this->requested_device->c_str() : nullptr
    );

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

    // Initialize metatensor system object
    auto options = MetatomicSystemOptions{
        this->type_mapping,
        mta_data->max_cutoff,
        mta_data->check_consistency,
        /* requires_grad */ false,
    };
    this->system_adaptor = std::make_unique<MetatomicSystemAdaptor>(lmp, options);

    // We ask LAMMPS for a full neighbor lists because we need to know about
    // ALL pairs, even if options->full_list() is false. We will then filter
    // the pairs to only include each pair once where needed.
    auto request = neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);
    request->set_cutoff(mta_data->max_cutoff);

    auto mincut = mta_data->max_cutoff + neighbor->skin;
    if (comm->get_comm_cutoff() < mincut) {
        if (comm->me == 0) {
            error->warning(FLERR,
                "Increasing communication cutoff to {:.8} for fix metatomic",
                mincut
            );
        }
        comm->cutghostuser = mincut;
    }

    // Translate from the metatomic neighbor lists requests to LAMMPS neighbor
    // lists requests.
    auto requested_nl = mta_data->model->run_method("requested_neighbor_lists");
    for (const auto& ivalue: requested_nl.toList()) {
        auto options = ivalue.get().toCustomClass<metatomic_torch::NeighborListOptionsHolder>();
        auto cutoff = options->engine_cutoff(mta_data->evaluation_options->length_unit());
        assert(cutoff <= mta_data->max_cutoff);

        this->system_adaptor->add_nl_request(cutoff, options);
    }

    // HACK: Explicitly set the binsize for the neighbor list if there is no
    // pair_style that would set it instead.
    //
    // Otherwise, the default binsize of box[0] is used, which crashes kokkos
    // for large-ish boxes (~40A), and slow down the simulation for non-kokkos.
    if (strcmp(force->pair_style, "none") == 0) {
        neighbor->binsize_user = 0.5 * mta_data->max_cutoff;
        neighbor->binsizeflag = 1;
    }
    // END HACK
}

void FixMetatomic::pick_device(c10::Device& device, const char* requested) {
    torch::optional<std::string> requested_string;
    torch::DeviceType device_type;

    if (requested != nullptr) {
        requested_string = std::string(requested);
    } else {
        requested_string = torch::nullopt;
    }

    try {
        device_type = metatomic_torch::pick_device(
            this->mta_data->capabilities->supported_devices,
            requested_string
        );
    } catch (const c10::Error& e) {
        error->one(FLERR, "fix metatomic: {}", e.what());
    }

    if (device_type == torch::DeviceType::CUDA) {
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
        auto device_index = local_rank % torch::cuda::device_count();
        device = torch::Device(device_type, static_cast<torch::DeviceIndex>(device_index));
    } else {
        device = torch::Device(device_type);
    }
}

void FixMetatomic::init_list(int id, NeighList *ptr) {
    mta_list = ptr;
}

void FixMetatomic::initial_integrate(int /*vflag*/) {
    // This function performs ML-driven position and momentum updates
    // It uses a trained model to predict new positions and momenta at each timestep

    double** x = atom->x;
    double** v = atom->v;
    double** f = atom->f;
    double* rmass = atom->rmass;

    int nlocal = atom->nlocal;

    double *mass = atom->mass;
    int *type = atom->type;
    int *mask = atom->mask;
    if (igroup == atom->firstgroup) {
        nlocal = atom->nfirst;
    }

    // Apply velocity corrections from forces added after post_force
    // This handles stochastic forces from Langevin thermostats by applying only the
    // incremental force (f_current - f_snapshot) to velocities. The remaining half-step
    // is done in final_integrate(); this is the first O in an OBABO integrator.
    double dtf = 0.5 * update->dt * force->ftm2v;
    double m_i;
    for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
            // Apply any force added by other fixes to velocities
            // rmass is per-atom mass (if used), otherwise use type-based mass
            m_i = rmass ? rmass[i] : mass[type[i]];
            v[i][0] += f[i][0] * dtf / m_i;
            v[i][1] += f[i][1] * dtf / m_i;
            v[i][2] += f[i][2] * dtf / m_i;
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
        dtype,
        mta_data->device
    );

    // add the required additional inputs
    this->system_adaptor->add_masses(system, 1.0);
    this->system_adaptor->add_momenta(system, this->momentum_conversion_factor);

    // Configure selected atoms for evaluation
    // Only run the calculation for atoms in the current group
    mta_data->selected_atoms_values.resize_({group->count(igroup), 2});
    mta_data->selected_atoms_values.index_put_({torch::indexing::Slice(), 0}, 0);
    int64_t idx = 0;
    for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
            mta_data->selected_atoms_values.index_put_({idx, 1}, i);
            idx++;
        }
    }

    auto selected_atoms = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        std::vector<std::string>{"system", "atom"}, mta_data->selected_atoms_values
    );
    mta_data->evaluation_options->set_selected_atoms(selected_atoms);

    // Call the ML model to predict new positions and momenta
    torch::IValue result_ivalue;
    try {
        result_ivalue = mta_data->model->forward({
            std::vector<metatomic_torch::System>{system},
            mta_data->evaluation_options,
            mta_data->check_consistency
        });
    } catch (const std::exception& e) {
        error->all(FLERR, "error evaluating the torch model: {}", e.what());
    }

    // Extract results from the model output
    auto result = result_ivalue.toGenericDict();

    // Extract predicted positions
    auto positions_map = result.at("positions").toCustomClass<metatensor_torch::TensorMapHolder>();
    auto positions_block = metatensor_torch::TensorMapHolder::block_by_id(positions_map, 0);
    auto positions = positions_block->values().squeeze(-1).to(torch::kCPU).to(torch::kFloat64);
    auto positions_samples = positions_block->samples()->values().to(torch::kCPU).contiguous();
    assert(positions_block->samples()->size() == 2);
    assert(positions_block->samples()->names()[0] == "system");
    assert(positions_block->samples()->names()[1] == "atom");

    // Extract predicted momenta
    auto momenta_map = result.at("momenta").toCustomClass<metatensor_torch::TensorMapHolder>();
    auto momenta_block = metatensor_torch::TensorMapHolder::block_by_id(momenta_map, 0);
    auto momenta = momenta_block->values().squeeze(-1).to(torch::kCPU).to(torch::kFloat64);

    // we use the positions samples to map back to LAMMPS atoms, so we need to
    // check that the samples are the same for momenta
    assert(*momenta_block->samples() == *positions_block->samples());

    // Convert momenta back from model units to LAMMPS velocity units
    // This reverses the unit conversion applied before the model call
    momenta = momenta / this->momentum_conversion_factor;

    // Get old center of mass (and its velocity) before updating positions and velocities
    std::array<double, 3> com_old = {0.0, 0.0, 0.0};
    std::array<double, 3> com_velocity_old = {0.0, 0.0, 0.0};
    double total_mass = 0.0;
    for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
            double m_i = rmass ? rmass[i] : mass[type[i]];
            com_old[0] += x[i][0] * m_i;
            com_old[1] += x[i][1] * m_i;
            com_old[2] += x[i][2] * m_i;
            com_velocity_old[0] += v[i][0] * m_i;
            com_velocity_old[1] += v[i][1] * m_i;
            com_velocity_old[2] += v[i][2] * m_i;
            total_mass += m_i;
        }
    }
    if (total_mass > 0.0) {
        com_old[0] /= total_mass;
        com_old[1] /= total_mass;
        com_old[2] /= total_mass;
        com_velocity_old[0] /= total_mass;
        com_velocity_old[1] /= total_mass;
        com_velocity_old[2] /= total_mass;
    }

    // Apply ML predictions to LAMMPS atoms
    auto positions_accessor = positions.accessor<double, 2>();
    auto momenta_accessor = momenta.accessor<double, 2>();
    auto positions_samples_accessor = positions_samples.accessor<int32_t, 2>();
    auto& mta_to_lmp = this->system_adaptor->mta_to_lmp;
    for (int64_t i = 0; i < positions.size(0); i++) {
        auto atom_i = mta_to_lmp[positions_samples_accessor[i][1]];
        assert(atom_i < nlocal);
        if (mask[atom_i] & groupbit) {
            double m_i = rmass ? rmass[atom_i] : mass[type[atom_i]];

            // Update positions with ML predictions
            x[atom_i][0] = positions_accessor[i][0];
            x[atom_i][1] = positions_accessor[i][1];
            x[atom_i][2] = positions_accessor[i][2];

            // Update velocities from predicted momenta
            v[atom_i][0] = momenta_accessor[i][0] / m_i;
            v[atom_i][1] = momenta_accessor[i][1] / m_i;
            v[atom_i][2] = momenta_accessor[i][2] / m_i;
        }
    }

    std::array<double, 3> com_new = {0.0, 0.0, 0.0};
    std::array<double, 3> com_velocity_new = {0.0, 0.0, 0.0};
    for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
            double m_i = rmass ? rmass[i] : mass[type[i]];
            com_new[0] += x[i][0] * m_i;
            com_new[1] += x[i][1] * m_i;
            com_new[2] += x[i][2] * m_i;
            com_velocity_new[0] += v[i][0] * m_i;
            com_velocity_new[1] += v[i][1] * m_i;
            com_velocity_new[2] += v[i][2] * m_i;
        }
    }
    if (total_mass > 0.0) {
        com_new[0] /= total_mass;
        com_new[1] /= total_mass;
        com_new[2] /= total_mass;
        com_velocity_new[0] /= total_mass;
        com_velocity_new[1] /= total_mass;
        com_velocity_new[2] /= total_mass;
    }

    // Adjust positions and velocities to preserve center of mass motion, namely
    // conservation of momentum of the center of mass and uniform linear motion of the
    // center of mass.
    for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
            // Update positions with ML predictions
            x[i][0] = x[i][0] - com_new[0] + com_old[0] + com_velocity_old[0] * update->dt;
            x[i][1] = x[i][1] - com_new[1] + com_old[1] + com_velocity_old[1] * update->dt;
            x[i][2] = x[i][2] - com_new[2] + com_old[2] + com_velocity_old[2] * update->dt;
            v[i][0] = v[i][0] - com_velocity_new[0] + com_velocity_old[0];
            v[i][1] = v[i][1] - com_velocity_new[1] + com_velocity_old[1];
            v[i][2] = v[i][2] - com_velocity_new[2] + com_velocity_old[2];
        }
    }
}

void FixMetatomic::post_force(int /*vflag*/) {
    // Set the forces that comes from pair_style, bond_style, etc. to zero.
    //
    // This way we can isolate any forces that are added after this point (e.g.
    // Langevin thermostat forces) and add them during final_integrate().
    //
    // Crucially, this means that fix metatomic needs to be the first fix in the
    // post_force() sequence, i.e., the user must have it before any other fix
    // that adds forces in the input script.

    double **f = atom->f;
    int *mask = atom->mask;

    int nlocal = atom->nlocal;
    if (igroup == atom->firstgroup) {
        nlocal = atom->nfirst;
    }

    for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
            f[i][0] = 0.0;
            f[i][1] = 0.0;
            f[i][2] = 0.0;
        }
    }
}

void FixMetatomic::final_integrate() {
    // Apply velocity corrections from forces that were added after post_force
    //
    // This handles stochastic forces from Langevin thermostats:
    // - initial_integrate: ML model updates positions and velocities
    // - post_force: we snapshot forces (includes pair, bond, and Langevin
    //   forces)
    // - Between post_force and final_integrate: additional forces may be added
    // - final_integrate: we apply only the force difference as a velocity
    //   correction
    //
    // This ensures Langevin forces properly affect the dynamics while allowing
    // the ML model to handle the deterministic evolution. The first half-step
    // is done in initial_integrate(); this is the second O in an OBABO integrator.

    double dtf = 0.5 * update->dt * force->ftm2v;

    double** v = atom->v;
    double** f = atom->f;
    double* rmass = atom->rmass;
    double* mass = atom->mass;
    int* type = atom->type;
    double m_i;

    int nlocal = atom->nlocal;
    int* mask = atom->mask;
    if (igroup == atom->firstgroup) nlocal = atom->nfirst;

    for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
            // Apply any force added by other fixes to velocities
            m_i = rmass ? rmass[i] : mass[type[i]];
            v[i][0] += f[i][0] * dtf / m_i;
            v[i][1] += f[i][1] * dtf / m_i;
            v[i][2] += f[i][2] * dtf / m_i;
        }
    }
}
