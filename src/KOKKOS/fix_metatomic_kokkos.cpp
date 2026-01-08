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
   Fix metatomic/kk: Kokkos version of ML-driven position and momentum prediction

   This is the Kokkos-enabled version of fix metatomic. It uses Kokkos views
   for data access and the MetatomicSystemAdaptorKokkos for efficient data
   transfer between LAMMPS and the ML model.
------------------------------------------------------------------------- */

#include "fix_metatomic_kokkos.h"

#include "error.h"
#include "neigh_request.h"
#include "atom_masks.h"
#include "group.h"
#include "force.h"
#include "update.h"
#include "neighbor_kokkos.h"

#include "atom_kokkos.h"
#include "metatomic_system_kokkos.h"
#include "metatomic_types.h"

#include <algorithm>
#include <cctype>

using namespace LAMMPS_NS;
using namespace FixConst;

// LAMMPS uses `LAMMPS_NS::tagint` and `int` for tags and neighbor lists, respectively.
// For the moment, we require both to be int32_t for this interface
static_assert(std::is_same_v<LAMMPS_NS::tagint, int32_t>, "Error: LAMMPS_NS::tagint must be int32_t to compile fix metatomic/kk");
static_assert(std::is_same_v<int, int32_t>, "Error: int must be int32_t to compile fix metatomic/kk");

template<typename T, class DeviceType>
using UnmanagedView = Kokkos::View<T, Kokkos::LayoutRight, DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixMetatomicKokkos<DeviceType>::FixMetatomicKokkos(LAMMPS *lmp, int narg, char **arg): FixMetatomic(lmp, narg, arg) {
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = X_MASK | V_MASK | F_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixMetatomicKokkos<DeviceType>::~FixMetatomicKokkos() {}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMetatomicKokkos<DeviceType>::init() {
    FixMetatomic::init();

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
        /* requires_grad */ false,
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

    // Sync mass data to device
    atomKK->k_mass.modify_host();
    atomKK->k_mass.sync<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMetatomicKokkos<DeviceType>::pick_device(c10::Device& device, const char* requested) {
    // Pick device based on Kokkos execution space
    device = KokkosDeviceToTorch<DeviceType>::convert();

    if (requested != nullptr) {
        auto requested_str = std::string(requested);
        std::transform(requested_str.begin(), requested_str.end(), requested_str.begin(), ::tolower);
        if (c10::DeviceTypeName(device.type(), /*lower_case=*/true) != requested_str) {
            error->all(FLERR,
                "requested device '{}' does not match the device being used by kokkos '{}', "
                "use the non-kokkos version of this fix to use a different "
                "device for the model and LAMMPS",
                requested, device.str()
            );
        }
    }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMetatomicKokkos<DeviceType>::initial_integrate(int /*vflag*/) {
    // ML-driven position and momentum updates using Kokkos
    // This is the main integration step where the ML model predicts new positions and momenta

    // Get views to atom data on device
    auto x = atomKK->k_x.view<DeviceType>();
    auto v = atomKK->k_v.view<DeviceType>();
    auto f = atomKK->k_f.view<DeviceType>();
    auto rmass = atomKK->k_rmass.view<DeviceType>();
    auto mass = atomKK->k_mass.view<DeviceType>();
    auto type = atomKK->k_type.view<DeviceType>();
    auto mask = atomKK->k_mask.view<DeviceType>();

    auto groupbit = this->groupbit;

    // Sync data to execution space and immediately claim ownership
    // This prevents output->write() from causing data corruption on next timestep
    atomKK->sync(execution_space, datamask_read);
    atomKK->modified(execution_space, datamask_modify);

    int nlocal = atomKK->nlocal;
    if (igroup == atomKK->firstgroup) {
        nlocal = atomKK->nfirst;
    }

    // Apply velocity corrections from forces added after `post_force`.
    //
    // This handles stochastic forces from Langevin thermostats by applying only
    // the incremental force to velocities. The remaining half-step is done in
    // final_integrate(); this is the first O in an OBABO integrator.
    double dtf = 0.5 * update->dt * force->ftm2v;
    bool use_rmass = rmass.data() != nullptr;
    Kokkos::parallel_for(
        nlocal,
        KOKKOS_LAMBDA(int i) {
            if (mask[i] & groupbit) {
                double mass_i = use_rmass ? rmass[i] : mass[type[i]];
                double dtfm = dtf / mass_i;

                v(i, 0) += f(i, 0) * dtfm;
                v(i, 1) += f(i, 1) * dtfm;
                v(i, 2) += f(i, 2) * dtfm;
            }
        }
    );

    // Determine dtype for the model
    auto dtype = torch::kFloat64;
    if (mta_data->capabilities->dtype() == "float64") {
        dtype = torch::kFloat64;
    } else if (mta_data->capabilities->dtype() == "float32") {
        dtype = torch::kFloat32;
    } else {
        error->all(FLERR, "the model requested an unsupported dtype '{}'", mta_data->capabilities->dtype());
    }

    // The v (and x) data is about to be read by torch to create the metatomic System,
    // we need to synchronize to make sure kokkos kernels finished their execution
    Kokkos::fence();

    // Transform from LAMMPS to metatomic System using Kokkos adaptor
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
    auto selected_atoms_kk = UnmanagedView<int32_t**, DeviceType>(
        mta_data->selected_atoms_values.data_ptr<int32_t>(),
        mta_data->selected_atoms_values.size(0),
        2
    );
    auto d_atom_i = Kokkos::View<int32_t, LMPDeviceType>("atom_i");
    Kokkos::deep_copy(d_atom_i, 0);
    Kokkos::parallel_for(
        nlocal,
        KOKKOS_LAMBDA(size_t i) {
            if (mask[i] & groupbit) {
                selected_atoms_kk(i, 0) = 0; // system index
                auto atom_i = Kokkos::atomic_fetch_add(&d_atom_i(), 1);
                selected_atoms_kk(i, 1) = atom_i;
            }
        }
    );

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

    // Extract predicted positions (keep on device)
    auto positions_map = result.at("positions").toCustomClass<metatensor_torch::TensorMapHolder>();
    auto positions_block = metatensor_torch::TensorMapHolder::block_by_id(positions_map, 0);
    auto positions = positions_block->values().squeeze(-1).to(mta_data->device).to(torch::kFloat64).contiguous();
    auto positions_samples = positions_block->samples()->values().contiguous();
    assert(positions_block->samples()->size() == 2);
    assert(positions_block->samples()->names()[0] == "system");
    assert(positions_block->samples()->names()[1] == "atom");


    // Extract predicted momenta (keep on device)
    auto momenta_map = result.at("momenta").toCustomClass<metatensor_torch::TensorMapHolder>();
    auto momenta_block = metatensor_torch::TensorMapHolder::block_by_id(momenta_map, 0);
    auto momenta = momenta_block->values().squeeze(-1).to(mta_data->device).to(torch::kFloat64);

    // we use the positions samples to map back to LAMMPS atoms, so we need to
    // check that the samples are the same for momenta
    assert(*momenta_block->samples() == *positions_block->samples());

    // Convert momenta back from model units to LAMMPS velocity units
    momenta = momenta / this->momentum_conversion_factor;
    momenta = momenta.contiguous();

    // The torch data is about to be read by kokkos to update the LAMMPS variables,
    // we need to synchronize to make sure the torch operations finished their execution
    if (mta_data->device.type() == torch::kCUDA) {
        torch::cuda::synchronize();
    }

    // Wrap torch tensors with UnmanagedView for device access
    auto positions_kk = UnmanagedView<double**, DeviceType>(
        positions.template data_ptr<double>(),
        positions.size(0), 3
    );
    auto positions_samples_kk = UnmanagedView<int32_t**, DeviceType>(
        positions_samples.data_ptr<int32_t>(),
        positions_samples.size(0),
        positions_samples.size(1)
    );

    auto momenta_kk = UnmanagedView<double**, DeviceType>(
        momenta.template data_ptr<double>(),
        momenta.size(0), 3
    );

    auto system_adaptor_kk = dynamic_cast<MetatomicSystemAdaptorKokkos<DeviceType>*>(this->system_adaptor.get());
    assert(system_adaptor_kk != nullptr);
    auto mta_to_lmp_kk = UnmanagedView<int32_t*, DeviceType>(
        system_adaptor_kk->mta_to_lmp_tensor.template data_ptr<int32_t>(),
        system_adaptor_kk->mta_to_lmp_tensor.size(0)
    );

    // Prepare masses view for device access
    // Copy masses to device if needed
    typename AT::t_kkfloat_1d masses_kk;
    if (rmass.data()) {
        masses_kk = rmass;
    } else {
        // Create a per-atom mass array from type-based masses
        masses_kk = typename AT::t_kkfloat_1d("fix_metatomic:masses", nlocal);
        Kokkos::parallel_for(nlocal,
            KOKKOS_LAMBDA(int i) { masses_kk[i] = mass[type[i]]; }
        );
    }

    // Helper reducers to avoid repetitive parallel_reduce calls
    auto reduce_mass = [&](void) -> double {
        double out = 0.0;
        Kokkos::parallel_reduce(
            nlocal,
            KOKKOS_LAMBDA(int i, double &sum) {
                if (mask[i] & groupbit) {
                    sum += masses_kk[i];
                }
            },
            out
        );
        return out;
    };

    auto reduce_pos_component = [&](int comp) -> double {
        double out = 0.0;
        Kokkos::parallel_reduce(
            nlocal,
            KOKKOS_LAMBDA(int i, double &sum) {
                if (mask[i] & groupbit) {
                    double mi = masses_kk[i];
                    sum += mi * x(i, comp);
                }
            },
            out
        );
        return out;
    };

    auto reduce_vel_component = [&](int comp) -> double {
        double out = 0.0;
        Kokkos::parallel_reduce(
            nlocal,
            KOKKOS_LAMBDA(int i, double &sum) {
                if (mask[i] & groupbit) {
                    double mi = masses_kk[i];
                    sum += mi * v(i, comp);
                }
            },
            out
        );
        return out;
    };

    // Compute total mass and mass-weighted sums before update
    double total_mass = reduce_mass();

    double com_old_x_sum = reduce_pos_component(0);
    double com_old_y_sum = reduce_pos_component(1);
    double com_old_z_sum = reduce_pos_component(2);

    double com_vel_old_x_sum = reduce_vel_component(0);
    double com_vel_old_y_sum = reduce_vel_component(1);
    double com_vel_old_z_sum = reduce_vel_component(2);

    Kokkos::parallel_for(
        positions.size(0),
        KOKKOS_LAMBDA(int i) {
            auto atom_i = mta_to_lmp_kk[positions_samples_kk(i, 1)];
            assert(atom_i < nlocal);
            if (mask[atom_i] & groupbit) {
                // Update positions with ML predictions
                x(atom_i, 0) = positions_kk(i, 0);
                x(atom_i, 1) = positions_kk(i, 1);
                x(atom_i, 2) = positions_kk(i, 2);

                // Update velocities from predicted momenta
                double mass_i = masses_kk[atom_i];
                v(atom_i, 0) = momenta_kk(i, 0) / mass_i;
                v(atom_i, 1) = momenta_kk(i, 1) / mass_i;
                v(atom_i, 2) = momenta_kk(i, 2) / mass_i;
            }
        }
    );

    // Compute mass-weighted sums after update
    double com_new_x_sum = reduce_pos_component(0);
    double com_new_y_sum = reduce_pos_component(1);
    double com_new_z_sum = reduce_pos_component(2);

    double com_vel_new_x_sum = reduce_vel_component(0);
    double com_vel_new_y_sum = reduce_vel_component(1);
    double com_vel_new_z_sum = reduce_vel_component(2);

    // Normalize to get COM positions and velocities (if mass > 0)
    std::array<double, 3> com_old;
    std::array<double, 3> com_velocity_old;
    std::array<double, 3> com_new;
    std::array<double, 3> com_velocity_new;

    if (total_mass > 0.0) {
        com_old[0] = com_old_x_sum / total_mass;
        com_old[1] = com_old_y_sum / total_mass;
        com_old[2] = com_old_z_sum / total_mass;

        com_velocity_old[0] = com_vel_old_x_sum / total_mass;
        com_velocity_old[1] = com_vel_old_y_sum / total_mass;
        com_velocity_old[2] = com_vel_old_z_sum / total_mass;

        com_new[0] = com_new_x_sum / total_mass;
        com_new[1] = com_new_y_sum / total_mass;
        com_new[2] = com_new_z_sum / total_mass;

        com_velocity_new[0] = com_vel_new_x_sum / total_mass;
        com_velocity_new[1] = com_vel_new_y_sum / total_mass;
        com_velocity_new[2] = com_vel_new_z_sum / total_mass;
    }

    // Compute adjustments to preserve COM motion:
    // pos_adjust = com_old - com_new + com_velocity_old * dt
    // vel_adjust = com_velocity_old - com_velocity_new
    double pos_adj0 = 0.0, pos_adj1 = 0.0, pos_adj2 = 0.0;
    double vel_adj0 = 0.0, vel_adj1 = 0.0, vel_adj2 = 0.0;

    if (total_mass > 0.0) {
        pos_adj0 = com_old[0] - com_new[0] + com_velocity_old[0] * update->dt;
        pos_adj1 = com_old[1] - com_new[1] + com_velocity_old[1] * update->dt;
        pos_adj2 = com_old[2] - com_new[2] + com_velocity_old[2] * update->dt;

        vel_adj0 = com_velocity_old[0] - com_velocity_new[0];
        vel_adj1 = com_velocity_old[1] - com_velocity_new[1];
        vel_adj2 = com_velocity_old[2] - com_velocity_new[2];
    }

    // Apply COM adjustments in a second pass
    Kokkos::parallel_for(
        nlocal,
        KOKKOS_LAMBDA(int i) {
            if (mask[i] & groupbit) {
                x(i, 0) += pos_adj0;
                x(i, 1) += pos_adj1;
                x(i, 2) += pos_adj2;

                v(i, 0) += vel_adj0;
                v(i, 1) += vel_adj1;
                v(i, 2) += vel_adj2;
            }
        }
    );
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMetatomicKokkos<DeviceType>::post_force(int /*vflag*/) {
    // Set the forces that comes from pair_style, bond_style, etc. to zero.
    //
    // This way we can isolate any forces that are added after this point (e.g.
    // Langevin thermostat forces) and add them during final_integrate().
    //
    // Crucially, this means that fix metatomic needs to be the first fix in the
    // post_force() sequence, i.e., the user must have it before any other fix
    // that adds forces in the input script.
    atomKK->sync(execution_space, F_MASK | MASK_MASK);
    atomKK->modified(execution_space, F_MASK);
    auto f = atomKK->k_f.template view<DeviceType>();
    auto mask = atomKK->k_mask.template view<DeviceType>();

    auto groupbit = this->groupbit;

    int nlocal = atomKK->nlocal;
    if (igroup == atomKK->firstgroup) {
        nlocal = atomKK->nfirst;
    }

    Kokkos::parallel_for(
        nlocal,
        KOKKOS_LAMBDA(int i) {
            if (mask[i] & groupbit) {
                f(i, 0) = 0.0;
                f(i, 1) = 0.0;
                f(i, 2) = 0.0;
            }
        }
    );
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMetatomicKokkos<DeviceType>::final_integrate() {
    // Apply velocity corrections from forces added after post_force.
    // This handles stochastic forces from Langevin thermostats by applying only
    // the incremental force to velocities. The first half-step is done in
    // initial_integrate(); this is the second O in an OBABO integrator.

    auto v = atomKK->k_v.template view<DeviceType>();
    auto f = atomKK->k_f.template view<DeviceType>();
    auto rmass = atomKK->k_rmass.template view<DeviceType>();
    auto mass = atomKK->k_mass.template view<DeviceType>();
    auto type = atomKK->k_type.template view<DeviceType>();
    auto mask = atomKK->k_mask.template view<DeviceType>();

    // Sync data and mark velocities as modified
    atomKK->sync(execution_space, V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK);
    atomKK->modified(execution_space, V_MASK);

    auto groupbit = this->groupbit;

    int nlocal = atomKK->nlocal;
    if (igroup == atomKK->firstgroup) {
        nlocal = atomKK->nfirst;
    }

    double dtf = 0.5 * update->dt * force->ftm2v;
    bool use_rmass = rmass.data() != nullptr;
    Kokkos::parallel_for(
        nlocal,
        KOKKOS_LAMBDA(int i) {
            if (mask[i] & groupbit) {
                double mass_i = use_rmass ? rmass[i] : mass[type[i]];
                double dtfm = dtf / mass_i;

                // Apply any force added by other fixes to velocities
                v(i, 0) += f(i, 0) * dtfm;
                v(i, 1) += f(i, 1) * dtfm;
                v(i, 2) += f(i, 2) * dtfm;
            }
        }
    );
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixMetatomicKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixMetatomicKokkos<LMPHostType>;
#endif
}
