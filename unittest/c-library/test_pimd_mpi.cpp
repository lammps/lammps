// unit tests for path-integral dynamics through the library interface

#define LAMMPS_LIB_MPI 1
#include "exceptions.h"
#include "fix.h"
#include "lammps.h"
#include "library.h"
#include "lmptype.h"
#include "modify.h"
#include "platform.h"
#include "random_mars.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#include "gtest/gtest.h"

#include "../testing/test_mpi_main.h"

using ::LAMMPS_NS::tagint;

namespace {

enum class Ownership { DEFAULT, ASYMMETRIC, MIGRATING };

struct NPTState {
    std::array<double, 36> atoms{};
    std::array<double, 6> box{};
    double barostat_velocity = 0.0;
};

struct PressureSelectorState {
    NPTState npt;
    double total_enthalpy = 0.0;
};

void *open_two_bead_partition(MPI_Comm communicator, const char *partition)
{
    const char *args[] = {"LAMMPS_test", "-screen", "none", "-log",    "none", "-partition",
                          partition,     "-in",     "none", "-nocite", nullptr};
    char **argv        = (char **)args;
    int argc           = (sizeof(args) / sizeof(char *)) - 1;
    return lammps_open(argc, argv, communicator, nullptr);
}

void *open_uneven_two_bead_partition(MPI_Comm communicator)
{
    const char *args[] = {"LAMMPS_test", "-screen", "none", "-log", "none",    "-partition",
                          "1x1",         "1x3",     "-in",  "none", "-nocite", nullptr};
    char **argv        = (char **)args;
    int argc           = (sizeof(args) / sizeof(char *)) - 1;
    return lammps_open(argc, argv, communicator, nullptr);
}

void create_two_bead_test_system(void *lmp, int world_size, Ownership ownership,
                                 const char *unit_style = "lj", int dimension = 3)
{
    if (ownership == Ownership::MIGRATING) {
        lammps_command(lmp, "variable x1 world 9.9 5.0");
        lammps_command(lmp, "variable x2 world 15.0 6.0");
        lammps_command(lmp, "variable y2 world 10.0 14.0");
    } else if (ownership == Ownership::ASYMMETRIC) {
        lammps_command(lmp, "variable x1 world 5.0 5.0");
        lammps_command(lmp, "variable x2 world 15.0 6.0");
        lammps_command(lmp, "variable y2 world 10.0 14.0");
    } else {
        lammps_command(lmp, "variable x1 world 10.0 10.0");
        lammps_command(lmp, "variable x2 world 19.0 10.0");
        lammps_command(lmp, "variable y2 world 10.0 17.0");
    }
    const std::string units_command = "units " + std::string(unit_style);
    lammps_command(lmp, units_command.c_str());
    if (dimension == 2) lammps_command(lmp, "dimension 2");
    lammps_command(lmp, "atom_style atomic");
    lammps_command(lmp, "atom_modify map array");
    lammps_command(lmp, "boundary p p p");
    if (ownership != Ownership::DEFAULT && world_size > 1) lammps_command(lmp, "processors 2 1 1");
    if (dimension == 2)
        lammps_command(lmp, "region box block 0.0 20.0 0.0 20.0 -0.5 0.5");
    else
        lammps_command(lmp, "region box block 0.0 20.0 0.0 20.0 0.0 20.0");
    lammps_command(lmp, "create_box 2 box");
    if (dimension == 2) {
        lammps_command(lmp, "create_atoms 1 single ${x1} 10.0 0.0");
        lammps_command(lmp, "create_atoms 1 single ${x2} ${y2} 0.0");
    } else {
        lammps_command(lmp, "create_atoms 1 single ${x1} 10.0 2.0");
        lammps_command(lmp, "create_atoms 1 single ${x2} ${y2} 2.0");
    }
    lammps_command(lmp, "mass * 1.0");
    lammps_command(lmp, "pair_style lj/cut 2.5");
    lammps_command(lmp, "pair_coeff * * 0.0 1.0");
    lammps_command(lmp, "velocity all set 0.0 0.0 0.0");
}

void add_two_bead_fix(void *lmp, Ownership ownership, int thermostat_seed)
{
    if (thermostat_seed > 0) {
        const std::string fix_command = "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                                        "integrator obabo thermostat PILE_L " +
                                        std::to_string(thermostat_seed) +
                                        " tau 1.0 temp 1.0 fixcom no";
        lammps_command(lmp, fix_command.c_str());
    } else if (ownership == Ownership::MIGRATING)
        lammps_command(lmp, "fix fpimd all pimd/langevin method pimd ensemble nve "
                            "integrator obabo thermostat PILE_L 1234 tau 1.0 temp 0.01 sp 100.0 "
                            "fixcom no");
    else
        lammps_command(lmp, "fix fpimd all pimd/langevin method pimd ensemble nve "
                            "integrator obabo thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");
}

void add_two_bead_integrator_fix(void *lmp, const std::string &method,
                                 const std::string &integrator,
                                 const std::string &nonfinite_trace_prefix = {})
{
    std::string command = "fix fpimd all pimd/langevin method " + method +
                          " ensemble nvt integrator " + integrator +
                          " thermostat PILE_L 2468 tau 1.0 temp 1.0 fixcom no";
    if (!nonfinite_trace_prefix.empty()) command += " nonfinite_trace " + nonfinite_trace_prefix;
    lammps_command(lmp, command.c_str());
}

void set_two_bead_initial_velocities(void *lmp, Ownership ownership)
{
    if (ownership == Ownership::MIGRATING) {
        lammps_command(lmp, "variable first_vx world 20.0 0.0");
        lammps_command(lmp, "variable first_vy world 0.0 0.0");
        lammps_command(lmp, "variable second_vx world 0.0 0.0");
        lammps_command(lmp, "variable second_vy world 0.0 0.0");
    } else {
        lammps_command(lmp, "variable first_vx world 0.10 -0.15");
        lammps_command(lmp, "variable first_vy world -0.05 0.20");
        lammps_command(lmp, "variable second_vx world -0.25 0.30");
        lammps_command(lmp, "variable second_vy world 0.15 -0.10");
    }
    lammps_command(lmp, "group first id 1");
    lammps_command(lmp, "group second id 2");
    lammps_command(lmp, "velocity first set ${first_vx} ${first_vy} 0.0");
    lammps_command(lmp, "velocity second set ${second_vx} ${second_vy} 0.0");
}

void configure_two_bead_dynamics(void *lmp, Ownership ownership, int thermostat_seed)
{
    set_two_bead_initial_velocities(lmp, ownership);
    if (ownership == Ownership::MIGRATING) {
        lammps_command(lmp, "neigh_modify every 1 delay 0 check no");
        lammps_command(lmp, "timestep 0.01");
    } else {
        lammps_command(lmp, "timestep 0.00001");
    }
    lammps_command(lmp, "thermo 0");
    add_two_bead_fix(lmp, ownership, thermostat_seed);
}

void add_two_bead_npt_fix(void *lmp, int thermostat_seed)
{
    const std::string fix_command = "fix fpimd all pimd/langevin method nmpimd ensemble npt "
                                    "integrator obabo thermostat PILE_L " +
                                    std::to_string(thermostat_seed) +
                                    " tau 1.0 temp 1.0 iso 0.5 barostat BZP taup 1.0 fixcom no";
    lammps_command(lmp, fix_command.c_str());
}

void configure_two_bead_npt(void *lmp, int thermostat_seed)
{
    set_two_bead_initial_velocities(lmp, Ownership::DEFAULT);
    lammps_command(lmp, "timestep 0.00001");
    lammps_command(lmp, "thermo 0");
    add_two_bead_npt_fix(lmp, thermostat_seed);
}

std::array<double, 18> extract_two_atom_state(void *lmp)
{
    std::array<double, 18> result{};
    auto *nlocal      = (int *)lammps_extract_global(lmp, "nlocal");
    auto *tags        = (tagint *)lammps_extract_atom(lmp, "id");
    auto **positions  = (double **)lammps_extract_atom(lmp, "x");
    auto **velocities = (double **)lammps_extract_atom(lmp, "v");
    auto **forces     = (double **)lammps_extract_atom(lmp, "f");
    EXPECT_NE(nlocal, nullptr);
    EXPECT_NE(tags, nullptr);
    EXPECT_NE(positions, nullptr);
    EXPECT_NE(velocities, nullptr);
    EXPECT_NE(forces, nullptr);
    if (nlocal && tags && positions && velocities && forces) {
        for (int i = 0; i < *nlocal; ++i) {
            EXPECT_GE(tags[i], 1);
            EXPECT_LE(tags[i], 2);
            if (tags[i] < 1 || tags[i] > 2) continue;
            const int atom = tags[i] - 1;
            for (int dimension = 0; dimension < 3; ++dimension) {
                result[3 * atom + dimension]      = positions[i][dimension];
                result[6 + 3 * atom + dimension]  = velocities[i][dimension];
                result[12 + 3 * atom + dimension] = forces[i][dimension];
            }
        }
    }
    return result;
}

std::array<double, 3> extract_nvt_global_outputs(void *lmp)
{
    std::array<double, 3> result{};
    for (int i = 0; i < 3; ++i) {
        auto *value =
            (double *)lammps_extract_fix(lmp, "fpimd", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, i, 0);
        EXPECT_NE(value, nullptr);
        if (value) {
            result[i] = *value;
            lammps_free(value);
        }
    }

    auto *lammps = (LAMMPS_NS::LAMMPS *)lmp;
    auto *fix    = lammps->modify->get_fix_by_id("fpimd");
    EXPECT_NE(fix, nullptr);
    if (fix) {
        EXPECT_EQ(fix->extvector, -1);
        EXPECT_NE(fix->extlist, nullptr);
        if (fix->extlist) {
            EXPECT_EQ(fix->extlist[0], 1);
            EXPECT_EQ(fix->extlist[1], 0);
            EXPECT_EQ(fix->extlist[2], 1);
        }
    }
    return result;
}

std::array<double, 10> extract_langevin_nvt_global_outputs(void *lmp)
{
    std::array<double, 10> result{};
    for (int i = 0; i < 10; ++i) {
        auto *value =
            (double *)lammps_extract_fix(lmp, "fpimd", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, i, 0);
        EXPECT_NE(value, nullptr);
        if (value) {
            result[i] = *value;
            lammps_free(value);
        }
    }
    return result;
}

void expect_langevin_output_metadata(const char *fix_command, const std::array<int, 17> &expected,
                                     int size)
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return;

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_two_bead_test_system(lmp, 2, Ownership::DEFAULT);
    lammps_command(lmp, fix_command);

    auto *lammps = (LAMMPS_NS::LAMMPS *)lmp;
    auto *fix    = lammps->modify->get_fix_by_id("fpimd");
    EXPECT_NE(fix, nullptr);
    if (fix) {
        EXPECT_EQ(fix->size_vector, size);
        EXPECT_EQ(fix->extvector, -1);
        EXPECT_NE(fix->extlist, nullptr);
        if (fix->extlist) {
            for (int i = 0; i < size; ++i)
                EXPECT_EQ(fix->extlist[i], expected[i]) << i;
        }
    }
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
}

std::array<double, 18> run_two_bead_case(MPI_Comm communicator, const char *partition,
                                         int world_size, Ownership ownership, int run_steps,
                                         int thermostat_seed                    = 0,
                                         std::array<double, 10> *global_outputs = nullptr)
{
    std::array<double, 18> result{};
    void *lmp = open_two_bead_partition(communicator, partition);
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return result;

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), world_size);
    create_two_bead_test_system(lmp, world_size, ownership);
    if (ownership != Ownership::DEFAULT && world_size == 2) {
        int rank;
        MPI_Comm_rank(communicator, &rank);
        const std::array<int, 4> expected_nlocal = {1, 1, 2, 0};
        auto *nlocal                             = (int *)lammps_extract_global(lmp, "nlocal");
        EXPECT_NE(nlocal, nullptr);
        if (nlocal) EXPECT_EQ(*nlocal, expected_nlocal[rank]);
    }
    configure_two_bead_dynamics(lmp, ownership, thermostat_seed);
    const std::string run_command = "run " + std::to_string(run_steps);
    lammps_command(lmp, run_command.c_str());

    const int has_error = lammps_has_error(lmp);
    EXPECT_EQ(has_error, 0);
    if (!has_error) {
        auto *nlocal = (int *)lammps_extract_global(lmp, "nlocal");
        if (ownership == Ownership::MIGRATING && world_size == 2) {
            int rank;
            MPI_Comm_rank(communicator, &rank);
            const std::array<int, 4> expected_nlocal = {0, 2, 2, 0};
            EXPECT_NE(nlocal, nullptr);
            if (nlocal) EXPECT_EQ(*nlocal, expected_nlocal[rank]);
        }
        result = extract_two_atom_state(lmp);
        if (global_outputs) *global_outputs = extract_langevin_nvt_global_outputs(lmp);
        EXPECT_EQ(lammps_has_error(lmp), 0);
    }
    lammps_close(lmp);
    return result;
}

std::array<double, 18> run_two_bead_nvt_case(MPI_Comm communicator, const char *partition,
                                             int world_size, Ownership ownership, int run_steps,
                                             std::array<double, 3> *global_outputs = nullptr)
{
    std::array<double, 18> result{};
    void *lmp = open_two_bead_partition(communicator, partition);
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return result;

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), world_size);
    create_two_bead_test_system(lmp, world_size, ownership, "real");
    if (ownership != Ownership::DEFAULT && world_size == 2) {
        int rank;
        MPI_Comm_rank(communicator, &rank);
        const std::array<int, 4> expected_nlocal = {1, 1, 2, 0};
        auto *nlocal                             = (int *)lammps_extract_global(lmp, "nlocal");
        EXPECT_NE(nlocal, nullptr);
        if (nlocal) EXPECT_EQ(*nlocal, expected_nlocal[rank]);
    }

    set_two_bead_initial_velocities(lmp, ownership);
    if (ownership == Ownership::MIGRATING) {
        lammps_command(lmp, "neigh_modify every 1 delay 0 check no");
        lammps_command(lmp, "timestep 0.01");
        lammps_command(lmp,
                       "fix fpimd all pimd/nvt method pimd temp 0.01 fmass 1.0 sp 100.0 nhc 2");
    } else {
        lammps_command(lmp, "timestep 0.00001");
        lammps_command(lmp, "fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2");
    }
    lammps_command(lmp, "thermo 0");
    const std::string run_command = "run " + std::to_string(run_steps);
    lammps_command(lmp, run_command.c_str());

    const int has_error = lammps_has_error(lmp);
    EXPECT_EQ(has_error, 0);
    if (!has_error) {
        auto *nlocal = (int *)lammps_extract_global(lmp, "nlocal");
        if (ownership == Ownership::MIGRATING && world_size == 2) {
            int rank;
            MPI_Comm_rank(communicator, &rank);
            const std::array<int, 4> expected_nlocal = {0, 2, 2, 0};
            EXPECT_NE(nlocal, nullptr);
            if (nlocal) EXPECT_EQ(*nlocal, expected_nlocal[rank]);
        }
        result = extract_two_atom_state(lmp);
        if (global_outputs) *global_outputs = extract_nvt_global_outputs(lmp);
        EXPECT_EQ(lammps_has_error(lmp), 0);
    }
    lammps_close(lmp);
    return result;
}

std::array<double, 36> collect_two_bead_state(const std::array<double, 18> &local_state)
{
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    std::array<double, 36> state{};
    const int offset = 18 * (me / 2);
    for (std::size_t i = 0; i < local_state.size(); ++i)
        state[offset + i] = local_state[i];
    MPI_Allreduce(MPI_IN_PLACE, state.data(), static_cast<int>(state.size()), MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    return state;
}

std::array<double, 36>
run_two_bead_integrator_segments(int first_steps, int second_steps, bool restart = false,
                                 std::array<double, 10> *global_outputs = nullptr,
                                 const std::string &method              = "pimd",
                                 const std::string &integrator = "obabo", bool interacting = false,
                                 const std::string &nonfinite_trace_prefix = {})
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_two_bead_test_system(lmp, 2, Ownership::DEFAULT);
    if (interacting) {
        // Both bead geometries have nonzero LJ forces (separations 7 and 9).
        lammps_command(lmp, "pair_style lj/cut 9.5");
        lammps_command(lmp, "pair_coeff * * 0.1 6.0");
    }
    set_two_bead_initial_velocities(lmp, Ownership::DEFAULT);
    lammps_command(lmp, "timestep 0.00001");
    lammps_command(lmp, "thermo 0");
    add_two_bead_integrator_fix(lmp, method, integrator, nonfinite_trace_prefix);
    const std::string first_run = "run " + std::to_string(first_steps);
    lammps_command(lmp, first_run.c_str());
    if (second_steps > 0) {
        if (restart) {
            lammps_command(lmp, "variable restart_file world pimd_nvt_restart.0 "
                                "pimd_nvt_restart.1");
            lammps_command(lmp, "write_restart ${restart_file}");
            lammps_close(lmp);
            lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
            EXPECT_NE(lmp, nullptr);
            if (!lmp) return {};
            lammps_command(lmp, "variable restart_file world pimd_nvt_restart.0 "
                                "pimd_nvt_restart.1");
            lammps_command(lmp, "read_restart ${restart_file}");
            lammps_command(lmp, "thermo 0");
            add_two_bead_integrator_fix(lmp, method, integrator, nonfinite_trace_prefix);
        }
        const std::string second_run = "run " + std::to_string(second_steps);
        lammps_command(lmp, second_run.c_str());
    }

    const auto state = collect_two_bead_state(extract_two_atom_state(lmp));
    if (global_outputs) *global_outputs = extract_langevin_nvt_global_outputs(lmp);
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

std::array<double, 36> run_two_bead_damping_case(double tau, double scale, int seed,
                                                 int &floating_point_exceptions)
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_two_bead_test_system(lmp, 2, Ownership::DEFAULT);
    set_two_bead_initial_velocities(lmp, Ownership::DEFAULT);
    lammps_command(lmp, "timestep 0.00001");
    lammps_command(lmp, "thermo 0");
    const std::string fix_command =
        "fix fpimd all pimd/langevin method nmpimd ensemble nvt integrator obabo "
        "thermostat PILE_L " +
        std::to_string(seed) + " tau " + std::to_string(tau) + " scale " + std::to_string(scale) +
        " temp 1.0 fixcom no";
    lammps_command(lmp, fix_command.c_str());

    std::feclearexcept(FE_ALL_EXCEPT);
    lammps_command(lmp, "run 1");
    floating_point_exceptions = std::fetestexcept(FE_ALL_EXCEPT);

    const auto state = collect_two_bead_state(extract_two_atom_state(lmp));
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

NPTState extract_npt_state(void *lmp)
{
    NPTState state;
    state.atoms = collect_two_bead_state(extract_two_atom_state(lmp));

    double boxlo[3], boxhi[3], xy, yz, xz;
    int periodicity[3], boxflag;
    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, periodicity, &boxflag);
    for (int i = 0; i < 3; ++i) {
        state.box[i]     = boxlo[i];
        state.box[3 + i] = boxhi[i];
    }

    auto *barostat_velocity =
        (double *)lammps_extract_fix(lmp, "fpimd", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, 10, 0);
    EXPECT_NE(barostat_velocity, nullptr);
    if (barostat_velocity) {
        state.barostat_velocity = *barostat_velocity;
        lammps_free(barostat_velocity);
    }
    return state;
}

void remove_nvt_restart_files()
{
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (me == 0) {
        LAMMPS_NS::platform::unlink("pimd_nvt_restart.0");
        LAMMPS_NS::platform::unlink("pimd_nvt_restart.1");
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

std::array<double, 36> write_two_bead_nvt_restart()
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_two_bead_test_system(lmp, 2, Ownership::DEFAULT);
    configure_two_bead_dynamics(lmp, Ownership::DEFAULT, 2468);
    lammps_command(lmp, "run 3");
    const auto state = collect_two_bead_state(extract_two_atom_state(lmp));
    lammps_command(lmp, "variable restart_file world pimd_nvt_restart.0 "
                        "pimd_nvt_restart.1");
    lammps_command(lmp, "write_restart ${restart_file}");
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

std::array<double, 36> run_two_bead_nvt_restart(std::array<double, 36> &restored_state,
                                                std::array<double, 10> *global_outputs = nullptr)
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    lammps_command(lmp, "variable restart_file world pimd_nvt_restart.0 "
                        "pimd_nvt_restart.1");
    lammps_command(lmp, "read_restart ${restart_file}");
    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    lammps_command(lmp, "thermo 0");
    add_two_bead_fix(lmp, Ownership::DEFAULT, 2468);
    restored_state = collect_two_bead_state(extract_two_atom_state(lmp));
    lammps_command(lmp, "run 3");
    const auto state = collect_two_bead_state(extract_two_atom_state(lmp));
    if (global_outputs) *global_outputs = extract_langevin_nvt_global_outputs(lmp);
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

void remove_npt_restart_files()
{
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (me == 0) {
        LAMMPS_NS::platform::unlink("pimd_npt_restart.0");
        LAMMPS_NS::platform::unlink("pimd_npt_restart.1");
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

NPTState run_two_bead_npt_case(int run_steps)
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_two_bead_test_system(lmp, 2, Ownership::DEFAULT);
    configure_two_bead_npt(lmp, 8642);
    const std::string run_command = "run " + std::to_string(run_steps);
    lammps_command(lmp, run_command.c_str());
    const auto state = extract_npt_state(lmp);
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

NPTState write_two_bead_npt_restart()
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_two_bead_test_system(lmp, 2, Ownership::DEFAULT);
    configure_two_bead_npt(lmp, 8642);
    lammps_command(lmp, "run 3");
    const auto state = extract_npt_state(lmp);
    lammps_command(lmp, "variable restart_file world pimd_npt_restart.0 "
                        "pimd_npt_restart.1");
    lammps_command(lmp, "write_restart ${restart_file}");
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

NPTState run_two_bead_npt_restart(NPTState &restored_state)
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    lammps_command(lmp, "variable restart_file world pimd_npt_restart.0 "
                        "pimd_npt_restart.1");
    lammps_command(lmp, "read_restart ${restart_file}");
    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    lammps_command(lmp, "thermo 0");
    add_two_bead_npt_fix(lmp, 8642);
    restored_state = extract_npt_state(lmp);
    lammps_command(lmp, "run 3");
    const auto state = extract_npt_state(lmp);
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

std::array<double, 12> run_two_bead_anisotropic_nph(const std::array<double, 3> &pressure_targets)
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_two_bead_test_system(lmp, 2, Ownership::DEFAULT);
    set_two_bead_initial_velocities(lmp, Ownership::DEFAULT);
    lammps_command(lmp, "timestep 0.00001");
    lammps_command(lmp, "thermo 0");
    lammps_command(lmp, "fix drive all addforce 0.2 -0.3 0.4");
    const std::string fix_command =
        "fix fpimd all pimd/langevin method nmpimd ensemble nph integrator obabo "
        "thermostat PILE_L 8642 tau 1.0 temp 1.0 x " +
        std::to_string(pressure_targets[0]) + " y " + std::to_string(pressure_targets[1]) + " z " +
        std::to_string(pressure_targets[2]) + " barostat BZP taup 1.0 fixcom no";
    lammps_command(lmp, fix_command.c_str());
    lammps_command(lmp, "run 1");

    std::array<double, 3> local{};
    for (int dimension = 0; dimension < 3; ++dimension) {
        auto *value = (double *)lammps_extract_fix(lmp, "fpimd", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR,
                                                   10 + dimension, 0);
        EXPECT_NE(value, nullptr);
        if (value) {
            local[dimension] = *value;
            lammps_free(value);
        }
    }

    std::array<double, 12> gathered{};
    MPI_Allgather(local.data(), static_cast<int>(local.size()), MPI_DOUBLE, gathered.data(),
                  static_cast<int>(local.size()), MPI_DOUBLE, MPI_COMM_WORLD);
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return gathered;
}

PressureSelectorState run_two_bead_pressure_selector(const char *pressure_keywords)
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_two_bead_test_system(lmp, 2, Ownership::DEFAULT);
    set_two_bead_initial_velocities(lmp, Ownership::DEFAULT);
    lammps_command(lmp, "timestep 0.00001");
    lammps_command(lmp, "thermo 0");
    const std::string fix_command =
        "fix fpimd all pimd/langevin method nmpimd ensemble nph integrator obabo "
        "thermostat PILE_L 8642 tau 1.0 temp 1.0 " +
        std::string(pressure_keywords) + " barostat BZP taup 1.0 fixcom no";
    lammps_command(lmp, fix_command.c_str());
    lammps_command(lmp, "run 1");

    PressureSelectorState state;
    state.npt = extract_npt_state(lmp);
    auto *total_enthalpy =
        (double *)lammps_extract_fix(lmp, "fpimd", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, 16, 0);
    EXPECT_NE(total_enthalpy, nullptr);
    if (total_enthalpy) {
        state.total_enthalpy = *total_enthalpy;
        lammps_free(total_enthalpy);
    }
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

NPTState run_two_bead_baoab_bzp()
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_two_bead_test_system(lmp, 2, Ownership::DEFAULT);
    set_two_bead_initial_velocities(lmp, Ownership::DEFAULT);
    lammps_command(lmp, "timestep 0.00001");
    lammps_command(lmp, "thermo 0");
    const std::string fix_command =
        "fix fpimd all pimd/langevin method nmpimd ensemble nph integrator baoab "
        "thermostat PILE_L 8642 tau 1.0 temp 1.0 iso 0.5 barostat BZP taup 1.0 fixcom no";
    lammps_command(lmp, fix_command.c_str());
    lammps_command(lmp, "run 1");

    const auto state = extract_npt_state(lmp);
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

NPTState run_two_bead_zero_barostat_velocity(const char *pressure_keyword)
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_two_bead_test_system(lmp, 2, Ownership::DEFAULT);
    // The two-bead matched pressure is exactly 1/4096 for this zero-force, zero-velocity box.
    lammps_command(lmp, "change_box all x final 0.0 16.0 y final 0.0 16.0 z final 0.0 16.0 remap");
    lammps_command(lmp, "timestep 0.00001");
    lammps_command(lmp, "thermo 0");
    const std::string fix_command =
        "fix fpimd all pimd/langevin method nmpimd ensemble nph integrator baoab "
        "thermostat PILE_L 8642 tau 1.0 temp 1.0 " +
        std::string(pressure_keyword) + " 0.000244140625 barostat BZP taup 1.0 fixcom no";
    lammps_command(lmp, fix_command.c_str());
    lammps_command(lmp, "run 1");

    NPTState state{};
    const int has_error = lammps_has_error(lmp);
    EXPECT_EQ(has_error, 0);
    if (!has_error) state = extract_npt_state(lmp);
    lammps_close(lmp);
    return state;
}

void expect_pimd_fix_error_on_partition(const std::string &fix_command, const char *message,
                                        const char *unit_style, MPI_Comm communicator,
                                        const char *partition, int expected_world_size,
                                        int dimension = 3, const char *setup_command = nullptr)
{
    void *lmp = open_two_bead_partition(communicator, partition);
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return;

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), expected_world_size);
    create_two_bead_test_system(lmp, expected_world_size, Ownership::DEFAULT, unit_style,
                                dimension);
    if (setup_command) lammps_command(lmp, setup_command);
    lammps_set_show_error(lmp, 0);
    lammps_command(lmp, fix_command.c_str());

    const int has_error = lammps_has_error(lmp);
    EXPECT_EQ(has_error, 1);
    if (has_error) {
        char error_message[512];
        EXPECT_NE(lammps_get_last_error_message(lmp, error_message, sizeof(error_message)), 0);
        EXPECT_NE(std::string(error_message).find(message), std::string::npos) << error_message;
    }
    lammps_close(lmp);
}

void expect_pimd_fix_error(const std::string &fix_command, const char *message,
                           const char *unit_style)
{
    expect_pimd_fix_error_on_partition(fix_command, message, unit_style, MPI_COMM_WORLD, "2x2", 2);
}

void expect_pimd_fix_error_on_uneven_partition(const std::string &fix_command, const char *message,
                                               const char *unit_style)
{
    void *lmp = open_uneven_two_bead_partition(MPI_COMM_WORLD);
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return;

    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    const int world_size = lammps_extract_setting(lmp, "world_size");
    EXPECT_EQ(world_size, me == 0 ? 1 : 3);
    create_two_bead_test_system(lmp, world_size, Ownership::DEFAULT, unit_style);
    lammps_set_show_error(lmp, 0);
    lammps_command(lmp, fix_command.c_str());

    const int has_error = lammps_has_error(lmp);
    EXPECT_EQ(has_error, 1);
    if (has_error) {
        char error_message[512];
        EXPECT_NE(lammps_get_last_error_message(lmp, error_message, sizeof(error_message)), 0);
        EXPECT_NE(std::string(error_message).find(message), std::string::npos) << error_message;
    }
    lammps_close(lmp);
}

void expect_pimd_run_error_on_mismatched_atom_count(const std::string &fix_command,
                                                    const char *message, const char *unit_style,
                                                    MPI_Comm communicator)
{
    void *lmp = open_two_bead_partition(communicator, "2x1");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return;

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 1);
    create_two_bead_test_system(lmp, 1, Ownership::DEFAULT, unit_style);

    int me;
    MPI_Comm_rank(communicator, &me);
    if (me == 1) {
        lammps_command(lmp, "group removed id 2");
        lammps_command(lmp, "delete_atoms group removed");
    }

    lammps_set_show_error(lmp, 0);
    lammps_command(lmp, fix_command.c_str());
    if (!lammps_has_error(lmp)) lammps_command(lmp, "run 0");

    const int has_error = lammps_has_error(lmp);
    EXPECT_EQ(has_error, 1);
    if (has_error) {
        char error_message[512];
        EXPECT_NE(lammps_get_last_error_message(lmp, error_message, sizeof(error_message)), 0);
        EXPECT_NE(std::string(error_message).find(message), std::string::npos) << error_message;
    }
    lammps_close(lmp);
}

void expect_pimd_run_success_on_partition(const std::string &fix_command, const char *unit_style,
                                          MPI_Comm communicator,
                                          const char *setup_command = nullptr)
{
    void *lmp = open_two_bead_partition(communicator, "2x1");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return;

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 1);
    create_two_bead_test_system(lmp, 1, Ownership::DEFAULT, unit_style);
    if (setup_command) lammps_command(lmp, setup_command);
    lammps_set_show_error(lmp, 0);
    lammps_command(lmp, fix_command.c_str());
    if (!lammps_has_error(lmp)) lammps_command(lmp, "run 0");
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
}

void expect_pimd_run_error_after_partition_change(const std::string &fix_command,
                                                  const char *message, const char *unit_style,
                                                  const char *change_command, MPI_Comm communicator,
                                                  const char *setup_command = nullptr)
{
    void *lmp = open_two_bead_partition(communicator, "2x1");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return;

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 1);
    create_two_bead_test_system(lmp, 1, Ownership::DEFAULT, unit_style);
    if (setup_command) lammps_command(lmp, setup_command);

    int me;
    MPI_Comm_rank(communicator, &me);
    if (me == 1) lammps_command(lmp, change_command);

    lammps_set_show_error(lmp, 0);
    lammps_command(lmp, fix_command.c_str());
    if (!lammps_has_error(lmp)) lammps_command(lmp, "run 0");

    const int has_error = lammps_has_error(lmp);
    EXPECT_EQ(has_error, 1);
    if (has_error) {
        char error_message[512];
        EXPECT_NE(lammps_get_last_error_message(lmp, error_message, sizeof(error_message)), 0);
        EXPECT_NE(std::string(error_message).find(message), std::string::npos) << error_message;
    }
    lammps_close(lmp);
}

void expect_pimd_fix_success(const std::string &fix_command, const char *unit_style)
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return;

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_two_bead_test_system(lmp, 2, Ownership::DEFAULT, unit_style);
    lammps_set_show_error(lmp, 0);
    lammps_command(lmp, fix_command.c_str());
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
}

void check_decomposition_equivalence(Ownership ownership, int run_steps)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm reference_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me % 2, me, &reference_comm);
    const auto reference = run_two_bead_case(reference_comm, "2x1", 1, ownership, run_steps);
    MPI_Comm_free(&reference_comm);

    std::array<double, 36> reference_state{};
    if (me % 2 == 0) {
        const int offset = 18 * (me / 2);
        for (std::size_t i = 0; i < reference.size(); ++i)
            reference_state[offset + i] = reference[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, reference_state.data(), static_cast<int>(reference_state.size()),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    const auto parallel       = run_two_bead_case(MPI_COMM_WORLD, "2x2", 2, ownership, run_steps);
    const auto parallel_state = collect_two_bead_state(parallel);

    double initial_x_bead_zero = 10.0;
    double initial_x_bead_one  = 10.0;
    if (ownership == Ownership::ASYMMETRIC)
        initial_x_bead_zero = initial_x_bead_one = 5.0;
    else if (ownership == Ownership::MIGRATING) {
        initial_x_bead_zero = 9.9;
        initial_x_bead_one  = 5.0;
    }
    EXPECT_NE(reference_state[0], initial_x_bead_zero);
    EXPECT_NE(reference_state[15], 0.0);
    EXPECT_NE(reference_state[18], initial_x_bead_one);
    EXPECT_NE(reference_state[33], 0.0);
    for (std::size_t i = 0; i < reference_state.size(); ++i)
        EXPECT_NEAR(parallel_state[i], reference_state[i], 1.0e-11);
}

void check_nvt_decomposition_equivalence(Ownership ownership, int run_steps)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm reference_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me % 2, me, &reference_comm);
    const auto reference = run_two_bead_nvt_case(reference_comm, "2x1", 1, ownership, run_steps);
    MPI_Comm_free(&reference_comm);

    std::array<double, 36> reference_state{};
    if (me % 2 == 0) {
        const int offset = 18 * (me / 2);
        for (std::size_t i = 0; i < reference.size(); ++i)
            reference_state[offset + i] = reference[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, reference_state.data(), static_cast<int>(reference_state.size()),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    const auto parallel = run_two_bead_nvt_case(MPI_COMM_WORLD, "2x2", 2, ownership, run_steps);
    const auto parallel_state = collect_two_bead_state(parallel);

    for (std::size_t i = 0; i < reference_state.size(); ++i) {
        EXPECT_TRUE(std::isfinite(reference_state[i]));
        EXPECT_TRUE(std::isfinite(parallel_state[i]));
        EXPECT_NEAR(parallel_state[i], reference_state[i], 1.0e-11);
    }
}

void check_nvt_global_output_equivalence(Ownership ownership)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm reference_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me % 2, me, &reference_comm);
    std::array<double, 3> reference_local{};
    run_two_bead_nvt_case(reference_comm, "2x1", 1, ownership, 1, &reference_local);
    MPI_Comm_free(&reference_comm);

    std::array<double, 6> reference_by_bead{};
    if (me % 2 == 0) {
        const int offset = 3 * (me / 2);
        for (std::size_t i = 0; i < reference_local.size(); ++i)
            reference_by_bead[offset + i] = reference_local[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, reference_by_bead.data(),
                  static_cast<int>(reference_by_bead.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < 3; ++i) {
        EXPECT_TRUE(std::isfinite(reference_by_bead[i]));
        EXPECT_TRUE(std::isfinite(reference_by_bead[3 + i]));
        EXPECT_NEAR(reference_by_bead[3 + i], reference_by_bead[i], 1.0e-11);
    }

    std::array<double, 3> parallel_outputs{};
    run_two_bead_nvt_case(MPI_COMM_WORLD, "2x2", 2, ownership, 1, &parallel_outputs);
    for (int i = 0; i < 3; ++i) {
        EXPECT_TRUE(std::isfinite(parallel_outputs[i]));
        EXPECT_NEAR(parallel_outputs[i], reference_by_bead[i], 1.0e-11);
    }
}

void expect_nonfinite_trace(const std::string &prefix, bool inject_position_nan)
{
    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    ASSERT_NE(lmp, nullptr);
    create_two_bead_test_system(lmp, 2, Ownership::DEFAULT);
    set_two_bead_initial_velocities(lmp, Ownership::DEFAULT);
    lammps_command(lmp, "timestep 0.00001");
    lammps_command(lmp, "thermo 0");
    const std::string fix_command =
        "fix fpimd all pimd/langevin method nmpimd ensemble nvt integrator baoab "
        "thermostat PILE_L 2468 tau 1.0 temp 1.0 fixcom no nonfinite_trace " +
        prefix;
    lammps_command(lmp, fix_command.c_str());
    lammps_command(lmp, "run 0 post no");
    ASSERT_EQ(lammps_has_error(lmp), 0);

    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    auto *nlocal = (int *)lammps_extract_global(lmp, "nlocal");
    auto *tags   = (tagint *)lammps_extract_atom(lmp, "id");
    ASSERT_NE(nlocal, nullptr);
    ASSERT_NE(tags, nullptr);

    tagint local_tag       = std::numeric_limits<tagint>::max();
    int local_index        = -1;
    const int target_world = inject_position_nan ? 1 : 0;
    if (me / 2 == target_world) {
        for (int i = 0; i < *nlocal; ++i) {
            if (tags[i] < local_tag) {
                local_tag   = tags[i];
                local_index = i;
            }
        }
    }
    std::array<tagint, 4> candidate_tags{};
    MPI_Allgather(&local_tag, 1, MPI_LMP_TAGINT, candidate_tags.data(), 1, MPI_LMP_TAGINT,
                  MPI_COMM_WORLD);
    int owner = 2 * target_world;
    for (int rank = 2 * target_world; rank < 2 * target_world + 2; ++rank)
        if (candidate_tags[rank] < candidate_tags[owner]) owner = rank;
    const tagint expected_tag = candidate_tags[owner];
    ASSERT_NE(expected_tag, std::numeric_limits<tagint>::max());

    if (me == owner) {
        ASSERT_GE(local_index, 0);
        if (inject_position_nan) {
            auto **positions = (double **)lammps_extract_atom(lmp, "x");
            ASSERT_NE(positions, nullptr);
            positions[local_index][1] = std::numeric_limits<double>::quiet_NaN();
        } else {
            auto **velocities = (double **)lammps_extract_atom(lmp, "v");
            ASSERT_NE(velocities, nullptr);
            velocities[local_index][2] = std::numeric_limits<double>::infinity();
        }
    }

    lammps_set_show_error(lmp, 0);
    lammps_command(lmp, "run 1 pre no post no");
    ASSERT_EQ(lammps_has_error(lmp), 1);
    char error_message[1024];
    ASSERT_NE(lammps_get_last_error_message(lmp, error_message, sizeof(error_message)), 0);
    EXPECT_NE(std::string(error_message).find("nonfinite trace detected"), std::string::npos);
    EXPECT_NE(std::string(error_message).find("stage initial-entry"), std::string::npos);
    lammps_close(lmp);

    const int world_rank   = owner - 2 * target_world;
    const std::string path = prefix + ".pimd.step1.u" + std::to_string(owner) + ".w" +
                             std::to_string(target_world) + ".r" + std::to_string(world_rank) +
                             ".txt";
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        std::ifstream input(path);
        ASSERT_TRUE(input.good()) << path;
        std::ostringstream buffer;
        buffer << input.rdbuf();
        const std::string report = buffer.str();
        EXPECT_NE(report.find("schema=pimd-nonfinite-trace-v1"), std::string::npos);
        EXPECT_NE(report.find("stage=initial-entry"), std::string::npos);
        EXPECT_NE(report.find("atom_tag=" + std::to_string(expected_tag)), std::string::npos);
        EXPECT_NE(report.find(inject_position_nan ? "field=position" : "field=velocity"),
                  std::string::npos);
        EXPECT_NE(report.find(inject_position_nan ? "component=1" : "component=2"),
                  std::string::npos);
        EXPECT_NE(report.find("previous_available=1"), std::string::npos);
        LAMMPS_NS::platform::unlink(path);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace

TEST(PIMD, multirank_dynamics_equivalence)
{
    check_decomposition_equivalence(Ownership::DEFAULT, 5);
}

TEST(PIMD, multirank_asymmetric_ownership_equivalence)
{
    check_decomposition_equivalence(Ownership::ASYMMETRIC, 1);
}

TEST(PIMD, multirank_migrating_ownership_equivalence)
{
    check_decomposition_equivalence(Ownership::MIGRATING, 1);
}

TEST(PIMD, multirank_nvt_asymmetric_ownership_equivalence)
{
    check_nvt_decomposition_equivalence(Ownership::ASYMMETRIC, 1);
}

TEST(PIMD, multirank_nvt_migrating_ownership_equivalence)
{
    check_nvt_decomposition_equivalence(Ownership::MIGRATING, 1);
}

TEST(PIMD, multirank_nvt_asymmetric_global_output_equivalence)
{
    check_nvt_global_output_equivalence(Ownership::ASYMMETRIC);
}

TEST(PIMD, multirank_nvt_migrating_global_output_equivalence)
{
    check_nvt_global_output_equivalence(Ownership::MIGRATING);
}

TEST(PIMD, langevin_output_metadata)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    expect_langevin_output_metadata(
        "fix fpimd all pimd/langevin method nmpimd ensemble nvt thermostat PILE_L 8642",
        {1, 1, 1, 1, 1, 1, 1, 0, 0, 0}, 10);
    expect_langevin_output_metadata(
        "fix fpimd all pimd/langevin method nmpimd ensemble nph iso 0.5 barostat BZP taup 1.0",
        {1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1}, 15);
    expect_langevin_output_metadata(
        "fix fpimd all pimd/langevin method nmpimd ensemble nph x 0.5 y 0.5 z 0.5 "
        "barostat BZP taup 1.0",
        {1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}, 17);
}

TEST(PIMD, multirank_nvt_seed_replay)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    auto run_nvt = [](int seed) {
        const auto local = run_two_bead_case(MPI_COMM_WORLD, "2x2", 2, Ownership::DEFAULT, 3, seed);
        return collect_two_bead_state(local);
    };

    const auto first          = run_nvt(2468);
    const auto replay         = run_nvt(2468);
    const auto different_seed = run_nvt(9753);

    for (std::size_t i = 0; i < first.size(); ++i) {
        EXPECT_TRUE(std::isfinite(first[i]));
        EXPECT_TRUE(std::isfinite(replay[i]));
        EXPECT_TRUE(std::isfinite(different_seed[i]));
        EXPECT_NEAR(replay[i], first[i], 1.0e-14);
    }

    double max_velocity_difference = 0.0;
    for (int bead = 0; bead < 2; ++bead) {
        const int offset = 18 * bead + 6;
        for (int i = 0; i < 6; ++i)
            max_velocity_difference = std::max(
                max_velocity_difference, std::fabs(different_seed[offset + i] - first[offset + i]));
    }
    EXPECT_GT(max_velocity_difference, 1.0e-6);
}

TEST(PIMD, multirank_direct_pimd_run_and_restart_continuity)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    remove_nvt_restart_files();
    std::array<double, 10> continuous_outputs{}, segmented_outputs{}, restarted_outputs{};
    const auto continuous = run_two_bead_integrator_segments(6, 0, false, &continuous_outputs);
    const auto segmented  = run_two_bead_integrator_segments(3, 3, false, &segmented_outputs);
    const auto restarted  = run_two_bead_integrator_segments(3, 3, true, &restarted_outputs);
    remove_nvt_restart_files();
    for (std::size_t i = 0; i < continuous.size(); ++i) {
        EXPECT_NEAR(segmented[i], continuous[i], 1.0e-14);
        EXPECT_NEAR(restarted[i], continuous[i], 1.0e-14);
    }
    for (std::size_t i = 0; i < continuous_outputs.size(); ++i) {
        EXPECT_TRUE(std::isfinite(continuous_outputs[i]));
        EXPECT_NEAR(segmented_outputs[i], continuous_outputs[i], 1.0e-14) << i;
        EXPECT_NEAR(restarted_outputs[i], continuous_outputs[i], 1.0e-14) << i;
    }
}

TEST(PIMD, multirank_interacting_integrator_restart_seams)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    for (const auto &method : {"pimd", "nmpimd"}) {
        for (const auto &integrator : {"obabo", "baoab"}) {
            SCOPED_TRACE(std::string(method) + "/" + integrator);
            std::array<double, 10> reference_outputs{};
            const auto reference = run_two_bead_integrator_segments(8, 0, false, &reference_outputs,
                                                                    method, integrator, true);
            EXPECT_GT(std::abs(reference_outputs[2]), 1.0e-6);
            for (int seam : {1, 3, 4}) {
                SCOPED_TRACE(seam);
                for (bool restart : {false, true}) {
                    SCOPED_TRACE(restart);
                    remove_nvt_restart_files();
                    std::array<double, 10> outputs{};
                    const auto state = run_two_bead_integrator_segments(
                        seam, 8 - seam, restart, &outputs, method, integrator, true);
                    remove_nvt_restart_files();
                    for (std::size_t i = 0; i < state.size(); ++i) {
                        EXPECT_TRUE(std::isfinite(state[i]));
                        EXPECT_NEAR(state[i], reference[i], 1.0e-14) << i;
                    }
                    for (std::size_t i = 0; i < outputs.size(); ++i) {
                        EXPECT_TRUE(std::isfinite(outputs[i]));
                        const double tolerance =
                            1.0e-14 *
                            std::max({1.0, std::abs(outputs[i]), std::abs(reference_outputs[i])});
                        EXPECT_NEAR(outputs[i], reference_outputs[i], tolerance) << i;
                    }
                }
            }
        }
    }
}

TEST(PIMD, multirank_nonfinite_trace_finite_parity)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    const auto baseline =
        run_two_bead_integrator_segments(4, 0, false, nullptr, "nmpimd", "baoab", true);
    const auto traced = run_two_bead_integrator_segments(4, 0, false, nullptr, "nmpimd", "baoab",
                                                         true, "test_pimd_finite_trace");
    for (std::size_t i = 0; i < baseline.size(); ++i)
        EXPECT_DOUBLE_EQ(traced[i], baseline[i]) << i;
}

TEST(PIMD, multirank_nonfinite_trace_reports_nan_position)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);
    expect_nonfinite_trace("test_pimd_trace_nan_position", true);
}

TEST(PIMD, multirank_nonfinite_trace_reports_infinite_velocity)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);
    expect_nonfinite_trace("test_pimd_trace_infinite_velocity", false);
}

TEST(PIMD, multirank_restart_record_sizes)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);
    void *handle = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    ASSERT_NE(handle, nullptr);
    create_two_bead_test_system(handle, 2, Ownership::DEFAULT);
    configure_two_bead_dynamics(handle, Ownership::DEFAULT, 2468);
    auto *lmp = static_cast<LAMMPS_NS::LAMMPS *>(handle);
    auto *fix = lmp->modify->get_fix_by_id("fpimd");
    ASSERT_NE(fix, nullptr);

    // The public six-double legacy record must not be read past its end.
    double legacy[6]{};
    EXPECT_NO_THROW(fix->restart(reinterpret_cast<char *>(legacy), sizeof(legacy)));
    EXPECT_THROW(fix->restart(reinterpret_cast<char *>(legacy), sizeof(legacy) - 1),
                 LAMMPS_NS::LAMMPSException);
    EXPECT_THROW(fix->restart(reinterpret_cast<char *>(legacy), 0), LAMMPS_NS::LAMMPSException);

    std::array<double, 7 + 2 * LAMMPS_NS::RanMars::STATE_SIZE> state{};
    state[6] = 2;
    {
        LAMMPS_NS::RanMars random(lmp, 2468);
        random.gaussian();
        random.get_state(state.data() + 7);
        random.get_state(state.data() + 7 + LAMMPS_NS::RanMars::STATE_SIZE);
    }
    EXPECT_NO_THROW(fix->restart(reinterpret_cast<char *>(state.data()), sizeof(state)));
    // The maintained development fork used a zero marker with a full cache.
    state[7] = state[7 + LAMMPS_NS::RanMars::STATE_SIZE] = 0.0;
    EXPECT_NO_THROW(fix->restart(reinterpret_cast<char *>(state.data()), sizeof(state)));
    state[6] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(fix->restart(reinterpret_cast<char *>(state.data()), sizeof(state)),
                 LAMMPS_NS::LAMMPSException);
    state[6] = 1.5;
    EXPECT_THROW(fix->restart(reinterpret_cast<char *>(state.data()), sizeof(state)),
                 LAMMPS_NS::LAMMPSException);
    state[6] = 2;
    EXPECT_THROW(
        fix->restart(reinterpret_cast<char *>(state.data()), sizeof(state) - sizeof(double)),
        LAMMPS_NS::LAMMPSException);
    lammps_close(handle);
}

TEST(PIMD, multirank_nvt_restart_restores_rng)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    remove_nvt_restart_files();
    std::array<double, 10> continuous_outputs{}, first_outputs{}, replay_outputs{};
    const auto continuous = collect_two_bead_state(run_two_bead_case(
        MPI_COMM_WORLD, "2x2", 2, Ownership::DEFAULT, 6, 2468, &continuous_outputs));
    const auto boundary   = write_two_bead_nvt_restart();
    std::array<double, 36> first_restored{};
    const auto first = run_two_bead_nvt_restart(first_restored, &first_outputs);
    std::array<double, 36> replay_restored{};
    const auto replay = run_two_bead_nvt_restart(replay_restored, &replay_outputs);
    remove_nvt_restart_files();

    for (int bead = 0; bead < 2; ++bead) {
        const int offset = 18 * bead;
        for (int i = 0; i < 12; ++i) {
            EXPECT_NEAR(first_restored[offset + i], boundary[offset + i], 1.0e-14);
            EXPECT_NEAR(replay_restored[offset + i], boundary[offset + i], 1.0e-14);
        }
    }

    for (std::size_t i = 0; i < first.size(); ++i) {
        EXPECT_TRUE(std::isfinite(continuous[i]));
        EXPECT_TRUE(std::isfinite(first[i]));
        EXPECT_TRUE(std::isfinite(replay[i]));
        EXPECT_NEAR(replay[i], first[i], 1.0e-14);
    }

    for (std::size_t i = 0; i < first.size(); ++i)
        EXPECT_NEAR(first[i], continuous[i], 1.0e-14);
    for (std::size_t i = 0; i < continuous_outputs.size(); ++i) {
        EXPECT_TRUE(std::isfinite(continuous_outputs[i]));
        EXPECT_DOUBLE_EQ(replay_outputs[i], first_outputs[i]) << i;
        const double tolerance =
            1.0e-14 * std::max({1.0, std::abs(first_outputs[i]), std::abs(continuous_outputs[i])});
        EXPECT_NEAR(first_outputs[i], continuous_outputs[i], tolerance) << i;
    }
}

TEST(PIMD, multirank_pile_l_damping_sentinels)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    int zero_tau_exceptions, negative_tau_exceptions;
    const auto zero_tau     = run_two_bead_damping_case(0.0, 1.0, 8642, zero_tau_exceptions);
    const auto negative_tau = run_two_bead_damping_case(-1.0, 1.0, 8642, negative_tau_exceptions);

    EXPECT_EQ(zero_tau_exceptions & FE_DIVBYZERO, 0);
    EXPECT_EQ(negative_tau_exceptions & FE_DIVBYZERO, 0);
    for (std::size_t i = 0; i < zero_tau.size(); ++i) {
        EXPECT_TRUE(std::isfinite(zero_tau[i]));
        EXPECT_TRUE(std::isfinite(negative_tau[i]));
        EXPECT_NEAR(negative_tau[i], zero_tau[i], 1.0e-14);
    }

    int first_seed_exceptions, second_seed_exceptions;
    const auto first_seed  = run_two_bead_damping_case(1.0, 0.0, 2468, first_seed_exceptions);
    const auto second_seed = run_two_bead_damping_case(1.0, 0.0, 9753, second_seed_exceptions);

    EXPECT_EQ(first_seed_exceptions & FE_DIVBYZERO, 0);
    EXPECT_EQ(second_seed_exceptions & FE_DIVBYZERO, 0);
    for (std::size_t i = 0; i < first_seed.size(); ++i) {
        EXPECT_TRUE(std::isfinite(first_seed[i]));
        EXPECT_TRUE(std::isfinite(second_seed[i]));
    }

    double max_velocity_difference = 0.0;
    for (int i = 0; i < 6; ++i) {
        EXPECT_NEAR(second_seed[24 + i], first_seed[24 + i], 1.0e-14);
        max_velocity_difference =
            std::max(max_velocity_difference, std::fabs(second_seed[6 + i] - first_seed[6 + i]));
    }
    EXPECT_GT(max_velocity_difference, 1.0e-6);
}

TEST(PIMD, multirank_npt_restart_restores_barostat_and_rng)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    remove_npt_restart_files();
    const auto continuous = run_two_bead_npt_case(6);
    const auto boundary   = write_two_bead_npt_restart();
    NPTState first_restored;
    const auto first = run_two_bead_npt_restart(first_restored);
    NPTState replay_restored;
    const auto replay = run_two_bead_npt_restart(replay_restored);
    remove_npt_restart_files();

    for (int bead = 0; bead < 2; ++bead) {
        const int offset = 18 * bead;
        for (int i = 0; i < 12; ++i) {
            EXPECT_NEAR(first_restored.atoms[offset + i], boundary.atoms[offset + i], 1.0e-14);
            EXPECT_NEAR(replay_restored.atoms[offset + i], boundary.atoms[offset + i], 1.0e-14);
        }
    }
    for (std::size_t i = 0; i < boundary.box.size(); ++i) {
        EXPECT_NEAR(first_restored.box[i], boundary.box[i], 1.0e-14);
        EXPECT_NEAR(replay_restored.box[i], boundary.box[i], 1.0e-14);
    }
    EXPECT_NEAR(first_restored.barostat_velocity, boundary.barostat_velocity, 1.0e-14);
    EXPECT_NEAR(replay_restored.barostat_velocity, boundary.barostat_velocity, 1.0e-14);

    for (std::size_t i = 0; i < first.atoms.size(); ++i) {
        EXPECT_TRUE(std::isfinite(continuous.atoms[i]));
        EXPECT_TRUE(std::isfinite(first.atoms[i]));
        EXPECT_TRUE(std::isfinite(replay.atoms[i]));
        EXPECT_NEAR(replay.atoms[i], first.atoms[i], 1.0e-14);
    }
    for (std::size_t i = 0; i < first.box.size(); ++i) {
        EXPECT_TRUE(std::isfinite(continuous.box[i]));
        EXPECT_TRUE(std::isfinite(first.box[i]));
        EXPECT_TRUE(std::isfinite(replay.box[i]));
        EXPECT_NEAR(replay.box[i], first.box[i], 1.0e-14);
    }
    EXPECT_TRUE(std::isfinite(continuous.barostat_velocity));
    EXPECT_TRUE(std::isfinite(first.barostat_velocity));
    EXPECT_TRUE(std::isfinite(replay.barostat_velocity));
    EXPECT_NEAR(replay.barostat_velocity, first.barostat_velocity, 1.0e-14);

    const double volume = (first.box[3] - first.box[0]) * (first.box[4] - first.box[1]) *
                          (first.box[5] - first.box[2]);
    EXPECT_GT(std::fabs(volume - 8000.0), 1.0e-10);
    EXPECT_GT(std::fabs(first.barostat_velocity), 1.0e-8);

    for (std::size_t i = 0; i < first.atoms.size(); ++i)
        EXPECT_NEAR(first.atoms[i], continuous.atoms[i], 1.0e-14);
    for (std::size_t i = 0; i < first.box.size(); ++i)
        EXPECT_NEAR(first.box[i], continuous.box[i], 1.0e-14);
    EXPECT_NEAR(first.barostat_velocity, continuous.barostat_velocity, 1.0e-14);
}

TEST(PIMD, multirank_anisotropic_pressure_targets)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    const auto baseline = run_two_bead_anisotropic_nph({0.5, 0.5, 0.5});
    const auto shifted  = run_two_bead_anisotropic_nph({0.8, 0.5, 0.5});

    for (int rank = 0; rank < nprocs; ++rank) {
        for (int dimension = 0; dimension < 3; ++dimension) {
            EXPECT_TRUE(std::isfinite(baseline[3 * rank + dimension]));
            EXPECT_TRUE(std::isfinite(shifted[3 * rank + dimension]));
            EXPECT_NEAR(baseline[3 * rank + dimension], baseline[dimension], 1.0e-14);
            EXPECT_NEAR(shifted[3 * rank + dimension], shifted[dimension], 1.0e-14);
        }
    }

    EXPECT_GT(std::fabs(shifted[0] - baseline[0]), 1.0e-6);
    EXPECT_NEAR(shifted[1], baseline[1], 1.0e-8);
    EXPECT_NEAR(shifted[2], baseline[2], 1.0e-8);
}

TEST(PIMD, multirank_duplicate_pressure_selector)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    const auto single    = run_two_bead_pressure_selector("x 0.8");
    const auto duplicate = run_two_bead_pressure_selector("x 0.5 x 0.8");

    for (std::size_t i = 0; i < single.npt.atoms.size(); ++i) {
        EXPECT_TRUE(std::isfinite(single.npt.atoms[i]));
        EXPECT_TRUE(std::isfinite(duplicate.npt.atoms[i]));
        EXPECT_NEAR(duplicate.npt.atoms[i], single.npt.atoms[i], 1.0e-14);
    }
    for (std::size_t i = 0; i < single.npt.box.size(); ++i) {
        EXPECT_TRUE(std::isfinite(single.npt.box[i]));
        EXPECT_TRUE(std::isfinite(duplicate.npt.box[i]));
        EXPECT_NEAR(duplicate.npt.box[i], single.npt.box[i], 1.0e-14);
    }
    EXPECT_NEAR(duplicate.npt.barostat_velocity, single.npt.barostat_velocity, 1.0e-14);
    EXPECT_TRUE(std::isfinite(single.total_enthalpy));
    EXPECT_TRUE(std::isfinite(duplicate.total_enthalpy));
    EXPECT_DOUBLE_EQ(duplicate.total_enthalpy, single.total_enthalpy);
}

TEST(PIMD, multirank_baoab_pressure_contract)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    const auto bzp = run_two_bead_baoab_bzp();

    for (double value : bzp.atoms)
        EXPECT_TRUE(std::isfinite(value));
    for (double value : bzp.box)
        EXPECT_TRUE(std::isfinite(value));
    EXPECT_TRUE(std::isfinite(bzp.barostat_velocity));

    const auto volume = [](const NPTState &state) {
        return (state.box[3] - state.box[0]) * (state.box[4] - state.box[1]) *
               (state.box[5] - state.box[2]);
    };
    EXPECT_GT(std::fabs(volume(bzp) - 8000.0), 1.0e-10);
    EXPECT_GT(std::fabs(bzp.barostat_velocity), 1.0e-8);

    void *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
    ASSERT_NE(lmp, nullptr);
    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_two_bead_test_system(lmp, 2, Ownership::DEFAULT);
    lammps_set_show_error(lmp, 0);
    lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nph integrator baoab "
                        "thermostat PILE_L 8642 tau 1.0 temp 1.0 iso 0.5 barostat MTTK taup 1.0 "
                        "fixcom no");
    const int has_error = lammps_has_error(lmp);
    EXPECT_EQ(has_error, 1);
    if (has_error) {
        char error_message[512];
        EXPECT_NE(lammps_get_last_error_message(lmp, error_message, sizeof(error_message)), 0);
        EXPECT_NE(
            std::string(error_message).find("MTTK barostat is not supported by fix pimd/langevin"),
            std::string::npos);
    }
    lammps_close(lmp);
}

TEST(PIMD, multirank_zero_barostat_velocity_limit)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    const auto volume = [](const NPTState &state) {
        return (state.box[3] - state.box[0]) * (state.box[4] - state.box[1]) *
               (state.box[5] - state.box[2]);
    };

    for (const char *pressure_keyword : {"iso", "x"}) {
        SCOPED_TRACE(pressure_keyword);
        const auto state = run_two_bead_zero_barostat_velocity(pressure_keyword);
        for (double value : state.atoms)
            EXPECT_TRUE(std::isfinite(value));
        for (double value : state.box)
            EXPECT_TRUE(std::isfinite(value));
        EXPECT_TRUE(std::isfinite(state.barostat_velocity));
        EXPECT_NEAR(volume(state), 4096.0, 1.0e-12);
        EXPECT_NEAR(state.barostat_velocity, 0.0, 1.0e-14);
    }
}

TEST(PIMD, multirank_singular_parameter_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    struct InvalidParameter {
        const char *setting;
        const char *message;
    };

    for (const auto &parameter :
         {InvalidParameter{"temp 0.0", "Invalid temp value for fix pimd/langevin"},
          InvalidParameter{"fmass 0.0", "Invalid fmass value for fix pimd/langevin"},
          InvalidParameter{"sp 0.0", "Invalid sp value for fix pimd/langevin"}}) {
        SCOPED_TRACE(parameter.setting);
        const std::string fix_command =
            "fix fpimd all pimd/langevin method nmpimd ensemble nve integrator obabo "
            "thermostat PILE_L 8642 tau 1.0 temp 1.0 fixcom no " +
            std::string(parameter.setting);
        expect_pimd_fix_error(fix_command, parameter.message, "lj");
    }
}

TEST(PIMD, multirank_langevin_dimension_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    expect_pimd_fix_error_on_partition("fix fpimd all pimd/langevin method nmpimd ensemble nve",
                                       "Fix pimd/langevin requires a 3d system", "lj",
                                       MPI_COMM_WORLD, "2x2", 2, 2);
}

TEST(PIMD, multirank_nvt_dimension_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    expect_pimd_fix_error_on_partition(
        "fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2",
        "Fix pimd/nvt requires a 3d system", "real", MPI_COMM_WORLD, "2x2", 2, 2);
}

TEST(PIMD, multirank_langevin_uneven_partition_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    expect_pimd_fix_error_on_uneven_partition(
        "fix fpimd all pimd/langevin method nmpimd ensemble nve",
        "Fix pimd/langevin requires the same number of processors in every partition", "lj");
}

TEST(PIMD, multirank_nvt_uneven_partition_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    expect_pimd_fix_error_on_uneven_partition(
        "fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2",
        "Fix pimd/nvt requires the same number of processors in every partition", "real");
}

TEST(PIMD, multirank_uniform_partition_acceptance)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    expect_pimd_fix_success("fix fpimd all pimd/langevin method nmpimd ensemble nve", "lj");
    expect_pimd_fix_success("fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2",
                            "real");
}

TEST(PIMD, multirank_langevin_atom_count_rejection)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_run_error_on_mismatched_atom_count(
        "fix fpimd all pimd/langevin method nmpimd ensemble nve",
        "Fix pimd/langevin requires the same number of atoms in every partition", "lj",
        two_rank_comm);
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_nvt_atom_count_rejection)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_run_error_on_mismatched_atom_count(
        "fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2",
        "Fix pimd/nvt requires the same number of atoms in every partition", "real", two_rank_comm);
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_matching_atom_count_acceptance)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_run_success_on_partition("fix fpimd all pimd/langevin method nmpimd ensemble nve",
                                         "lj", two_rank_comm);
    expect_pimd_run_success_on_partition(
        "fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2", "real",
        two_rank_comm);
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_langevin_mass_rejection)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_run_error_after_partition_change(
        "fix fpimd all pimd/langevin method nmpimd ensemble nve",
        "Fix pimd/langevin requires the same atom masses in every partition", "lj", "mass 1 2.0",
        two_rank_comm);
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_nvt_mass_rejection)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_run_error_after_partition_change(
        "fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2",
        "Fix pimd/nvt requires the same atom masses in every partition", "real", "mass 1 2.0",
        two_rank_comm);
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_matching_mass_acceptance)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_run_success_on_partition("fix fpimd all pimd/langevin method nmpimd ensemble nve",
                                         "lj", two_rank_comm);
    expect_pimd_run_success_on_partition(
        "fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2", "real",
        two_rank_comm);
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_langevin_atom_type_rejection)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_run_error_after_partition_change(
        "fix fpimd all pimd/langevin method nmpimd ensemble nve",
        "Fix pimd/langevin requires the same atom types for every atom ID in every partition", "lj",
        "set atom 2 type 2", two_rank_comm);
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_nvt_atom_type_rejection)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_run_error_after_partition_change(
        "fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2",
        "Fix pimd/nvt requires the same atom types for every atom ID in every partition", "real",
        "set atom 2 type 2", two_rank_comm);
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_matching_atom_type_acceptance)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_run_success_on_partition("fix fpimd all pimd/langevin method nmpimd ensemble nve",
                                         "lj", two_rank_comm);
    expect_pimd_run_success_on_partition(
        "fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2", "real",
        two_rank_comm);
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_nvt_group_rejection)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_fix_error_on_partition(
        "fix fpimd selected pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2",
        "Fix pimd/nvt only supports group all", "real", two_rank_comm, "2x1", 1, 3,
        "group selected id 1");
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_nvt_bosonic_group_rejection)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_fix_error_on_partition(
        "fix fpimd selected pimd/nvt/bosonic method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2",
        "Fix pimd/nvt/bosonic only supports group all", "real", two_rank_comm, "2x1", 1, 3,
        "group selected id 1");
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_langevin_bosonic_group_rejection)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_fix_error_on_partition("fix fpimd selected pimd/langevin/bosonic ensemble nve",
                                       "Fix pimd/langevin/bosonic only supports group all", "lj",
                                       two_rank_comm, "2x1", 1, 3, "group selected id 1");
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_langevin_bosonic_esynch_acceptance)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_run_success_on_partition(
        "fix fpimd all pimd/langevin/bosonic ensemble nve esynch yes", "lj", two_rank_comm);
    expect_pimd_run_success_on_partition(
        "fix fpimd all pimd/langevin/bosonic ensemble nve esynch no", "lj", two_rank_comm);
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, langevin_esynch_rejection)
{
    expect_pimd_fix_error("fix fpimd all pimd/langevin ensemble nve esynch yes",
                          "Unknown keyword esynch for fix pimd/langevin", "lj");
}

TEST(PIMD, multirank_langevin_bosonic_missing_esynch_value)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_fix_error_on_partition(
        "fix fpimd all pimd/langevin/bosonic ensemble nve esynch",
        "Illegal fix pimd/langevin/bosonic esynch command: missing argument(s)", "lj",
        two_rank_comm, "2x1", 1);
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_langevin_group_membership_rejection)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_run_error_after_partition_change(
        "fix fpimd selected pimd/langevin method nmpimd ensemble nve",
        "Fix pimd/langevin requires the same group membership for every atom ID in every partition",
        "lj", "group selected id 2", two_rank_comm, "group selected id 1");
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_langevin_matching_group_acceptance)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_run_success_on_partition(
        "fix fpimd selected pimd/langevin method nmpimd ensemble nve", "lj", two_rank_comm,
        "group selected id 1");
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_langevin_cell_rejection)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    for (const char *mismatch :
         {"change_box all x final 0.0 21.0 remap", "change_box all boundary f p p"}) {
        SCOPED_TRACE(mismatch);
        expect_pimd_run_error_after_partition_change(
            "fix fpimd all pimd/langevin method nmpimd ensemble nve",
            "Fix pimd/langevin requires the same simulation cell in every partition", "lj",
            mismatch, two_rank_comm);
    }
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_nvt_cell_rejection)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    for (const char *mismatch :
         {"change_box all x final 0.0 21.0 remap", "change_box all boundary f p p"}) {
        SCOPED_TRACE(mismatch);
        expect_pimd_run_error_after_partition_change(
            "fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2",
            "Fix pimd/nvt requires the same simulation cell in every partition", "real", mismatch,
            two_rank_comm);
    }
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_matching_cell_acceptance)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);
    expect_pimd_run_success_on_partition("fix fpimd all pimd/langevin method nmpimd ensemble nve",
                                         "lj", two_rank_comm);
    expect_pimd_run_success_on_partition(
        "fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2", "real",
        two_rank_comm);
    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_scale_method_order_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    expect_pimd_fix_error(
        "fix fpimd all pimd/langevin scale 1.0 method pimd ensemble nve "
        "integrator obabo thermostat PILE_L 8642 tau 1.0 temp 1.0 fixcom no",
        "Scale parameter of the PILE_L thermostat is not supported with method pimd", "lj");
}

TEST(PIMD, multirank_fixcom_value_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    expect_pimd_fix_error("fix fpimd all pimd/langevin method nmpimd ensemble nve integrator obabo "
                          "thermostat PILE_L 8642 tau 1.0 temp 1.0 fixcom invalid",
                          "Unknown fixcom value invalid for fix pimd/langevin", "lj");
}

TEST(PIMD, multirank_thermostat_value_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    expect_pimd_fix_error("fix fpimd all pimd/langevin method nmpimd ensemble nve temp 1.0 "
                          "thermostat INVALID 8642",
                          "Unknown thermostat parameter for fix pimd/langevin", "lj");
}

TEST(PIMD, multirank_thermostat_seed_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    expect_pimd_fix_error("fix fpimd all pimd/langevin method nmpimd ensemble nvt integrator obabo "
                          "thermostat PILE_L -10 tau 1.0 temp 1.0 fixcom no",
                          "Invalid thermostat seed value for fix pimd/langevin", "lj");
}

TEST(PIMD, multirank_thermostat_seed_presence_rejection)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    MPI_Comm two_rank_comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / 2, me, &two_rank_comm);

    struct ThermostattedEnsemble {
        const char *name;
        const char *label;
        const char *pressure;
    };

    for (const auto &ensemble :
         {ThermostattedEnsemble{"nvt", "NVT", ""},
          ThermostattedEnsemble{"npt", "NPT", "iso 0.5 barostat BZP taup 1.0"}}) {
        SCOPED_TRACE(ensemble.name);
        const std::string fix_command =
            "fix fpimd all pimd/langevin method nmpimd ensemble " + std::string(ensemble.name) +
            " integrator obabo temp 1.0 tau 1.0 fixcom no " + ensemble.pressure;
        const std::string message =
            "Thermostat seed must be specified for fix pimd/langevin with " +
            std::string(ensemble.label) + " ensemble";
        expect_pimd_fix_error_on_partition(fix_command, message.c_str(), "lj", two_rank_comm, "2x1",
                                           1);
    }

    MPI_Comm_free(&two_rank_comm);
}

TEST(PIMD, multirank_unthermostatted_seed_omission)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    struct UnthermostattedEnsemble {
        const char *name;
        const char *pressure;
    };

    for (const auto &ensemble : {UnthermostattedEnsemble{"nve", ""},
                                 UnthermostattedEnsemble{"nph", "iso 0.5 barostat BZP taup 1.0"}}) {
        SCOPED_TRACE(ensemble.name);
        const std::string fix_command =
            "fix fpimd all pimd/langevin method nmpimd ensemble " + std::string(ensemble.name) +
            " integrator obabo temp 1.0 tau 1.0 fixcom no " + ensemble.pressure;
        expect_pimd_fix_success(fix_command, "lj");
    }
}

TEST(PIMD, multirank_langevin_missing_value_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    expect_pimd_fix_error("fix fpimd all pimd/langevin method nmpimd ensemble nve integrator obabo "
                          "thermostat PILE_L 8642 tau 1.0 temp 1.0 fixcom no temp",
                          "Illegal fix pimd/langevin temp command: missing argument(s)", "lj");
}

TEST(PIMD, multirank_langevin_missing_seed_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    expect_pimd_fix_error("fix fpimd all pimd/langevin method nmpimd ensemble nve temp 1.0 "
                          "thermostat PILE_L",
                          "Illegal fix pimd/langevin thermostat command: missing argument(s)",
                          "lj");
}

TEST(PIMD, multirank_nvt_missing_value_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    expect_pimd_fix_error("fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2 temp",
                          "Illegal fix pimd/nvt temp command: missing argument(s)", "real");
}

TEST(PIMD, multirank_nvt_singular_parameter_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    struct InvalidParameter {
        const char *setting;
        const char *message;
    };

    for (const auto &parameter :
         {InvalidParameter{"temp 0.0", "Invalid temp value for fix pimd/nvt"},
          InvalidParameter{"fmass 0.0", "Invalid fmass value 0 for fix pimd/nvt"},
          InvalidParameter{"sp 0.0", "Invalid sp value for fix pimd/nvt"}}) {
        SCOPED_TRACE(parameter.setting);
        const std::string fix_command =
            "fix fpimd all pimd/nvt method pimd temp 1.0 fmass 1.0 sp 1.0 nhc 2 " +
            std::string(parameter.setting);
        expect_pimd_fix_error(fix_command, parameter.message, "real");
    }
}

TEST(PIMD, multirank_per_atom_mass_rejection)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    for (const char *style :
         {"pimd/langevin", "pimd/langevin/bosonic", "pimd/nvt", "pimd/nvt/bosonic"}) {
        SCOPED_TRACE(style);
        auto *lmp = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
        lammps_command(lmp, "units real");
        lammps_command(lmp, "atom_style sphere");
        lammps_command(lmp, "atom_modify map array");
        lammps_command(lmp, "region box block 0 8 0 8 0 8");
        lammps_command(lmp, "create_box 1 box");
        lammps_command(lmp, "create_atoms 1 single 4 4 4 units box");
        lammps_set_show_error(lmp, 0);
        std::string command = "fix fpimd all " + std::string(style);
        if (std::string(style).find("langevin") != std::string::npos) command += " ensemble nve";
        lammps_command(lmp, command.c_str());
        const int has_error = lammps_has_error(lmp);
        EXPECT_EQ(has_error, 1);
        if (has_error) {
            char message[512];
            lammps_get_last_error_message(lmp, message, sizeof(message));
            EXPECT_NE(std::string(message).find("requires per-type atom masses"),
                      std::string::npos);
        }
        lammps_close(lmp);
    }
}

TEST(PIMD, multirank_inactive_barostat_restart_state)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    for (const char *ensemble : {"nve", "nvt"}) {
        SCOPED_TRACE(ensemble);
        auto *handle = open_two_bead_partition(MPI_COMM_WORLD, "2x2");
        create_two_bead_test_system(handle, 2, Ownership::DEFAULT);
        const std::string command = "fix fpimd all pimd/langevin ensemble " +
                                    std::string(ensemble) +
                                    " thermostat PILE_L 2468 temp 1.0 fixcom no";
        lammps_command(handle, command.c_str());
        lammps_command(handle, "run 0 post no");
        ASSERT_EQ(lammps_has_error(handle), 0);
        auto *lmp = static_cast<LAMMPS_NS::LAMMPS *>(handle);
        auto *fix = lmp->modify->get_fix_by_id("fpimd");
        ASSERT_NE(fix, nullptr);
        const int root = lammps_extract_setting(handle, "world_rank") == 0;
        FILE *file     = root ? std::tmpfile() : nullptr;
        int file_ok    = !root || file != nullptr;
        MPI_Allreduce(MPI_IN_PLACE, &file_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if (!file_ok) {
            if (file) std::fclose(file);
            lammps_close(handle);
        }
        ASSERT_TRUE(file_ok);
        fix->write_restart(file);
        if (root && file) {
            std::rewind(file);
            int bytes = 0;
            std::array<double, 6> barostat{};
            EXPECT_EQ(std::fread(&bytes, sizeof(bytes), 1, file), 1U);
            EXPECT_GE(bytes, int(sizeof(barostat)));
            EXPECT_EQ(std::fread(barostat.data(), sizeof(double), barostat.size(), file),
                      barostat.size());
            for (double value : barostat)
                EXPECT_DOUBLE_EQ(value, 0.0);
            std::fclose(file);
        }
        lammps_close(handle);
    }
}
