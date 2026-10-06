// unit tests for PLUMED coupling to PIMD through the library interface

#define LAMMPS_LIB_MPI 1
#include "exceptions.h"
#include "fix.h"
#include "lammps.h"
#include "library.h"
#include "lmptype.h"
#include "modify.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../testing/test_mpi_main.h"

using ::LAMMPS_NS::tagint;
using ::testing::HasSubstr;

namespace {

void *open_multirank_partition()
{
    const char *args[] = {"LAMMPS_test", "-screen", "none", "-log",    "none", "-partition",
                          "2x2",         "-in",     "none", "-nocite", nullptr};
    char **argv        = (char **)args;
    int argc           = (sizeof(args) / sizeof(char *)) - 1;
    return lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
}

void create_multirank_two_atom_system(void *lmp)
{
    lammps_command(lmp, "variable x2 world 19.0 10.0");
    lammps_command(lmp, "variable y2 world 10.0 17.0");
    lammps_command(lmp, "units lj");
    lammps_command(lmp, "atom_style atomic");
    lammps_command(lmp, "atom_modify map array");
    lammps_command(lmp, "boundary p p p");
    lammps_command(lmp, "region box block 0.0 20.0 0.0 20.0 0.0 20.0");
    lammps_command(lmp, "create_box 1 box");
    lammps_command(lmp, "create_atoms 1 single 10.0 10.0 2.0");
    lammps_command(lmp, "create_atoms 1 single ${x2} ${y2} 2.0");
    lammps_command(lmp, "mass 1 1.0");
    lammps_command(lmp, "pair_style lj/cut 2.5");
    lammps_command(lmp, "pair_coeff * * 0.0 1.0");
    lammps_command(lmp, "timestep 0.001");
    lammps_command(lmp, "velocity all set 0.0 0.0 0.0");
}

} // namespace
namespace {

double extract_pimd_vector(void *lmp, int index)
{
    auto *value =
        (double *)lammps_extract_fix(lmp, "fpimd", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, index, 0);
    EXPECT_NE(value, nullptr);
    if (!value) return 0.0;
    const double result = *value;
    lammps_free(value);
    return result;
}

struct CentroidNPTState {
    double p_cv              = 0.0;
    double pe_bead           = 0.0;
    double tote              = 0.0;
    double barostat_velocity = 0.0;
    double volume            = 0.0;
    double bias              = 0.0;
    std::array<double, 6> pressure{};
};

CentroidNPTState run_centroid_npt_volume_leg(const char *plumed_file, const char *plumed_log,
                                             double external_pressure)
{
    void *lmp = open_multirank_partition();
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_multirank_two_atom_system(lmp);
    lammps_command(lmp, "thermo 0");
    lammps_command(lmp, "compute plumed_pressure all pressure NULL virial");
    const std::string pimd_command =
        "fix fpimd all pimd/langevin method nmpimd ensemble npt integrator obabo "
        "thermostat PILE_L 2468 tau 1.0 temp 1.0 iso " +
        std::to_string(external_pressure) + " barostat BZP taup 1.0 fixcom no";
    lammps_command(lmp, pimd_command.c_str());
    const std::string plumed_command = "fix bias all plumed plumedfile " +
                                       std::string(plumed_file) + " outfile " + plumed_log +
                                       " path_integral centroid pimd_fix fpimd";
    lammps_command(lmp, plumed_command.c_str());
    lammps_command(lmp, "run 0 post no");

    CentroidNPTState state;
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (me == 0) state.p_cv = extract_pimd_vector(lmp, 9);
    MPI_Bcast(&state.p_cv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    state.pe_bead = extract_pimd_vector(lmp, 2);
    state.tote    = extract_pimd_vector(lmp, 3);
    auto *pressure =
        (double *)lammps_extract_compute(lmp, "plumed_pressure", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
    EXPECT_NE(pressure, nullptr);
    if (pressure)
        for (int i = 0; i < 6; ++i)
            state.pressure[i] = pressure[i];
    auto *bias =
        (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
    EXPECT_NE(bias, nullptr);
    if (bias) {
        state.bias = *bias;
        lammps_free(bias);
    }

    lammps_command(lmp, "run 1 post no");
    state.barostat_velocity = extract_pimd_vector(lmp, 10);
    double boxlo[3], boxhi[3], xy, yz, xz;
    int periodicity[3], boxflag;
    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, periodicity, &boxflag);
    state.volume = (boxhi[0] - boxlo[0]) * (boxhi[1] - boxlo[1]) * (boxhi[2] - boxlo[2]);

    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

CentroidNPTState run_nmpimd_bead_mode_volume_leg(const char *plumed_file, const char *plumed_log,
                                                 const char *mode, const char *ensemble,
                                                 double external_pressure)
{
    void *lmp = open_multirank_partition();
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    create_multirank_two_atom_system(lmp);
    lammps_command(lmp, "thermo 0");
    lammps_command(lmp, "compute plumed_pressure all pressure NULL virial");
    std::string pimd_command = "fix fpimd all pimd/langevin method nmpimd ensemble " +
                               std::string(ensemble) + " integrator obabo ";
    if (std::string(ensemble) == "npt") pimd_command += "thermostat PILE_L 2468 tau 1.0 ";
    pimd_command +=
        "temp 1.0 iso " + std::to_string(external_pressure) + " barostat BZP taup 1.0 fixcom no";
    lammps_command(lmp, pimd_command.c_str());
    const std::string plumed_command = "fix bias all plumed plumedfile " +
                                       std::string(plumed_file) + " outfile " + plumed_log +
                                       " path_integral " + mode + " pimd_fix fpimd";
    lammps_command(lmp, plumed_command.c_str());
    lammps_command(lmp, "run 0 post no");

    CentroidNPTState state;
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (me == 0) state.p_cv = extract_pimd_vector(lmp, 9);
    MPI_Bcast(&state.p_cv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    state.pe_bead = extract_pimd_vector(lmp, 2);
    state.tote    = extract_pimd_vector(lmp, 3);
    auto *pressure =
        (double *)lammps_extract_compute(lmp, "plumed_pressure", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
    EXPECT_NE(pressure, nullptr);
    if (pressure)
        for (int i = 0; i < 6; ++i)
            state.pressure[i] = pressure[i];
    auto *bias =
        (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
    EXPECT_NE(bias, nullptr);
    if (bias) {
        state.bias = *bias;
        lammps_free(bias);
    }

    lammps_command(lmp, "run 1 post no");
    state.barostat_velocity = extract_pimd_vector(lmp, 10);
    double boxlo[3], boxhi[3], xy, yz, xz;
    int periodicity[3], boxflag;
    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, periodicity, &boxflag);
    state.volume = (boxhi[0] - boxlo[0]) * (boxhi[1] - boxlo[1]) * (boxhi[2] - boxlo[2]);

    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

struct CentroidNVTState {
    double p_cv   = 0.0;
    double volume = 0.0;
    double bias   = 0.0;
    std::array<double, 6> pressure{};
};

CentroidNVTState run_centroid_nvt_volume_leg(const char *plumed_file, const char *plumed_log)
{
    void *lmp = open_multirank_partition();
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_multirank_two_atom_system(lmp);
    lammps_command(lmp, "thermo 0");
    lammps_command(lmp, "compute plumed_pressure all pressure NULL virial");
    lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                        "integrator baoab thermostat PILE_L 2468 tau 1.0 temp 1.0 fixcom no");
    const std::string plumed_command = "fix bias all plumed plumedfile " +
                                       std::string(plumed_file) + " outfile " + plumed_log +
                                       " path_integral centroid pimd_fix fpimd";
    lammps_command(lmp, plumed_command.c_str());
    lammps_command(lmp, "run 0 post no");

    CentroidNVTState state;
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (me == 0) state.p_cv = extract_pimd_vector(lmp, 9);
    MPI_Bcast(&state.p_cv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    auto *pressure =
        (double *)lammps_extract_compute(lmp, "plumed_pressure", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
    EXPECT_NE(pressure, nullptr);
    if (pressure)
        for (int i = 0; i < 6; ++i)
            state.pressure[i] = pressure[i];
    auto *bias =
        (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
    EXPECT_NE(bias, nullptr);
    if (bias) {
        state.bias = *bias;
        lammps_free(bias);
    }

    lammps_command(lmp, "run 1 post no");
    double boxlo[3], boxhi[3], xy, yz, xz;
    int periodicity[3], boxflag;
    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, periodicity, &boxflag);
    state.volume = (boxhi[0] - boxlo[0]) * (boxhi[1] - boxlo[1]) * (boxhi[2] - boxlo[2]);

    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

struct CentroidNVTContinuationState {
    std::array<double, 36> atoms{};
    std::array<double, 10> pimd{};
    std::array<double, 4> bias{};
    double volume = 0.0;
};

CentroidNVTContinuationState extract_centroid_nvt_continuation_state(void *lmp)
{
    CentroidNVTContinuationState state;
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

    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    const int world_size = lammps_extract_setting(lmp, "world_size");
    const int offset     = 18 * (me / world_size);
    if (nlocal && tags && positions && velocities && forces) {
        for (int i = 0; i < *nlocal; ++i) {
            EXPECT_GE(tags[i], 1);
            EXPECT_LE(tags[i], 2);
            if (tags[i] < 1 || tags[i] > 2) continue;
            const int atom = tags[i] - 1;
            for (int d = 0; d < 3; ++d) {
                state.atoms[offset + 3 * atom + d]      = positions[i][d];
                state.atoms[offset + 6 + 3 * atom + d]  = velocities[i][d];
                state.atoms[offset + 12 + 3 * atom + d] = forces[i][d];
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, state.atoms.data(), static_cast<int>(state.atoms.size()),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < 10; ++i)
        state.pimd[i] = extract_pimd_vector(lmp, i);

    double local_bias = 0.0;
    auto *bias =
        (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
    EXPECT_NE(bias, nullptr);
    if (bias) {
        local_bias = *bias;
        lammps_free(bias);
    }
    MPI_Allgather(&local_bias, 1, MPI_DOUBLE, state.bias.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);

    double boxlo[3], boxhi[3], xy, yz, xz;
    int periodicity[3], boxflag;
    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, periodicity, &boxflag);
    state.volume = (boxhi[0] - boxlo[0]) * (boxhi[1] - boxlo[1]) * (boxhi[2] - boxlo[2]);
    return state;
}

void add_centroid_nmpimd_nvt_fixes(void *lmp, const char *plumed_file, const char *plumed_log,
                                   const std::string &nonfinite_trace_prefix = {})
{
    std::string pimd_command = "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                               "integrator baoab thermostat PILE_L 2468 tau 1.0 temp 1.0 "
                               "fixcom no";
    if (!nonfinite_trace_prefix.empty())
        pimd_command += " nonfinite_trace " + nonfinite_trace_prefix;
    lammps_command(lmp, pimd_command.c_str());
    std::string plumed_command = "fix bias all plumed plumedfile " + std::string(plumed_file) +
                                 " outfile " + plumed_log +
                                 " path_integral centroid pimd_fix fpimd";
    if (!nonfinite_trace_prefix.empty())
        plumed_command += " nonfinite_trace " + nonfinite_trace_prefix;
    lammps_command(lmp, plumed_command.c_str());
}

void remove_centroid_nmpimd_nvt_restart_files()
{
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (me == 0) {
        std::remove("test_plumed_nmpimd_nvt_restart.0");
        std::remove("test_plumed_nmpimd_nvt_restart.1");
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

CentroidNVTContinuationState
run_centroid_nmpimd_nvt_segments(int first_steps, int second_steps, bool restart,
                                 const char *initial_plumed_file, const char *restart_plumed_file,
                                 const char *initial_log, const char *restart_log,
                                 const std::string &nonfinite_trace_prefix = {})
{
    void *lmp = open_multirank_partition();
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
    create_multirank_two_atom_system(lmp);
    lammps_command(lmp, "timestep 0.00001");
    lammps_command(lmp, "thermo 0");
    add_centroid_nmpimd_nvt_fixes(lmp, initial_plumed_file, initial_log, nonfinite_trace_prefix);
    const std::string first_run = "run " + std::to_string(first_steps);
    lammps_command(lmp, first_run.c_str());

    if (second_steps > 0) {
        if (restart) {
            lammps_command(lmp, "variable restart_file world test_plumed_nmpimd_nvt_restart.0 "
                                "test_plumed_nmpimd_nvt_restart.1");
            lammps_command(lmp, "write_restart ${restart_file}");
            lammps_close(lmp);

            lmp = open_multirank_partition();
            EXPECT_NE(lmp, nullptr);
            if (!lmp) return {};
            lammps_command(lmp, "variable restart_file world test_plumed_nmpimd_nvt_restart.0 "
                                "test_plumed_nmpimd_nvt_restart.1");
            lammps_command(lmp, "read_restart ${restart_file}");
            lammps_command(lmp, "thermo 0");
            add_centroid_nmpimd_nvt_fixes(lmp, restart_plumed_file, restart_log,
                                          nonfinite_trace_prefix);
        }
        const std::string second_run = "run " + std::to_string(second_steps);
        lammps_command(lmp, second_run.c_str());
    }

    const auto state = extract_centroid_nvt_continuation_state(lmp);
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

CentroidNVTContinuationState run_nmpimd_bead_mode_segments(
    const char *mode, const char *restart_prefix, int first_steps, int second_steps, bool restart,
    const char *initial_plumed_file, const char *restart_plumed_file, const char *initial_log,
    const char *restart_log, const std::string &nonfinite_trace_prefix = {})
{
    auto add_fixes = [&](void *lmp_instance, const char *plumed_file, const char *plumed_log) {
        std::string pimd_command = "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                                   "integrator baoab thermostat PILE_L 2468 tau 1.0 temp 1.0 "
                                   "fixcom no";
        if (!nonfinite_trace_prefix.empty())
            pimd_command += " nonfinite_trace " + nonfinite_trace_prefix;
        lammps_command(lmp_instance, pimd_command.c_str());
        std::string plumed_command = "fix bias all plumed plumedfile " + std::string(plumed_file) +
                                     " outfile " + plumed_log + " path_integral " + mode +
                                     " pimd_fix fpimd";
        if (!nonfinite_trace_prefix.empty())
            plumed_command += " nonfinite_trace " + nonfinite_trace_prefix;
        lammps_command(lmp_instance, plumed_command.c_str());
    };
    auto restart_variable = [&]() {
        return "variable restart_file world " + std::string(restart_prefix) + ".0 " +
               restart_prefix + ".1";
    };

    void *lmp = open_multirank_partition();
    EXPECT_NE(lmp, nullptr);
    if (!lmp) return {};

    create_multirank_two_atom_system(lmp);
    lammps_command(lmp, "timestep 0.00001");
    lammps_command(lmp, "thermo 0");
    add_fixes(lmp, initial_plumed_file, initial_log);
    lammps_command(lmp, ("run " + std::to_string(first_steps)).c_str());

    if (second_steps > 0) {
        if (restart) {
            const std::string variable = restart_variable();
            lammps_command(lmp, variable.c_str());
            lammps_command(lmp, "write_restart ${restart_file}");
            lammps_close(lmp);

            lmp = open_multirank_partition();
            EXPECT_NE(lmp, nullptr);
            if (!lmp) return {};
            lammps_command(lmp, variable.c_str());
            lammps_command(lmp, "read_restart ${restart_file}");
            lammps_command(lmp, "thermo 0");
            add_fixes(lmp, restart_plumed_file, restart_log);
        }
        std::string second_run = "run " + std::to_string(second_steps);
        if (!restart) second_run += " pre no";
        lammps_command(lmp, second_run.c_str());
    }

    const auto state = extract_centroid_nvt_continuation_state(lmp);
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    return state;
}

void expect_plumed_nonfinite_trace(const std::string &prefix)
{
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    const std::string plumed_file = prefix + ".dat";
    const std::string plumed_log  = prefix + ".log";
    if (me == 0) {
        std::ofstream input(plumed_file);
        input << "d: DISTANCE ATOMS=1,2 NOPBC\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    void *lmp = open_multirank_partition();
    ASSERT_NE(lmp, nullptr);
    create_multirank_two_atom_system(lmp);
    lammps_command(lmp, "thermo 0");
    lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                        "integrator baoab thermostat PILE_L 2468 tau 1.0 temp 1.0 fixcom no");
    const std::string plumed_command =
        "fix bias all plumed plumedfile " + plumed_file + " outfile " + plumed_log +
        " path_integral bead_density pimd_fix fpimd nonfinite_trace " + prefix;
    lammps_command(lmp, plumed_command.c_str());
    lammps_command(lmp, "run 0 post no");
    ASSERT_EQ(lammps_has_error(lmp), 0);

    auto *nlocal = (int *)lammps_extract_global(lmp, "nlocal");
    auto *tags   = (tagint *)lammps_extract_atom(lmp, "id");
    ASSERT_NE(nlocal, nullptr);
    ASSERT_NE(tags, nullptr);
    tagint local_tag           = std::numeric_limits<tagint>::max();
    int local_index            = -1;
    constexpr int target_world = 1;
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
        auto **positions = (double **)lammps_extract_atom(lmp, "x");
        ASSERT_NE(positions, nullptr);
        positions[local_index][0] = std::numeric_limits<double>::quiet_NaN();
    }

    auto *lammps     = static_cast<LAMMPS_NS::LAMMPS *>(lmp);
    auto *plumed_fix = lammps->modify->get_fix_by_id("bias");
    ASSERT_NE(plumed_fix, nullptr);
    int caught = 0;
    std::string message;
    try {
        plumed_fix->post_force(0);
    } catch (const LAMMPS_NS::LAMMPSAbortException &exception) {
        caught  = 1;
        message = exception.what();
    }
    EXPECT_EQ(caught, 1);
    EXPECT_NE(message.find("nonfinite trace detected"), std::string::npos);
    EXPECT_NE(message.find("stage pre-plumed"), std::string::npos);
    lammps_close(lmp);

    const int world_rank   = owner - 2 * target_world;
    const std::string path = prefix + ".plumed.step0.u" + std::to_string(owner) + ".w" +
                             std::to_string(target_world) + ".r" + std::to_string(world_rank) +
                             ".txt";
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        std::ifstream input(path);
        ASSERT_TRUE(input.good()) << path;
        std::ostringstream buffer;
        buffer << input.rdbuf();
        const std::string report = buffer.str();
        EXPECT_THAT(report, HasSubstr("schema=plumed-nonfinite-trace-v1"));
        EXPECT_THAT(report, HasSubstr("stage=pre-plumed"));
        EXPECT_THAT(report, HasSubstr("atom_tag=" + std::to_string(expected_tag)));
        EXPECT_THAT(report, HasSubstr("field=position"));
        EXPECT_THAT(report, HasSubstr("component=0"));
        std::remove(path.c_str());
        std::remove(plumed_file.c_str());
        std::remove((plumed_log + ".0").c_str());
        std::remove((plumed_log + ".1").c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void expect_plumed_force_delta_trace(const std::string &prefix)
{
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    const std::string plumed_file = prefix + ".dat";
    const std::string plumed_log  = prefix + ".log";
    if (me == 0) {
        std::ofstream input(plumed_file);
        input << "d: DISTANCE ATOMS=1,2 NOPBC\n"
              << "bias: RESTRAINT ARG=d AT=0.0 KAPPA=1e308\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    void *lmp = open_multirank_partition();
    ASSERT_NE(lmp, nullptr);
    create_multirank_two_atom_system(lmp);
    lammps_command(lmp, "thermo 0");

    auto *nlocal = (int *)lammps_extract_global(lmp, "nlocal");
    auto *tags   = (tagint *)lammps_extract_atom(lmp, "id");
    ASSERT_NE(nlocal, nullptr);
    ASSERT_NE(tags, nullptr);
    tagint local_tag           = std::numeric_limits<tagint>::max();
    constexpr int target_world = 0;
    if (me / 2 == target_world)
        for (int i = 0; i < *nlocal; ++i)
            local_tag = std::min(local_tag, tags[i]);
    std::array<tagint, 4> candidate_tags{};
    MPI_Allgather(&local_tag, 1, MPI_LMP_TAGINT, candidate_tags.data(), 1, MPI_LMP_TAGINT,
                  MPI_COMM_WORLD);
    int owner = 2 * target_world;
    for (int rank = 2 * target_world; rank < 2 * target_world + 2; ++rank)
        if (candidate_tags[rank] < candidate_tags[owner]) owner = rank;
    const tagint expected_tag = candidate_tags[owner];
    ASSERT_NE(expected_tag, std::numeric_limits<tagint>::max());

    lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                        "integrator baoab thermostat PILE_L 2468 tau 1.0 temp 1.0 fixcom no");
    const std::string plumed_command =
        "fix bias all plumed plumedfile " + plumed_file + " outfile " + plumed_log +
        " path_integral bead_density pimd_fix fpimd nonfinite_trace " + prefix;
    lammps_command(lmp, plumed_command.c_str());
    lammps_set_show_error(lmp, 0);
    lammps_command(lmp, "run 0 post no");
    ASSERT_EQ(lammps_has_error(lmp), 1);
    char error_message[1024];
    ASSERT_NE(lammps_get_last_error_message(lmp, error_message, sizeof(error_message)), 0);
    EXPECT_THAT(error_message, HasSubstr("plumed-force-delta"));
    EXPECT_THAT(error_message, HasSubstr("stage post-perform-pre-scale"));
    lammps_close(lmp);

    const int world_rank   = owner - 2 * target_world;
    const std::string path = prefix + ".plumed.step0.u" + std::to_string(owner) + ".w" +
                             std::to_string(target_world) + ".r" + std::to_string(world_rank) +
                             ".txt";
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        std::ifstream input(path);
        ASSERT_TRUE(input.good()) << path;
        std::ostringstream buffer;
        buffer << input.rdbuf();
        const std::string report = buffer.str();
        EXPECT_THAT(report, HasSubstr("stage=post-perform-pre-scale"));
        EXPECT_THAT(report, HasSubstr("atom_tag=" + std::to_string(expected_tag)));
        EXPECT_THAT(report, HasSubstr("field=plumed-force-delta"));
        EXPECT_THAT(report, HasSubstr("physical_force=0 0 0"));
        std::remove(path.c_str());
        std::remove(plumed_file.c_str());
        std::remove((plumed_log + ".0").c_str());
        std::remove((plumed_log + ".1").c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace

TEST(MPI, plumed_nonfinite_trace_finite_parity)
{
    auto compare = [](const CentroidNVTContinuationState &baseline,
                      const CentroidNVTContinuationState &traced) {
        for (std::size_t i = 0; i < baseline.atoms.size(); ++i)
            EXPECT_DOUBLE_EQ(traced.atoms[i], baseline.atoms[i]) << i;
        for (std::size_t i = 0; i < baseline.pimd.size(); ++i)
            EXPECT_DOUBLE_EQ(traced.pimd[i], baseline.pimd[i]) << i;
        for (std::size_t i = 0; i < baseline.bias.size(); ++i)
            EXPECT_DOUBLE_EQ(traced.bias[i], baseline.bias[i]) << i;
        EXPECT_DOUBLE_EQ(traced.volume, baseline.volume);
    };

    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    const char *centroid_file = "test_plumed_nonfinite_trace_centroid.dat";
    const char *mean_file     = "test_plumed_nonfinite_trace_mean.dat";
    const char *density_file  = "test_plumed_nonfinite_trace_density.dat";
    if (me == 0) {
        const std::string distance = "d: DISTANCE ATOMS=1,2 NOPBC\n";
        std::ofstream centroid(centroid_file);
        centroid << distance << "bias: RESTRAINT ARG=d AT=6.5 KAPPA=4.0\n";
        std::ofstream mean(mean_file);
        mean << distance << "mean: ENSEMBLE ARG=d\n"
             << "bias: RESTRAINT ARG=mean.d AT=6.5 KAPPA=4.0\n";
        std::ofstream density(density_file);
        density << distance << "bias: RESTRAINT ARG=d AT=6.5 KAPPA=4.0\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const auto centroid_baseline = run_centroid_nmpimd_nvt_segments(
        3, 0, false, centroid_file, centroid_file, "test_trace_centroid_base.log",
        "test_trace_centroid_base_restart.log");
    const auto centroid_traced = run_centroid_nmpimd_nvt_segments(
        3, 0, false, centroid_file, centroid_file, "test_trace_centroid_on.log",
        "test_trace_centroid_on_restart.log", "test_centroid_finite_trace");
    compare(centroid_baseline, centroid_traced);

    const std::array<std::pair<const char *, const char *>, 2> bead_modes = {
        std::pair{"bead_mean", mean_file}, std::pair{"bead_density", density_file}};
    for (const auto &[mode, file] : bead_modes) {
        const std::string stem = std::string("test_trace_") + mode;
        const auto baseline    = run_nmpimd_bead_mode_segments(
            mode, (stem + "_base_restart").c_str(), 3, 0, false, file, file,
            (stem + "_base.log").c_str(), (stem + "_base_restart.log").c_str());
        const auto traced = run_nmpimd_bead_mode_segments(
            mode, (stem + "_on_restart").c_str(), 3, 0, false, file, file,
            (stem + "_on.log").c_str(), (stem + "_on_restart.log").c_str(), stem + "_finite_trace");
        compare(baseline, traced);
    }

    if (me == 0) {
        std::remove(centroid_file);
        std::remove(mean_file);
        std::remove(density_file);
        for (const char *mode : {"centroid", "bead_mean", "bead_density"}) {
            const std::string stem = std::string("test_trace_") + mode;
            for (const char *leg : {"base", "on"}) {
                const std::string log = stem + "_" + leg + ".log";
                std::remove(log.c_str());
                std::remove((log + ".0").c_str());
                std::remove((log + ".1").c_str());
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

TEST(MPI, plumed_nonfinite_trace_reports_nan_position)
{
    expect_plumed_nonfinite_trace("test_plumed_trace_nan_position");
}

TEST(MPI, plumed_nonfinite_trace_reports_force_delta)
{
    expect_plumed_force_delta_trace("test_plumed_trace_force_delta");
}

TEST(MPI, plumed_pimd_input_contract)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    auto open_lammps = []() {
        const char *args[] = {"LAMMPS_test", "-screen", "none", "-log",    "none", "-partition",
                              "4x1",         "-in",     "none", "-nocite", nullptr};
        char **argv        = (char **)args;
        int argc           = (sizeof(args) / sizeof(char *)) - 1;
        void *lmp          = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
        lammps_command(lmp, "units lj");
        lammps_command(lmp, "atom_style atomic");
        lammps_command(lmp, "atom_modify map array");
        lammps_command(lmp, "boundary p p p");
        lammps_command(lmp, "region box block 0.0 20.0 0.0 20.0 0.0 20.0");
        lammps_command(lmp, "create_box 1 box");
        lammps_command(lmp, "create_atoms 1 single 10.0 10.0 2.0");
        lammps_command(lmp, "create_atoms 1 single 19.0 10.0 2.0");
        lammps_command(lmp, "mass 1 1.0");
        lammps_command(lmp, "pair_style lj/cut 2.5");
        lammps_command(lmp, "pair_coeff * * 0.0 1.0");
        lammps_command(lmp, "timestep 0.001");
        lammps_command(lmp, "velocity all set 0.0 0.0 0.0");
        lammps_set_show_error(lmp, 0);
        return lmp;
    };
    auto expect_error = [](void *lmp, const char *command, const char *expected_message) {
        lammps_command(lmp, command);
        const int has_error = lammps_has_error(lmp);
        EXPECT_EQ(has_error, 1);
        if (has_error) {
            char error_message[512];
            EXPECT_NE(lammps_get_last_error_message(lmp, error_message, sizeof(error_message)), 0);
            EXPECT_THAT(error_message, HasSubstr(expected_message));
        }
    };

    const std::array<std::array<const char *, 2>, 3> invalid_commands = {{
        {"fix guard all plumed path_integral invalid",
         "Unknown fix plumed path_integral value: invalid"},
        {"fix guard all plumed path_integral centroid",
         "Fix plumed path_integral mode requires the pimd_fix keyword"},
        {"fix guard all plumed pimd_fix fpimd",
         "Fix plumed pimd_fix requires a path_integral mode"},
    }};
    for (const auto &test_case : invalid_commands) {
        void *lmp = open_lammps();
        expect_error(lmp, test_case[0], test_case[1]);
        lammps_close(lmp);
    }

    void *lmp = open_lammps();
    lammps_command(lmp, "fix fpimd all pimd/langevin method pimd ensemble nve "
                        "integrator obabo thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");
    lammps_command(lmp, "fix guard all plumed path_integral centroid pimd_fix fpimd");
    EXPECT_EQ(lammps_has_error(lmp), 0);
    expect_error(lmp, "run 0 post no",
                 "Fix plumed path_integral centroid requires method pimd with ensemble nvt or "
                 "method nmpimd with ensemble nvt, nph, or npt");
    lammps_close(lmp);

    lmp = open_lammps();
    lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                        "integrator obabo thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");
    lammps_command(lmp, "fix guard all plumed path_integral centroid pimd_fix fpimd");
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_command(lmp, "run 0 post no");
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);

    lmp = open_lammps();
    lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nve "
                        "integrator obabo temp 1.0 fixcom no");
    lammps_command(lmp, "fix guard all plumed path_integral centroid pimd_fix fpimd");
    EXPECT_EQ(lammps_has_error(lmp), 0);
    expect_error(lmp, "run 0 post no",
                 "Fix plumed path_integral centroid requires method pimd with ensemble nvt or "
                 "method nmpimd with ensemble nvt, nph, or npt");
    lammps_close(lmp);

    for (const char *mode : {"bead_mean", "bead_density"}) {
        const std::array<const char *, 3> valid_commands = {
            "fix fpimd all pimd/langevin method nmpimd ensemble nvt integrator obabo "
            "thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no",
            "fix fpimd all pimd/langevin method nmpimd ensemble nph integrator obabo "
            "temp 1.0 iso 0.0 barostat BZP taup 1.0 fixcom no",
            "fix fpimd all pimd/langevin method nmpimd ensemble npt integrator obabo "
            "thermostat PILE_L 1234 tau 1.0 temp 1.0 iso 0.0 barostat BZP taup 1.0 "
            "fixcom no"};
        for (const char *pimd_command : valid_commands) {
            lmp = open_lammps();
            lammps_command(lmp, pimd_command);
            const std::string command =
                "fix guard all plumed path_integral " + std::string(mode) + " pimd_fix fpimd";
            lammps_command(lmp, command.c_str());
            EXPECT_EQ(lammps_has_error(lmp), 0);
            lammps_command(lmp, "run 0 post no");
            EXPECT_EQ(lammps_has_error(lmp), 0);
            lammps_close(lmp);
        }

        lmp = open_lammps();
        lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nve "
                            "integrator obabo temp 1.0 fixcom no");
        const std::string command =
            "fix guard all plumed path_integral " + std::string(mode) + " pimd_fix fpimd";
        lammps_command(lmp, command.c_str());
        EXPECT_EQ(lammps_has_error(lmp), 0);
        expect_error(lmp, "run 0 post no",
                     "Fix plumed path_integral bead modes require method pimd with ensemble nvt or "
                     "method nmpimd with ensemble nvt, nph, or npt");
        lammps_close(lmp);
    }

    for (const char *mode : {"centroid", "bead_mean", "bead_density"}) {
        lmp = open_lammps();
        lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                            "integrator obabo thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");
        lammps_command(lmp, "fix early all addforce 0.1 0.0 0.0");
        const std::string valid_command =
            "fix guard all plumed path_integral " + std::string(mode) + " pimd_fix fpimd";
        lammps_command(lmp, valid_command.c_str());
        EXPECT_EQ(lammps_has_error(lmp), 0);
        lammps_command(lmp, "run 0 post no");
        EXPECT_EQ(lammps_has_error(lmp), 0);
        lammps_close(lmp);

        lmp = open_lammps();
        lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                            "integrator obabo thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");
        lammps_command(lmp, valid_command.c_str());
        lammps_command(lmp, "fix late all addforce 0.1 0.0 0.0");
        expect_error(lmp, "run 0 post no",
                     "Fix plumed with NMPIMD path_integral modes must be defined after fix "
                     "addforce because it has a post-force callback");
        lammps_close(lmp);
    }
}

TEST(MPI, plumed_nmpimd_bead_modes_force)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);
    // Tables contain -grad(B); the beta/P Hamiltonian requires -P*grad(B).
    constexpr int beads = 4;

    const char *mean_zero_file    = "test_plumed_nmpimd_mean_zero.dat";
    const char *mean_bias_file    = "test_plumed_nmpimd_mean_bias.dat";
    const char *density_zero_file = "test_plumed_nmpimd_density_zero.dat";
    const char *density_bias_file = "test_plumed_nmpimd_density_bias.dat";
    const char *mean_zero_log     = "test_plumed_nmpimd_mean_zero.log";
    const char *mean_bias_log     = "test_plumed_nmpimd_mean_bias.log";
    const char *density_zero_log  = "test_plumed_nmpimd_density_zero.log";
    const char *density_bias_log  = "test_plumed_nmpimd_density_bias.log";
    if (me == 0) {
        std::ofstream mean_zero(mean_zero_file);
        mean_zero << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                  << "mean: ENSEMBLE ARG=d\n";
        std::ofstream mean_bias(mean_bias_file);
        mean_bias << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                  << "mean: ENSEMBLE ARG=d\n"
                  << "bias: RESTRAINT ARG=mean.d AT=6.5 KAPPA=4.0\n";
        std::ofstream density_zero(density_zero_file);
        density_zero << "d: DISTANCE ATOMS=1,2 NOPBC\n";
        std::ofstream density_bias(density_bias_file);
        density_bias << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                     << "bias: RESTRAINT ARG=d AT=6.5 KAPPA=4.0\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const char *args[] = {"LAMMPS_test", "-screen", "none", "-log",    "none", "-partition",
                          "4x1",         "-in",     "none", "-nocite", nullptr};
    char **argv        = (char **)args;
    int argc           = (sizeof(args) / sizeof(char *)) - 1;
    void *lmp          = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
    ASSERT_NE(lmp, nullptr);

    lammps_command(lmp, "variable x2 world 19.0 10.0 3.0 10.0");
    lammps_command(lmp, "variable y2 world 10.0 18.0 10.0 4.0");
    lammps_command(lmp, "units lj");
    lammps_command(lmp, "atom_style atomic");
    lammps_command(lmp, "atom_modify map array");
    lammps_command(lmp, "boundary p p p");
    lammps_command(lmp, "region box block 0.0 20.0 0.0 20.0 0.0 20.0");
    lammps_command(lmp, "create_box 1 box");
    lammps_command(lmp, "create_atoms 1 single 10.0 10.0 2.0");
    lammps_command(lmp, "create_atoms 1 single ${x2} ${y2} 2.0");
    lammps_command(lmp, "mass 1 1.0");
    lammps_command(lmp, "pair_style lj/cut 2.5");
    lammps_command(lmp, "pair_coeff * * 0.0 1.0");
    lammps_command(lmp, "timestep 0.001");
    lammps_command(lmp, "velocity all set 0.0 0.0 0.0");
    lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                        "integrator baoab thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");

    const double sqrt_two                                       = std::sqrt(2.0);
    const std::array<std::array<double, 6>, 4> mean_force_delta = {
        {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
         {sqrt_two, 0.0, 0.0, -sqrt_two, 0.0, 0.0},
         {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
         {0.0, -sqrt_two, 0.0, 0.0, sqrt_two, 0.0}}};
    const std::array<std::array<double, 6>, 4> density_force_delta = {
        {{1.0, 1.0, 0.0, -1.0, -1.0, 0.0},
         {3.0 / sqrt_two, 0.0, 0.0, -3.0 / sqrt_two, 0.0, 0.0},
         {1.0, -1.0, 0.0, -1.0, 1.0, 0.0},
         {0.0, -1.0 / sqrt_two, 0.0, 0.0, 1.0 / sqrt_two, 0.0}}};

    auto extract_forces = [&]() {
        std::array<double, 6> result{};
        auto **forces = (double **)lammps_extract_atom(lmp, "f");
        EXPECT_NE(forces, nullptr);
        if (!forces) return result;
        for (int atom = 0; atom < 2; ++atom) {
            const tagint atom_id = atom + 1;
            const int index      = lammps_map_atom(lmp, &atom_id);
            EXPECT_GE(index, 0);
            if (index < 0) continue;
            for (int dimension = 0; dimension < 3; ++dimension)
                result[3 * atom + dimension] = forces[index][dimension];
        }
        return result;
    };
    auto check_mode = [&](const std::string &zero_command, const std::string &bias_command,
                          double expected_bias, const std::array<double, 6> &expected_force_delta) {
        lammps_command(lmp, zero_command.c_str());
        lammps_command(lmp, "run 0 post no");
        const auto zero_forces = extract_forces();
        lammps_command(lmp, "unfix zero");

        lammps_command(lmp, bias_command.c_str());
        lammps_command(lmp, "run 0 post no");
        const auto biased_forces = extract_forces();
        auto *bias =
            (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
        ASSERT_NE(bias, nullptr);
        EXPECT_NEAR(*bias, me == 0 ? expected_bias : 0.0, 1.0e-12);
        lammps_free(bias);
        for (std::size_t i = 0; i < expected_force_delta.size(); ++i)
            EXPECT_NEAR(biased_forces[i] - zero_forces[i], beads * expected_force_delta[i], 1.0e-12)
                << i;
        lammps_command(lmp, "unfix bias");
    };

    check_mode("fix zero all plumed plumedfile " + std::string(mean_zero_file) + " outfile " +
                   mean_zero_log + " path_integral bead_mean pimd_fix fpimd",
               "fix bias all plumed plumedfile " + std::string(mean_bias_file) + " outfile " +
                   mean_bias_log + " path_integral bead_mean pimd_fix fpimd",
               2.0, mean_force_delta[me]);
    check_mode("fix zero all plumed plumedfile " + std::string(density_zero_file) + " outfile " +
                   density_zero_log + " path_integral bead_density pimd_fix fpimd",
               "fix bias all plumed plumedfile " + std::string(density_bias_file) + " outfile " +
                   density_bias_log + " path_integral bead_density pimd_fix fpimd",
               4.5, density_force_delta[me]);

    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        for (const char *file :
             {mean_zero_file, mean_bias_file, density_zero_file, density_bias_file})
            std::remove(file);
        for (int bead = 0; bead < nprocs; ++bead) {
            for (const char *log :
                 {mean_zero_log, mean_bias_log, density_zero_log, density_bias_log})
                std::remove((std::string(log) + "." + std::to_string(bead)).c_str());
        }
    }
}

TEST(MPI, plumed_nmpimd_multirank_bead_modes)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);
    // Tables contain -grad(B); the beta/P Hamiltonian requires -P*grad(B).
    constexpr int beads = 2;

    const char *mean_zero_file    = "test_plumed_nmpimd_multirank_mean_zero.dat";
    const char *mean_bias_file    = "test_plumed_nmpimd_multirank_mean_bias.dat";
    const char *density_zero_file = "test_plumed_nmpimd_multirank_density_zero.dat";
    const char *density_bias_file = "test_plumed_nmpimd_multirank_density_bias.dat";
    const char *mean_zero_log     = "test_plumed_nmpimd_multirank_mean_zero.log";
    const char *mean_bias_log     = "test_plumed_nmpimd_multirank_mean_bias.log";
    const char *density_zero_log  = "test_plumed_nmpimd_multirank_density_zero.log";
    const char *density_bias_log  = "test_plumed_nmpimd_multirank_density_bias.log";
    if (me == 0) {
        std::ofstream mean_zero(mean_zero_file);
        mean_zero << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                  << "mean: ENSEMBLE ARG=d\n";
        std::ofstream mean_bias(mean_bias_file);
        mean_bias << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                  << "mean: ENSEMBLE ARG=d\n"
                  << "bias: RESTRAINT ARG=mean.d AT=7.5 KAPPA=4.0\n";
        std::ofstream density_zero(density_zero_file);
        density_zero << "d: DISTANCE ATOMS=1,2 NOPBC\n";
        std::ofstream density_bias(density_bias_file);
        density_bias << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                     << "bias: RESTRAINT ARG=d AT=7.5 KAPPA=4.0\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const int bead                                              = me / 2;
    const double inverse_sqrt_two                               = 1.0 / std::sqrt(2.0);
    const std::array<std::array<double, 6>, 2> mean_force_delta = {
        {{inverse_sqrt_two, inverse_sqrt_two, 0.0, -inverse_sqrt_two, -inverse_sqrt_two, 0.0},
         {inverse_sqrt_two, -inverse_sqrt_two, 0.0, -inverse_sqrt_two, inverse_sqrt_two, 0.0}}};
    const std::array<std::array<double, 6>, 2> density_force_delta = {
        {{3.0 * inverse_sqrt_two, -inverse_sqrt_two, 0.0, -3.0 * inverse_sqrt_two, inverse_sqrt_two,
          0.0},
         {3.0 * inverse_sqrt_two, inverse_sqrt_two, 0.0, -3.0 * inverse_sqrt_two, -inverse_sqrt_two,
          0.0}}};

    auto run_case = [&](const char *plumed_file, const char *plumed_log, const char *mode) {
        std::array<double, 7> result{};
        void *lmp = open_multirank_partition();
        EXPECT_NE(lmp, nullptr);
        if (!lmp) return result;

        create_multirank_two_atom_system(lmp);
        auto *nlocal = (int *)lammps_extract_global(lmp, "nlocal");
        EXPECT_NE(nlocal, nullptr);
        int zero_atom_ranks = nlocal && *nlocal == 0;
        MPI_Allreduce(MPI_IN_PLACE, &zero_atom_ranks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        EXPECT_GT(zero_atom_ranks, 0);
        lammps_command(lmp, "group first id 1");
        lammps_command(lmp, "group second id 2");
        lammps_command(lmp, "compute first_force first reduce sum fx fy fz");
        lammps_command(lmp, "compute second_force second reduce sum fx fy fz");
        lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                            "integrator baoab thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");
        const std::string fix_command = "fix bias all plumed plumedfile " +
                                        std::string(plumed_file) + " outfile " + plumed_log +
                                        " path_integral " + mode + " pimd_fix fpimd";
        lammps_command(lmp, fix_command.c_str());
        lammps_command(lmp, "run 0 post no");

        auto *first =
            (double *)lammps_extract_compute(lmp, "first_force", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
        auto *second = (double *)lammps_extract_compute(lmp, "second_force", LMP_STYLE_GLOBAL,
                                                        LMP_TYPE_VECTOR);
        EXPECT_NE(first, nullptr);
        EXPECT_NE(second, nullptr);
        if (first && second) {
            for (int dimension = 0; dimension < 3; ++dimension) {
                result[dimension]     = first[dimension];
                result[3 + dimension] = second[dimension];
            }
        }
        auto *bias =
            (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
        EXPECT_NE(bias, nullptr);
        if (bias) {
            result[6] = *bias;
            lammps_free(bias);
        }
        EXPECT_EQ(lammps_has_error(lmp), 0);
        lammps_close(lmp);
        return result;
    };

    const auto mean_zero = run_case(mean_zero_file, mean_zero_log, "bead_mean");
    const auto mean_bias = run_case(mean_bias_file, mean_bias_log, "bead_mean");
    EXPECT_NEAR(mean_zero[6], 0.0, 1.0e-12);
    EXPECT_NEAR(mean_bias[6], bead == 0 ? 0.5 : 0.0, 1.0e-12);
    for (std::size_t i = 0; i < mean_force_delta[bead].size(); ++i)
        EXPECT_NEAR(mean_bias[i] - mean_zero[i], beads * mean_force_delta[bead][i], 1.0e-12) << i;

    const auto density_zero = run_case(density_zero_file, density_zero_log, "bead_density");
    const auto density_bias = run_case(density_bias_file, density_bias_log, "bead_density");
    EXPECT_NEAR(density_zero[6], 0.0, 1.0e-12);
    EXPECT_NEAR(density_bias[6], bead == 0 ? 2.5 : 0.0, 1.0e-12);
    for (std::size_t i = 0; i < density_force_delta[bead].size(); ++i)
        EXPECT_NEAR(density_bias[i] - density_zero[i], beads * density_force_delta[bead][i],
                    1.0e-12)
            << i;

    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        for (const char *file :
             {mean_zero_file, mean_bias_file, density_zero_file, density_bias_file})
            std::remove(file);
        for (int bead_index = 0; bead_index < 2; ++bead_index) {
            for (const char *log :
                 {mean_zero_log, mean_bias_log, density_zero_log, density_bias_log})
                std::remove((std::string(log) + "." + std::to_string(bead_index)).c_str());
        }
    }
}

TEST(MPI, plumed_nmpimd_bead_modes_pbc)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);
    // Tables contain -grad(B); the beta/P Hamiltonian requires -P*grad(B).
    constexpr int beads = 2;

    const char *mean_zero_file    = "test_plumed_nmpimd_pbc_mean_zero.dat";
    const char *mean_bias_file    = "test_plumed_nmpimd_pbc_mean_bias.dat";
    const char *density_zero_file = "test_plumed_nmpimd_pbc_density_zero.dat";
    const char *density_bias_file = "test_plumed_nmpimd_pbc_density_bias.dat";
    if (me == 0) {
        std::ofstream mean_zero(mean_zero_file);
        mean_zero << "d: DISTANCE ATOMS=1,2\n"
                  << "mean: ENSEMBLE ARG=d\n";
        std::ofstream mean_bias(mean_bias_file);
        mean_bias << "d: DISTANCE ATOMS=1,2\n"
                  << "mean: ENSEMBLE ARG=d\n"
                  << "bias: RESTRAINT ARG=mean.d AT=8.0 KAPPA=4.0\n";
        std::ofstream density_zero(density_zero_file);
        density_zero << "d: DISTANCE ATOMS=1,2\n";
        std::ofstream density_bias(density_bias_file);
        density_bias << "d: DISTANCE ATOMS=1,2\n"
                     << "bias: RESTRAINT ARG=d AT=8.0 KAPPA=4.0\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const int bead                                                  = me / 2;
    const std::array<std::array<double, 6>, 2> expected_force_delta = {
        {{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
         {{2.0 * std::sqrt(2.0), 0.0, 0.0, -2.0 * std::sqrt(2.0), 0.0, 0.0}}}};
    auto run_case = [&](const char *plumed_file, const std::string &plumed_log, const char *mode) {
        std::array<double, 7> result{};
        void *lmp = open_multirank_partition();
        EXPECT_NE(lmp, nullptr);
        if (!lmp) return result;

        create_multirank_two_atom_system(lmp);
        lammps_command(lmp, "variable wrapped_x world 19.0 1.0");
        lammps_command(lmp, "set atom 2 x ${wrapped_x} y 10.0 z 2.0");
        lammps_command(lmp, "group first id 1");
        lammps_command(lmp, "group second id 2");
        lammps_command(lmp, "compute first_force first reduce sum fx fy fz");
        lammps_command(lmp, "compute second_force second reduce sum fx fy fz");
        lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                            "integrator baoab thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");
        const std::string fix_command = "fix bias all plumed plumedfile " +
                                        std::string(plumed_file) + " outfile " + plumed_log +
                                        " path_integral " + mode + " pimd_fix fpimd";
        lammps_command(lmp, fix_command.c_str());
        lammps_command(lmp, "run 0 post no");

        auto *first =
            (double *)lammps_extract_compute(lmp, "first_force", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
        auto *second = (double *)lammps_extract_compute(lmp, "second_force", LMP_STYLE_GLOBAL,
                                                        LMP_TYPE_VECTOR);
        EXPECT_NE(first, nullptr);
        EXPECT_NE(second, nullptr);
        if (first && second) {
            for (int dimension = 0; dimension < 3; ++dimension) {
                result[dimension]     = first[dimension];
                result[3 + dimension] = second[dimension];
            }
        }
        auto *bias =
            (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
        EXPECT_NE(bias, nullptr);
        if (bias) {
            result[6] = *bias;
            lammps_free(bias);
        }
        EXPECT_EQ(lammps_has_error(lmp), 0);
        lammps_close(lmp);
        return result;
    };

    for (const char *mode : {"bead_mean", "bead_density"}) {
        const bool mean          = std::string(mode) == "bead_mean";
        const char *zero_file    = mean ? mean_zero_file : density_zero_file;
        const char *bias_file    = mean ? mean_bias_file : density_bias_file;
        const std::string prefix = "test_plumed_nmpimd_pbc_" + std::string(mode);
        const auto zero          = run_case(zero_file, prefix + "_zero.log", mode);
        const auto biased        = run_case(bias_file, prefix + "_bias.log", mode);
        EXPECT_NEAR(zero[6], 0.0, 1.0e-12);
        EXPECT_NEAR(biased[6], bead == 0 ? 2.0 : 0.0, 1.0e-12);
        for (std::size_t i = 0; i < expected_force_delta[bead].size(); ++i)
            EXPECT_NEAR(biased[i] - zero[i], beads * expected_force_delta[bead][i], 1.0e-12) << i;

        MPI_Barrier(MPI_COMM_WORLD);
        if (me == 0) {
            for (int bead_index = 0; bead_index < 2; ++bead_index) {
                std::remove((prefix + "_zero.log." + std::to_string(bead_index)).c_str());
                std::remove((prefix + "_bias.log." + std::to_string(bead_index)).c_str());
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        for (const char *file :
             {mean_zero_file, mean_bias_file, density_zero_file, density_bias_file})
            std::remove(file);
    }
}

TEST(MPI, plumed_centroid_multirank_force_modes)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    const char *plumed_file = "test_plumed_nmpimd_centroid_force.dat";
    const char *plumed_log  = "test_plumed_nmpimd_centroid_force.log";
    const char *direct_log  = "test_plumed_pimd_centroid_multirank_force.log";
    if (me == 0) {
        std::ofstream input(plumed_file);
        input << "d: DISTANCE ATOMS=1,2 NOPBC\n"
              << "bias: RESTRAINT ARG=d AT=0.0 KAPPA=4.0\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const std::array<std::string, 2> normal_mode_commands = {
        "fix fpimd all pimd/langevin method nmpimd ensemble nvt integrator baoab "
        "thermostat PILE_L 1357 tau 1.0 temp 1.0 fixcom no",
        "fix fpimd all pimd/langevin method nmpimd ensemble nph integrator obabo "
        "iso 0.0 barostat BZP taup 1.0 fixcom no"};
    for (std::size_t mode = 0; mode < normal_mode_commands.size(); ++mode) {
        void *lmp = open_multirank_partition();
        ASSERT_NE(lmp, nullptr);
        EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
        create_multirank_two_atom_system(lmp);
        lammps_command(lmp, "group first id 1");
        lammps_command(lmp, "group second id 2");
        lammps_command(lmp, "compute first_force first reduce sum fx fy fz");
        lammps_command(lmp, "compute second_force second reduce sum fx fy fz");
        lammps_command(lmp, normal_mode_commands[mode].c_str());
        const std::string fix_command =
            "fix bias all plumed plumedfile " + std::string(plumed_file) + " outfile " +
            plumed_log + "." + std::to_string(mode) + " path_integral centroid pimd_fix fpimd";
        lammps_command(lmp, fix_command.c_str());
        lammps_command(lmp, "run 0 post no");

        auto *first =
            (double *)lammps_extract_compute(lmp, "first_force", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
        auto *second = (double *)lammps_extract_compute(lmp, "second_force", LMP_STYLE_GLOBAL,
                                                        LMP_TYPE_VECTOR);
        ASSERT_NE(first, nullptr);
        ASSERT_NE(second, nullptr);
        const double mode_scale = me / 2 == 0 ? std::sqrt(2.0) : 0.0;
        EXPECT_NEAR(first[0], 18.0 * mode_scale, 1.0e-12);
        EXPECT_NEAR(first[1], 14.0 * mode_scale, 1.0e-12);
        EXPECT_NEAR(first[2], 0.0, 1.0e-12);
        EXPECT_NEAR(second[0], -18.0 * mode_scale, 1.0e-12);
        EXPECT_NEAR(second[1], -14.0 * mode_scale, 1.0e-12);
        EXPECT_NEAR(second[2], 0.0, 1.0e-12);

        auto *bias =
            (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
        ASSERT_NE(bias, nullptr);
        EXPECT_NEAR(*bias, me / 2 == 0 ? 65.0 : 0.0, 1.0e-12);
        lammps_free(bias);
        EXPECT_EQ(lammps_has_error(lmp), 0);
        lammps_close(lmp);
    }

    void *lmp = open_multirank_partition();
    ASSERT_NE(lmp, nullptr);
    create_multirank_two_atom_system(lmp);
    lammps_command(lmp, "set atom 2 x 19.0 y 10.0 z 2.0");
    lammps_command(lmp, "group first id 1");
    lammps_command(lmp, "group second id 2");
    lammps_command(lmp, "compute first_force first reduce sum fx fy fz");
    lammps_command(lmp, "compute second_force second reduce sum fx fy fz");
    lammps_command(lmp, "fix fpimd all pimd/langevin method pimd ensemble nvt "
                        "integrator obabo thermostat PILE_L 2468 tau 1.0 temp 1.0 fixcom no");
    const std::string direct_command = "fix bias all plumed plumedfile " +
                                       std::string(plumed_file) + " outfile " + direct_log +
                                       " path_integral centroid pimd_fix fpimd";
    lammps_command(lmp, direct_command.c_str());
    lammps_command(lmp, "run 0 post no");

    auto *first =
        (double *)lammps_extract_compute(lmp, "first_force", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
    auto *second =
        (double *)lammps_extract_compute(lmp, "second_force", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
    ASSERT_NE(first, nullptr);
    ASSERT_NE(second, nullptr);
    EXPECT_NEAR(first[0], 36.0, 1.0e-12);
    EXPECT_NEAR(first[1], 0.0, 1.0e-12);
    EXPECT_NEAR(first[2], 0.0, 1.0e-12);
    EXPECT_NEAR(second[0], -36.0, 1.0e-12);
    EXPECT_NEAR(second[1], 0.0, 1.0e-12);
    EXPECT_NEAR(second[2], 0.0, 1.0e-12);
    auto *bias =
        (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
    ASSERT_NE(bias, nullptr);
    EXPECT_NEAR(*bias, me / 2 == 0 ? 162.0 : 0.0, 1.0e-12);
    lammps_free(bias);
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);

    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        std::remove(plumed_file);
        std::remove((std::string(plumed_log) + ".0").c_str());
        std::remove((std::string(plumed_log) + ".1").c_str());
        std::remove(direct_log);
    }
}

TEST(MPI, plumed_nmpimd_centroid_npt_volume)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);
    // B(V)=-2V gives physical pressure +2 and dynamical energy P*B.
    constexpr double bias_pressure = 2.0;

    const char *zero_file = "test_plumed_nmpimd_volume_zero.dat";
    const char *bias_file = "test_plumed_nmpimd_volume_bias.dat";
    const char *zero_log  = "test_plumed_nmpimd_volume_zero.log";
    const char *bias_log  = "test_plumed_nmpimd_volume_bias.log";
    if (me == 0) {
        std::ofstream zero(zero_file);
        zero << "v: VOLUME\n";
        std::ofstream biased(bias_file);
        biased << "v: VOLUME\n"
               << "bias: RESTRAINT ARG=v AT=0.0 KAPPA=0.0 SLOPE=-2.0\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const auto zero   = run_centroid_npt_volume_leg(zero_file, zero_log, 0.0);
    const auto biased = run_centroid_npt_volume_leg(bias_file, bias_log, bias_pressure);

    EXPECT_NEAR(biased.p_cv - zero.p_cv, bias_pressure, 1.0e-12);
    EXPECT_NEAR(biased.pe_bead - zero.pe_bead, me / 2 == 0 ? -32000.0 : 0.0, 1.0e-10);
    EXPECT_NEAR(biased.tote - zero.tote, -32000.0, 1.0e-10);
    for (int i = 0; i < 3; ++i)
        EXPECT_NEAR(biased.pressure[i] - zero.pressure[i], me / 2 == 0 ? 2.0 * bias_pressure : 0.0,
                    1.0e-12);
    for (int i = 3; i < 6; ++i)
        EXPECT_NEAR(biased.pressure[i] - zero.pressure[i], 0.0, 1.0e-12);
    EXPECT_NEAR(zero.bias, 0.0, 1.0e-12);
    EXPECT_NEAR(biased.bias, me / 2 == 0 ? -16000.0 : 0.0, 1.0e-10);
    EXPECT_NEAR(biased.barostat_velocity, zero.barostat_velocity, 1.0e-12);
    EXPECT_NEAR(biased.volume, zero.volume, 1.0e-10);
    EXPECT_GT(std::fabs(zero.barostat_velocity), 1.0e-8);
    EXPECT_GT(std::fabs(zero.volume - 8000.0), 1.0e-10);

    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        std::remove(zero_file);
        std::remove(bias_file);
        std::remove(zero_log);
        std::remove(bias_log);
    }
}

TEST(MPI, plumed_nmpimd_bead_modes_bzp_volume)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);
    // B(V)=-2V gives physical pressure +2 and dynamical energy P*B.
    constexpr double bias_pressure = 2.0;

    const char *mean_zero_file    = "test_plumed_nmpimd_bzp_mean_zero.dat";
    const char *mean_bias_file    = "test_plumed_nmpimd_bzp_mean_bias.dat";
    const char *density_zero_file = "test_plumed_nmpimd_bzp_density_zero.dat";
    const char *density_bias_file = "test_plumed_nmpimd_bzp_density_bias.dat";
    if (me == 0) {
        std::ofstream mean_zero(mean_zero_file);
        mean_zero << "v: VOLUME\n"
                  << "mean: ENSEMBLE ARG=v\n";
        std::ofstream mean_bias(mean_bias_file);
        mean_bias << "v: VOLUME\n"
                  << "mean: ENSEMBLE ARG=v\n"
                  << "bias: RESTRAINT ARG=mean.v AT=0.0 KAPPA=0.0 SLOPE=-2.0\n";
        std::ofstream density_zero(density_zero_file);
        density_zero << "v: VOLUME\n";
        std::ofstream density_bias(density_bias_file);
        density_bias << "v: VOLUME\n"
                     << "bias: RESTRAINT ARG=v AT=0.0 KAPPA=0.0 SLOPE=-2.0\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const int bead = me / 2;
    for (const char *ensemble : {"nph", "npt"}) {
        for (const char *mode : {"bead_mean", "bead_density"}) {
            const bool mean       = std::string(mode) == "bead_mean";
            const char *zero_file = mean ? mean_zero_file : density_zero_file;
            const char *bias_file = mean ? mean_bias_file : density_bias_file;
            const std::string prefix =
                "test_plumed_nmpimd_bzp_" + std::string(ensemble) + "_" + mode;
            const auto zero = run_nmpimd_bead_mode_volume_leg(
                zero_file, (prefix + "_zero.log").c_str(), mode, ensemble, 0.0);
            const auto biased = run_nmpimd_bead_mode_volume_leg(
                bias_file, (prefix + "_bias.log").c_str(), mode, ensemble, bias_pressure);

            EXPECT_NEAR(biased.p_cv - zero.p_cv, bias_pressure, 1.0e-12);
            EXPECT_NEAR(biased.pe_bead - zero.pe_bead, bead == 0 ? -32000.0 : 0.0, 1.0e-10);
            EXPECT_NEAR(biased.tote - zero.tote, -32000.0, 1.0e-10);
            for (int i = 0; i < 3; ++i)
                EXPECT_NEAR(biased.pressure[i] - zero.pressure[i], bias_pressure, 1.0e-12);
            for (int i = 3; i < 6; ++i)
                EXPECT_NEAR(biased.pressure[i] - zero.pressure[i], 0.0, 1.0e-12);
            EXPECT_NEAR(zero.bias, 0.0, 1.0e-12);
            EXPECT_NEAR(biased.bias, bead == 0 ? -16000.0 : 0.0, 1.0e-10);
            EXPECT_NEAR(biased.barostat_velocity, zero.barostat_velocity, 1.0e-12);
            EXPECT_NEAR(biased.volume, zero.volume, 1.0e-10);
            EXPECT_GT(std::fabs(zero.barostat_velocity), 1.0e-8);
            EXPECT_GT(std::fabs(zero.volume - 8000.0), 1.0e-10);

            MPI_Barrier(MPI_COMM_WORLD);
            if (me == 0) {
                for (int bead_index = 0; bead_index < 2; ++bead_index) {
                    std::remove((prefix + "_zero.log." + std::to_string(bead_index)).c_str());
                    std::remove((prefix + "_bias.log." + std::to_string(bead_index)).c_str());
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        for (const char *file :
             {mean_zero_file, mean_bias_file, density_zero_file, density_bias_file})
            std::remove(file);
    }
}

TEST(MPI, plumed_nmpimd_centroid_nvt_virial)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);
    // B(V)=-2V gives physical pressure +2 and dynamical energy P*B.
    constexpr double bias_pressure = 2.0;

    const char *zero_file = "test_plumed_nmpimd_nvt_volume_zero.dat";
    const char *bias_file = "test_plumed_nmpimd_nvt_volume_bias.dat";
    const char *zero_log  = "test_plumed_nmpimd_nvt_volume_zero.log";
    const char *bias_log  = "test_plumed_nmpimd_nvt_volume_bias.log";
    if (me == 0) {
        std::ofstream zero(zero_file);
        zero << "v: VOLUME\n";
        std::ofstream biased(bias_file);
        biased << "v: VOLUME\n"
               << "bias: RESTRAINT ARG=v AT=0.0 KAPPA=0.0 SLOPE=-2.0\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const auto zero   = run_centroid_nvt_volume_leg(zero_file, zero_log);
    const auto biased = run_centroid_nvt_volume_leg(bias_file, bias_log);

    EXPECT_NEAR(biased.p_cv - zero.p_cv, bias_pressure, 1.0e-12);
    for (int i = 0; i < 3; ++i)
        EXPECT_NEAR(biased.pressure[i] - zero.pressure[i], me / 2 == 0 ? 2.0 * bias_pressure : 0.0,
                    1.0e-12);
    for (int i = 3; i < 6; ++i)
        EXPECT_NEAR(biased.pressure[i] - zero.pressure[i], 0.0, 1.0e-12);
    EXPECT_NEAR(zero.bias, 0.0, 1.0e-12);
    EXPECT_NEAR(biased.bias, me / 2 == 0 ? -16000.0 : 0.0, 1.0e-10);
    EXPECT_EQ(zero.volume, 8000.0);
    EXPECT_EQ(biased.volume, zero.volume);

    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        std::remove(zero_file);
        std::remove(bias_file);
        std::remove(zero_log);
        std::remove(bias_log);
    }
}

TEST(MPI, plumed_nmpimd_centroid_nvt_wrapped_four_bead)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    const char *plumed_file = "test_plumed_nmpimd_centroid_nvt_wrapped.dat";
    const char *plumed_log  = "test_plumed_nmpimd_centroid_nvt_wrapped.log";
    if (me == 0) {
        std::ofstream input(plumed_file);
        input << "d: DISTANCE ATOMS=1,2 NOPBC\n"
              << "bias: RESTRAINT ARG=d AT=0.0 KAPPA=4.0\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const char *args[] = {"LAMMPS_test", "-screen", "none", "-log",    "none", "-partition",
                          "4x1",         "-in",     "none", "-nocite", nullptr};
    char **argv        = (char **)args;
    int argc           = (sizeof(args) / sizeof(char *)) - 1;
    void *lmp          = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
    ASSERT_NE(lmp, nullptr);
    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 1);

    lammps_command(lmp, "variable x2 world 19.0 1.0 3.0 9.0");
    lammps_command(lmp, "variable ix world 0 1 0 0");
    lammps_command(lmp, "units lj");
    lammps_command(lmp, "atom_style atomic");
    lammps_command(lmp, "atom_modify map array");
    lammps_command(lmp, "boundary p p p");
    lammps_command(lmp, "region box block 0.0 20.0 0.0 20.0 0.0 20.0");
    lammps_command(lmp, "create_box 1 box");
    lammps_command(lmp, "create_atoms 1 single 10.0 10.0 2.0");
    lammps_command(lmp, "create_atoms 1 single ${x2} 10.0 2.0");
    lammps_command(lmp, "set atom 2 image ${ix} 0 0");
    lammps_command(lmp, "mass 1 1.0");
    lammps_command(lmp, "pair_style lj/cut 2.5");
    lammps_command(lmp, "pair_coeff * * 0.0 1.0");
    lammps_command(lmp, "timestep 0.001");
    lammps_command(lmp, "velocity all set 0.0 0.0 0.0");
    lammps_command(lmp, "group first id 1");
    lammps_command(lmp, "group second id 2");
    lammps_command(lmp, "compute first_force first reduce sum fx fy fz");
    lammps_command(lmp, "compute second_force second reduce sum fx fy fz");
    lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                        "integrator baoab thermostat PILE_L 2468 tau 1.0 temp 1.0 fixcom no");
    const std::string fix_command = "fix bias all plumed plumedfile " + std::string(plumed_file) +
                                    " outfile " + plumed_log +
                                    " path_integral centroid pimd_fix fpimd";
    lammps_command(lmp, fix_command.c_str());
    lammps_command(lmp, "run 0 post no");

    auto *first =
        (double *)lammps_extract_compute(lmp, "first_force", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
    auto *second =
        (double *)lammps_extract_compute(lmp, "second_force", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
    ASSERT_NE(first, nullptr);
    ASSERT_NE(second, nullptr);
    const double mode_scale = me == 0 ? std::sqrt(4.0) : 0.0;
    EXPECT_NEAR(first[0], 12.0 * mode_scale, 1.0e-12);
    EXPECT_NEAR(first[1], 0.0, 1.0e-12);
    EXPECT_NEAR(first[2], 0.0, 1.0e-12);
    EXPECT_NEAR(second[0], -12.0 * mode_scale, 1.0e-12);
    EXPECT_NEAR(second[1], 0.0, 1.0e-12);
    EXPECT_NEAR(second[2], 0.0, 1.0e-12);

    auto *bias =
        (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
    ASSERT_NE(bias, nullptr);
    EXPECT_NEAR(*bias, me == 0 ? 18.0 : 0.0, 1.0e-12);
    lammps_free(bias);

    double boxlo[3], boxhi[3], xy, yz, xz;
    int periodicity[3], boxflag;
    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, periodicity, &boxflag);
    const double volume_before =
        (boxhi[0] - boxlo[0]) * (boxhi[1] - boxlo[1]) * (boxhi[2] - boxlo[2]);
    lammps_command(lmp, "run 1 post no");
    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, periodicity, &boxflag);
    const double volume_after =
        (boxhi[0] - boxlo[0]) * (boxhi[1] - boxlo[1]) * (boxhi[2] - boxlo[2]);
    EXPECT_EQ(volume_before, 8000.0);
    EXPECT_EQ(volume_after, volume_before);
    EXPECT_EQ(lammps_has_error(lmp), 0);
    lammps_close(lmp);

    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        std::remove(plumed_file);
        std::remove(plumed_log);
    }
}

TEST(MPI, plumed_nmpimd_centroid_nvt_restart_continuity)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    const char *initial_file = "test_plumed_nmpimd_nvt_restart_initial.dat";
    const char *restart_file = "test_plumed_nmpimd_nvt_restart_continue.dat";
    const std::string graph  = "d: DISTANCE ATOMS=1,2 NOPBC\n"
                               "bias: RESTRAINT ARG=d AT=0.0 KAPPA=0.01\n";
    if (me == 0) {
        std::ofstream initial(initial_file);
        initial << graph;
        std::ofstream restarted(restart_file);
        restarted << "RESTART\n" << graph;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    remove_centroid_nmpimd_nvt_restart_files();
    const auto continuous = run_centroid_nmpimd_nvt_segments(
        6, 0, false, initial_file, restart_file, "test_plumed_nmpimd_nvt_continuous.log", "");
    const auto segmented = run_centroid_nmpimd_nvt_segments(
        3, 3, false, initial_file, restart_file, "test_plumed_nmpimd_nvt_segmented.log", "");
    const auto restarted = run_centroid_nmpimd_nvt_segments(
        3, 3, true, initial_file, restart_file, "test_plumed_nmpimd_nvt_restart_first.log",
        "test_plumed_nmpimd_nvt_restart_second.log");
    remove_centroid_nmpimd_nvt_restart_files();

    for (std::size_t i = 0; i < continuous.atoms.size(); ++i) {
        EXPECT_TRUE(std::isfinite(continuous.atoms[i]));
        EXPECT_NEAR(segmented.atoms[i], continuous.atoms[i], 1.0e-14) << i;
        EXPECT_NEAR(restarted.atoms[i], continuous.atoms[i], 1.0e-14) << i;
    }
    for (std::size_t i = 0; i < continuous.pimd.size(); ++i) {
        EXPECT_TRUE(std::isfinite(continuous.pimd[i]));
        const double tolerance =
            1.0e-14 * std::max({1.0, std::abs(continuous.pimd[i]), std::abs(restarted.pimd[i])});
        EXPECT_NEAR(segmented.pimd[i], continuous.pimd[i], tolerance) << i;
        EXPECT_NEAR(restarted.pimd[i], continuous.pimd[i], tolerance) << i;
    }
    for (std::size_t i = 0; i < continuous.bias.size(); ++i) {
        EXPECT_TRUE(std::isfinite(continuous.bias[i]));
        EXPECT_NEAR(segmented.bias[i], continuous.bias[i], 1.0e-14) << i;
        EXPECT_NEAR(restarted.bias[i], continuous.bias[i], 1.0e-14) << i;
    }
    EXPECT_GT(std::abs(continuous.bias[0]), 0.0);
    EXPECT_NEAR(continuous.bias[0], continuous.bias[1], 1.0e-14);
    EXPECT_NEAR(continuous.bias[2], 0.0, 1.0e-14);
    EXPECT_NEAR(continuous.bias[3], 0.0, 1.0e-14);
    EXPECT_EQ(continuous.volume, 8000.0);
    EXPECT_EQ(segmented.volume, continuous.volume);
    EXPECT_EQ(restarted.volume, continuous.volume);

    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        std::remove(initial_file);
        std::remove(restart_file);
        for (const char *log :
             {"test_plumed_nmpimd_nvt_continuous.log", "test_plumed_nmpimd_nvt_segmented.log",
              "test_plumed_nmpimd_nvt_restart_first.log",
              "test_plumed_nmpimd_nvt_restart_second.log"})
            std::remove(log);
    }
}

TEST(MPI, plumed_nmpimd_opes_langevin_restart_continuity)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    auto counter = [](const std::string &file) {
        std::ifstream input(file);
        std::string line;
        while (std::getline(input, line)) {
            if (line.rfind("#! SET counter", 0) != 0) continue;
            std::istringstream row(line);
            std::string marker, set, key;
            int value = -1;
            row >> marker >> set >> key >> value;
            return value;
        }
        return -1;
    };
    for (const char *mode : {"centroid", "bead_mean", "bead_density"}) {
        SCOPED_TRACE(mode);
        for (int seam : {3, 4}) {
            SCOPED_TRACE(seam);
            const std::string prefix =
                "test_nmpimd_opes_rng_" + std::string(mode) + "_" + std::to_string(seam);
            const bool mean          = std::string(mode) == "bead_mean";
            const bool density       = std::string(mode) == "bead_density";
            const std::string binary = prefix + "_binary";
            auto write_input         = [&](const std::string &leg, bool restart) {
                const std::string stem = prefix + "_" + leg;
                std::ofstream input(stem + ".dat");
                if (restart) input << "RESTART\n";
                input << "d: DISTANCE ATOMS=1,2 NOPBC\n";
                if (mean) input << "mean: ENSEMBLE ARG=d\n";
                input << "bias: OPES_METAD ARG=" << (mean ? "mean.d" : "d")
                      << " PACE=2 BARRIER=4 TEMP=1 SIGMA=0.5 FIXED_SIGMA FILE=" << stem
                      << ".kernels FMT=%24.17g STATE_WFILE=" << stem
                      << ".state STATE_WSTRIDE=" << seam;
                if (density) input << " WALKERS_MPI";
                if (restart) input << " STATE_RFILE=" << prefix << "_i.state";
                input << (restart ? " RESTART=YES\n" : " RESTART=NO\n");
            };
            if (me == 0) {
                write_input("c", false);
                write_input("s", false);
                write_input("i", false);
                write_input("r", true);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            const auto continuous = run_nmpimd_bead_mode_segments(
                mode, binary.c_str(), 2 * seam, 0, false, (prefix + "_c.dat").c_str(), "",
                (prefix + "_c.log").c_str(), "");
            const auto segmented = run_nmpimd_bead_mode_segments(
                mode, binary.c_str(), seam, seam, false, (prefix + "_s.dat").c_str(), "",
                (prefix + "_s.log").c_str(), "");
            const auto restarted = run_nmpimd_bead_mode_segments(
                mode, binary.c_str(), seam, seam, true, (prefix + "_i.dat").c_str(),
                (prefix + "_r.dat").c_str(), (prefix + "_i.log").c_str(),
                (prefix + "_r.log").c_str());
            // OPES deposits at end_of_step.  At a deposition-aligned restart seam, the
            // restarted run recomputes the current force from the updated state, while an
            // uninterrupted run retains the pre-deposition force until the next force step.
            // Non-deposition seams remain near-bitwise; the aligned-seam limit is still four
            // orders tighter than the duplicated-update error this test is designed to catch.
            const double restart_tolerance = seam % 2 == 0 ? 1.0e-10 : 1.0e-13;
            for (std::size_t i = 0; i < continuous.atoms.size(); ++i) {
                EXPECT_TRUE(std::isfinite(continuous.atoms[i]));
                EXPECT_TRUE(std::isfinite(segmented.atoms[i]));
                EXPECT_TRUE(std::isfinite(restarted.atoms[i]));
                EXPECT_NEAR(segmented.atoms[i], continuous.atoms[i], 1.0e-14) << i;
                EXPECT_NEAR(restarted.atoms[i], segmented.atoms[i], restart_tolerance) << i;
            }
            for (std::size_t i = 0; i < continuous.pimd.size(); ++i) {
                EXPECT_TRUE(std::isfinite(continuous.pimd[i]));
                EXPECT_TRUE(std::isfinite(segmented.pimd[i]));
                EXPECT_TRUE(std::isfinite(restarted.pimd[i]));
                const double scale =
                    std::max({1.0, std::abs(continuous.pimd[i]), std::abs(segmented.pimd[i]),
                              std::abs(restarted.pimd[i])});
                EXPECT_NEAR(segmented.pimd[i], continuous.pimd[i], 1.0e-14 * scale) << i;
                EXPECT_NEAR(restarted.pimd[i], segmented.pimd[i], restart_tolerance * scale) << i;
            }
            for (std::size_t i = 0; i < continuous.bias.size(); ++i) {
                EXPECT_TRUE(std::isfinite(continuous.bias[i]));
                EXPECT_TRUE(std::isfinite(segmented.bias[i]));
                EXPECT_TRUE(std::isfinite(restarted.bias[i]));
                EXPECT_NEAR(segmented.bias[i], continuous.bias[i], 1.0e-14) << i;
                EXPECT_NEAR(restarted.bias[i], segmented.bias[i], 1.0e-14) << i;
            }
            EXPECT_GT(std::abs(continuous.bias[0]), 0.0);
            MPI_Barrier(MPI_COMM_WORLD);
            if (me == 0) {
                // Every two steps, with two walkers only for shared density.
                const int walkers = density ? 2 : 1;
                for (int replica = 0; replica < (mean ? 2 : 1); ++replica) {
                    const std::string suffix = mean ? "." + std::to_string(replica) : "";
                    EXPECT_EQ(counter(prefix + "_c.state" + suffix), 1 + seam * walkers);
                    EXPECT_EQ(counter(prefix + "_s.state" + suffix), 1 + seam * walkers);
                    EXPECT_EQ(counter(prefix + "_i.state" + suffix), 1 + (seam / 2) * walkers);
                    EXPECT_EQ(counter(prefix + "_r.state" + suffix), 1 + seam * walkers);
                }
                if (!::testing::Test::HasFailure()) {
                    for (const char *leg : {"c", "s", "i", "r"}) {
                        for (const char *extension : {".dat", ".log", ".state", ".kernels"}) {
                            const std::string file = prefix + "_" + leg + extension;
                            for (const std::string &backup : {"", "bck.last."}) {
                                std::remove((backup + file).c_str());
                                for (int replica = 0; replica < 4; ++replica)
                                    std::remove(
                                        (backup + file + "." + std::to_string(replica)).c_str());
                            }
                        }
                    }
                    std::remove((binary + ".0").c_str());
                    std::remove((binary + ".1").c_str());
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}

TEST(MPI, plumed_nmpimd_bead_modes_restart_continuity)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    for (const char *mode : {"bead_mean", "bead_density"}) {
        const bool mean                  = std::string(mode) == "bead_mean";
        const std::string prefix         = "test_plumed_nmpimd_restart_" + std::string(mode);
        const std::string initial_file   = prefix + "_initial.dat";
        const std::string restart_file   = prefix + "_continue.dat";
        const std::string restart_prefix = prefix + "_binary";
        std::string graph                = "d: DISTANCE ATOMS=1,2 NOPBC\n";
        if (mean) graph += "mean: ENSEMBLE ARG=d\n";
        graph += mean ? "bias: RESTRAINT ARG=mean.d AT=7.5 KAPPA=0.01\n"
                      : "bias: RESTRAINT ARG=d AT=7.5 KAPPA=0.01\n";
        if (me == 0) {
            std::ofstream initial(initial_file);
            initial << graph;
            std::ofstream restarted(restart_file);
            restarted << "RESTART\n" << graph;
            std::remove((restart_prefix + ".0").c_str());
            std::remove((restart_prefix + ".1").c_str());
        }
        MPI_Barrier(MPI_COMM_WORLD);

        const std::array<std::string, 4> logs = {
            prefix + "_continuous.log", prefix + "_segmented.log", prefix + "_restart_first.log",
            prefix + "_restart_second.log"};
        const auto continuous = run_nmpimd_bead_mode_segments(
            mode, restart_prefix.c_str(), 6, 0, false, initial_file.c_str(), restart_file.c_str(),
            logs[0].c_str(), "");
        const auto segmented = run_nmpimd_bead_mode_segments(
            mode, restart_prefix.c_str(), 3, 3, false, initial_file.c_str(), restart_file.c_str(),
            logs[1].c_str(), "");
        const auto restarted = run_nmpimd_bead_mode_segments(
            mode, restart_prefix.c_str(), 3, 3, true, initial_file.c_str(), restart_file.c_str(),
            logs[2].c_str(), logs[3].c_str());

        for (std::size_t i = 0; i < continuous.atoms.size(); ++i) {
            EXPECT_TRUE(std::isfinite(continuous.atoms[i]));
            EXPECT_NEAR(segmented.atoms[i], continuous.atoms[i], 1.0e-14) << mode << " " << i;
            EXPECT_NEAR(restarted.atoms[i], continuous.atoms[i], 1.0e-14) << mode << " " << i;
        }
        for (std::size_t i = 0; i < continuous.pimd.size(); ++i) {
            EXPECT_TRUE(std::isfinite(continuous.pimd[i]));
            const double tolerance = 1.0e-14 * std::max({1.0, std::abs(continuous.pimd[i]),
                                                         std::abs(restarted.pimd[i])});
            EXPECT_NEAR(segmented.pimd[i], continuous.pimd[i], tolerance) << mode << " " << i;
            EXPECT_NEAR(restarted.pimd[i], continuous.pimd[i], tolerance) << mode << " " << i;
        }
        for (std::size_t i = 0; i < continuous.bias.size(); ++i) {
            EXPECT_TRUE(std::isfinite(continuous.bias[i]));
            EXPECT_NEAR(segmented.bias[i], continuous.bias[i], 1.0e-14) << mode << " " << i;
            EXPECT_NEAR(restarted.bias[i], continuous.bias[i], 1.0e-14) << mode << " " << i;
        }
        EXPECT_GT(std::abs(continuous.bias[0]), 0.0);
        EXPECT_NEAR(continuous.bias[0], continuous.bias[1], 1.0e-14);
        EXPECT_NEAR(continuous.bias[2], 0.0, 1.0e-14);
        EXPECT_NEAR(continuous.bias[3], 0.0, 1.0e-14);
        EXPECT_EQ(continuous.volume, 8000.0);
        EXPECT_EQ(segmented.volume, continuous.volume);
        EXPECT_EQ(restarted.volume, continuous.volume);

        MPI_Barrier(MPI_COMM_WORLD);
        if (me == 0) {
            std::remove(initial_file.c_str());
            std::remove(restart_file.c_str());
            std::remove((restart_prefix + ".0").c_str());
            std::remove((restart_prefix + ".1").c_str());
            for (const auto &log : logs)
                for (int bead = 0; bead < 2; ++bead)
                    std::remove((log + "." + std::to_string(bead)).c_str());
        }
    }
}

TEST(MPI, plumed_pimd_single_bead_limit)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    const char *plumed_file = "test_plumed_pimd_single_bead.dat";
    const std::string regular_log =
        "test_plumed_pimd_single_bead_regular." + std::to_string(me) + ".log";
    const std::string centroid_log =
        "test_plumed_pimd_single_bead_centroid." + std::to_string(me) + ".log";
    const std::string normal_mode_log =
        "test_plumed_nmpimd_single_bead_centroid." + std::to_string(me) + ".log";
    if (me == 0) {
        std::ofstream restraint(plumed_file);
        restraint << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                  << "bias: RESTRAINT ARG=d AT=6.0 KAPPA=4.0\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const char *args[] = {"LAMMPS_test", "-screen", "none", "-log", "none", "-nocite", nullptr};
    char **argv        = (char **)args;
    int argc           = (sizeof(args) / sizeof(char *)) - 1;
    void *lmp          = lammps_open(argc, argv, MPI_COMM_SELF, nullptr);
    ASSERT_NE(lmp, nullptr);
    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 1);

    lammps_command(lmp, "units lj");
    lammps_command(lmp, "atom_style atomic");
    lammps_command(lmp, "atom_modify map array");
    lammps_command(lmp, "boundary p p p");
    lammps_command(lmp, "region box block 0.0 20.0 0.0 20.0 0.0 20.0");
    lammps_command(lmp, "create_box 1 box");
    lammps_command(lmp, "create_atoms 1 single 10.0 10.0 2.0");
    lammps_command(lmp, "create_atoms 1 single 19.0 10.0 2.0");
    lammps_command(lmp, "mass 1 1.0");
    lammps_command(lmp, "pair_style lj/cut 2.5");
    lammps_command(lmp, "pair_coeff * * 0.0 1.0");
    lammps_command(lmp, "timestep 0.001");
    lammps_command(lmp, "velocity all set 0.0 0.0 0.0");

    auto run_case = [&](const std::string &fix_command) {
        std::array<double, 7> result{};
        lammps_command(lmp, fix_command.c_str());
        lammps_command(lmp, "run 0 post no");
        // Stock fix plumed exposes its first calculated bias on the next setup.
        lammps_command(lmp, "run 0 post no");
        auto **forces = (double **)lammps_extract_atom(lmp, "f");
        EXPECT_NE(forces, nullptr);
        if (forces) {
            for (int atom = 0; atom < 2; ++atom) {
                const tagint atom_id = atom + 1;
                const int index      = lammps_map_atom(lmp, &atom_id);
                EXPECT_GE(index, 0);
                if (index < 0) continue;
                for (int dimension = 0; dimension < 3; ++dimension)
                    result[3 * atom + dimension] = forces[index][dimension];
            }
        }
        auto *bias =
            (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
        EXPECT_NE(bias, nullptr);
        if (bias) {
            result[6] = *bias;
            lammps_free(bias);
        }
        lammps_command(lmp, "unfix bias");
        return result;
    };

    const auto regular = run_case("fix bias all plumed plumedfile " + std::string(plumed_file) +
                                  " outfile " + regular_log);
    lammps_command(lmp, "fix fpimd all pimd/langevin method pimd ensemble nvt "
                        "integrator obabo thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");
    const auto centroid =
        run_case("fix bias all plumed plumedfile " + std::string(plumed_file) + " outfile " +
                 centroid_log + " path_integral centroid pimd_fix fpimd");
    lammps_command(lmp, "unfix fpimd");
    lammps_command(lmp, "fix fpimd all pimd/langevin method nmpimd ensemble nvt "
                        "integrator baoab thermostat PILE_L 5678 tau 1.0 temp 1.0 fixcom no");
    const auto normal_mode =
        run_case("fix bias all plumed plumedfile " + std::string(plumed_file) + " outfile " +
                 normal_mode_log + " path_integral centroid pimd_fix fpimd");

    EXPECT_NEAR(regular[6], 18.0, 1.0e-12);
    EXPECT_NEAR(centroid[6], 18.0, 1.0e-12);
    EXPECT_NEAR(normal_mode[6], 18.0, 1.0e-12);
    for (std::size_t i = 0; i < regular.size(); ++i) {
        EXPECT_NEAR(centroid[i], regular[i], 1.0e-12);
        EXPECT_NEAR(normal_mode[i], regular[i], 1.0e-12);
    }

    lammps_close(lmp);
    std::remove(regular_log.c_str());
    std::remove(centroid_log.c_str());
    std::remove(normal_mode_log.c_str());
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) std::remove(plumed_file);
};

TEST(MPI, plumed_pimd_cyclic_bead_permutation)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    const char *centroid_file        = "test_plumed_pimd_cyclic_centroid.dat";
    const char *bead_mean_file       = "test_plumed_pimd_cyclic_bead_mean.dat";
    const char *path_spread_file     = "test_plumed_pimd_cyclic_path_spread.dat";
    const char *centroid_spread_file = "test_plumed_pimd_cyclic_centroid_spread.dat";
    const char *bead_density_file    = "test_plumed_pimd_cyclic_bead_density.dat";
    const char *centroid_log         = "test_plumed_pimd_cyclic_centroid.log";
    const char *bead_mean_log        = "test_plumed_pimd_cyclic_bead_mean.log";
    const char *path_spread_log      = "test_plumed_pimd_cyclic_path_spread.log";
    const char *centroid_spread_log  = "test_plumed_pimd_cyclic_centroid_spread.log";
    const char *bead_density_log     = "test_plumed_pimd_cyclic_bead_density.log";
    if (me == 0) {
        std::ofstream centroid(centroid_file);
        centroid << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                 << "bias: RESTRAINT ARG=d AT=0.0 KAPPA=4.0\n";
        std::ofstream bead_mean(bead_mean_file);
        bead_mean << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                  << "mean: ENSEMBLE ARG=d\n"
                  << "bias: RESTRAINT ARG=mean.d AT=6.5 KAPPA=4.0\n";
        std::ofstream path_spread(path_spread_file);
        path_spread << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                    << "d2: CUSTOM ARG=d FUNC=x*x PERIODIC=NO\n"
                    << "mean: ENSEMBLE ARG=d\n"
                    << "mean2: ENSEMBLE ARG=d2\n"
                    << "spread2: CUSTOM ARG=mean.d,mean2.d2 FUNC=y-x*x PERIODIC=NO\n"
                    << "bias: RESTRAINT ARG=spread2 AT=0.0 KAPPA=4.0\n";
        std::ofstream centroid_spread(centroid_spread_file);
        centroid_spread << "dvec: DISTANCE ATOMS=1,2 COMPONENTS NOPBC\n"
                        << "s: CUSTOM ARG=dvec.x,dvec.y,dvec.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO\n"
                        << "s2: CUSTOM ARG=s FUNC=x*x PERIODIC=NO\n"
                        << "mean: ENSEMBLE ARG=s\n"
                        << "mean2: ENSEMBLE ARG=s2\n"
                        << "spread2: CUSTOM ARG=mean.s,mean2.s2 FUNC=y-x*x PERIODIC=NO\n"
                        << "centroid: ENSEMBLE ARG=dvec.x,dvec.y,dvec.z\n"
                        << "sc: CUSTOM ARG=centroid.dvec.x,centroid.dvec.y,centroid.dvec.z "
                           "FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO\n"
                        << "coupled: CUSTOM ARG=sc,spread2 FUNC=x*y PERIODIC=NO\n"
                        << "bias: BIASVALUE ARG=coupled\n";
        std::ofstream bead_density(bead_density_file);
        bead_density << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                     << "bias: RESTRAINT ARG=d AT=6.5 KAPPA=4.0\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const std::array<std::array<double, 2>, 4> positions = {
        {{{19.0, 10.0}}, {{10.0, 18.0}}, {{3.0, 10.0}}, {{10.0, 4.0}}}};
    auto run_case = [&](const char *plumed_file, const char *plumed_log,
                        const char *path_integral_mode, int shift) {
        const char *args[] = {"LAMMPS_test", "-screen", "none", "-log",    "none", "-partition",
                              "4x1",         "-in",     "none", "-nocite", nullptr};
        char **argv        = (char **)args;
        int argc           = (sizeof(args) / sizeof(char *)) - 1;
        void *lmp          = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
        EXPECT_NE(lmp, nullptr);
        if (!lmp) return 0.0;

        EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 1);
        const auto &position = positions[(me + shift) % nprocs];
        lammps_command(lmp, "units lj");
        lammps_command(lmp, "atom_style atomic");
        lammps_command(lmp, "atom_modify map array");
        lammps_command(lmp, "boundary p p p");
        lammps_command(lmp, "region box block 0.0 20.0 0.0 20.0 0.0 20.0");
        lammps_command(lmp, "create_box 1 box");
        lammps_command(lmp, "create_atoms 1 single 10.0 10.0 2.0");
        const std::string create_second = "create_atoms 1 single " + std::to_string(position[0]) +
                                          " " + std::to_string(position[1]) + " 2.0";
        lammps_command(lmp, create_second.c_str());
        lammps_command(lmp, "mass 1 1.0");
        lammps_command(lmp, "pair_style lj/cut 2.5");
        lammps_command(lmp, "pair_coeff * * 0.0 1.0");
        lammps_command(lmp, "timestep 0.001");
        lammps_command(lmp, "velocity all set 0.0 0.0 0.0");
        lammps_command(lmp, "fix fpimd all pimd/langevin method pimd ensemble nvt "
                            "integrator obabo thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");
        const std::string fix_command = "fix bias all plumed plumedfile " +
                                        std::string(plumed_file) + " outfile " + plumed_log + "." +
                                        std::to_string(shift) + " path_integral " +
                                        path_integral_mode + " pimd_fix fpimd";
        lammps_command(lmp, fix_command.c_str());
        lammps_command(lmp, "run 0 post no");

        double result = 0.0;
        auto *bias =
            (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
        EXPECT_NE(bias, nullptr);
        if (bias) {
            result = *bias;
            lammps_free(bias);
        }
        lammps_close(lmp);
        return result;
    };

    const double centroid_original    = run_case(centroid_file, centroid_log, "centroid", 0);
    const double centroid_shifted     = run_case(centroid_file, centroid_log, "centroid", 1);
    const double bead_mean_original   = run_case(bead_mean_file, bead_mean_log, "bead_mean", 0);
    const double bead_mean_shifted    = run_case(bead_mean_file, bead_mean_log, "bead_mean", 1);
    const double path_spread_original = run_case(path_spread_file, path_spread_log, "bead_mean", 0);
    const double path_spread_shifted  = run_case(path_spread_file, path_spread_log, "bead_mean", 1);
    const double centroid_spread_original =
        run_case(centroid_spread_file, centroid_spread_log, "bead_mean", 0);
    const double centroid_spread_shifted =
        run_case(centroid_spread_file, centroid_spread_log, "bead_mean", 1);
    const double bead_density_original =
        run_case(bead_density_file, bead_density_log, "bead_density", 0);
    const double bead_density_shifted =
        run_case(bead_density_file, bead_density_log, "bead_density", 1);
    EXPECT_NEAR(centroid_original, centroid_shifted, 1.0e-12);
    EXPECT_NEAR(bead_mean_original, bead_mean_shifted, 1.0e-12);
    EXPECT_NEAR(path_spread_original, path_spread_shifted, 1.0e-12);
    EXPECT_NEAR(centroid_spread_original, centroid_spread_shifted, 1.0e-12);
    EXPECT_NEAR(bead_density_original, bead_density_shifted, 1.0e-12);
    EXPECT_NEAR(centroid_original, me == 0 ? 1.0 : 0.0, 1.0e-12);
    EXPECT_NEAR(bead_mean_original, me == 0 ? 2.0 : 0.0, 1.0e-12);
    EXPECT_NEAR(path_spread_original, me == 0 ? 3.125 : 0.0, 1.0e-12);
    EXPECT_NEAR(centroid_spread_original, me == 0 ? 1.25 / std::sqrt(2.0) : 0.0, 1.0e-12);
    EXPECT_NEAR(bead_density_original, me == 0 ? 4.5 : 0.0, 1.0e-12);

    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        std::remove(centroid_file);
        std::remove(bead_mean_file);
        std::remove(path_spread_file);
        std::remove(centroid_spread_file);
        std::remove(bead_density_file);
        std::remove((std::string(centroid_log) + ".0").c_str());
        std::remove((std::string(centroid_log) + ".1").c_str());
    }
    for (int shift = 0; shift < 2; ++shift)
        std::remove(
            (std::string(bead_mean_log) + "." + std::to_string(shift) + "." + std::to_string(me))
                .c_str());
    for (int shift = 0; shift < 2; ++shift)
        std::remove(
            (std::string(path_spread_log) + "." + std::to_string(shift) + "." + std::to_string(me))
                .c_str());
    for (int shift = 0; shift < 2; ++shift)
        std::remove((std::string(centroid_spread_log) + "." + std::to_string(shift) + "." +
                     std::to_string(me))
                        .c_str());
    for (int shift = 0; shift < 2; ++shift)
        std::remove(
            (std::string(bead_density_log) + "." + std::to_string(shift) + "." + std::to_string(me))
                .c_str());
};

TEST(MPI, plumed_pimd_bias_modes)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);
    // Tables contain -grad(B); the beta/P Hamiltonian requires -P*grad(B).
    constexpr int beads = 4;

    const char *bead_zero_file             = "test_plumed_pimd_bead_zero.dat";
    const char *bead_restraint_file        = "test_plumed_pimd_bead_restraint.dat";
    const char *path_spread_zero_file      = "test_plumed_pimd_path_spread_zero.dat";
    const char *path_spread_restraint_file = "test_plumed_pimd_path_spread_restraint.dat";
    const char *centroid_spread_zero_file  = "test_plumed_pimd_centroid_spread_zero.dat";
    const char *centroid_spread_bias_file  = "test_plumed_pimd_centroid_spread_bias.dat";
    const char *centroid_zero_file         = "test_plumed_pimd_centroid_zero.dat";
    const char *centroid_restraint_file    = "test_plumed_pimd_centroid_restraint.dat";
    const char *density_zero_file          = "test_plumed_pimd_density_zero.dat";
    const char *density_restraint_file     = "test_plumed_pimd_density_restraint.dat";
    const char *bead_zero_log              = "test_plumed_pimd_bead_zero.log";
    const char *bead_restraint_log         = "test_plumed_pimd_bead_restraint.log";
    const char *path_spread_zero_log       = "test_plumed_pimd_path_spread_zero.log";
    const char *path_spread_restraint_log  = "test_plumed_pimd_path_spread_restraint.log";
    const char *centroid_spread_zero_log   = "test_plumed_pimd_centroid_spread_zero.log";
    const char *centroid_spread_bias_log   = "test_plumed_pimd_centroid_spread_bias.log";
    const char *centroid_zero_log          = "test_plumed_pimd_centroid_zero.log";
    const char *centroid_restraint_log     = "test_plumed_pimd_centroid_restraint.log";
    const char *density_zero_log           = "test_plumed_pimd_density_zero.log";
    const char *density_restraint_log      = "test_plumed_pimd_density_restraint.log";
    if (me == 0) {
        std::ofstream bead_zero(bead_zero_file);
        bead_zero << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                  << "mean: ENSEMBLE ARG=d\n";
        std::ofstream bead_restraint(bead_restraint_file);
        bead_restraint << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                       << "mean: ENSEMBLE ARG=d\n"
                       << "bias: RESTRAINT ARG=mean.d AT=6.5 KAPPA=4.0\n";
        std::ofstream path_spread_zero(path_spread_zero_file);
        path_spread_zero << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                         << "d2: CUSTOM ARG=d FUNC=x*x PERIODIC=NO\n"
                         << "mean: ENSEMBLE ARG=d\n"
                         << "mean2: ENSEMBLE ARG=d2\n"
                         << "spread2: CUSTOM ARG=mean.d,mean2.d2 FUNC=y-x*x PERIODIC=NO\n";
        std::ofstream path_spread_restraint(path_spread_restraint_file);
        path_spread_restraint << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                              << "d2: CUSTOM ARG=d FUNC=x*x PERIODIC=NO\n"
                              << "mean: ENSEMBLE ARG=d\n"
                              << "mean2: ENSEMBLE ARG=d2\n"
                              << "spread2: CUSTOM ARG=mean.d,mean2.d2 FUNC=y-x*x PERIODIC=NO\n"
                              << "bias: RESTRAINT ARG=spread2 AT=0.0 KAPPA=4.0\n";
        const std::string centroid_spread_graph =
            "dvec: DISTANCE ATOMS=1,2 COMPONENTS NOPBC\n"
            "s: CUSTOM ARG=dvec.x,dvec.y,dvec.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO\n"
            "s2: CUSTOM ARG=s FUNC=x*x PERIODIC=NO\n"
            "mean: ENSEMBLE ARG=s\n"
            "mean2: ENSEMBLE ARG=s2\n"
            "spread2: CUSTOM ARG=mean.s,mean2.s2 FUNC=y-x*x PERIODIC=NO\n"
            "centroid: ENSEMBLE ARG=dvec.x,dvec.y,dvec.z\n"
            "sc: CUSTOM ARG=centroid.dvec.x,centroid.dvec.y,centroid.dvec.z "
            "FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO\n"
            "coupled: CUSTOM ARG=sc,spread2 FUNC=x*y PERIODIC=NO\n";
        std::ofstream centroid_spread_zero(centroid_spread_zero_file);
        centroid_spread_zero << centroid_spread_graph;
        std::ofstream centroid_spread_bias(centroid_spread_bias_file);
        centroid_spread_bias << centroid_spread_graph << "bias: BIASVALUE ARG=coupled\n";
        std::ofstream centroid_zero(centroid_zero_file);
        centroid_zero << "d: DISTANCE ATOMS=1,2 NOPBC\n";
        std::ofstream centroid_restraint(centroid_restraint_file);
        centroid_restraint << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                           << "bias: RESTRAINT ARG=d AT=0.0 KAPPA=4.0\n";
        std::ofstream density_zero(density_zero_file);
        density_zero << "d: DISTANCE ATOMS=1,2 NOPBC\n";
        std::ofstream density_restraint(density_restraint_file);
        density_restraint << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                          << "bias: RESTRAINT ARG=d AT=6.5 KAPPA=4.0\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const char *args[] = {"LAMMPS_test", "-screen", "none", "-log",    "none", "-partition",
                          "4x1",         "-in",     "none", "-nocite", nullptr};
    char **argv        = (char **)args;
    int argc           = (sizeof(args) / sizeof(char *)) - 1;
    void *lmp          = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
    ASSERT_NE(lmp, nullptr);

    lammps_command(lmp, "variable x2 world 19.0 10.0 3.0 10.0");
    lammps_command(lmp, "variable y2 world 10.0 18.0 10.0 4.0");
    lammps_command(lmp, "units lj");
    lammps_command(lmp, "atom_style atomic");
    lammps_command(lmp, "atom_modify map array");
    lammps_command(lmp, "boundary p p p");
    lammps_command(lmp, "region box block 0.0 20.0 0.0 20.0 0.0 20.0");
    lammps_command(lmp, "create_box 1 box");
    lammps_command(lmp, "create_atoms 1 single 10.0 10.0 2.0");
    lammps_command(lmp, "create_atoms 1 single ${x2} ${y2} 2.0");
    lammps_command(lmp, "mass 1 1.0");
    lammps_command(lmp, "pair_style lj/cut 2.5");
    lammps_command(lmp, "pair_coeff * * 0.0 1.0");
    lammps_command(lmp, "timestep 0.001");
    lammps_command(lmp, "velocity all set 0.0 0.0 0.0");
    // Keep nonzero springs modest so force differences do not lose precision.
    lammps_command(lmp,
                   "fix fpimd all pimd/langevin method pimd ensemble nvt "
                   "integrator obabo thermostat PILE_L 1234 tau 1.0 temp 1.0 sp 10.0 fixcom no");
    const std::array<std::array<double, 6>, 4> bead_force_delta = {
        {{1.0, 0.0, 0.0, -1.0, 0.0, 0.0},
         {0.0, 1.0, 0.0, 0.0, -1.0, 0.0},
         {-1.0, 0.0, 0.0, 1.0, 0.0, 0.0},
         {0.0, -1.0, 0.0, 0.0, 1.0, 0.0}}};
    const std::array<std::array<double, 6>, 4> path_spread_force_delta = {
        {{3.75, 0.0, 0.0, -3.75, 0.0, 0.0},
         {0.0, 1.25, 0.0, 0.0, -1.25, 0.0},
         {1.25, 0.0, 0.0, -1.25, 0.0, 0.0},
         {0.0, 3.75, 0.0, 0.0, -3.75, 0.0}}};
    const std::array<std::array<double, 6>, 4> centroid_spread_force_delta = {
        {{0.7513009550107068, 0.2209708691207961, 0.0, -0.7513009550107068, -0.2209708691207961,
          0.0},
         {0.2209708691207961, 0.3977475644174330, 0.0, -0.2209708691207961, -0.3977475644174330,
          0.0},
         {0.3977475644174330, 0.2209708691207961, 0.0, -0.3977475644174330, -0.2209708691207961,
          0.0},
         {0.2209708691207961, 0.7513009550107068, 0.0, -0.2209708691207961, -0.7513009550107068,
          0.0}}};
    const std::array<double, 6> centroid_force_delta = {0.5, 0.5, 0.0, -0.5, -0.5, 0.0};
    const std::array<std::array<double, 6>, 4> density_force_delta = {
        {{2.5, 0.0, 0.0, -2.5, 0.0, 0.0},
         {0.0, 1.5, 0.0, 0.0, -1.5, 0.0},
         {-0.5, 0.0, 0.0, 0.5, 0.0, 0.0},
         {0.0, 0.5, 0.0, 0.0, -0.5, 0.0}}};
    auto extract_forces = [&]() {
        std::array<double, 6> result{};
        auto **forces = (double **)lammps_extract_atom(lmp, "f");
        EXPECT_NE(forces, nullptr);
        if (!forces) return result;

        for (int atom = 0; atom < 2; ++atom) {
            const tagint atom_id = atom + 1;
            const int index      = lammps_map_atom(lmp, &atom_id);
            EXPECT_GE(index, 0);
            if (index < 0) continue;
            for (int dimension = 0; dimension < 3; ++dimension)
                result[3 * atom + dimension] = forces[index][dimension];
        }
        return result;
    };
    auto check_mode = [&](const char *zero_command, const char *bias_command, double expected_bias,
                          const std::array<double, 6> &expected_force_delta,
                          double force_tolerance = 1.0e-12) {
        lammps_command(lmp, zero_command);
        lammps_command(lmp, "run 0 post no");
        const auto zero_forces = extract_forces();
        lammps_command(lmp, "unfix zero");

        lammps_command(lmp, bias_command);
        lammps_command(lmp, "run 0 post no");
        const auto biased_forces = extract_forces();
        auto *bias =
            (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
        ASSERT_NE(bias, nullptr);
        EXPECT_NEAR(*bias, me == 0 ? expected_bias : 0.0, 1.0e-12);
        lammps_free(bias);

        for (std::size_t i = 0; i < expected_force_delta.size(); ++i)
            EXPECT_NEAR(biased_forces[i] - zero_forces[i], beads * expected_force_delta[i],
                        force_tolerance);
        lammps_command(lmp, "unfix bias");
    };

    check_mode("fix zero all plumed plumedfile test_plumed_pimd_bead_zero.dat "
               "outfile test_plumed_pimd_bead_zero.log path_integral bead_mean pimd_fix fpimd",
               "fix bias all plumed plumedfile test_plumed_pimd_bead_restraint.dat "
               "outfile test_plumed_pimd_bead_restraint.log path_integral bead_mean pimd_fix "
               "fpimd",
               2.0, bead_force_delta[me]);
    check_mode("fix zero all plumed plumedfile test_plumed_pimd_path_spread_zero.dat "
               "outfile test_plumed_pimd_path_spread_zero.log path_integral bead_mean "
               "pimd_fix fpimd",
               "fix bias all plumed plumedfile test_plumed_pimd_path_spread_restraint.dat "
               "outfile test_plumed_pimd_path_spread_restraint.log path_integral bead_mean "
               "pimd_fix fpimd",
               3.125, path_spread_force_delta[me]);
    check_mode("fix zero all plumed plumedfile test_plumed_pimd_centroid_spread_zero.dat "
               "outfile test_plumed_pimd_centroid_spread_zero.log path_integral bead_mean "
               "pimd_fix fpimd",
               "fix bias all plumed plumedfile test_plumed_pimd_centroid_spread_bias.dat "
               "outfile test_plumed_pimd_centroid_spread_bias.log path_integral bead_mean "
               "pimd_fix fpimd",
               1.25 / std::sqrt(2.0), centroid_spread_force_delta[me], 2.0e-12);
    check_mode("fix zero all plumed plumedfile test_plumed_pimd_centroid_zero.dat "
               "outfile test_plumed_pimd_centroid_zero.log path_integral centroid pimd_fix fpimd",
               "fix bias all plumed plumedfile test_plumed_pimd_centroid_restraint.dat "
               "outfile test_plumed_pimd_centroid_restraint.log path_integral centroid pimd_fix "
               "fpimd",
               1.0, centroid_force_delta);
    check_mode("fix zero all plumed plumedfile test_plumed_pimd_density_zero.dat "
               "outfile test_plumed_pimd_density_zero.log path_integral bead_density pimd_fix "
               "fpimd",
               "fix bias all plumed plumedfile test_plumed_pimd_density_restraint.dat "
               "outfile test_plumed_pimd_density_restraint.log path_integral bead_density "
               "pimd_fix fpimd",
               4.5, density_force_delta[me]);

    lammps_close(lmp);
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        std::remove(bead_zero_file);
        std::remove(bead_restraint_file);
        std::remove(path_spread_zero_file);
        std::remove(path_spread_restraint_file);
        std::remove(centroid_spread_zero_file);
        std::remove(centroid_spread_bias_file);
        std::remove(centroid_zero_file);
        std::remove(centroid_restraint_file);
        std::remove(density_zero_file);
        std::remove(density_restraint_file);
        for (int i = 0; i < nprocs; ++i) {
            std::remove((std::string(bead_zero_log) + "." + std::to_string(i)).c_str());
            std::remove((std::string(bead_restraint_log) + "." + std::to_string(i)).c_str());
            std::remove((std::string(path_spread_zero_log) + "." + std::to_string(i)).c_str());
            std::remove((std::string(path_spread_restraint_log) + "." + std::to_string(i)).c_str());
            std::remove((std::string(centroid_spread_zero_log) + "." + std::to_string(i)).c_str());
            std::remove((std::string(centroid_spread_bias_log) + "." + std::to_string(i)).c_str());
        }
        std::remove(centroid_zero_log);
        std::remove(centroid_restraint_log);
    }
    std::remove((std::string(density_zero_log) + "." + std::to_string(me)).c_str());
    std::remove((std::string(density_restraint_log) + "." + std::to_string(me)).c_str());
};

TEST(MPI, plumed_pimd_centroid_bead_density_dynamics_equivalence)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    const char *plumed_file  = "test_plumed_pimd_matched_dynamics.dat";
    const char *centroid_log = "test_plumed_pimd_matched_centroid.log";
    const char *density_log  = "test_plumed_pimd_matched_density.log";
    if (me == 0) {
        std::ofstream input(plumed_file);
        input << "d: DISTANCE ATOMS=1,2 COMPONENTS NOPBC\n"
              << "bias: BIASVALUE ARG=d.x\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    auto run_case = [&](const char *mode, const char *plumed_log) {
        std::array<double, 13> result{};
        const char *args[] = {"LAMMPS_test", "-screen", "none", "-log",    "none", "-partition",
                              "4x1",         "-in",     "none", "-nocite", nullptr};
        char **argv        = (char **)args;
        int argc           = (sizeof(args) / sizeof(char *)) - 1;
        void *lmp          = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
        EXPECT_NE(lmp, nullptr);
        if (!lmp) return result;

        lammps_command(lmp, "variable x2 world 19.0 18.0 17.0 16.0");
        lammps_command(lmp, "units lj");
        lammps_command(lmp, "atom_style atomic");
        lammps_command(lmp, "atom_modify map array");
        lammps_command(lmp, "boundary p p p");
        lammps_command(lmp, "region box block 0.0 20.0 0.0 20.0 0.0 20.0");
        lammps_command(lmp, "create_box 1 box");
        lammps_command(lmp, "create_atoms 1 single 10.0 10.0 2.0");
        lammps_command(lmp, "create_atoms 1 single ${x2} 10.0 2.0");
        lammps_command(lmp, "mass 1 1.0");
        lammps_command(lmp, "pair_style lj/cut 2.5");
        lammps_command(lmp, "pair_coeff * * 0.0 1.0");
        lammps_command(lmp, "timestep 0.001");
        lammps_command(lmp, "velocity all set 0.0 0.0 0.0");
        lammps_command(lmp, "fix fpimd all pimd/langevin method pimd ensemble nvt "
                            "integrator obabo thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");
        const std::string fix_command = "fix bias all plumed plumedfile " +
                                        std::string(plumed_file) + " outfile " + plumed_log +
                                        " path_integral " + mode + " pimd_fix fpimd";
        lammps_command(lmp, fix_command.c_str());
        lammps_command(lmp, "run 5 post no");

        auto **positions  = (double **)lammps_extract_atom(lmp, "x");
        auto **velocities = (double **)lammps_extract_atom(lmp, "v");
        EXPECT_NE(positions, nullptr);
        EXPECT_NE(velocities, nullptr);
        if (positions && velocities) {
            for (int atom = 0; atom < 2; ++atom) {
                const tagint atom_id = atom + 1;
                const int index      = lammps_map_atom(lmp, &atom_id);
                EXPECT_GE(index, 0);
                if (index < 0) continue;
                for (int dimension = 0; dimension < 3; ++dimension) {
                    result[3 * atom + dimension]     = positions[index][dimension];
                    result[6 + 3 * atom + dimension] = velocities[index][dimension];
                }
            }
        }
        auto *bias =
            (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
        EXPECT_NE(bias, nullptr);
        if (bias) {
            result[12] = *bias;
            lammps_free(bias);
        }
        lammps_close(lmp);
        return result;
    };

    const auto centroid = run_case("centroid", centroid_log);
    const auto density  = run_case("bead_density", density_log);
    for (std::size_t i = 0; i < centroid.size(); ++i)
        EXPECT_NEAR(density[i], centroid[i], 1.0e-12);
    if (me == 0) {
        EXPECT_GT(centroid[12], 0.0);
        EXPECT_GT(density[12], 0.0);
    } else {
        EXPECT_NEAR(centroid[12], 0.0, 1.0e-12);
        EXPECT_NEAR(density[12], 0.0, 1.0e-12);
    }
    double max_speed = 0.0;
    for (std::size_t i = 6; i < 12; ++i)
        max_speed = std::max(max_speed, std::fabs(centroid[i]));
    EXPECT_GT(max_speed, 1.0e-8);

    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        std::remove(plumed_file);
        std::remove(centroid_log);
    }
    std::remove((std::string(density_log) + "." + std::to_string(me)).c_str());
};

TEST(MPI, plumed_pimd_bead_density_shared_history_restart)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    const char *initial_file      = "test_plumed_pimd_density_metad.dat";
    const char *restart_file      = "test_plumed_pimd_density_metad_restart.dat";
    const char *hills_file        = "test_plumed_pimd_density_metad_hills";
    const char *initial_log       = "test_plumed_pimd_density_metad.log";
    const char *restart_log       = "test_plumed_pimd_density_metad_restart.log";
    const char *opes_initial_file = "test_plumed_pimd_density_opes.dat";
    const char *opes_restart_file = "test_plumed_pimd_density_opes_restart.dat";
    const char *kernels_initial   = "test_plumed_pimd_density_opes_kernels";
    const char *kernels_restart   = "test_plumed_pimd_density_opes_kernels_restart";
    const char *state_initial     = "test_plumed_pimd_density_opes_state";
    const char *state_restart     = "test_plumed_pimd_density_opes_state_restart";
    const char *colvar_initial    = "test_plumed_pimd_density_opes_colvar";
    const char *colvar_restart    = "test_plumed_pimd_density_opes_colvar_restart";
    const char *opes_initial_log  = "test_plumed_pimd_density_opes.log";
    const char *opes_restart_log  = "test_plumed_pimd_density_opes_restart.log";

    if (me == 0) {
        std::remove(initial_file);
        std::remove(restart_file);
        std::remove(hills_file);
        for (int walker = 0; walker < nprocs; ++walker)
            std::remove((std::string(hills_file) + "." + std::to_string(walker)).c_str());
        for (const char *file : {opes_initial_file, opes_restart_file, kernels_initial,
                                 kernels_restart, state_initial, state_restart})
            std::remove(file);
        std::remove((std::string("bck.last.") + state_initial).c_str());
        std::remove((std::string("bck.last.") + state_restart).c_str());
        for (int walker = 0; walker < nprocs; ++walker) {
            for (const char *file :
                 {kernels_initial, kernels_restart, state_initial, state_restart})
                std::remove((std::string(file) + "." + std::to_string(walker)).c_str());
            std::remove((std::string(colvar_initial) + "." + std::to_string(walker)).c_str());
            std::remove((std::string(colvar_restart) + "." + std::to_string(walker)).c_str());
        }

        std::ofstream initial(initial_file);
        initial << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                << "bias: METAD ARG=d SIGMA=0.5 HEIGHT=0.25 PACE=2 FILE=" << hills_file
                << " WALKERS_MPI RESTART=NO\n";
        std::ofstream restart(restart_file);
        restart << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                << "bias: METAD ARG=d SIGMA=0.5 HEIGHT=0.25 PACE=2 FILE=" << hills_file
                << " WALKERS_MPI RESTART=YES\n";
        std::ofstream opes_initial(opes_initial_file);
        opes_initial << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                     << "bias: OPES_METAD ARG=d PACE=2 BARRIER=4 TEMP=1 SIGMA=0.5 FILE="
                     << kernels_initial << " STATE_WFILE=" << state_initial
                     << " STATE_WSTRIDE=2 WALKERS_MPI RESTART=NO\n"
                     << "PRINT ARG=d,bias.bias,bias.rct,bias.zed,bias.neff,bias.nker FILE="
                     << colvar_initial << " STRIDE=1\n";
        std::ofstream opes_restart(opes_restart_file);
        opes_restart << "RESTART\n"
                     << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                     << "bias: OPES_METAD ARG=d PACE=2 BARRIER=4 TEMP=1 SIGMA=0.5 FILE="
                     << kernels_restart << " STATE_RFILE=" << state_initial
                     << " STATE_WFILE=" << state_restart
                     << " STATE_WSTRIDE=2 WALKERS_MPI RESTART=YES\n"
                     << "PRINT ARG=d,bias.bias,bias.rct,bias.zed,bias.neff,bias.nker FILE="
                     << colvar_restart << " STRIDE=1 RESTART=NO\n";
    }
    std::remove((std::string(initial_log) + "." + std::to_string(me)).c_str());
    std::remove((std::string(restart_log) + "." + std::to_string(me)).c_str());
    std::remove((std::string(opes_initial_log) + "." + std::to_string(me)).c_str());
    std::remove((std::string(opes_restart_log) + "." + std::to_string(me)).c_str());
    MPI_Barrier(MPI_COMM_WORLD);

    auto run_case = [&](const char *plumed_file, const char *plumed_log, int step) {
        const char *args[] = {"LAMMPS_test", "-screen", "none", "-log",    "none", "-partition",
                              "4x1",         "-in",     "none", "-nocite", nullptr};
        char **argv        = (char **)args;
        int argc           = (sizeof(args) / sizeof(char *)) - 1;
        void *lmp          = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
        EXPECT_NE(lmp, nullptr);
        if (!lmp) return 0.0;
        lammps_set_show_error(lmp, 0);

        auto has_error = [&](const char *stage) {
            if (!lammps_has_error(lmp)) return false;
            char error_message[512];
            lammps_get_last_error_message(lmp, error_message, sizeof(error_message));
            ADD_FAILURE() << stage << ": " << error_message;
            return true;
        };

        lammps_command(lmp, "variable x2 world 19.0 18.0 17.0 16.0");
        lammps_command(lmp, "units lj");
        lammps_command(lmp, "atom_style atomic");
        lammps_command(lmp, "atom_modify map array");
        lammps_command(lmp, "boundary p p p");
        lammps_command(lmp, "region box block 0.0 20.0 0.0 20.0 0.0 20.0");
        lammps_command(lmp, "create_box 1 box");
        lammps_command(lmp, "create_atoms 1 single 10.0 10.0 2.0");
        lammps_command(lmp, "create_atoms 1 single ${x2} 10.0 2.0");
        lammps_command(lmp, "mass 1 1.0");
        lammps_command(lmp, "pair_style lj/cut 2.5");
        lammps_command(lmp, "pair_coeff * * 0.0 1.0");
        lammps_command(lmp, "timestep 0.001");
        lammps_command(lmp, "velocity all set 0.0 0.0 0.0");
        lammps_command(lmp, ("reset_timestep " + std::to_string(step)).c_str());
        lammps_command(lmp, "fix fpimd all pimd/langevin method pimd ensemble nvt "
                            "integrator obabo thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");
        const std::string fix_command = "fix bias all plumed plumedfile " +
                                        std::string(plumed_file) + " outfile " + plumed_log +
                                        " path_integral bead_density pimd_fix fpimd";
        lammps_command(lmp, fix_command.c_str());
        if (has_error("fix plumed")) {
            lammps_close(lmp);
            return 0.0;
        }
        lammps_command(lmp, "run 1 post no");
        if (has_error("run")) {
            lammps_close(lmp);
            return 0.0;
        }

        double result = 0.0;
        auto *bias =
            (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
        EXPECT_NE(bias, nullptr);
        if (bias) {
            result = *bias;
            lammps_free(bias);
        }
        lammps_close(lmp);
        return result;
    };

    auto read_records = [](const char *file) {
        std::vector<std::vector<double>> records;
        std::ifstream input(file);
        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream row(line);
            std::vector<double> record;
            double value;
            while (row >> value)
                record.push_back(value);
            if (!record.empty()) records.push_back(record);
        }
        return records;
    };
    auto read_state_counter = [](const char *file) {
        std::ifstream input(file);
        std::string line;
        while (std::getline(input, line)) {
            if (line.rfind("#! SET counter", 0) != 0) continue;
            std::istringstream row(line);
            std::string marker, set, key;
            int counter = -1;
            row >> marker >> set >> key >> counter;
            return counter;
        }
        return -1;
    };

    const double initial_bias = run_case(initial_file, initial_log, 1);
    EXPECT_NEAR(initial_bias, 0.0, 1.0e-12);
    MPI_Barrier(MPI_COMM_WORLD);

    const auto initial_hills = read_records(hills_file);
    ASSERT_EQ(initial_hills.size(), 4);
    for (const auto &record : initial_hills) {
        ASSERT_EQ(record.size(), 5);
        EXPECT_NEAR(record[3], 0.25, 1.0e-12);
    }
    for (int walker = 0; walker < nprocs; ++walker) {
        std::ifstream split_hills(std::string(hills_file) + "." + std::to_string(walker));
        EXPECT_FALSE(split_hills.good());
    }

    const double restart_bias = run_case(restart_file, restart_log, 3);
    if (me == 0)
        EXPECT_GT(restart_bias, 0.0);
    else
        EXPECT_NEAR(restart_bias, 0.0, 1.0e-12);
    MPI_Barrier(MPI_COMM_WORLD);

    const auto restarted_hills = read_records(hills_file);
    ASSERT_EQ(restarted_hills.size(), 8);
    for (std::size_t i = 0; i < restarted_hills.size(); ++i) {
        ASSERT_EQ(restarted_hills[i].size(), 5);
        EXPECT_NEAR(restarted_hills[i][3], 0.25, 1.0e-12);
        if (i < initial_hills.size()) EXPECT_EQ(restarted_hills[i], initial_hills[i]);
    }
    EXPECT_GT(restarted_hills[4][0], restarted_hills[3][0]);

    const double initial_opes_bias = run_case(opes_initial_file, opes_initial_log, 1);
    EXPECT_TRUE(std::isfinite(initial_opes_bias));
    if (me == 0)
        EXPECT_LT(initial_opes_bias, 0.0);
    else
        EXPECT_NEAR(initial_opes_bias, 0.0, 1.0e-12);
    MPI_Barrier(MPI_COMM_WORLD);

    const auto initial_kernels = read_records(kernels_initial);
    const auto initial_state   = read_records(state_initial);
    ASSERT_EQ(initial_kernels.size(), 4);
    ASSERT_EQ(initial_state.size(), 4);
    EXPECT_EQ(read_state_counter(state_initial), 5);
    for (const auto &record : initial_kernels)
        ASSERT_EQ(record.size(), 5);
    for (const auto &record : initial_state)
        ASSERT_EQ(record.size(), 4);
    for (int walker = 0; walker < nprocs; ++walker) {
        std::ifstream split_kernels(std::string(kernels_initial) + "." + std::to_string(walker));
        std::ifstream split_state(std::string(state_initial) + "." + std::to_string(walker));
        EXPECT_FALSE(split_kernels.good());
        EXPECT_FALSE(split_state.good());
        const std::string colvar = std::string(colvar_initial) + "." + std::to_string(walker);
        const auto records       = read_records(colvar.c_str());
        ASSERT_EQ(records.size(), 2);
        ASSERT_EQ(records.back().size(), 7);
        EXPECT_NEAR(records.back()[6], 4.0, 1.0e-12);
    }

    const double restart_opes_bias = run_case(opes_restart_file, opes_restart_log, 3);
    EXPECT_TRUE(std::isfinite(restart_opes_bias));
    if (me == 0)
        EXPECT_NE(restart_opes_bias, 0.0);
    else
        EXPECT_NEAR(restart_opes_bias, 0.0, 1.0e-12);
    MPI_Barrier(MPI_COMM_WORLD);

    EXPECT_EQ(read_records(state_initial), initial_state);
    const auto restarted_kernels = read_records(kernels_restart);
    const auto restarted_state   = read_records(state_restart);
    ASSERT_EQ(restarted_kernels.size(), 4);
    ASSERT_EQ(restarted_state.size(), 4);
    EXPECT_EQ(read_state_counter(state_restart), 9);
    for (const auto &record : restarted_kernels)
        ASSERT_EQ(record.size(), 5);
    for (const auto &record : restarted_state)
        ASSERT_EQ(record.size(), 4);
    for (int walker = 0; walker < nprocs; ++walker) {
        for (const char *file : {kernels_restart, state_restart}) {
            std::ifstream split_file(std::string(file) + "." + std::to_string(walker));
            EXPECT_FALSE(split_file.good());
        }
        const std::string initial_colvar =
            std::string(colvar_initial) + "." + std::to_string(walker);
        const std::string restart_colvar =
            std::string(colvar_restart) + "." + std::to_string(walker);
        const auto initial_records = read_records(initial_colvar.c_str());
        const auto restart_records = read_records(restart_colvar.c_str());
        ASSERT_EQ(initial_records.size(), 2);
        ASSERT_EQ(restart_records.size(), 2);
        ASSERT_EQ(initial_records.back().size(), 7);
        ASSERT_EQ(restart_records.front().size(), 7);
        EXPECT_NEAR(restart_records.front()[3], initial_records.back()[3], 1.0e-12);
        EXPECT_NEAR(restart_records.front()[4], initial_records.back()[4], 1.0e-12);
        EXPECT_NEAR(restart_records.front()[5], initial_records.back()[5], 1.0e-12);
        EXPECT_NEAR(restart_records.front()[6], 4.0, 1.0e-12);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        std::remove(initial_file);
        std::remove(restart_file);
        std::remove(hills_file);
        for (const char *file : {opes_initial_file, opes_restart_file, kernels_initial,
                                 kernels_restart, state_initial, state_restart})
            std::remove(file);
        std::remove((std::string("bck.last.") + state_initial).c_str());
        std::remove((std::string("bck.last.") + state_restart).c_str());
        for (int walker = 0; walker < nprocs; ++walker) {
            std::remove((std::string(colvar_initial) + "." + std::to_string(walker)).c_str());
            std::remove((std::string(colvar_restart) + "." + std::to_string(walker)).c_str());
        }
    }
    std::remove((std::string(initial_log) + "." + std::to_string(me)).c_str());
    std::remove((std::string(restart_log) + "." + std::to_string(me)).c_str());
    std::remove((std::string(opes_initial_log) + "." + std::to_string(me)).c_str());
    std::remove((std::string(opes_restart_log) + "." + std::to_string(me)).c_str());
    std::remove((std::string("bck.0.") + opes_restart_log + "." + std::to_string(me)).c_str());
};

TEST(MPI, plumed_pimd_centroid_spread_shared_opes_state_restart)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);

    const char *initial_file    = "test_plumed_pimd_centroid_spread_opes.dat";
    const char *restart_file    = "test_plumed_pimd_centroid_spread_opes_restart.dat";
    const char *kernels_initial = "test_plumed_pimd_centroid_spread_opes_kernels";
    const char *kernels_restart = "test_plumed_pimd_centroid_spread_opes_kernels_restart";
    const char *state_initial   = "test_plumed_pimd_centroid_spread_opes_state";
    const char *state_restart   = "test_plumed_pimd_centroid_spread_opes_state_restart";
    const char *colvar_initial  = "test_plumed_pimd_centroid_spread_opes_colvar";
    const char *colvar_restart  = "test_plumed_pimd_centroid_spread_opes_colvar_restart";
    const char *initial_log     = "test_plumed_pimd_centroid_spread_opes.log";
    const char *restart_log     = "test_plumed_pimd_centroid_spread_opes_restart.log";
    const std::string cv_graph =
        "dvec: DISTANCE ATOMS=1,2 COMPONENTS NOPBC\n"
        "s: CUSTOM ARG=dvec.x,dvec.y,dvec.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO\n"
        "s2: CUSTOM ARG=s FUNC=x*x PERIODIC=NO\n"
        "mean: ENSEMBLE ARG=s\n"
        "mean2: ENSEMBLE ARG=s2\n"
        "spread2: CUSTOM ARG=mean.s,mean2.s2 FUNC=y-x*x PERIODIC=NO\n"
        "centroid: ENSEMBLE ARG=dvec.x,dvec.y,dvec.z\n"
        "sc: CUSTOM ARG=centroid.dvec.x,centroid.dvec.y,centroid.dvec.z "
        "FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO\n";

    if (me == 0) {
        for (const char *file : {initial_file, restart_file, kernels_initial, kernels_restart,
                                 state_initial, state_restart})
            std::remove(file);
        std::remove((std::string("bck.last.") + state_initial).c_str());
        std::remove((std::string("bck.last.") + state_restart).c_str());
        for (int walker = 0; walker < nprocs; ++walker) {
            for (const char *file :
                 {kernels_initial, kernels_restart, state_initial, state_restart})
                std::remove((std::string(file) + "." + std::to_string(walker)).c_str());
            std::remove((std::string(colvar_initial) + "." + std::to_string(walker)).c_str());
            std::remove((std::string(colvar_restart) + "." + std::to_string(walker)).c_str());
        }

        std::ofstream initial(initial_file);
        initial << cv_graph
                << "bias: OPES_METAD ARG=sc,spread2 PACE=2 BARRIER=4 TEMP=1 "
                   "SIGMA=0.5,0.25 FILE="
                << kernels_initial << " STATE_WFILE=" << state_initial
                << " STATE_WSTRIDE=2 WALKERS_MPI RESTART=NO\n"
                << "PRINT ARG=sc,spread2,bias.bias,bias.rct,bias.zed,bias.neff,bias.nker FILE="
                << colvar_initial << " STRIDE=1\n";
        std::ofstream restart(restart_file);
        restart << "RESTART\n"
                << cv_graph
                << "bias: OPES_METAD ARG=sc,spread2 PACE=2 BARRIER=4 TEMP=1 "
                   "SIGMA=0.5,0.25 FILE="
                << kernels_restart << " STATE_RFILE=" << state_initial
                << " STATE_WFILE=" << state_restart << " STATE_WSTRIDE=2 WALKERS_MPI RESTART=YES\n"
                << "PRINT ARG=sc,spread2,bias.bias,bias.rct,bias.zed,bias.neff,bias.nker FILE="
                << colvar_restart << " STRIDE=1 RESTART=NO\n";
    }
    std::remove((std::string(initial_log) + "." + std::to_string(me)).c_str());
    std::remove((std::string(restart_log) + "." + std::to_string(me)).c_str());
    MPI_Barrier(MPI_COMM_WORLD);

    auto run_case = [&](const char *plumed_file, const char *plumed_log, int step) {
        const char *args[] = {"LAMMPS_test", "-screen", "none", "-log",    "none", "-partition",
                              "4x1",         "-in",     "none", "-nocite", nullptr};
        char **argv        = (char **)args;
        int argc           = (sizeof(args) / sizeof(char *)) - 1;
        void *lmp          = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
        EXPECT_NE(lmp, nullptr);
        if (!lmp) return 0.0;
        lammps_set_show_error(lmp, 0);

        auto has_error = [&](const char *stage) {
            if (!lammps_has_error(lmp)) return false;
            char error_message[512];
            lammps_get_last_error_message(lmp, error_message, sizeof(error_message));
            ADD_FAILURE() << stage << ": " << error_message;
            return true;
        };

        lammps_command(lmp, "variable x2 world 19.0 18.0 17.0 16.0");
        lammps_command(lmp, "units lj");
        lammps_command(lmp, "atom_style atomic");
        lammps_command(lmp, "atom_modify map array");
        lammps_command(lmp, "boundary p p p");
        lammps_command(lmp, "region box block 0.0 20.0 0.0 20.0 0.0 20.0");
        lammps_command(lmp, "create_box 1 box");
        lammps_command(lmp, "create_atoms 1 single 10.0 10.0 2.0");
        lammps_command(lmp, "create_atoms 1 single ${x2} 10.0 2.0");
        lammps_command(lmp, "mass 1 1.0");
        lammps_command(lmp, "pair_style lj/cut 2.5");
        lammps_command(lmp, "pair_coeff * * 0.0 1.0");
        lammps_command(lmp, "timestep 0.001");
        lammps_command(lmp, "velocity all set 0.0 0.0 0.0");
        lammps_command(lmp, ("reset_timestep " + std::to_string(step)).c_str());
        lammps_command(lmp, "fix fpimd all pimd/langevin method pimd ensemble nvt "
                            "integrator obabo thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");
        const std::string fix_command = "fix bias all plumed plumedfile " +
                                        std::string(plumed_file) + " outfile " + plumed_log +
                                        " path_integral bead_mean pimd_fix fpimd";
        lammps_command(lmp, fix_command.c_str());
        if (has_error("fix plumed")) {
            lammps_close(lmp);
            return 0.0;
        }
        lammps_command(lmp, "run 1 post no");
        if (has_error("run")) {
            lammps_close(lmp);
            return 0.0;
        }

        double result = 0.0;
        auto *bias =
            (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
        EXPECT_NE(bias, nullptr);
        if (bias) {
            result = *bias;
            lammps_free(bias);
        }
        lammps_close(lmp);
        return result;
    };

    auto read_records = [](const char *file) {
        std::vector<std::vector<double>> records;
        std::ifstream input(file);
        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream row(line);
            std::vector<double> record;
            double value;
            while (row >> value)
                record.push_back(value);
            if (!record.empty()) records.push_back(record);
        }
        return records;
    };
    auto read_fields = [](const char *file) {
        std::vector<std::string> fields;
        std::ifstream input(file);
        std::string line;
        while (std::getline(input, line)) {
            if (line.rfind("#! FIELDS", 0) != 0) continue;
            std::istringstream row(line);
            std::string marker, label, field;
            row >> marker >> label;
            while (row >> field)
                fields.push_back(field);
            break;
        }
        return fields;
    };
    auto read_state_counter = [](const char *file) {
        std::ifstream input(file);
        std::string line;
        while (std::getline(input, line)) {
            if (line.rfind("#! SET counter", 0) != 0) continue;
            std::istringstream row(line);
            std::string marker, set, key;
            int counter = -1;
            row >> marker >> set >> key >> counter;
            return counter;
        }
        return -1;
    };

    const double initial_bias = run_case(initial_file, initial_log, 1);
    EXPECT_TRUE(std::isfinite(initial_bias));
    if (me == 0)
        EXPECT_LT(initial_bias, 0.0);
    else
        EXPECT_NEAR(initial_bias, 0.0, 1.0e-12);
    MPI_Barrier(MPI_COMM_WORLD);

    const auto initial_kernels = read_records(kernels_initial);
    const auto initial_state   = read_records(state_initial);
    ASSERT_EQ(initial_kernels.size(), 4);
    ASSERT_EQ(initial_state.size(), 1);
    EXPECT_EQ(read_state_counter(state_initial), 5);
    for (const auto &record : initial_kernels)
        ASSERT_EQ(record.size(), 7);
    for (const auto &record : initial_state)
        ASSERT_EQ(record.size(), 6);
    const auto kernel_fields = read_fields(kernels_initial);
    const auto state_fields  = read_fields(state_initial);
    EXPECT_NE(std::find(kernel_fields.begin(), kernel_fields.end(), "sc"), kernel_fields.end());
    EXPECT_NE(std::find(kernel_fields.begin(), kernel_fields.end(), "spread2"),
              kernel_fields.end());
    EXPECT_NE(std::find(state_fields.begin(), state_fields.end(), "sc"), state_fields.end());
    EXPECT_NE(std::find(state_fields.begin(), state_fields.end(), "spread2"), state_fields.end());
    for (int walker = 0; walker < nprocs; ++walker) {
        for (const char *file : {kernels_initial, state_initial}) {
            std::ifstream split_file(std::string(file) + "." + std::to_string(walker));
            EXPECT_FALSE(split_file.good());
        }
        const std::string colvar = std::string(colvar_initial) + "." + std::to_string(walker);
        const auto records       = read_records(colvar.c_str());
        ASSERT_EQ(records.size(), 2);
        ASSERT_EQ(records.back().size(), 8);
        EXPECT_NEAR(records.front()[1], 7.5, 1.0e-12);
        EXPECT_NEAR(records.front()[2], 1.25, 1.0e-12);
        EXPECT_NEAR(records.back()[7], 1.0, 1.0e-12);
    }

    const double restart_bias = run_case(restart_file, restart_log, 3);
    EXPECT_TRUE(std::isfinite(restart_bias));
    if (me == 0)
        EXPECT_NE(restart_bias, 0.0);
    else
        EXPECT_NEAR(restart_bias, 0.0, 1.0e-12);
    MPI_Barrier(MPI_COMM_WORLD);

    EXPECT_EQ(read_records(state_initial), initial_state);
    const auto restarted_kernels = read_records(kernels_restart);
    const auto restarted_state   = read_records(state_restart);
    ASSERT_EQ(restarted_kernels.size(), 4);
    ASSERT_EQ(restarted_state.size(), 1);
    EXPECT_EQ(read_state_counter(state_restart), 9);
    for (const auto &record : restarted_kernels)
        ASSERT_EQ(record.size(), 7);
    for (const auto &record : restarted_state)
        ASSERT_EQ(record.size(), 6);
    for (int walker = 0; walker < nprocs; ++walker) {
        for (const char *file : {kernels_restart, state_restart}) {
            std::ifstream split_file(std::string(file) + "." + std::to_string(walker));
            EXPECT_FALSE(split_file.good());
        }
        const std::string initial_colvar =
            std::string(colvar_initial) + "." + std::to_string(walker);
        const std::string restart_colvar =
            std::string(colvar_restart) + "." + std::to_string(walker);
        const auto initial_records = read_records(initial_colvar.c_str());
        const auto restart_records = read_records(restart_colvar.c_str());
        ASSERT_EQ(initial_records.size(), 2);
        ASSERT_EQ(restart_records.size(), 2);
        ASSERT_EQ(initial_records.back().size(), 8);
        ASSERT_EQ(restart_records.front().size(), 8);
        EXPECT_NEAR(restart_records.front()[4], initial_records.back()[4], 1.0e-12);
        EXPECT_NEAR(restart_records.front()[5], initial_records.back()[5], 1.0e-12);
        EXPECT_NEAR(restart_records.front()[6], initial_records.back()[6], 1.0e-12);
        EXPECT_NEAR(restart_records.front()[7], 1.0, 1.0e-12);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        for (const char *file : {initial_file, restart_file, kernels_initial, kernels_restart,
                                 state_initial, state_restart})
            std::remove(file);
        std::remove((std::string("bck.last.") + state_initial).c_str());
        std::remove((std::string("bck.last.") + state_restart).c_str());
        for (int walker = 0; walker < nprocs; ++walker) {
            std::remove((std::string(colvar_initial) + "." + std::to_string(walker)).c_str());
            std::remove((std::string(colvar_restart) + "." + std::to_string(walker)).c_str());
        }
    }
    std::remove((std::string(initial_log) + "." + std::to_string(me)).c_str());
    std::remove((std::string(restart_log) + "." + std::to_string(me)).c_str());
    std::remove((std::string("bck.0.") + restart_log + "." + std::to_string(me)).c_str());
};

TEST(MPI, plumed_pimd_multirank_bead_modes)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    ASSERT_EQ(nprocs, 4);
    // Tables contain -grad(B); the beta/P Hamiltonian requires -P*grad(B).
    constexpr int beads = 2;

    const char *zero_file              = "test_plumed_pimd_multirank_zero.dat";
    const char *restraint_file         = "test_plumed_pimd_multirank_restraint.dat";
    const char *density_zero_file      = "test_plumed_pimd_multirank_density_zero.dat";
    const char *density_restraint_file = "test_plumed_pimd_multirank_density_restraint.dat";
    const char *density_metad_file     = "test_plumed_pimd_multirank_density_metad.dat";
    const char *density_restart_file   = "test_plumed_pimd_multirank_density_metad_restart.dat";
    const char *density_hills_file     = "test_plumed_pimd_multirank_density_metad_hills";
    const char *zero_log               = "test_plumed_pimd_multirank_zero.log";
    const char *restraint_log          = "test_plumed_pimd_multirank_restraint.log";
    const char *density_zero_log       = "test_plumed_pimd_multirank_density_zero.log";
    const char *density_restraint_log  = "test_plumed_pimd_multirank_density_restraint.log";
    const char *density_metad_log      = "test_plumed_pimd_multirank_density_metad.log";
    const char *density_restart_log    = "test_plumed_pimd_multirank_density_metad_restart.log";
    if (me == 0) {
        std::remove(density_hills_file);
        for (int bead_index = 0; bead_index < 2; ++bead_index)
            std::remove(
                (std::string(density_hills_file) + "." + std::to_string(bead_index)).c_str());
        std::ofstream zero(zero_file);
        zero << "d: DISTANCE ATOMS=1,2 NOPBC\n"
             << "mean: ENSEMBLE ARG=d\n";
        std::ofstream restraint(restraint_file);
        restraint << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                  << "mean: ENSEMBLE ARG=d\n"
                  << "bias: RESTRAINT ARG=mean.d AT=7.5 KAPPA=4.0\n";
        std::ofstream density_zero(density_zero_file);
        density_zero << "d: DISTANCE ATOMS=1,2 NOPBC\n";
        std::ofstream density_restraint(density_restraint_file);
        density_restraint << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                          << "bias: RESTRAINT ARG=d AT=7.5 KAPPA=4.0\n";
        std::ofstream density_metad(density_metad_file);
        density_metad << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                      << "bias: METAD ARG=d SIGMA=0.5 HEIGHT=0.5 PACE=2 FILE=" << density_hills_file
                      << " WALKERS_MPI RESTART=NO\n";
        std::ofstream density_restart(density_restart_file);
        density_restart << "d: DISTANCE ATOMS=1,2 NOPBC\n"
                        << "bias: METAD ARG=d SIGMA=0.5 HEIGHT=0.5 PACE=2 FILE="
                        << density_hills_file << " WALKERS_MPI RESTART=YES\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const int bead = me / 2;
    ASSERT_GE(bead, 0);
    ASSERT_LT(bead, 2);
    auto run_case = [&](const char *plumed_file, const char *plumed_log, const char *mode,
                        int step) {
        std::array<double, 7> result{};
        void *lmp = open_multirank_partition();
        EXPECT_NE(lmp, nullptr);
        if (!lmp) return result;

        EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 2);
        create_multirank_two_atom_system(lmp);
        auto *nlocal = (int *)lammps_extract_global(lmp, "nlocal");
        EXPECT_NE(nlocal, nullptr);
        int zero_atom_ranks = nlocal && *nlocal == 0;
        MPI_Allreduce(MPI_IN_PLACE, &zero_atom_ranks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        EXPECT_GT(zero_atom_ranks, 0);
        lammps_command(lmp, "group first id 1");
        lammps_command(lmp, "group second id 2");
        lammps_command(lmp, "compute first_force first reduce sum fx fy fz");
        lammps_command(lmp, "compute second_force second reduce sum fx fy fz");
        if (step >= 0) lammps_command(lmp, ("reset_timestep " + std::to_string(step)).c_str());
        lammps_command(lmp, "fix fpimd all pimd/langevin method pimd ensemble nvt "
                            "integrator obabo thermostat PILE_L 1234 tau 1.0 temp 1.0 fixcom no");
        const std::string fix_command = "fix bias all plumed plumedfile " +
                                        std::string(plumed_file) + " outfile " + plumed_log +
                                        " path_integral " + mode + " pimd_fix fpimd";
        lammps_command(lmp, fix_command.c_str());
        lammps_command(lmp, step >= 0 ? "run 1 post no" : "run 0");

        auto *first =
            (double *)lammps_extract_compute(lmp, "first_force", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
        auto *second = (double *)lammps_extract_compute(lmp, "second_force", LMP_STYLE_GLOBAL,
                                                        LMP_TYPE_VECTOR);
        EXPECT_NE(first, nullptr);
        EXPECT_NE(second, nullptr);
        if (first && second) {
            for (int dimension = 0; dimension < 3; ++dimension) {
                result[dimension]     = first[dimension];
                result[3 + dimension] = second[dimension];
            }
        }
        auto *bias =
            (double *)lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
        EXPECT_NE(bias, nullptr);
        if (bias) {
            result[6] = *bias;
            lammps_free(bias);
        }
        lammps_close(lmp);
        return result;
    };

    const auto zero_result   = run_case(zero_file, zero_log, "bead_mean", -1);
    const auto biased_result = run_case(restraint_file, restraint_log, "bead_mean", -1);
    EXPECT_NEAR(zero_result[6], 0.0, 1.0e-12);
    EXPECT_NEAR(biased_result[6], bead == 0 ? 0.5 : 0.0, 1.0e-12);

    const std::array<std::array<double, 6>, 2> expected_force_delta = {
        {{1.0, 0.0, 0.0, -1.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0, -1.0, 0.0}}};
    for (std::size_t i = 0; i < expected_force_delta[bead].size(); ++i)
        EXPECT_NEAR(biased_result[i] - zero_result[i], beads * expected_force_delta[bead][i],
                    1.0e-12);

    const auto density_zero_result =
        run_case(density_zero_file, density_zero_log, "bead_density", -1);
    const auto density_biased_result =
        run_case(density_restraint_file, density_restraint_log, "bead_density", -1);
    EXPECT_NEAR(density_zero_result[6], 0.0, 1.0e-12);
    EXPECT_NEAR(density_biased_result[6], bead == 0 ? 2.5 : 0.0, 1.0e-12);

    const std::array<std::array<double, 6>, 2> expected_density_force_delta = {
        {{3.0, 0.0, 0.0, -3.0, 0.0, 0.0}, {0.0, -1.0, 0.0, 0.0, 1.0, 0.0}}};
    for (std::size_t i = 0; i < expected_density_force_delta[bead].size(); ++i)
        EXPECT_NEAR(density_biased_result[i] - density_zero_result[i],
                    beads * expected_density_force_delta[bead][i], 1.0e-12);

    auto read_hills = [&]() {
        std::vector<std::vector<double>> records;
        std::ifstream input(density_hills_file);
        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream row(line);
            std::vector<double> record;
            double value;
            while (row >> value)
                record.push_back(value);
            if (!record.empty()) records.push_back(record);
        }
        return records;
    };

    const auto initial_metad_result =
        run_case(density_metad_file, density_metad_log, "bead_density", 1);
    EXPECT_NEAR(initial_metad_result[6], 0.0, 1.0e-12);
    MPI_Barrier(MPI_COMM_WORLD);

    const auto initial_hills = read_hills();
    ASSERT_EQ(initial_hills.size(), 2);
    for (const auto &record : initial_hills) {
        ASSERT_EQ(record.size(), 5);
        EXPECT_NEAR(record[3], 0.5, 1.0e-12);
    }
    for (int bead_index = 0; bead_index < 2; ++bead_index) {
        std::ifstream split_hills(std::string(density_hills_file) + "." +
                                  std::to_string(bead_index));
        EXPECT_FALSE(split_hills.good());
    }

    const auto restarted_metad_result =
        run_case(density_restart_file, density_restart_log, "bead_density", 3);
    if (bead == 0)
        EXPECT_GT(restarted_metad_result[6], 0.0);
    else
        EXPECT_NEAR(restarted_metad_result[6], 0.0, 1.0e-12);
    MPI_Barrier(MPI_COMM_WORLD);

    const auto restarted_hills = read_hills();
    ASSERT_EQ(restarted_hills.size(), 4);
    for (std::size_t i = 0; i < restarted_hills.size(); ++i) {
        ASSERT_EQ(restarted_hills[i].size(), 5);
        EXPECT_NEAR(restarted_hills[i][3], 0.5, 1.0e-12);
        if (i < initial_hills.size()) EXPECT_EQ(restarted_hills[i], initial_hills[i]);
    }
    EXPECT_GT(restarted_hills[2][0], restarted_hills[1][0]);

    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        std::remove(zero_file);
        std::remove(restraint_file);
        std::remove(density_zero_file);
        std::remove(density_restraint_file);
        std::remove(density_metad_file);
        std::remove(density_restart_file);
        std::remove(density_hills_file);
        for (int bead_index = 0; bead_index < 2; ++bead_index) {
            std::remove((std::string(zero_log) + "." + std::to_string(bead_index)).c_str());
            std::remove((std::string(restraint_log) + "." + std::to_string(bead_index)).c_str());
            std::remove((std::string(density_zero_log) + "." + std::to_string(bead_index)).c_str());
            std::remove(
                (std::string(density_restraint_log) + "." + std::to_string(bead_index)).c_str());
            std::remove(
                (std::string(density_metad_log) + "." + std::to_string(bead_index)).c_str());
            std::remove(
                (std::string(density_restart_log) + "." + std::to_string(bead_index)).c_str());
        }
    }
};
