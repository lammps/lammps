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

// unit tests for the output data of compute and fix styles
//
// Each YAML file sets up a test system (input_file plus optional
// input_coeffs) and defines - via post_commands - a compute or fix with
// the ID "test" (plus any helpers like groups, chunks, variables, or
// other computes it consumes).  The driver runs for a fixed number of
// steps and then collects whatever output the style provides, based on
// its output flags: global scalar (global_scalar), global vector
// (global_vector), global array (global_array), per-atom vector or array
// (peratom_data, first column is the atom tag), and local vector or
// array (local_data).  The same code doubles as the reference generator.
//
// Styles whose output requires per-atom energy or virial tallies during
// the run (e.g. compute pe/atom, stress/atom, or heat/flux) are supported
// by scheduling the tallies for the final step of the run: every compute
// with the peatomflag or pressatomflag set - the tested compute or any
// helper compute it consumes - gets an addstep() request, which makes the
// integrator enable the per-atom energy/virial accumulation on that step.

#include "error_stats.h"
#include "test_config.h"
#include "test_main.h"
#include "yaml_writer.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "atom.h"
#include "compute.h"
#include "fix.h"
#include "info.h"
#include "utils.h"
#include "input.h"
#include "modify.h"
#include "update.h"

#include <exception>
#include <iostream>
#include <vector>

using ::testing::HasSubstr;
using ::testing::StartsWith;

using namespace LAMMPS_NS;

// the "kokkos_omp_full" and "kokkos_serial_full" test cases select "newton off"
// through the newton_pair and newton_bond index variables on the command line.
// several YAML files redefine those variables in their pre_commands (a
// convention taken over from the GPU package tests), which discards the command
// line setting, so the override has to be re-applied after the pre_commands
// have been processed and before the input template is read.

static bool kokkos_full_neigh = false;

static void enforce_kokkos_full_neigh(LAMMPS *lmp)
{
    if (!kokkos_full_neigh) return;
    lmp->input->one("variable newton_pair delete");
    lmp->input->one("variable newton_pair index off");
    lmp->input->one("variable newton_bond delete");
    lmp->input->one("variable newton_bond index off");
}

// styles that require "newton on" or a half neighbor list cannot run in the
// full neighbor list configuration of the KOKKOS package.  those are
// documented restrictions of the style, so the corresponding test case is
// skipped instead of failed when the setup stops with such an error.

static bool full_neigh_unsupported(const std::string &errmsg)
{
    // a style that refuses the neighbor list or the newton setting asked for
    // skips rather than fails.  that applies to the full neighbor list cases
    // and to any accelerator settings supplied through LAMMPS_KOKKOS_ARGS: a
    // run of the suite under a different "package kokkos" profile is meant to
    // report what a style does support, not to fail on what it does not
    const char *extra = std::getenv("LAMMPS_KOKKOS_ARGS");
    if (!kokkos_full_neigh && (!extra || (extra[0] == '\0'))) return false;
    return (LAMMPS_NS::utils::strmatch(errmsg, "newton") ||
            LAMMPS_NS::utils::strmatch(errmsg, "half neighbor list"));
}

// fixed number of MD steps before collecting the output, so that
// history-dependent styles (msd, vacf, ave/time, ave/atom) have data
static constexpr int RUN_STEPS = 10;

static void cleanup_lammps(LAMMPS *&lmp)
{
    delete lmp;
    lmp = nullptr;
}

// normalized snapshot of all output data a compute or fix provides

struct OutputData {
    bool has_scalar = false;
    double scalar   = 0.0;
    std::vector<double> vector;
    std::vector<std::vector<double>> array;
    std::vector<std::vector<double>> peratom;
    std::vector<std::vector<double>> local;
};

static void collect_peratom(Atom *atom, double *vec, double **arr, int ncols, OutputData &data)
{
    const int natoms = (int)atom->natoms;
    for (int i = 1; i <= natoms; ++i) {
        const int j = atom->map(i);
        std::vector<double> row;
        row.push_back(i);
        if (ncols == 0)
            row.push_back(vec[j]);
        else
            for (int c = 0; c < ncols; ++c)
                row.push_back(arr[j][c]);
        data.peratom.push_back(row);
    }
}

static OutputData collect_compute(Compute *icompute, Atom *atom)
{
    OutputData data;
    if (icompute->scalar_flag) {
        data.has_scalar = true;
        data.scalar     = icompute->compute_scalar();
    }
    if (icompute->vector_flag) {
        icompute->compute_vector();
        for (int i = 0; i < icompute->size_vector; ++i)
            data.vector.push_back(icompute->vector[i]);
    }
    if (icompute->array_flag) {
        // for variable-size arrays (e.g. per-chunk computes) the row count
        // is only valid after the array was computed
        icompute->compute_array();
        for (int i = 0; i < icompute->size_array_rows; ++i) {
            std::vector<double> row;
            for (int j = 0; j < icompute->size_array_cols; ++j)
                row.push_back(icompute->array[i][j]);
            data.array.push_back(row);
        }
    }
    if (icompute->peratom_flag) {
        icompute->compute_peratom();
        collect_peratom(atom, icompute->vector_atom, icompute->array_atom,
                        icompute->size_peratom_cols, data);
    }
    if (icompute->local_flag) {
        icompute->compute_local();
        const int ncols = icompute->size_local_cols;
        for (int i = 0; i < icompute->size_local_rows; ++i) {
            std::vector<double> row;
            if (ncols == 0)
                row.push_back(icompute->vector_local[i]);
            else
                for (int j = 0; j < ncols; ++j)
                    row.push_back(icompute->array_local[i][j]);
            data.local.push_back(row);
        }
    }
    return data;
}

static OutputData collect_fix(Fix *ifix, Atom *atom)
{
    OutputData data;
    if (ifix->scalar_flag) {
        data.has_scalar = true;
        data.scalar     = ifix->compute_scalar();
    }
    if (ifix->vector_flag) {
        for (int i = 0; i < ifix->size_vector; ++i)
            data.vector.push_back(ifix->compute_vector(i));
    }
    if (ifix->array_flag) {
        // for variable-size arrays (e.g. fix ave/chunk) the row count is
        // determined by the fix during the run
        for (int i = 0; i < ifix->size_array_rows; ++i) {
            std::vector<double> row;
            for (int j = 0; j < ifix->size_array_cols; ++j)
                row.push_back(ifix->compute_array(i, j));
            data.array.push_back(row);
        }
    }
    if (ifix->peratom_flag) {
        collect_peratom(atom, ifix->vector_atom, ifix->array_atom, ifix->size_peratom_cols, data);
    }
    return data;
}

static LAMMPS *init_lammps(LAMMPS::argv &args, const TestConfig &cfg)
{
    LAMMPS *lmp = new LAMMPS(args, MPI_COMM_WORLD);

    // check if prerequisite styles are available
    Info *info = new Info(lmp);
    int nfail  = 0;
    for (const auto &prerequisite : cfg.prerequisites) {
        if (!info->has_style(prerequisite.first, prerequisite.second)) ++nfail;
    }
    delete info;
    if (nfail > 0) {
        cleanup_lammps(lmp);
        return nullptr;
    }

    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    command("variable input_dir index " + INPUT_FOLDER);
    for (const auto &pre_command : cfg.pre_commands)
        command(pre_command);
    enforce_kokkos_full_neigh(lmp);

    std::string input_file = platform::path_join(INPUT_FOLDER, cfg.input_file);
    lmp->input->file(input_file.c_str());

    // optional force field setup from a coeffs file (molecular templates);
    // systems that configure the force field in post_commands omit it

    if (!cfg.input_coeffs.empty()) {
        std::string coeffs_file = platform::path_join(INPUT_FOLDER, cfg.input_coeffs);
        lmp->input->file(coeffs_file.c_str());
    }

    // the compute or fix under test (ID "test") and any helper commands

    for (const auto &post_command : cfg.post_commands)
        command(post_command);

    // time integration so that history-dependent output has data to report;
    // the templates provide initial velocities.  skipped when the test
    // defines its own time-integrating fix (e.g. fix rigid)
    bool has_integrator = false;
    for (const auto *ifix : lmp->modify->get_fix_list())
        if (ifix->time_integrate) has_integrator = true;
    if (!has_integrator) command("fix output_nve all nve");

    // schedule per-atom energy/virial tallies for the final step of the
    // run, so that computes consuming them (pe/atom, stress/atom,
    // centroid/stress/atom, heat/flux, ...) can be collected after the
    // run.  this covers the tested compute and any helper computes.
    const bigint laststep = lmp->update->ntimestep + RUN_STEPS;
    for (auto *icompute : lmp->modify->get_compute_list())
        if (icompute->peatomflag || icompute->pressatomflag) icompute->addstep(laststep);
    command(fmt::format("timestep {}", (cfg.timestep > 0.0) ? cfg.timestep : 0.25));
    command("thermo 5");
    command(fmt::format("run {} post no", RUN_STEPS));
    return lmp;
}

// collect the output of the compute or fix with ID "test".
// returns false when neither exists.

static bool collect_output(LAMMPS *lmp, OutputData &data)
{
    auto *icompute = lmp->modify->get_compute_by_id("test");
    if (icompute) {
        data = collect_compute(icompute, lmp->atom);
        return true;
    }
    auto *ifix = lmp->modify->get_fix_by_id("test");
    if (ifix) {
        data = collect_fix(ifix, lmp->atom);
        return true;
    }
    return false;
}

static void compare_rows(const std::string &name,
                         const std::vector<std::vector<double>> &reference,
                         const std::vector<std::vector<double>> &current, double epsilon,
                         ErrorStats &stats)
{
    SCOPED_TRACE(name);
    ASSERT_EQ(reference.size(), current.size());
    for (std::size_t i = 0; i < reference.size(); ++i) {
        ASSERT_EQ(reference[i].size(), current[i].size());
        for (std::size_t j = 0; j < reference[i].size(); ++j) {
            EXPECT_FP_LE_WITH_EPS(current[i][j], reference[i][j], epsilon);
        }
    }
}

// append the words of the LAMMPS_KOKKOS_ARGS environment variable to the
// command line of the KOKKOS test cases.  this lets the whole suite be re-run
// with the "package kokkos" settings a GPU would choose --
// LAMMPS_KOKKOS_ARGS="-pk kokkos comm device sort device atom/map device gpu/aware on"
// -- which is what the host/device transfer checking of a build configured with
// -D KOKKOS_DEBUG_SYNC=on needs in order to see anything

static void append_kokkos_env_args(LAMMPS_NS::LAMMPS::argv &args)
{
    const char *extra = std::getenv("LAMMPS_KOKKOS_ARGS");
    if (!extra || (extra[0] == '\0')) return;
    auto words = LAMMPS_NS::utils::split_words(extra);
    args.insert(args.end(), words.begin(), words.end());
}

static void run_output_test(LAMMPS::argv &args, double epsilon, bool kokkos)
{
    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config);
    } catch (std::exception &e) {
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        if (full_neigh_unsupported(e.what())) GTEST_SKIP() << e.what();
        FAIL() << e.what();
    }
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    if (!lmp) {
        std::cerr << "One or more prerequisite styles are not available "
                     "in this LAMMPS configuration:\n";
        for (auto &prerequisite : test_config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        GTEST_SKIP();
    }

    EXPECT_THAT(output, StartsWith("LAMMPS ("));

    // init_lammps() always runs the system for RUN_STEPS steps, so the
    // timing summary of that run has to be part of the output
    if (kokkos) EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    ASSERT_EQ(lmp->atom->natoms, lmp->atom->nlocal);

    OutputData data;
    if (!collect_output(lmp, data)) {
        cleanup_lammps(lmp);
        FAIL() << "no compute or fix with ID 'test' defined";
    }

    ErrorStats stats;

    if (data.has_scalar) EXPECT_FP_LE_WITH_EPS(data.scalar, test_config.global_scalar, epsilon);

    {
        SCOPED_TRACE("global vector");
        ASSERT_EQ(test_config.global_vector.size(), data.vector.size());
        for (std::size_t i = 0; i < data.vector.size(); ++i)
            EXPECT_FP_LE_WITH_EPS(data.vector[i], test_config.global_vector[i], epsilon);
    }

    compare_rows("global array", test_config.global_array, data.array, epsilon, stats);
    compare_rows("per-atom data", test_config.peratom_data, data.peratom, epsilon, stats);
    compare_rows("local data", test_config.local_data, data.local, epsilon, stats);

    if (print_stats) std::cerr << "output stats:" << stats << std::endl;

    cleanup_lammps(lmp);
}

TEST(OutputStyle, plain)
{
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"OutputStyle", "-log", "none", "-echo", "screen", "-nocite"};

    run_output_test(args, test_config.epsilon, false);
}

// precision of the KOKKOS package as selected with -D KOKKOS_PREC at compile time
static std::string kokkos_precision()
{
    if (Info::has_accelerator_feature("KOKKOS", "precision", "mixed")) return "mixed";
    if (Info::has_accelerator_feature("KOKKOS", "precision", "single")) return "single";
    return "double";
}

// the KOKKOS package accumulates in a different order and - depending on how it
// was compiled - with reduced precision, so the tolerance has to be relaxed
static double kokkos_epsilon()
{
    double epsilon                 = 5.0 * test_config.epsilon;
    const std::string kk_precision = kokkos_precision();
    if (kk_precision == "mixed")
        epsilon *= 2.0e9;
    else if (kk_precision == "single")
        epsilon *= 1.0e10;
    return epsilon;
}

// the KOKKOS tests below use the same prerequisites as the plain test, i.e. no
// "/kk" suffix is appended to them.  Unlike the other force style test drivers
// this one has no single tested style category, and a compute or fix without a
// KOKKOS variant is still worth running inside a KOKKOS enabled run: it
// exercises the automatic synchronization between the host and device copies
// of the atom data.

TEST(OutputStyle, kokkos_omp)
{
    if (!Info::has_package("KOKKOS")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    // skip entries may also be qualified by the KOKKOS package precision,
    // e.g. "kokkos_omp_single" skips only single precision KOKKOS builds
    if (test_config.skip_tests.count(std::string(test_info_->name()) + "_" + kokkos_precision()))
        GTEST_SKIP();
    // skip entries qualified with "_devicerng" apply only to builds where the
    // KOKKOS styles use the device random number generator
    if (Info::has_accelerator_feature("KOKKOS", "rng", "device") &&
        test_config.skip_tests.count(std::string(test_info_->name()) + "_devicerng"))
        GTEST_SKIP();
    // this test requires the OpenMP backend of KOKKOS
    if (!Info::has_accelerator_feature("KOKKOS", "api", "openmp"))
        GTEST_SKIP() << "KOKKOS OpenMP backend not enabled";
    // if KOKKOS has GPU support enabled, it *must* be used. We cannot test OpenMP only.
    if (Info::has_accelerator_feature("KOKKOS", "api", "cuda") ||
        Info::has_accelerator_feature("KOKKOS", "api", "hip") ||
        Info::has_accelerator_feature("KOKKOS", "api", "sycl")) {
        GTEST_SKIP() << "Cannot test KOKKOS/OpenMP with GPU support enabled";
    }

    LAMMPS::argv args = {"OutputStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k",          "on",   "t",    "4",     "-sf",    "kk"};

    append_kokkos_env_args(args);
    run_output_test(args, kokkos_epsilon(), true);
}

TEST(OutputStyle, kokkos_omp_full)
{
    if (!Info::has_package("KOKKOS")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    // skip entries may also be qualified by the KOKKOS package precision,
    // e.g. "kokkos_omp_full_single" skips only single precision KOKKOS builds
    if (test_config.skip_tests.count(std::string(test_info_->name()) + "_" + kokkos_precision()))
        GTEST_SKIP();
    // skip entries qualified with "_devicerng" apply only to builds where the
    // KOKKOS styles use the device random number generator
    if (Info::has_accelerator_feature("KOKKOS", "rng", "device") &&
        test_config.skip_tests.count(std::string(test_info_->name()) + "_devicerng"))
        GTEST_SKIP();
    // a style that cannot be tested with KOKKOS at all cannot be tested
    // with a full neighbor list either, so the plain "kokkos_omp"
    // skip entries apply here as well
    if (test_config.skip_tests.count("kokkos_omp")) GTEST_SKIP();
    if (test_config.skip_tests.count("kokkos_omp_" + kokkos_precision()))
        GTEST_SKIP();
    if (Info::has_accelerator_feature("KOKKOS", "rng", "device") &&
        test_config.skip_tests.count("kokkos_omp_devicerng"))
        GTEST_SKIP();
    // this test requires the OpenMP backend of KOKKOS
    if (!Info::has_accelerator_feature("KOKKOS", "api", "openmp"))
        GTEST_SKIP() << "KOKKOS OpenMP backend not enabled";
    // if KOKKOS has GPU support enabled, it *must* be used. We cannot test OpenMP only.
    if (Info::has_accelerator_feature("KOKKOS", "api", "cuda") ||
        Info::has_accelerator_feature("KOKKOS", "api", "hip") ||
        Info::has_accelerator_feature("KOKKOS", "api", "sycl")) {
        GTEST_SKIP() << "Cannot test KOKKOS/OpenMP with GPU support enabled";
    }

    // exercise the NEIGHFLAG == FULL kernels of the KOKKOS package.  those are
    // what the GPU backends select by default, but they are never reached in a
    // CPU only test build, which always uses a half neighbor list with newton
    // on.  the KOKKOS package requires "newton off" with "neigh full", so the
    // newton settings of the input template must be overridden as well: an
    // index style variable defined with -var on the command line takes
    // precedence over the "variable ... index" definition inside the template
    LAMMPS::argv args = {"OutputStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k", "on", "t", "4", "-sf", "kk",
                         "-pk", "kokkos", "neigh", "full", "newton", "off",
                         "-var", "newton_pair", "off", "-var", "newton_bond", "off"};

    kokkos_full_neigh = true;
    append_kokkos_env_args(args);
    run_output_test(args, kokkos_epsilon(), true);
    kokkos_full_neigh = false;
}

TEST(OutputStyle, kokkos_serial)
{
    if (!Info::has_package("KOKKOS")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    // skip entries may also be qualified by the KOKKOS package precision,
    // e.g. "kokkos_serial_single" skips only single precision KOKKOS builds
    if (test_config.skip_tests.count(std::string(test_info_->name()) + "_" + kokkos_precision()))
        GTEST_SKIP();
    // skip entries qualified with "_devicerng" apply only to builds where the
    // KOKKOS styles use the device random number generator
    if (Info::has_accelerator_feature("KOKKOS", "rng", "device") &&
        test_config.skip_tests.count(std::string(test_info_->name()) + "_devicerng"))
        GTEST_SKIP();
    // this test requires the KOKKOS package compiled with only the Serial backend: when the
    // OpenMP (or a GPU) backend is enabled, the host execution space is not Serial
    if (!Info::has_accelerator_feature("KOKKOS", "api", "serial"))
        GTEST_SKIP() << "KOKKOS Serial backend not enabled";
    if (Info::has_accelerator_feature("KOKKOS", "api", "openmp") ||
        Info::has_accelerator_feature("KOKKOS", "api", "pthreads"))
        GTEST_SKIP() << "Cannot test KOKKOS/Serial with threading support enabled";
    if (Info::has_accelerator_feature("KOKKOS", "api", "cuda") ||
        Info::has_accelerator_feature("KOKKOS", "api", "hip") ||
        Info::has_accelerator_feature("KOKKOS", "api", "sycl")) {
        GTEST_SKIP() << "Cannot test KOKKOS/Serial with GPU support enabled";
    }

    LAMMPS::argv args = {"OutputStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k",          "on",   "t",    "1",     "-sf",    "kk"};

    append_kokkos_env_args(args);
    run_output_test(args, kokkos_epsilon(), true);
}

TEST(OutputStyle, kokkos_serial_full)
{
    if (!Info::has_package("KOKKOS")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    // skip entries may also be qualified by the KOKKOS package precision,
    // e.g. "kokkos_serial_full_single" skips only single precision KOKKOS builds
    if (test_config.skip_tests.count(std::string(test_info_->name()) + "_" + kokkos_precision()))
        GTEST_SKIP();
    // skip entries qualified with "_devicerng" apply only to builds where the
    // KOKKOS styles use the device random number generator
    if (Info::has_accelerator_feature("KOKKOS", "rng", "device") &&
        test_config.skip_tests.count(std::string(test_info_->name()) + "_devicerng"))
        GTEST_SKIP();
    // a style that cannot be tested with KOKKOS at all cannot be tested
    // with a full neighbor list either, so the plain "kokkos_serial"
    // skip entries apply here as well
    if (test_config.skip_tests.count("kokkos_serial")) GTEST_SKIP();
    if (test_config.skip_tests.count("kokkos_serial_" + kokkos_precision()))
        GTEST_SKIP();
    if (Info::has_accelerator_feature("KOKKOS", "rng", "device") &&
        test_config.skip_tests.count("kokkos_serial_devicerng"))
        GTEST_SKIP();
    // this test requires the KOKKOS package compiled with only the Serial backend: when the
    // OpenMP (or a GPU) backend is enabled, the host execution space is not Serial
    if (!Info::has_accelerator_feature("KOKKOS", "api", "serial"))
        GTEST_SKIP() << "KOKKOS Serial backend not enabled";
    if (Info::has_accelerator_feature("KOKKOS", "api", "openmp") ||
        Info::has_accelerator_feature("KOKKOS", "api", "pthreads"))
        GTEST_SKIP() << "Cannot test KOKKOS/Serial with threading support enabled";
    if (Info::has_accelerator_feature("KOKKOS", "api", "cuda") ||
        Info::has_accelerator_feature("KOKKOS", "api", "hip") ||
        Info::has_accelerator_feature("KOKKOS", "api", "sycl")) {
        GTEST_SKIP() << "Cannot test KOKKOS/Serial with GPU support enabled";
    }

    // exercise the NEIGHFLAG == FULL kernels of the KOKKOS package.  those are
    // what the GPU backends select by default, but they are never reached in a
    // CPU only test build, which always uses a half neighbor list with newton
    // on.  the KOKKOS package requires "newton off" with "neigh full", so the
    // newton settings of the input template must be overridden as well: an
    // index style variable defined with -var on the command line takes
    // precedence over the "variable ... index" definition inside the template
    LAMMPS::argv args = {"OutputStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k", "on", "t", "1", "-sf", "kk",
                         "-pk", "kokkos", "neigh", "full", "newton", "off",
                         "-var", "newton_pair", "off", "-var", "newton_bond", "off"};

    kokkos_full_neigh = true;
    append_kokkos_env_args(args);
    run_output_test(args, kokkos_epsilon(), true);
    kokkos_full_neigh = false;
}

TEST(OutputStyle, kokkos_gpu)
{
    if (!Info::has_package("KOKKOS")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    // skip entries may also be qualified by the KOKKOS package precision,
    // e.g. "kokkos_gpu_single" skips only single precision KOKKOS builds
    if (test_config.skip_tests.count(std::string(test_info_->name()) + "_" + kokkos_precision()))
        GTEST_SKIP();
    // skip entries qualified with "_devicerng" apply only to builds where the
    // KOKKOS styles use the device random number generator
    if (Info::has_accelerator_feature("KOKKOS", "rng", "device") &&
        test_config.skip_tests.count(std::string(test_info_->name()) + "_devicerng"))
        GTEST_SKIP();
    // this test requires a GPU backend of the KOKKOS package
    if (!Info::has_accelerator_feature("KOKKOS", "api", "cuda") &&
        !Info::has_accelerator_feature("KOKKOS", "api", "hip") &&
        !Info::has_accelerator_feature("KOKKOS", "api", "sycl"))
        GTEST_SKIP() << "KOKKOS GPU backend not enabled";
    // transparently skip when no compatible GPU device is present
    if (!Info::has_kokkos_gpu_device()) GTEST_SKIP() << "No compatible GPU device available";

    // use a half neighbor list so the GPU kernels run with the input's default
    // "newton on"; with the default "neigh full" the KOKKOS package requires
    // newton off, which the force-style input templates do not use
    LAMMPS::argv args = {"OutputStyle", "-log",   "none",  "-echo", "screen", "-nocite",
                         "-k",          "on",     "g",     "1",     "-sf",    "kk",
                         "-pk",         "kokkos", "neigh", "half",  "newton", "on"};

    append_kokkos_env_args(args);
    run_output_test(args, kokkos_epsilon(), true);
}

static void write_rows(YamlWriter &writer, const std::string &key,
                       const std::vector<std::vector<double>> &rows)
{
    if (rows.empty()) return;
    std::string block;
    for (const auto &row : rows) {
        std::string line;
        for (const auto &value : row)
            line += fmt::format(" {:23.16e}", value);
        block += line.substr(1) + "\n";
    }
    writer.emit_block(key, block);
}

void generate_yaml_file(const char *outfile, const TestConfig &config)
{
    // initialize system geometry
    LAMMPS::argv args = {"OutputStyle", "-log", "none", "-echo", "screen", "-nocite"};
    LAMMPS *lmp       = nullptr;
    try {
        lmp = init_lammps(args, config);
    } catch (std::exception &e) {
        // must abort instead of writing a reference file without data
        std::cerr << "ERROR: system setup failed: " << e.what() << "\n";
        exit(1);
    }
    if (!lmp) {
        std::cerr << "One or more prerequisite styles are not available "
                     "in this LAMMPS configuration:\n";
        for (auto prerequisite : config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        return;
    }

    OutputData data;
    if (!collect_output(lmp, data)) {
        std::cerr << "ERROR: no compute or fix with ID 'test' defined\n";
        cleanup_lammps(lmp);
        exit(1);
    }

    YamlWriter writer(outfile);

    // write yaml header
    write_yaml_header(&writer, &test_config, lmp->version);

    writer.emit("natoms", (int)lmp->atom->natoms);

    if (data.has_scalar) writer.emit("global_scalar", data.scalar);

    if (!data.vector.empty()) {
        std::string block = std::to_string(data.vector.size());
        for (const auto &value : data.vector)
            block += fmt::format(" {:23.16e}", value);
        writer.emit_block("global_vector", block);
    }

    write_rows(writer, "global_array", data.array);
    write_rows(writer, "peratom_data", data.peratom);
    write_rows(writer, "local_data", data.local);

    cleanup_lammps(lmp);
}
