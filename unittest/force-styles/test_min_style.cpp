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

// unit tests for minimizer styles (min_style command and min_modify options)
//
// Each YAML file sets up the molecular test system (input_file plus
// input_coeffs, like the fix-timestep tests), selects the minimizer via
// post_commands (min_style, min_modify, and optional fixes like box/relax),
// and stores reference data taken after a minimization with a fixed
// iteration budget: the potential energy before (init_energy) and after
// (run_energy) the minimization, the total force norm (global_scalar), and
// the box dimensions and tilt factors (global_vector).  Using a fixed
// number of iterations instead of a convergence tolerance keeps the
// reference data deterministic.
//
// The line search and step acceptance logic of the minimizers branches on
// floating-point comparisons, so differences in the last bits of the
// forces (from a different compiler or math library) can send the descent
// along a different but equally valid path.  Only observables that are
// stable against such path changes are compared to reference data: the
// potential energy (all paths descend into the same basin, so the final
// energies agree far more closely than any two different minimizers do)
// and the box dimensions.  The force norm is checked as an upper bound
// only and per-atom positions are not compared at all.

#include "error_stats.h"
#include "test_config.h"
#include "test_main.h"
#include "yaml_writer.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "info.h"
#include "input.h"
#include "modify.h"

#include <cmath>
#include <exception>
#include <stdexcept>
#include <iostream>

using ::testing::HasSubstr;
using ::testing::StartsWith;

using namespace LAMMPS_NS;

// fixed iteration budget: enough to exercise the algorithms and the line
// searches while keeping the tests fast and the reference deterministic
static constexpr int MAX_ITER = 100;
static constexpr int MAX_EVAL = 10000;

static void cleanup_lammps(LAMMPS *&lmp)
{
    delete lmp;
    lmp = nullptr;
}

static double potential_energy(LAMMPS *lmp)
{
    auto *pe = lmp->modify->get_compute_by_id("thermo_pe");
    return pe ? pe->compute_scalar() : 0.0;
}

// total force norm computed independently of the Min classes.
// the tests run in serial, so all atoms are local.

static double force_norm(LAMMPS *lmp)
{
    double **f = lmp->atom->f;
    double sumsq = 0.0;
    for (int i = 0; i < lmp->atom->nlocal; ++i)
        sumsq += f[i][0] * f[i][0] + f[i][1] * f[i][1] + f[i][2] * f[i][2];
    return sqrt(sumsq);
}

static void box_dims(LAMMPS *lmp, double dims[6])
{
    dims[0] = lmp->domain->xprd;
    dims[1] = lmp->domain->yprd;
    dims[2] = lmp->domain->zprd;
    dims[3] = lmp->domain->xy;
    dims[4] = lmp->domain->xz;
    dims[5] = lmp->domain->yz;
}

static LAMMPS *init_lammps(LAMMPS::argv &args, const TestConfig &cfg)
{
    LAMMPS *lmp = new LAMMPS(args, MPI_COMM_WORLD);

    // check if prerequisite styles are available
    Info *info = new Info(lmp);
    int nfail  = 0;
    for (const auto &prerequisite : cfg.prerequisites) {
        std::string style = prerequisite.second;

        // this is a test for minimizer styles, so if the suffixed
        // version is not available, there is no reason to test.
        if (prerequisite.first == "minimize") {
            if (lmp->suffix_enable) {
                style += "/";
                style += lmp->suffix;
            }
        }

        if (!info->has_style(prerequisite.first, style)) ++nfail;
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

    std::string input_file = platform::path_join(INPUT_FOLDER, cfg.input_file);
    lmp->input->file(input_file.c_str());

    // set up the molecular system force field from the coeffs file
    // indicated by the YAML file (the input template only defines the geometry)

    // a missing coeffs file is a malformed YAML file, not a missing
    // prerequisite: fail loudly instead of returning nullptr, which the
    // callers would report as a (silently passing) skipped test
    if (cfg.input_coeffs.empty()) {
        cleanup_lammps(lmp);
        throw std::runtime_error("no 'input_coeffs' file given in the YAML file");
    }
    std::string coeffs_file = platform::path_join(INPUT_FOLDER, cfg.input_coeffs);
    lmp->input->file(coeffs_file.c_str());

    // the minimizer selection (min_style, min_modify) and optional
    // fixes (e.g. box/relax) come from the post_commands block

    for (const auto &post_command : cfg.post_commands)
        command(post_command);

    // the timestep matters for the damped-dynamics minimizers (fire,
    // abcfire, quickmin); the default assumes the real units templates
    command(fmt::format("timestep {}", (cfg.timestep > 0.0) ? cfg.timestep : 0.25));
    command("run 0 post no");
    return lmp;
}

static void run_minimize(LAMMPS *lmp)
{
    lmp->input->one(fmt::format("minimize 0.0 0.0 {} {}", MAX_ITER, MAX_EVAL));
    // re-evaluate forces and energy at the final geometry so that the
    // potential energy compute and the force arrays are consistent
    lmp->input->one("run 0 post no");
}

static void run_min_test(LAMMPS::argv &args, double epsilon)
{
    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config);
    } catch (std::exception &e) {
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
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

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms, nlocal);

    ErrorStats stats;

    const double init_pe = potential_energy(lmp);
    EXPECT_FP_LE_WITH_EPS(init_pe, test_config.init_energy, epsilon);

    ::testing::internal::CaptureStdout();
    try {
        run_minimize(lmp);
    } catch (std::exception &e) {
        std::string mnout = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << mnout;
        cleanup_lammps(lmp);
        FAIL() << e.what();
    }
    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    // a minimization from a non-minimum start must lower the energy
    const double min_pe = potential_energy(lmp);
    EXPECT_LT(min_pe, init_pe);

    EXPECT_FP_LE_WITH_EPS(min_pe, test_config.run_energy, epsilon);

    // the exact force norm after a fixed iteration budget is not portable
    // across platforms (see comment at the top), so only require that the
    // minimizer reduces the forces comparably to the reference
    EXPECT_LE(force_norm(lmp), 10.0 * test_config.global_scalar);

    // box dimensions and tilt factors (changed by fix box/relax)
    if (test_config.global_vector.size() == 6) {
        double dims[6];
        box_dims(lmp, dims);
        SCOPED_TRACE("box dimensions after minimization");
        for (int i = 0; i < 6; ++i)
            EXPECT_FP_LE_WITH_EPS(dims[i], test_config.global_vector[i], epsilon);
    }

    if (print_stats) std::cerr << "min_style stats:" << stats << std::endl;

    cleanup_lammps(lmp);
}

TEST(MinStyle, plain)
{
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"MinStyle", "-log", "none", "-echo", "screen", "-nocite"};

    run_min_test(args, test_config.epsilon);
}

// precision of the KOKKOS package as selected with -D KOKKOS_PREC at compile time
static std::string kokkos_precision()
{
    if (Info::has_accelerator_feature("KOKKOS", "precision", "mixed")) return "mixed";
    if (Info::has_accelerator_feature("KOKKOS", "precision", "single")) return "single";
    return "double";
}

// the KOKKOS minimizers descend along a different (but equally valid) path
// than their plain counterparts, so the tolerances have to be relaxed the
// same way as for the other KOKKOS force style tests

static double kokkos_epsilon()
{
    // relax error a bit for KOKKOS package
    double epsilon = 10.0 * test_config.epsilon;
    // relax error a lot for reduced precision KOKKOS builds
    const std::string kk_precision = kokkos_precision();
    if (kk_precision == "mixed")
        epsilon *= 2.0e9;
    else if (kk_precision == "single")
        epsilon *= 1.0e10;
    return epsilon;
}

TEST(MinStyle, kokkos_omp)
{
    if (!Info::has_package("KOKKOS")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    // skip entries may also be qualified by the KOKKOS package precision,
    // e.g. "kokkos_omp_single" skips only single precision KOKKOS builds
    if (test_config.skip_tests.count(std::string(test_info_->name()) + "_" + kokkos_precision()))
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

    LAMMPS::argv args = {"MinStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k",       "on",   "t",    "4",     "-sf",    "kk"};

    run_min_test(args, kokkos_epsilon());
}

TEST(MinStyle, kokkos_omp_full)
{
    if (!Info::has_package("KOKKOS")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    // skip entries may also be qualified by the KOKKOS package precision,
    // e.g. "kokkos_omp_full_single" skips only single precision KOKKOS builds
    if (test_config.skip_tests.count(std::string(test_info_->name()) + "_" + kokkos_precision()))
        GTEST_SKIP();
    // a style that cannot be tested with KOKKOS at all cannot be tested
    // with a full neighbor list either, so the plain "kokkos_omp"
    // skip entries apply here as well
    if (test_config.skip_tests.count("kokkos_omp")) GTEST_SKIP();
    if (test_config.skip_tests.count("kokkos_omp_" + kokkos_precision()))
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
    LAMMPS::argv args = {"MinStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k", "on", "t", "4", "-sf", "kk",
                         "-pk", "kokkos", "neigh", "full", "newton", "off",
                         "-var", "newton_pair", "off", "-var", "newton_bond", "off"};

    run_min_test(args, kokkos_epsilon());
}

TEST(MinStyle, kokkos_serial)
{
    if (!Info::has_package("KOKKOS")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    // skip entries may also be qualified by the KOKKOS package precision,
    // e.g. "kokkos_serial_single" skips only single precision KOKKOS builds
    if (test_config.skip_tests.count(std::string(test_info_->name()) + "_" + kokkos_precision()))
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

    LAMMPS::argv args = {"MinStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k",       "on",   "t",    "1",     "-sf",    "kk"};

    run_min_test(args, kokkos_epsilon());
}

TEST(MinStyle, kokkos_serial_full)
{
    if (!Info::has_package("KOKKOS")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    // skip entries may also be qualified by the KOKKOS package precision,
    // e.g. "kokkos_serial_full_single" skips only single precision KOKKOS builds
    if (test_config.skip_tests.count(std::string(test_info_->name()) + "_" + kokkos_precision()))
        GTEST_SKIP();
    // a style that cannot be tested with KOKKOS at all cannot be tested
    // with a full neighbor list either, so the plain "kokkos_serial"
    // skip entries apply here as well
    if (test_config.skip_tests.count("kokkos_serial")) GTEST_SKIP();
    if (test_config.skip_tests.count("kokkos_serial_" + kokkos_precision()))
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
    LAMMPS::argv args = {"MinStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k", "on", "t", "1", "-sf", "kk",
                         "-pk", "kokkos", "neigh", "full", "newton", "off",
                         "-var", "newton_pair", "off", "-var", "newton_bond", "off"};

    run_min_test(args, kokkos_epsilon());
}

TEST(MinStyle, kokkos_gpu)
{
    if (!Info::has_package("KOKKOS")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    // skip entries may also be qualified by the KOKKOS package precision,
    // e.g. "kokkos_gpu_single" skips only single precision KOKKOS builds
    if (test_config.skip_tests.count(std::string(test_info_->name()) + "_" + kokkos_precision()))
        GTEST_SKIP();
    // this test requires a GPU backend of the KOKKOS package
    if (!Info::has_accelerator_feature("KOKKOS", "api", "cuda") &&
        !Info::has_accelerator_feature("KOKKOS", "api", "hip") &&
        !Info::has_accelerator_feature("KOKKOS", "api", "sycl"))
        GTEST_SKIP() << "KOKKOS GPU backend not enabled";
    // transparently skip when no compatible GPU device is present
    if (!Info::has_kokkos_gpu_device())
        GTEST_SKIP() << "No compatible GPU device available";

    // use a half neighbor list so the GPU kernels run with the input's default
    // "newton on"; with the default "neigh full" the KOKKOS package requires
    // newton off, which the force-style input templates do not use
    LAMMPS::argv args = {"MinStyle", "-log",   "none", "-echo", "screen", "-nocite", "-k", "on",
                         "g",        "1",      "-sf",  "kk",    "-pk",    "kokkos",  "neigh",
                         "half",     "newton", "on"};

    run_min_test(args, kokkos_epsilon());
}

void generate_yaml_file(const char *outfile, const TestConfig &config)
{
    // initialize system geometry
    LAMMPS::argv args = {"MinStyle", "-log", "none", "-echo", "screen", "-nocite"};
    LAMMPS *lmp       = nullptr;
    try {
        lmp = init_lammps(args, config);
    } catch (std::exception &e) {
        FAIL() << e.what();
    }
    if (!lmp) {
        std::cerr << "One or more prerequisite styles are not available "
                     "in this LAMMPS configuration:\n";
        for (auto prerequisite : config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        return;
    }

    const int natoms     = lmp->atom->natoms;
    const double init_pe = potential_energy(lmp);

    run_minimize(lmp);

    std::string block;
    YamlWriter writer(outfile);

    // write yaml header
    write_yaml_header(&writer, &test_config, lmp->version);

    writer.emit("natoms", natoms);
    writer.emit("init_energy", init_pe);
    writer.emit("run_energy", potential_energy(lmp));
    writer.emit("global_scalar", force_norm(lmp));

    // box dimensions and tilt factors
    double dims[6];
    box_dims(lmp, dims);
    block = "6";
    for (double dim : dims)
        block += fmt::format(" {:23.16e}", dim);
    writer.emit_block("global_vector", block);

    cleanup_lammps(lmp);
}
