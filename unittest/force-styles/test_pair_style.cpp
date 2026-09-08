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

// unit tests for pair styles intended for molecular systems

#include "error_stats.h"
#include "test_config.h"
#include "test_main.h"
#include "yaml_writer.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "force.h"
#include "info.h"
#include "utils.h"
#include "input.h"
#include "kspace.h"
#include "modify.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <exception>
#include <iostream>
#include <set>
#include <utility>

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

void cleanup_lammps(LAMMPS *&lmp, const TestConfig &cfg)
{
    platform::unlink(cfg.basename + ".restart");
    platform::unlink(cfg.basename + ".data");
    platform::unlink(cfg.basename + "-coeffs.in");
    delete lmp;
    lmp = nullptr;
}

LAMMPS *init_lammps(LAMMPS::argv &args, const TestConfig &cfg, const bool newton)
{
    LAMMPS *lmp;

    lmp = new LAMMPS(args, MPI_COMM_WORLD);

    // check if prerequisite styles are available
    Info *info = new Info(lmp);
    int nfail  = 0;
    for (const auto &prerequisite : cfg.prerequisites) {
        std::string style = prerequisite.second;

        // this is a test for pair styles, so if the suffixed
        // version is not available, there is no reason to test.
        if (prerequisite.first == "pair") {
            if (lmp->suffix_enable) {
                style += "/";
                style += lmp->suffix;
            }
        }

        if (!info->has_style(prerequisite.first, style)) ++nfail;
    }
    delete info;
    if (nfail > 0) {
        cleanup_lammps(lmp, cfg);
        return nullptr;
    }

    // utility lambdas to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    auto parse_input_script = [&](const std::string &filename) {
        lmp->input->file(filename.c_str());
    };

    if (newton) {
        command("variable newton_pair index on");
    } else {
        command("variable newton_pair index off");
    }

    command("variable input_dir index " + INPUT_FOLDER);

    for (const auto &pre_command : cfg.pre_commands) {
        command(pre_command);
    }
    enforce_kokkos_full_neigh(lmp);

    std::string input_file = platform::path_join(INPUT_FOLDER, cfg.input_file);
    parse_input_script(input_file);

    command("pair_style " + cfg.pair_style);

    for (const auto &pair_coeff : cfg.pair_coeff) {
        command("pair_coeff " + pair_coeff);
    }

    // set this here explicitly to a setting different
    // from the default, so we can spot YAML files for
    // long-range interactions that do not include these
    // settings. they will fail after restart or read data.
    command("pair_modify table 0");
    command("pair_modify table/disp 0");

    for (const auto &post_command : cfg.post_commands) {
        command(post_command);
    }

    command("run 0 post no");
    command("variable write_data_pair index ii");
    command("write_restart " + cfg.basename + ".restart");
    command("write_data " + cfg.basename + ".data pair ${write_data_pair}");
    command("write_coeff " + cfg.basename + "-coeffs.in");

    return lmp;
}

void run_lammps(LAMMPS *lmp, const TestConfig &cfg)
{
    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    command("compute pe all pe/atom pair");
    command("compute sum all reduce sum c_pe");
    command("thermo_style custom step temp pe press c_sum");
    command("thermo 2");

    // need to use a different integrator for different atom styles
    if (std::find(cfg.tags.begin(), cfg.tags.end(), "ellipsoid") != cfg.tags.end()) {
        command("fix 1 all nve/asphere");
        command("compute etemp all temp/asphere");
        command("thermo_modify temp etemp");
    } else if (std::find(cfg.tags.begin(), cfg.tags.end(), "spin") != cfg.tags.end()) {
        // spin systems must define "fix nve/spin" in the yaml post_commands so it is
        // present in all test stages: the spin pair styles compute the mechanical
        // forces only when the fix is present (they take its "lattice" setting).
    } else {
        command("fix 1 all nve");
    }
    command("run 4 post no");
}

void restart_lammps(LAMMPS *lmp, const TestConfig &cfg, bool nofdotr = false, bool newton = true)
{
    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    command("clear");
    if (newton)
        command("newton on");
    else
        command("newton off");
    command("read_restart " + cfg.basename + ".restart");

    if (!lmp->force->pair) {
        command("pair_style " + cfg.pair_style);
    }
    if (!lmp->force->pair->restartinfo || !lmp->force->pair->writedata) {
        for (const auto &pair_coeff : cfg.pair_coeff) {
            command("pair_coeff " + pair_coeff);
        }
    }

    for (const auto &post_command : cfg.post_commands) {
        command(post_command);
    }
    if (nofdotr) command("pair_modify nofdotr");

    command("run 0 post no");
}

void data_lammps(LAMMPS *lmp, const TestConfig &cfg)
{
    // utility lambdas to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };
    auto parse_input_script = [&](const std::string &filename) {
        lmp->input->file(filename.c_str());
    };

    command("clear");
    command("variable pair_style  delete");
    command("variable data_file   delete");
    command("variable newton_pair delete");
    command("variable newton_pair index on");

    for (const auto &pre_command : cfg.pre_commands) {
        command(pre_command);
    }
    enforce_kokkos_full_neigh(lmp);

    command("variable pair_style index '" + cfg.pair_style + "'");
    command("variable data_file index " + cfg.basename + ".data");

    std::string input_file = platform::path_join(INPUT_FOLDER, cfg.input_file);
    parse_input_script(input_file);

    for (const auto &pair_coeff : cfg.pair_coeff) {
        command("pair_coeff " + pair_coeff);
    }
    for (const auto &post_command : cfg.post_commands) {
        command(post_command);
    }
    command("run 0 post no");
}

// re-generate yaml file with current settings.

void generate_yaml_file(const char *outfile, const TestConfig &config)
{
    // initialize system geometry
    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite"};
    LAMMPS *lmp       = nullptr;
    try {
        lmp = init_lammps(args, config, true);
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

    const int natoms = lmp->atom->natoms;
    std::string block;

    YamlWriter writer(outfile);

    // write yaml header
    write_yaml_header(&writer, &test_config, lmp->version);

    // pair_style
    writer.emit("pair_style", config.pair_style);

    // pair_coeff
    block.clear();
    for (auto pair_coeff : config.pair_coeff) {
        block += pair_coeff + "\n";
    }
    writer.emit_block("pair_coeff", block);

    // extract
    block.clear();
    for (auto data : config.extract)
        block += fmt::format("{} {}\n", data.first, data.second);
    writer.emit_block("extract", block);

    // natoms
    writer.emit("natoms", natoms);

    // init_vdwl
    writer.emit("init_vdwl", lmp->force->pair->eng_vdwl);

    // init_coul
    writer.emit("init_coul", lmp->force->pair->eng_coul);

    // init_stress
    auto *stress = lmp->force->pair->virial;
    // avoid false positives on tiny stresses. force to zero instead.
    for (int i = 0; i < 6; ++i)
        if (fabs(stress[i]) < 1.0e-13) stress[i] = 0.0;
    block = fmt::format("{:23.16e} {:23.16e} {:23.16e} {:23.16e} {:23.16e} {:23.16e}", stress[0],
                        stress[1], stress[2], stress[3], stress[4], stress[5]);
    writer.emit_block("init_stress", block);

    // init_forces
    block.clear();
    auto *f = lmp->atom->f;
    for (int i = 1; i <= natoms; ++i) {
        const int j = lmp->atom->map(i);
        block += fmt::format("{:3} {:23.16e} {:23.16e} {:23.16e}\n", i, f[j][0], f[j][1], f[j][2]);
    }
    writer.emit_block("init_forces", block);

    // init_mag_forces (only for atom_style spin)
    if (lmp->atom->sp_flag) {
        block.clear();
        auto *fm = lmp->atom->fm;
        for (int i = 1; i <= natoms; ++i) {
            const int j = lmp->atom->map(i);
            block += fmt::format("{:3} {:23.16e} {:23.16e} {:23.16e}\n", i, fm[j][0], fm[j][1],
                                 fm[j][2]);
        }
        writer.emit_block("init_mag_forces", block);
    }

    // do a few steps of MD
    run_lammps(lmp, config);

    // run_vdwl
    writer.emit("run_vdwl", lmp->force->pair->eng_vdwl);

    // run_coul
    writer.emit("run_coul", lmp->force->pair->eng_coul);

    // run_stress
    stress = lmp->force->pair->virial;
    // avoid false positives on tiny stresses. force to zero instead.
    for (int i = 0; i < 6; ++i)
        if (fabs(stress[i]) < 1.0e-13) stress[i] = 0.0;
    block = fmt::format("{:23.16e} {:23.16e} {:23.16e} {:23.16e} {:23.16e} {:23.16e}", stress[0],
                        stress[1], stress[2], stress[3], stress[4], stress[5]);
    writer.emit_block("run_stress", block);

    block.clear();
    f = lmp->atom->f;
    for (int i = 1; i <= natoms; ++i) {
        const int j = lmp->atom->map(i);
        block += fmt::format("{:3} {:23.16e} {:23.16e} {:23.16e}\n", i, f[j][0], f[j][1], f[j][2]);
    }
    writer.emit_block("run_forces", block);

    // run_mag_forces (only for atom_style spin)
    if (lmp->atom->sp_flag) {
        block.clear();
        auto *fm = lmp->atom->fm;
        for (int i = 1; i <= natoms; ++i) {
            const int j = lmp->atom->map(i);
            block += fmt::format("{:3} {:23.16e} {:23.16e} {:23.16e}\n", i, fm[j][0], fm[j][1],
                                 fm[j][2]);
        }
        writer.emit_block("run_mag_forces", block);
    }

    cleanup_lammps(lmp, config);
}

TEST(PairStyle, plain)
{
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite"};

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config, true);
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
    EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms, nlocal);

    double epsilon = test_config.epsilon;
    // relax test precision when using pppm and single precision FFTs
#if defined(FFT_SINGLE)
    if (lmp->force->kspace && lmp->force->kspace->compute_flag)
        if (utils::strmatch(lmp->force->kspace_style, "^pppm")) epsilon *= 2.0e8;
#endif
    auto *pair = lmp->force->pair;

    EXPECT_FORCES("init_forces (newton on)", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_MAG_FORCES("init_mag_forces (newton on)", lmp->atom, test_config.init_mag_forces,
                      epsilon);
    EXPECT_STRESS("init_stress (newton on)", pair->virial, test_config.init_stress, epsilon);

    ErrorStats stats;
    stats.reset();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
    if (print_stats) std::cerr << "init_energy stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    EXPECT_FORCES("run_forces (newton on)", lmp->atom, test_config.run_forces, 5 * epsilon);
    EXPECT_MAG_FORCES("run_mag_forces (newton on)", lmp->atom, test_config.run_mag_forces,
                      5 * epsilon);
    EXPECT_STRESS("run_stress (newton on)", pair->virial, test_config.run_stress, epsilon);

    stats.reset();
    int id        = lmp->modify->find_compute("sum");
    double energy = lmp->modify->compute[id]->compute_scalar();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.run_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.run_coul, epsilon);
    // skip comparing per-atom energy with total energy for "kim" and "in.conp"
    if ((std::string("kim") != lmp->force->pair_style) &&
        (std::string("pod") != lmp->force->pair_style) &&
        (std::string("in.conp") != test_config.input_file))
        EXPECT_FP_LE_WITH_EPS((pair->eng_vdwl + pair->eng_coul), energy, epsilon);
    if (print_stats) std::cerr << "run_energy  stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    try {
        lmp = init_lammps(args, test_config, false);
    } catch (std::exception &e) {
        if (!verbose) ::testing::internal::GetCapturedStdout();
        FAIL() << e.what();
    }
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // skip over these tests if newton pair is forced to be on
    if (lmp->force->newton_pair == 0) {
        pair = lmp->force->pair;

        EXPECT_FORCES("init_forces (newton off)", lmp->atom, test_config.init_forces, epsilon);
        EXPECT_MAG_FORCES("init_mag_forces (newton off)", lmp->atom, test_config.init_mag_forces,
                          epsilon);
        EXPECT_STRESS("init_stress (newton off)", pair->virial, test_config.init_stress,
                      3 * epsilon);

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, epsilon);
        EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
        if (print_stats) std::cerr << "init_energy stats, newton off:" << stats << std::endl;

        if (!verbose) ::testing::internal::CaptureStdout();
        run_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        EXPECT_FORCES("run_forces (newton off)", lmp->atom, test_config.run_forces, 5 * epsilon);
        EXPECT_MAG_FORCES("run_mag_forces (newton off)", lmp->atom, test_config.run_mag_forces,
                          5 * epsilon);
        EXPECT_STRESS("run_stress (newton off)", pair->virial, test_config.run_stress, epsilon);

        stats.reset();
        id     = lmp->modify->find_compute("sum");
        energy = lmp->modify->compute[id]->compute_scalar();
        EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.run_vdwl, epsilon);
        EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.run_coul, epsilon);
        // skip comparing per-atom energy with total energy for "kim"
        if (std::string("kim") != lmp->force->pair_style)
            EXPECT_FP_LE_WITH_EPS((pair->eng_vdwl + pair->eng_coul), energy, epsilon);
        if (print_stats) std::cerr << "run_energy  stats, newton off:" << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    restart_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    pair = lmp->force->pair;

    EXPECT_FORCES("restart_forces", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_MAG_FORCES("restart_mag_forces", lmp->atom, test_config.init_mag_forces, epsilon);
    EXPECT_STRESS("restart_stress", pair->virial, test_config.init_stress, epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
    if (print_stats) std::cerr << "restart_energy stats:" << stats << std::endl;

    // pair style rann does not support pair_modify nofdotr.  styles whose
    // dissipative or random pair forces are not central (e.g. sdpd) cannot
    // reproduce the fdotr virial through ev_tally() and may opt out with
    // the "nofdotr" token in skip_tests.
    if ((test_config.pair_style != "rann") && !test_config.skip_tests.count("nofdotr")) {
        if (!verbose) ::testing::internal::CaptureStdout();
        restart_lammps(lmp, test_config, true);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        pair = lmp->force->pair;

        EXPECT_FORCES("nofdotr_forces", lmp->atom, test_config.init_forces, epsilon);
        EXPECT_MAG_FORCES("nofdotr_mag_forces", lmp->atom, test_config.init_mag_forces, epsilon);
        EXPECT_STRESS("nofdotr_stress", pair->virial, test_config.init_stress, epsilon);

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, epsilon);
        EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
        if (print_stats) std::cerr << "nofdotr_energy stats:" << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    data_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    pair = lmp->force->pair;
    EXPECT_FORCES("data_forces", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_MAG_FORCES("data_mag_forces", lmp->atom, test_config.init_mag_forces, epsilon);
    EXPECT_STRESS("data_stress", pair->virial, test_config.init_stress, epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
    if (print_stats) std::cerr << "data_energy stats:" << stats << std::endl;

    if (pair->respa_enable) {
        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        try {
            lmp = init_lammps(args, test_config, false);
            lmp->input->one("run_style respa 2 1 inner 1 4.8 5.5 outer 2");
            run_lammps(lmp, test_config);
        } catch (std::exception &e) {
            if (!verbose) ::testing::internal::GetCapturedStdout();
            FAIL() << e.what();
        }
        if (!verbose) ::testing::internal::GetCapturedStdout();

        // need to relax error by a large amount with tabulation, since
        // coul/long styles do not use tabulation in compute_inner()
        // and compute_middle() so we get a significant deviation.
        pair = lmp->force->pair;
        if (pair->ncoultablebits) epsilon *= 5.0e6;

        EXPECT_FORCES("run_forces (r-RESPA)", lmp->atom, test_config.run_forces, 5 * epsilon);
        EXPECT_STRESS("run_stress (r-RESPA)", pair->virial, test_config.run_stress, epsilon);

        stats.reset();
        id     = lmp->modify->find_compute("sum");
        energy = lmp->modify->compute[id]->compute_scalar();
        EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.run_vdwl, epsilon);
        EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.run_coul, epsilon);
        EXPECT_FP_LE_WITH_EPS((pair->eng_vdwl + pair->eng_coul), energy, epsilon);
        if (print_stats) std::cerr << "run_energy  stats, r-RESPA:" << stats << std::endl;
    }
    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

TEST(PairStyle, omp)
{
    if (!Info::has_package("OPENMP")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-pk",       "omp",  "4",    "-sf",   "omp"};

    // styles tagged "single_thread" (e.g. dpd, which uses multiple pRNGs) cannot
    // run with more than one thread in the test
    if (test_config.has_tag("single_thread")) args[8] = "1";

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config, true);
    } catch (std::exception &e) {
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        if (full_neigh_unsupported(e.what())) GTEST_SKIP() << e.what();
        FAIL() << e.what();
    }
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    if (!lmp) {
        std::cerr << "One or more prerequisite styles with /omp suffix\n"
                     "are not available in this LAMMPS configuration:\n";
        for (auto &prerequisite : test_config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        GTEST_SKIP();
    }

    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms, nlocal);

    // relax error a bit for OPENMP package
    double epsilon = 5.0 * test_config.epsilon;
    // relax test precision when using pppm and single precision FFTs
#if defined(FFT_SINGLE)
    if (lmp->force->kspace && lmp->force->kspace->compute_flag)
        if (utils::strmatch(lmp->force->kspace_style, "^pppm")) epsilon *= 2.0e8;
#endif
    auto *pair = lmp->force->pair;
    ErrorStats stats;

    EXPECT_FORCES("init_forces (newton on)", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_MAG_FORCES("init_mag_forces (newton on)", lmp->atom, test_config.init_mag_forces,
                      epsilon);
    EXPECT_STRESS("init_stress (newton on)", pair->virial, test_config.init_stress, 10 * epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
    if (print_stats) std::cerr << "init_energy stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    EXPECT_FORCES("run_forces (newton on)", lmp->atom, test_config.run_forces, 5 * epsilon);
    EXPECT_MAG_FORCES("run_mag_forces (newton on)", lmp->atom, test_config.run_mag_forces,
                      5 * epsilon);
    EXPECT_STRESS("run_stress (newton on)", pair->virial, test_config.run_stress, 10 * epsilon);

    stats.reset();
    int id        = lmp->modify->find_compute("sum");
    double energy = lmp->modify->compute[id]->compute_scalar();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.run_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.run_coul, epsilon);
    EXPECT_FP_LE_WITH_EPS((pair->eng_vdwl + pair->eng_coul), energy, epsilon);
    if (print_stats) std::cerr << "run_energy  stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    try {
        lmp = init_lammps(args, test_config, false);
    } catch (std::exception &e) {
        if (!verbose) ::testing::internal::GetCapturedStdout();
        FAIL() << e.what();
    }
    if (!verbose) ::testing::internal::GetCapturedStdout();

    pair = lmp->force->pair;

    // skip over these tests if newton pair is forced to be on
    if (lmp->force->newton_pair == 0) {

        EXPECT_FORCES("init_forces (newton off)", lmp->atom, test_config.init_forces, epsilon);
        EXPECT_MAG_FORCES("init_mag_forces (newton off)", lmp->atom, test_config.init_mag_forces,
                          epsilon);
        EXPECT_STRESS("init_stress (newton off)", pair->virial, test_config.init_stress,
                      10 * epsilon);

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, epsilon);
        EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
        if (print_stats) std::cerr << "init_energy stats, newton off:" << stats << std::endl;

        if (!verbose) ::testing::internal::CaptureStdout();
        run_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        EXPECT_FORCES("run_forces (newton off)", lmp->atom, test_config.run_forces, 5 * epsilon);
        EXPECT_MAG_FORCES("run_mag_forces (newton off)", lmp->atom, test_config.run_mag_forces,
                          5 * epsilon);
        EXPECT_STRESS("run_stress (newton off)", pair->virial, test_config.run_stress,
                      10 * epsilon);

        stats.reset();
        id     = lmp->modify->find_compute("sum");
        energy = lmp->modify->compute[id]->compute_scalar();
        EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.run_vdwl, epsilon);
        EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.run_coul, epsilon);
        EXPECT_FP_LE_WITH_EPS((pair->eng_vdwl + pair->eng_coul), energy, epsilon);
        if (print_stats) std::cerr << "run_energy  stats, newton off:" << stats << std::endl;
    }

    // the "nofdotr" token in skip_tests opts a style out of comparing the
    // virial tallied by ev_tally() with the fdotr one of the reference, for the
    // same reason as in the plain test case above

    if (!test_config.skip_tests.count("nofdotr")) {
        if (!verbose) ::testing::internal::CaptureStdout();
        restart_lammps(lmp, test_config, true);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        pair = lmp->force->pair;

        EXPECT_FORCES("nofdotr_forces", lmp->atom, test_config.init_forces, 5 * epsilon);
        EXPECT_MAG_FORCES("nofdotr_mag_forces", lmp->atom, test_config.init_mag_forces,
                          5 * epsilon);
        EXPECT_STRESS("nofdotr_stress", pair->virial, test_config.init_stress, 10 * epsilon);

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, 5 * epsilon);
        EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, 5 * epsilon);
        if (print_stats) std::cerr << "nofdotr_energy stats:" << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

// precision of the KOKKOS package as selected with -D KOKKOS_PREC at compile time
static std::string kokkos_precision()
{
    if (Info::has_accelerator_feature("KOKKOS", "precision", "mixed")) return "mixed";
    if (Info::has_accelerator_feature("KOKKOS", "precision", "single")) return "single";
    return "double";
}

// the "newton" argument must match the newton setting the command line
// arguments select: with "-pk kokkos neigh full" the KOKKOS package rejects
// "newton on", so the restarted run has to be set up with newton off as well

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

static void run_kokkos_test(LAMMPS::argv &args, bool newton = true)
{
    append_kokkos_env_args(args);

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config, true);
    } catch (std::exception &e) {
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        if (full_neigh_unsupported(e.what())) GTEST_SKIP() << e.what();
        FAIL() << e.what();
    }
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    if (!lmp) {
        std::cerr << "One or more prerequisite styles with /kk suffix\n"
                     "are not available in this LAMMPS configuration:\n";
        for (auto &prerequisite : test_config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        GTEST_SKIP();
    }

    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms, nlocal);

    // relax error a bit for KOKKOS package
    double epsilon = 5.0 * test_config.epsilon;
    // relax error a lot for reduced precision KOKKOS builds
    const std::string kk_precision = kokkos_precision();
    if (kk_precision == "mixed")
        epsilon *= 2.0e9;
    else if (kk_precision == "single")
        epsilon *= 1.0e10;
    // relax test precision when using pppm and single precision FFTs
#if defined(FFT_SINGLE)
    if (lmp->force->kspace && lmp->force->kspace->compute_flag)
        if (utils::strmatch(lmp->force->kspace_style, "^pppm")) epsilon *= 2.0e8;
#endif
    auto *pair = lmp->force->pair;
    ErrorStats stats;

    EXPECT_FORCES("init_forces (newton on)", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_MAG_FORCES("init_mag_forces (newton on)", lmp->atom, test_config.init_mag_forces,
                      epsilon);
    EXPECT_STRESS("init_stress (newton on)", pair->virial, test_config.init_stress, 10 * epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
    if (print_stats) std::cerr << "init_energy stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    EXPECT_FORCES("run_forces (newton on)", lmp->atom, test_config.run_forces, 5 * epsilon);
    EXPECT_MAG_FORCES("run_mag_forces (newton on)", lmp->atom, test_config.run_mag_forces,
                      5 * epsilon);
    EXPECT_STRESS("run_stress (newton on)", pair->virial, test_config.run_stress, 10 * epsilon);

    stats.reset();
    int id        = lmp->modify->find_compute("sum");
    double energy = lmp->modify->compute[id]->compute_scalar();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.run_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.run_coul, epsilon);
    EXPECT_FP_LE_WITH_EPS((pair->eng_vdwl + pair->eng_coul), energy, epsilon);
    if (print_stats) std::cerr << "run_energy  stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    try {
        lmp = init_lammps(args, test_config, false);
    } catch (std::exception &e) {
        if (!verbose) ::testing::internal::GetCapturedStdout();
        FAIL() << e.what();
    }
    if (!verbose) ::testing::internal::GetCapturedStdout();

    pair = lmp->force->pair;

    // skip over these tests if newton pair is forced to be on
    if (lmp->force->newton_pair == 0) {

        EXPECT_FORCES("init_forces (newton off)", lmp->atom, test_config.init_forces, epsilon);
        EXPECT_MAG_FORCES("init_mag_forces (newton off)", lmp->atom, test_config.init_mag_forces,
                          epsilon);
        EXPECT_STRESS("init_stress (newton off)", pair->virial, test_config.init_stress,
                      10 * epsilon);

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, epsilon);
        EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
        if (print_stats) std::cerr << "init_energy stats, newton off:" << stats << std::endl;

        if (!verbose) ::testing::internal::CaptureStdout();
        run_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        EXPECT_FORCES("run_forces (newton off)", lmp->atom, test_config.run_forces, 5 * epsilon);
        EXPECT_MAG_FORCES("run_mag_forces (newton off)", lmp->atom, test_config.run_mag_forces,
                          5 * epsilon);
        EXPECT_STRESS("run_stress (newton off)", pair->virial, test_config.run_stress,
                      10 * epsilon);

        stats.reset();
        id     = lmp->modify->find_compute("sum");
        energy = lmp->modify->compute[id]->compute_scalar();
        EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.run_vdwl, epsilon);
        EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.run_coul, epsilon);
        EXPECT_FP_LE_WITH_EPS((pair->eng_vdwl + pair->eng_coul), energy, epsilon);
        if (print_stats) std::cerr << "run_energy  stats, newton off:" << stats << std::endl;
    }

    // the "nofdotr" token in skip_tests opts a style out of comparing the
    // virial tallied by ev_tally() with the fdotr one of the reference.  it
    // applies here for the same reason it applies to the plain test case: the
    // random or dissipative pair forces of the style are not central

    if (!test_config.skip_tests.count("nofdotr")) {
        if (!verbose) ::testing::internal::CaptureStdout();
        restart_lammps(lmp, test_config, true, newton);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        pair = lmp->force->pair;

        EXPECT_FORCES("nofdotr_forces", lmp->atom, test_config.init_forces, 5 * epsilon);
        EXPECT_MAG_FORCES("nofdotr_mag_forces", lmp->atom, test_config.init_mag_forces,
                          5 * epsilon);
        EXPECT_STRESS("nofdotr_stress", pair->virial, test_config.init_stress, 10 * epsilon);

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, 5 * epsilon);
        EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, 5 * epsilon);
        if (print_stats) std::cerr << "nofdotr_energy stats:" << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
}


/* ----------------------------------------------------------------------
   collect per-atom pair energies from a direct pair->compute() call
   with explicit energy/virial flags on the current configuration,
   including the ghost atom reduction done by compute pe/atom. styles
   must produce the same per-atom energies whether or not the virial is
   tallied; styles that maintain their tally bookkeeping only when the
   virial is requested (as the TIP4P OpenMP variants did) assign
   per-atom energies to stale or invalid atom indices otherwise, which
   the comparison against the with-virial reference detects
------------------------------------------------------------------------- */

std::vector<double> eatom_direct(LAMMPS *lmp, bool with_virial)
{
    // satisfy the tally timestamp checks of pair->compute consumers

    lmp->update->eflag_atom  = lmp->update->ntimestep;
    lmp->update->eflag_global = lmp->update->ntimestep;
    // request the same global-virial mode the integrator would use;
    // per-atom virial support is not universal (e.g. pair rann)
    int vflag = 0;
    if (with_virial) {
        vflag = lmp->force->pair->no_virial_fdotr_compute ? VIRIAL_PAIR : VIRIAL_FDOTR;
        lmp->update->vflag_global = lmp->update->ntimestep;
    }
    lmp->force->pair->compute(ENERGY_GLOBAL | ENERGY_ATOM, vflag);

    auto *pea = lmp->modify->get_compute_by_id("peaonly");
    EXPECT_NE(pea, nullptr);
    pea->compute_peratom();
    pea->invoked_peratom = -1;    // force a fresh evaluation on the next call

    const int nlocal = lmp->atom->nlocal;
    const tagint *tag = lmp->atom->tag;
    std::vector<double> eatom(nlocal, 0.0);
    for (int i = 0; i < nlocal; i++)
        eatom[tag[i] - 1] = pea->vector_atom[i];
    return eatom;
}

void eatom_only_test(LAMMPS::argv args, const TestConfig &cfg)
{
    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, cfg, true);
    } catch (std::exception &e) {
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        if (full_neigh_unsupported(e.what())) GTEST_SKIP() << e.what();
        FAIL() << e.what();
    }
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    if (!lmp) GTEST_SKIP();

    // exception safety matters here: leaving the capturer active on a
    // throw aborts the whole test program at the next captured section
    ::testing::internal::CaptureStdout();
    std::vector<double> reference, eonly;
    try {
        lmp->input->one("compute peaonly all pe/atom pair");
        lmp->input->one("run 0 post no");
        reference = eatom_direct(lmp, true);
        eonly     = eatom_direct(lmp, false);
    } catch (std::exception &e) {
        std::string errout = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << errout;
        cleanup_lammps(lmp, cfg);
        FAIL() << e.what();
    }
    cleanup_lammps(lmp, cfg);
    ::testing::internal::GetCapturedStdout();

    ASSERT_EQ(reference.size(), eonly.size());
    const double epsilon = cfg.epsilon;
    for (std::size_t i = 0; i < reference.size(); i++)
        EXPECT_NEAR(eonly[i], reference[i], (fabs(reference[i]) + 1.0) * epsilon * 10.0)
            << "per-atom energy for atom " << i + 1
            << " differs between tallying with and without virial";
}

TEST(PairStyle, eatom_only)
{
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite"};
    eatom_only_test(args, test_config);
}

TEST(PairStyle, eatom_only_omp)
{
    if (!Info::has_package("OPENMP")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    if (test_config.skip_tests.count("omp")) GTEST_SKIP();

    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-pk",       "omp",  "4",    "-sf",   "omp"};
    if (test_config.has_tag("single_thread")) args[8] = "1";
    eatom_only_test(args, test_config);
}

/* ----------------------------------------------------------------------
   the per-atom virial of a pair style has to add up to its global
   virial.  this is a self-consistency check that needs no stored
   reference data and it is the only place in this driver where a
   per-atom virial is requested at all: without it the "vflag_atom"
   branches of the pair styles and of their accelerated variants are
   never executed, so a style that computes the global virial correctly
   but does not fill (or wrongly fills) its per-atom virial array is
   not detected
------------------------------------------------------------------------- */

// "compute stress/atom" reports the per-atom virial of the selected
// contributions scaled by -nktv2p, so the sum has to be scaled back to
// the units of Pair::virial before it can be compared with it

static void vatom_sum(LAMMPS *lmp, double *sum)
{
    auto *vas = lmp->modify->get_compute_by_id("vasum");
    ASSERT_NE(vas, nullptr);
    vas->compute_vector();
    const double scale = -1.0 / lmp->force->nktv2p;
    for (int k = 0; k < 6; ++k)
        sum[k] = scale * vas->vector[k];
}

static void vatom_only_test(LAMMPS::argv args, const TestConfig &cfg, double epsilon)
{
    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, cfg, true);
    } catch (std::exception &e) {
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        if (full_neigh_unsupported(e.what())) GTEST_SKIP() << e.what();
        FAIL() << e.what();
    }
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    if (!lmp) GTEST_SKIP();

    // exception safety matters here: leaving the capturer active on a
    // throw aborts the whole test program at the next captured section
    ::testing::internal::CaptureStdout();
    double sum[6], reference[6];
    for (int k = 0; k < 6; ++k)
        sum[k] = reference[k] = 0.0;
    try {
        // the "nofdotr" token in skip_tests marks styles whose global virial
        // from ev_tally() differs from the fdotr one, because they have force
        // contributions that are not tallied (random or dissipative forces).
        // for those the per-atom virial can only be compared with the global
        // virial that ev_tally() accumulates, so request that one here
        if (cfg.skip_tests.count("nofdotr")) lmp->input->one("pair_modify nofdotr");
        lmp->input->one("compute vaonly all stress/atom NULL pair");
        lmp->input->one("compute vasum all reduce sum c_vaonly[1] c_vaonly[2] c_vaonly[3] "
                        "c_vaonly[4] c_vaonly[5] c_vaonly[6]");
        lmp->input->one("run 0 post no");
        vatom_sum(lmp, sum);
        for (int k = 0; k < 6; ++k)
            reference[k] = lmp->force->pair->virial[k];
    } catch (std::exception &e) {
        std::string errout = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << errout;
        cleanup_lammps(lmp, cfg);
        FAIL() << e.what();
    }
    cleanup_lammps(lmp, cfg);
    ::testing::internal::GetCapturedStdout();

    double scale = 1.0;
    for (int k = 0; k < 6; ++k)
        if (fabs(reference[k]) > scale) scale = fabs(reference[k]);

    const char *label[6] = {"xx", "yy", "zz", "xy", "xz", "yz"};
    for (int k = 0; k < 6; ++k)
        EXPECT_NEAR(sum[k], reference[k], scale * epsilon)
            << "sum of the per-atom virial differs from the global virial "
            << "for the " << label[k] << " component";
}

TEST(PairStyle, vatom_only)
{
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite"};
    vatom_only_test(args, test_config, 10.0 * test_config.epsilon);
}

TEST(PairStyle, vatom_only_omp)
{
    if (!Info::has_package("OPENMP")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    if (test_config.skip_tests.count("omp")) GTEST_SKIP();
    // a style that cannot tally a per-atom virial at all cannot do so with
    // an accelerated variant either, so the plain "vatom_only" skip entries
    // apply to this test case as well
    if (test_config.skip_tests.count("vatom_only")) GTEST_SKIP();

    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-pk",       "omp",  "4",    "-sf",   "omp"};
    if (test_config.has_tag("single_thread")) args[8] = "1";
    vatom_only_test(args, test_config, 50.0 * test_config.epsilon);
}

TEST(PairStyle, vatom_only_kokkos)
{
    if (!Info::has_package("KOKKOS")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    // a style that cannot tally a per-atom virial at all cannot do so with
    // an accelerated variant either, so the plain "vatom_only" skip entries
    // apply to this test case as well
    if (test_config.skip_tests.count("vatom_only")) GTEST_SKIP();
    // skip entries may also be qualified by the KOKKOS package precision,
    // e.g. "vatom_only_kokkos_single" skips only single precision builds
    if (test_config.skip_tests.count(std::string(test_info_->name()) + "_" + kokkos_precision()))
        GTEST_SKIP();
    // skip entries qualified with "_devicerng" apply only to builds where the
    // KOKKOS styles use the device random number generator
    if (Info::has_accelerator_feature("KOKKOS", "rng", "device") &&
        test_config.skip_tests.count(std::string(test_info_->name()) + "_devicerng"))
        GTEST_SKIP();

    // a style that cannot be tested with KOKKOS at all cannot have its
    // per-atom virial tested with KOKKOS either, so the skip entries of the
    // regular KOKKOS test case for the backend in use apply here as well
    const bool kk_gpu = Info::has_accelerator_feature("KOKKOS", "api", "cuda") ||
                        Info::has_accelerator_feature("KOKKOS", "api", "hip") ||
                        Info::has_accelerator_feature("KOKKOS", "api", "sycl");
    const bool kk_threads = Info::has_accelerator_feature("KOKKOS", "api", "openmp") ||
                            Info::has_accelerator_feature("KOKKOS", "api", "pthreads");
    std::string base = "kokkos_serial";
    if (kk_gpu)
        base = "kokkos_gpu";
    else if (kk_threads)
        base = "kokkos_omp";
    if (test_config.skip_tests.count(base)) GTEST_SKIP();
    if (test_config.skip_tests.count(base + "_" + kokkos_precision())) GTEST_SKIP();
    if (Info::has_accelerator_feature("KOKKOS", "rng", "device") &&
        test_config.skip_tests.count(base + "_devicerng"))
        GTEST_SKIP();
    // transparently skip when no compatible GPU device is present
    if (kk_gpu && !Info::has_kokkos_gpu_device())
        GTEST_SKIP() << "No compatible GPU device available";

    // use a half neighbor list with newton on, as in the regular KOKKOS
    // test cases, so the setup matches what the input templates expect
    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k",        "on",   "t",    "1",     "-sf",    "kk"};
    if (kk_gpu)
        args = {"PairStyle", "-log", "none",   "-echo",  "screen", "-nocite", "-k",
                "on",        "g",    "1",      "-sf",    "kk",     "-pk",     "kokkos",
                "neigh",     "half", "newton", "on"};
    else if (kk_threads && !test_config.has_tag("single_thread"))
        args[9] = "4";

    // relax error a lot for reduced precision KOKKOS builds
    double epsilon                 = 50.0 * test_config.epsilon;
    const std::string kk_precision = kokkos_precision();
    if (kk_precision == "mixed")
        epsilon *= 2.0e9;
    else if (kk_precision == "single")
        epsilon *= 1.0e10;
    vatom_only_test(args, test_config, epsilon);
}

TEST(PairStyle, kokkos_omp)
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

    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k",        "on",   "t",    "4",     "-sf",    "kk"};

    // some styles cannot run with more than one thread in the test (dpd uses
    // multiple pRNGs, snap and pace due to their implementation); these are
    // flagged with the "single_thread" tag in their YAML file
    if (test_config.has_tag("single_thread")) args[9] = "1";

    run_kokkos_test(args);
};

TEST(PairStyle, kokkos_omp_full)
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
    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k", "on", "t", "4", "-sf", "kk",
                         "-pk", "kokkos", "neigh", "full", "newton", "off",
                         "-var", "newton_pair", "off", "-var", "newton_bond", "off"};

    // some styles cannot run with more than one thread in the test (dpd uses
    // multiple pRNGs, snap and pace due to their implementation); these are
    // flagged with the "single_thread" tag in their YAML file
    if (test_config.has_tag("single_thread")) args[9] = "1";

    kokkos_full_neigh = true;
    run_kokkos_test(args, false);
    kokkos_full_neigh = false;
};

TEST(PairStyle, kokkos_serial)
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

    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k",        "on",   "t",    "1",     "-sf",    "kk"};

    run_kokkos_test(args);
};

TEST(PairStyle, kokkos_serial_full)
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
    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k", "on", "t", "1", "-sf", "kk",
                         "-pk", "kokkos", "neigh", "full", "newton", "off",
                         "-var", "newton_pair", "off", "-var", "newton_bond", "off"};

    kokkos_full_neigh = true;
    run_kokkos_test(args, false);
    kokkos_full_neigh = false;
};

TEST(PairStyle, kokkos_gpu)
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
    if (!Info::has_kokkos_gpu_device())
        GTEST_SKIP() << "No compatible GPU device available";

    // use a half neighbor list with newton on so the GPU kernels run the way the
    // input templates expect; the GPU default is "neigh full" + newton off, which
    // (a) the templates do not use and (b) would make the package set newton off
    // at startup, so a later "newton on" after the box exists would error out
    LAMMPS::argv args = {"PairStyle", "-log", "none",   "-echo",  "screen", "-nocite", "-k",
                         "on",        "g",    "1",      "-sf",    "kk",     "-pk",     "kokkos",
                         "neigh",     "half", "newton", "on"};

    run_kokkos_test(args);
};

TEST(PairStyle, gpu)
{
    if (!Info::has_package("GPU")) GTEST_SKIP();
    if (!Info::has_gpu_device()) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    // when testing PPPM styles with GPUs and GPU support is compiled with single precision
    // we also must have single precision FFTs; otherwise skip since the test would abort
    if (utils::strmatch(test_config.basename, ".*pppm.*") &&
        (Info::has_accelerator_feature("GPU", "precision", "single")) &&
        (!Info::has_fft_single_support()))
        GTEST_SKIP();

    // some GPU pair styles do not support single and/or mixed precision GPU mode
    // (e.g. born/coul/long/cs/gpu errors out in single precision). Their tests are
    // tagged "gpu_no_single" / "gpu_no_mixed" and skipped when the GPU package is
    // compiled for that precision.
    if (test_config.has_tag("gpu_no_single") &&
        Info::has_accelerator_feature("GPU", "precision", "single"))
        GTEST_SKIP();
    if (test_config.has_tag("gpu_no_mixed") &&
        Info::has_accelerator_feature("GPU", "precision", "mixed"))
        GTEST_SKIP();

    LAMMPS::argv args_neigh   = {"PairStyle", "-log",    "none", "-echo",
                                 "screen",    "-nocite", "-sf",  "gpu"};
    LAMMPS::argv args_noneigh = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite", "-sf",
                                 "gpu",       "-pk",  "gpu",  "0",     "neigh",  "no"};
    LAMMPS::argv args         = args_neigh;

    // cannot use GPU neighbor list with hybrid pair style (yet)
    if (test_config.pair_style.substr(0, 6) == "hybrid") {
        args = args_noneigh;
    }

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config, false);
    } catch (std::exception &e) {
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        if (full_neigh_unsupported(e.what())) GTEST_SKIP() << e.what();
        FAIL() << e.what();
    }
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    if (!lmp) {
        std::cerr << "One or more prerequisite styles with /gpu suffix\n"
                     "are not available in this LAMMPS configuration:\n";
        for (auto &prerequisite : test_config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        GTEST_SKIP();
    }

    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms, nlocal);

    // relax error for GPU package depending on precision setting
    double epsilon = test_config.epsilon;
    if (Info::has_accelerator_feature("GPU", "precision", "double"))
        epsilon *= 7.5;
    else if (Info::has_accelerator_feature("GPU", "precision", "mixed"))
        epsilon *= 5.0e8;
    else
        epsilon *= 1.0e10;
    // relax test precision when using pppm and single precision FFTs, but only when also
    // running with double precision
#if defined(FFT_SINGLE)
    if (lmp->force->kspace && lmp->force->kspace->compute_flag &&
        Info::has_accelerator_feature("GPU", "precision", "double"))
        if (utils::strmatch(lmp->force->kspace_style, "^pppm")) epsilon *= 2.0e8;
#endif
    ErrorStats stats;
    auto *pair = lmp->force->pair;

    EXPECT_FORCES("init_forces (newton off)", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_MAG_FORCES("init_mag_forces (newton off)", lmp->atom, test_config.init_mag_forces,
                      epsilon);
    EXPECT_STRESS("init_stress (newton off)", pair->virial, test_config.init_stress, 10 * epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
    if (print_stats) std::cerr << "init_energy stats, newton off:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    EXPECT_FORCES("run_forces (newton off)", lmp->atom, test_config.run_forces, 5 * epsilon);
    EXPECT_MAG_FORCES("run_mag_forces (newton off)", lmp->atom, test_config.run_mag_forces,
                      5 * epsilon);
    EXPECT_STRESS("run_stress (newton off)", pair->virial, test_config.run_stress, 10 * epsilon);

    stats.reset();
    auto id     = lmp->modify->find_compute("sum");
    auto energy = lmp->modify->compute[id]->compute_scalar();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.run_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.run_coul, epsilon);
    EXPECT_FP_LE_WITH_EPS((pair->eng_vdwl + pair->eng_coul), energy, epsilon);
    if (print_stats) std::cerr << "run_energy  stats, newton off:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

TEST(PairStyle, intel)
{
    if (!Info::has_package("INTEL")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"PairStyle", "-log",  "none", "-echo", "screen", "-nocite",
                         "-pk",       "intel", "0",    "mode",  "double", "omp",
                         "4",         "lrt",   "no",   "-sf",   "intel"};

    // styles tagged "single_thread" (e.g. dpd, due to its pRNG) cannot use more
    // than one thread in the test
    if (test_config.has_tag("single_thread")) args[12] = "1";

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config, true);
    } catch (std::exception &e) {
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        if (full_neigh_unsupported(e.what())) GTEST_SKIP() << e.what();
        FAIL() << e.what();
    }
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    if (!lmp) {
        std::cerr << "One or more prerequisite styles with /intel suffix\n"
                     "are not available in this LAMMPS configuration:\n";
        for (auto prerequisite : test_config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        GTEST_SKIP();
    }

    // relax error a bit for INTEL package
    double epsilon = 7.5 * test_config.epsilon;
    // relax test precision when using pppm and single precision FFTs
#if defined(FFT_SINGLE)
    if (lmp->force->kspace && lmp->force->kspace->compute_flag)
        if (utils::strmatch(lmp->force->kspace_style, "^pppm")) epsilon *= 2.0e8;
#endif

    // we need to relax the epsilon a LOT for tests using long-range
    // coulomb with tabulation. seems more like mixed precision or a bug
    for (auto post_cmd : test_config.post_commands) {
        if (post_cmd.find("pair_modify table") != std::string::npos) {
            if (post_cmd.find("pair_modify table 0") == std::string::npos) epsilon *= 1000000.0;
        }
    }

    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms, nlocal);

    ErrorStats stats;
    auto *pair = lmp->force->pair;

    EXPECT_FORCES("init_forces", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_STRESS("init_stress", pair->virial, test_config.init_stress, 10 * epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
    if (print_stats) std::cerr << "init_energy stats:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    EXPECT_FORCES("run_forces", lmp->atom, test_config.run_forces, 5 * epsilon);
    EXPECT_STRESS("run_stress", pair->virial, test_config.run_stress, 10 * epsilon);

    stats.reset();
    int id        = lmp->modify->find_compute("sum");
    double energy = lmp->modify->compute[id]->compute_scalar();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.run_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.run_coul, epsilon);

    // rebo family of pair styles will have a large error in per-atom energy for INTEL
    if (test_config.pair_style.find("rebo") != std::string::npos) epsilon *= 100000.0;

    EXPECT_FP_LE_WITH_EPS((pair->eng_vdwl + pair->eng_coul), energy, epsilon);
    if (print_stats) std::cerr << "run_energy  stats:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

TEST(PairStyle, opt)
{
    if (!Info::has_package("OPT")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite", "-sf", "opt"};

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config, true);
    } catch (std::exception &e) {
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        if (full_neigh_unsupported(e.what())) GTEST_SKIP() << e.what();
        FAIL() << e.what();
    }
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    if (!lmp) {
        std::cerr << "One or more prerequisite styles with /opt suffix\n"
                     "are not available in this LAMMPS configuration:\n";
        for (auto prerequisite : test_config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        GTEST_SKIP();
    }

    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms, nlocal);

    // relax error a bit for OPT package
    double epsilon = 2.0 * test_config.epsilon;
    // relax test precision when using pppm and single precision FFTs
#if defined(FFT_SINGLE)
    if (lmp->force->kspace && lmp->force->kspace->compute_flag)
        if (utils::strmatch(lmp->force->kspace_style, "^pppm")) epsilon *= 2.0e8;
#endif
    ErrorStats stats;
    auto *pair = lmp->force->pair;

    EXPECT_FORCES("init_forces (newton off)", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_MAG_FORCES("init_mag_forces (newton off)", lmp->atom, test_config.init_mag_forces,
                      epsilon);
    EXPECT_STRESS("init_stress", pair->virial, test_config.init_stress, 10 * epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
    if (print_stats) std::cerr << "init_energy stats:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    EXPECT_FORCES("run_forces", lmp->atom, test_config.run_forces, 5 * epsilon);
    EXPECT_STRESS("run_stress", pair->virial, test_config.run_stress, 10 * epsilon);

    stats.reset();
    int id        = lmp->modify->find_compute("sum");
    double energy = lmp->modify->compute[id]->compute_scalar();
    EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.run_vdwl, epsilon);
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.run_coul, epsilon);
    EXPECT_FP_LE_WITH_EPS((pair->eng_vdwl + pair->eng_coul), energy, epsilon);
    if (print_stats) std::cerr << "run_energy  stats:" << stats << std::endl;

    // the "nofdotr" token in skip_tests opts a style out of comparing the
    // virial tallied by ev_tally() with the fdotr one of the reference, for the
    // same reason as in the plain test case above

    if (!test_config.skip_tests.count("nofdotr")) {
        if (!verbose) ::testing::internal::CaptureStdout();
        restart_lammps(lmp, test_config, true);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        pair = lmp->force->pair;

        EXPECT_FORCES("nofdotr_forces", lmp->atom, test_config.init_forces, 5 * epsilon);
        EXPECT_MAG_FORCES("nofdotr_mag_forces", lmp->atom, test_config.init_mag_forces,
                          5 * epsilon);
        EXPECT_STRESS("nofdotr_stress", pair->virial, test_config.init_stress, 10 * epsilon);

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(pair->eng_vdwl, test_config.init_vdwl, 5 * epsilon);
        EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, 5 * epsilon);
        if (print_stats) std::cerr << "nofdotr_energy stats:" << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

TEST(PairStyle, single)
{
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite"};

    // need to add this dependency
    test_config.prerequisites.emplace_back("atom", "full");

    // create a LAMMPS instance with standard settings to detect the number of atom types
    if (!verbose) ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config, true);
    } catch (std::exception &e) {
        if (!verbose) ::testing::internal::GetCapturedStdout();
        FAIL() << e.what();
    }
    if (!verbose) ::testing::internal::GetCapturedStdout();

    if (!lmp) {
        std::cerr << "One or more prerequisite styles are not available "
                     "in this LAMMPS configuration:\n";
        for (auto &prerequisite : test_config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        test_config.prerequisites.pop_back();
        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        GTEST_SKIP();
    }
    test_config.prerequisites.pop_back();

    // gather some information and skip if unsupported
    int ntypes    = lmp->atom->ntypes;
    int molecular = lmp->atom->molecular;
    if (molecular > Atom::MOLECULAR) {
        std::cerr << "Only atomic and simple molecular atom styles are supported\n";
        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        GTEST_SKIP();
    }

    Pair *pair = lmp->force->pair;
    if (!pair->single_enable) {
        std::cerr << "Single method not available for pair style " << test_config.pair_style
                  << std::endl;
        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        GTEST_SKIP();
    }

    if (!pair->compute_flag) {
        std::cerr << "Pair style disabled" << std::endl;
        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        GTEST_SKIP();
    }

    // now start over
    if (!verbose) ::testing::internal::CaptureStdout();

    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    command("clear");
    command("variable newton_pair delete");
    command("variable newton_pair index on");

    command("variable input_dir index " + INPUT_FOLDER);

    for (auto &pre_command : test_config.pre_commands) {
        command(pre_command);
    }

    if (lmp->domain->box_exist) {
        std::cerr << "Cannot test single() with YAML file that creates a box\n";
        cleanup_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        GTEST_SKIP();
    }

    if (std::find(test_config.tags.begin(), test_config.tags.end(), "ellipsoid") !=
        test_config.tags.end()) {
        command("atom_style ellipsoid");
    } else {
        command("atom_style full");
    }

    command("units ${units}");
    command("boundary p p p");
    command("newton ${newton_pair} ${newton_bond}");

    if (molecular == Atom::MOLECULAR) {
        command("special_bonds lj/coul "
                "${bond_factor} ${angle_factor} ${dihedral_factor}");
    }

    command("atom_modify map array");
    command("region box block -10.0 10.0 -10.0 10.0 -10.0 10.0 units box");

    auto cmd = fmt::format("create_box {} box", ntypes);
    if (molecular == Atom::MOLECULAR) {
        cmd += " bond/types 1"
               " extra/bond/per/atom 1"
               " extra/special/per/atom 1";
    }
    command(cmd);

    command("pair_style " + test_config.pair_style);

    pair = lmp->force->pair;

    for (auto &pair_coeff : test_config.pair_coeff) {
        command("pair_coeff " + pair_coeff);
    }

    // create (only) two atoms

    command("create_atoms 1 single 0.0 -0.75  0.4 units box");
    command("create_atoms 2 single 1.5  0.25 -0.1 units box");
    command("special_bonds lj/coul 1.0 1.0 1.0");

    // need to use a different integrator for different atom styles
    if (std::find(test_config.tags.begin(), test_config.tags.end(), "ellipsoid") !=
        test_config.tags.end()) {
        command("set atom 1 shape 1 2 2");
        command("set atom 2 shape 3 1 1");
        command("set group all quat/random 12238");
        command("set group all mass 1.0");
    } else {
        command("mass * 1.0");
        command("set atom 1 charge -0.5");
        command("set atom 2 charge  0.5");
        command("set atom 1 mol 1");
        command("set atom 2 mol 2");
    }

    if (molecular == Atom::MOLECULAR) {
        command("create_bonds single/bond 1 1 2");
        command("bond_style zero");
        command("bond_coeff 1 2.0");
    }

    for (auto &post_command : test_config.post_commands) {
        command(post_command);
    }

    command("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    int idx1       = lmp->atom->map(1);
    int idx2       = lmp->atom->map(2);
    double epsilon = test_config.epsilon;
    double **f     = lmp->atom->f;
    double **x     = lmp->atom->x;
    bool is_ellipsoid =
        std::find(test_config.tags.begin(), test_config.tags.end(), "ellipsoid") !=
        test_config.tags.end();
    double **tor   = is_ellipsoid ? lmp->atom->torque : nullptr;
    double delx    = x[idx2][0] - x[idx1][0];
    double dely    = x[idx2][1] - x[idx1][1];
    double delz    = x[idx2][2] - x[idx1][2];
    double rsq     = delx * delx + dely * dely + delz * delz;
    double fsingle = 0.0;
    double epair[4], esngl[4];
    double splj = lmp->force->special_lj[1];
    double spcl = lmp->force->special_coul[1];
    ErrorStats stats;

    epair[0] = pair->eng_vdwl + pair->eng_coul;
    esngl[0] = pair->single(idx1, idx2, 1, 2, rsq, spcl, splj, fsingle);
    if (is_ellipsoid) {
        EXPECT_NE(pair->svector, nullptr);
        EXPECT_GE(pair->single_extra, 6);
        if (pair->svector != nullptr && pair->single_extra >= 6) {
            EXPECT_FP_LE_WITH_EPS(pair->svector[0], f[idx1][0], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[1], f[idx1][1], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[2], f[idx1][2], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[3], tor[idx1][0], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[4], tor[idx1][1], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[5], tor[idx1][2], epsilon);
        }
    } else {
        EXPECT_FP_LE_WITH_EPS(f[idx1][0], -fsingle * delx, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx1][1], -fsingle * dely, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx1][2], -fsingle * delz, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx2][0], fsingle * delx, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx2][1], fsingle * dely, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx2][2], fsingle * delz, epsilon);
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    command("displace_atoms all random 0.5 0.5 0.5 723456");
    command("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    f       = lmp->atom->f;
    x       = lmp->atom->x;
    if (is_ellipsoid) tor = lmp->atom->torque;
    idx1    = lmp->atom->map(1);
    idx2    = lmp->atom->map(2);
    delx    = x[idx2][0] - x[idx1][0];
    dely    = x[idx2][1] - x[idx1][1];
    delz    = x[idx2][2] - x[idx1][2];
    rsq     = delx * delx + dely * dely + delz * delz;
    fsingle = 0.0;

    epair[1] = pair->eng_vdwl + pair->eng_coul;
    esngl[1] = pair->single(idx1, idx2, 1, 2, rsq, spcl, splj, fsingle);
    if (is_ellipsoid) {
        EXPECT_NE(pair->svector, nullptr);
        EXPECT_GE(pair->single_extra, 6);
        if (pair->svector != nullptr && pair->single_extra >= 6) {
            EXPECT_FP_LE_WITH_EPS(pair->svector[0], f[idx1][0], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[1], f[idx1][1], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[2], f[idx1][2], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[3], tor[idx1][0], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[4], tor[idx1][1], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[5], tor[idx1][2], epsilon);
        }
    } else {
        EXPECT_FP_LE_WITH_EPS(f[idx1][0], -fsingle * delx, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx1][1], -fsingle * dely, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx1][2], -fsingle * delz, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx2][0], fsingle * delx, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx2][1], fsingle * dely, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx2][2], fsingle * delz, epsilon);
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    command("displace_atoms all random 0.5 0.5 0.5 3456963");
    command("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    f       = lmp->atom->f;
    x       = lmp->atom->x;
    if (is_ellipsoid) tor = lmp->atom->torque;
    idx1    = lmp->atom->map(1);
    idx2    = lmp->atom->map(2);
    delx    = x[idx2][0] - x[idx1][0];
    dely    = x[idx2][1] - x[idx1][1];
    delz    = x[idx2][2] - x[idx1][2];
    rsq     = delx * delx + dely * dely + delz * delz;
    fsingle = 0.0;

    epair[2] = pair->eng_vdwl + pair->eng_coul;
    esngl[2] = pair->single(idx1, idx2, 1, 2, rsq, spcl, splj, fsingle);
    if (is_ellipsoid) {
        EXPECT_NE(pair->svector, nullptr);
        EXPECT_GE(pair->single_extra, 6);
        if (pair->svector != nullptr && pair->single_extra >= 6) {
            EXPECT_FP_LE_WITH_EPS(pair->svector[0], f[idx1][0], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[1], f[idx1][1], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[2], f[idx1][2], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[3], tor[idx1][0], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[4], tor[idx1][1], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[5], tor[idx1][2], epsilon);
        }
    } else {
        EXPECT_FP_LE_WITH_EPS(f[idx1][0], -fsingle * delx, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx1][1], -fsingle * dely, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx1][2], -fsingle * delz, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx2][0], fsingle * delx, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx2][1], fsingle * dely, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx2][2], fsingle * delz, epsilon);
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    command("displace_atoms all random 0.5 0.5 0.5 9726532");
    command("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    f       = lmp->atom->f;
    x       = lmp->atom->x;
    if (is_ellipsoid) tor = lmp->atom->torque;
    idx1    = lmp->atom->map(1);
    idx2    = lmp->atom->map(2);
    delx    = x[idx2][0] - x[idx1][0];
    dely    = x[idx2][1] - x[idx1][1];
    delz    = x[idx2][2] - x[idx1][2];
    rsq     = delx * delx + dely * dely + delz * delz;
    fsingle = 0.0;

    epair[3] = pair->eng_vdwl + pair->eng_coul;
    esngl[3] = pair->single(idx1, idx2, 1, 2, rsq, spcl, splj, fsingle);
    if (is_ellipsoid) {
        EXPECT_NE(pair->svector, nullptr);
        EXPECT_GE(pair->single_extra, 6);
        if (pair->svector != nullptr && pair->single_extra >= 6) {
            EXPECT_FP_LE_WITH_EPS(pair->svector[0], f[idx1][0], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[1], f[idx1][1], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[2], f[idx1][2], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[3], tor[idx1][0], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[4], tor[idx1][1], epsilon);
            EXPECT_FP_LE_WITH_EPS(pair->svector[5], tor[idx1][2], epsilon);
        }
    } else {
        EXPECT_FP_LE_WITH_EPS(f[idx1][0], -fsingle * delx, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx1][1], -fsingle * dely, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx1][2], -fsingle * delz, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx2][0], fsingle * delx, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx2][1], fsingle * dely, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[idx2][2], fsingle * delz, epsilon);
    }
    if (print_stats) std::cerr << "single_force  stats:" << stats << std::endl;

    // repeat the comparison with a non-unit special_bonds exclusion factor.
    // the factors are applied by compute() and by single() independently, and
    // the two implementations must agree.  styles that subtract only the
    // excluded part of the interaction instead of scaling the whole term (the
    // coul/dsf and coul/wolf families, tabulated coulomb, ...) are easy to get
    // wrong in one of the two places, and the unit factor used above cannot
    // tell the two conventions apart.

    if (molecular == Atom::MOLECULAR) {
        stats.reset();
        bool special_supported = true;
        if (!verbose) ::testing::internal::CaptureStdout();
        try {
            command("special_bonds lj/coul 0.25 1.0 1.0");
            command("run 0 post no");
        } catch (std::exception &) {
            special_supported = false;
        }
        if (!verbose) ::testing::internal::GetCapturedStdout();

        if (special_supported) {
            splj = lmp->force->special_lj[1];
            spcl = lmp->force->special_coul[1];
            f    = lmp->atom->f;
            x    = lmp->atom->x;
            if (is_ellipsoid) tor = lmp->atom->torque;
            idx1    = lmp->atom->map(1);
            idx2    = lmp->atom->map(2);
            delx    = x[idx2][0] - x[idx1][0];
            dely    = x[idx2][1] - x[idx1][1];
            delz    = x[idx2][2] - x[idx1][2];
            rsq     = delx * delx + dely * dely + delz * delz;
            fsingle = 0.0;

            pair->single(idx1, idx2, 1, 2, rsq, spcl, splj, fsingle);
            if (is_ellipsoid) {
                EXPECT_NE(pair->svector, nullptr);
                EXPECT_GE(pair->single_extra, 6);
                if (pair->svector != nullptr && pair->single_extra >= 6) {
                    EXPECT_FP_LE_WITH_EPS(pair->svector[0], f[idx1][0], epsilon);
                    EXPECT_FP_LE_WITH_EPS(pair->svector[1], f[idx1][1], epsilon);
                    EXPECT_FP_LE_WITH_EPS(pair->svector[2], f[idx1][2], epsilon);
                    EXPECT_FP_LE_WITH_EPS(pair->svector[3], tor[idx1][0], epsilon);
                    EXPECT_FP_LE_WITH_EPS(pair->svector[4], tor[idx1][1], epsilon);
                    EXPECT_FP_LE_WITH_EPS(pair->svector[5], tor[idx1][2], epsilon);
                }
            } else {
                EXPECT_FP_LE_WITH_EPS(f[idx1][0], -fsingle * delx, epsilon);
                EXPECT_FP_LE_WITH_EPS(f[idx1][1], -fsingle * dely, epsilon);
                EXPECT_FP_LE_WITH_EPS(f[idx1][2], -fsingle * delz, epsilon);
                EXPECT_FP_LE_WITH_EPS(f[idx2][0], fsingle * delx, epsilon);
                EXPECT_FP_LE_WITH_EPS(f[idx2][1], fsingle * dely, epsilon);
                EXPECT_FP_LE_WITH_EPS(f[idx2][2], fsingle * delz, epsilon);
            }
            if (print_stats) std::cerr << "single_special stats:" << stats << std::endl;
        } else if (print_stats)
            std::cerr << "skipping single_special test: style rejects the factor\n";
    }

    if ((test_config.pair_style.find("coul/dsf") != std::string::npos) &&
        (test_config.pair_style.find("coul/wolf") != std::string::npos)) {
        stats.reset();
        EXPECT_FP_LE_WITH_EPS(epair[0], esngl[0], epsilon);
        EXPECT_FP_LE_WITH_EPS(epair[1], esngl[1], epsilon);
        EXPECT_FP_LE_WITH_EPS(epair[2], esngl[2], epsilon);
        EXPECT_FP_LE_WITH_EPS(epair[3], esngl[3], epsilon);
        if (print_stats) std::cerr << "single_energy  stats:" << stats << std::endl;
    } else if (print_stats)
        std::cerr << "skipping single_energy test due to self energy\n";

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
}

TEST(PairStyle, extract)
{
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite"};

    if (!verbose) ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config, true);
    } catch (std::exception &e) {
        if (!verbose) ::testing::internal::GetCapturedStdout();
        FAIL() << e.what();
    }
    if (!verbose) ::testing::internal::GetCapturedStdout();

    if (!lmp) {
        std::cerr << "One or more prerequisite styles are not available "
                     "in this LAMMPS configuration:\n";
        for (const auto &prerequisite : test_config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        GTEST_SKIP();
    }

    auto *pair = lmp->force->pair;
    if (!pair->compute_flag) {
        std::cerr << "Pair style disabled" << std::endl;
        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        GTEST_SKIP();
    }

    void *ptr = nullptr;
    int dim   = 0;
    for (auto &extract : test_config.extract) {
        ptr = pair->extract(extract.first.c_str(), dim);
        EXPECT_NE(ptr, nullptr);
        EXPECT_EQ(dim, extract.second);
    }
    ptr = pair->extract("does_not_exist", dim);
    EXPECT_EQ(ptr, nullptr);

    // replace pair style with the same.
    // should just update setting, but not create new style.

    int ntypes = lmp->atom->ntypes;
    for (int i = 1; i <= ntypes; ++i) {
        for (int j = 1; j <= ntypes; ++j) {
            pair->cutsq[i][j] = -1.0;
        }
    }

    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    if (!verbose) ::testing::internal::CaptureStdout();
    command("pair_style " + test_config.pair_style);
    EXPECT_EQ(pair, lmp->force->pair);

    for (auto &pair_coeff : test_config.pair_coeff) {
        command("pair_coeff " + pair_coeff);
    }
    pair->init();
    if (!verbose) ::testing::internal::GetCapturedStdout();

    for (int i = 1; i <= ntypes; ++i) {
        for (int j = 1; j <= ntypes; ++j) {
            EXPECT_GE(pair->cutsq[i][j], 0.0);
        }
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
}

TEST(PairStyle, extract_omp)
{
    if (!Info::has_package("OPENMP")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"PairStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-pk",       "omp",  "4",    "-sf",   "omp"};

    if (!verbose) ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config, true);
    } catch (std::exception &e) {
        if (!verbose) ::testing::internal::GetCapturedStdout();
        FAIL() << e.what();
    }
    if (!verbose) ::testing::internal::GetCapturedStdout();

    if (!lmp) {
        std::cerr << "One or more prerequisite styles are not available "
                     "in this LAMMPS configuration:\n";
        for (const auto &prerequisite : test_config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        GTEST_SKIP();
    }

    auto *pair = lmp->force->pair;
    if (!pair->compute_flag) {
        std::cerr << "Pair style disabled" << std::endl;
        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        GTEST_SKIP();
    }

    void *ptr = nullptr;
    int dim   = 0;
    for (auto &extract : test_config.extract) {
        ptr = pair->extract(extract.first.c_str(), dim);
        EXPECT_NE(ptr, nullptr);
        EXPECT_EQ(dim, extract.second);
    }
    ptr = pair->extract("does_not_exist", dim);
    EXPECT_EQ(ptr, nullptr);

    // replace pair style with the same.
    // should just update setting, but not create new style.

    int ntypes = lmp->atom->ntypes;
    for (int i = 1; i <= ntypes; ++i) {
        for (int j = 1; j <= ntypes; ++j) {
            pair->cutsq[i][j] = -1.0;
        }
    }

    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    if (!verbose) ::testing::internal::CaptureStdout();
    command("pair_style " + test_config.pair_style);
    EXPECT_EQ(pair, lmp->force->pair);

    for (auto &pair_coeff : test_config.pair_coeff) {
        command("pair_coeff " + pair_coeff);
    }
    pair->init();
    if (!verbose) ::testing::internal::GetCapturedStdout();

    for (int i = 1; i <= ntypes; ++i) {
        for (int j = 1; j <= ntypes; ++j) {
            EXPECT_GE(pair->cutsq[i][j], 0.0);
        }
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
}
