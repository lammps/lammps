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
#include "exceptions.h"
#include "force.h"
#include "info.h"
#include "input.h"
#include "kspace.h"
#include "lammps.h"
#include "neighbor.h"
#include "pair.h"
#include "platform.h"
#include "utils.h"

#include <cstdio>
#include <cstdlib>
#include <mpi.h>

#include <string>
#include <utility>

using ::testing::HasSubstr;
using ::testing::StartsWith;

using namespace LAMMPS_NS;

void cleanup_lammps(LAMMPS *lmp, const TestConfig &cfg)
{
    platform::unlink(cfg.basename + ".restart");
    platform::unlink(cfg.basename + ".data");
    platform::unlink(cfg.basename + "-coeffs.in");
    delete lmp;
}

LAMMPS *init_lammps(LAMMPS::argv &args, const TestConfig &cfg, const bool newton = true,
                    const bool etypes = false)
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
        try {
            lmp->input->one(line);
        } catch (LAMMPSAbortException &ae) {
            fprintf(stderr, "LAMMPS Error: %s\n", ae.what());
            exit(2);
        } catch (LAMMPSException &e) {
            fprintf(stderr, "LAMMPS Error: %s\n", e.what());
            exit(3);
        } catch (fmt::format_error &fe) {
            fprintf(stderr, "fmt::format_error: %s\n", fe.what());
            exit(4);
        } catch (std::exception &e) {
            fprintf(stderr, "General exception: %s\n", e.what());
            exit(5);
        }
    };

    auto parse_input_script = [&](const std::string &filename) {
        lmp->input->file(filename.c_str());
    };

    if (newton) {
        command("variable newton_pair index on");
    } else {
        command("variable newton_pair index off");
    }

    if (etypes) {
        command("variable etypes index on");
    } else {
        command("variable etypes index off");
    }

    command("variable input_dir index " + INPUT_FOLDER);

    for (const auto &pre_command : cfg.pre_commands) {
        command(pre_command);
    }

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
    command("write_data " + cfg.basename + ".data pair ${write_data_pair} nofix");
    command("write_coeff " + cfg.basename + "-coeffs.in");

    return lmp;
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
    LAMMPS::argv args = {"FixElectrode", "-log", "none", "-echo", "screen", "-nocite"};

    LAMMPS *lmp = init_lammps(args, config);
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

    // init_coul
    writer.emit("init_coul", lmp->force->pair->eng_coul);

    // init_charges
    block.clear();
    double *q = lmp->atom->q;
    for (int i = 1; i <= natoms; ++i) {
        const int j = lmp->atom->map(i);
        block += fmt::format("{:3} {:23.16e}\n", i, q[j]);
    }
    writer.emit_block("init_charges", block);

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

    cleanup_lammps(lmp, config);
}

TEST(FixElectrode, plain)
{
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"FixElectrode", "-log", "none", "-echo", "screen", "-nocite"};

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = init_lammps(args, test_config, true);

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

    int const nlist_etypes_off = lmp->neighbor->nlist;

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

    EXPECT_CHARGES("init_charges (newton on)", lmp->atom, test_config.init_charges, epsilon);
    EXPECT_FORCES("init_forces (newton on)", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_STRESS("init_stress (newton on)", pair->virial, test_config.init_stress, epsilon);

    ErrorStats stats;
    stats.reset();
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
    if (print_stats) std::cerr << "init_energy stats, newton on: " << stats << std::endl;

    // etypes on
    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    lmp = init_lammps(args, test_config, true, true);
    if (!verbose) ::testing::internal::GetCapturedStdout();
    // skip over these tests if etypes keyword is not used
    if (lmp->neighbor->nlist > nlist_etypes_off) {
        pair = lmp->force->pair;
        EXPECT_CHARGES("init_charges (etypes on)", lmp->atom, test_config.init_charges, epsilon);
        EXPECT_FORCES("init_forces (etypes on)", lmp->atom, test_config.init_forces, epsilon);
        EXPECT_STRESS("init_stress (etypes on)", pair->virial, test_config.init_stress,
                      3 * epsilon);

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
        if (print_stats) std::cerr << "init_energy stats, etypes on:" << stats << std::endl;
    }

    // newton off
    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    lmp = init_lammps(args, test_config, false);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // skip over these tests if newton pair is forced to be on
    if (lmp->force->newton_pair == 0) {
        pair = lmp->force->pair;

        EXPECT_CHARGES("init_charges (newton off)", lmp->atom, test_config.init_charges, epsilon);
        EXPECT_FORCES("init_forces (newton off)", lmp->atom, test_config.init_forces, epsilon);
        EXPECT_STRESS("init_stress (newton off)", pair->virial, test_config.init_stress,
                      3 * epsilon);

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
        if (print_stats) std::cerr << "init_energy stats, newton off:" << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    restart_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    pair = lmp->force->pair;

    EXPECT_CHARGES("restart_charges", lmp->atom, test_config.init_charges, epsilon);
    EXPECT_FORCES("restart_forces", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_STRESS("restart_stress", pair->virial, test_config.init_stress, epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
    if (print_stats) std::cerr << "restart_energy stats:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    restart_lammps(lmp, test_config, true);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    pair = lmp->force->pair;
    EXPECT_CHARGES("nofdotr_charges", lmp->atom, test_config.init_charges, epsilon);
    EXPECT_FORCES("nofdotr_forces", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_STRESS("nofdotr_stress", pair->virial, test_config.init_stress, epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
    if (print_stats) std::cerr << "nofdotr_energy stats:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    data_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    pair = lmp->force->pair;
    EXPECT_CHARGES("data_charges", lmp->atom, test_config.init_charges, epsilon);
    EXPECT_FORCES("data_forces", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_STRESS("data_stress", pair->virial, test_config.init_stress, epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
    if (print_stats) std::cerr << "data_energy stats:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

TEST(FixElectrode, intel)
{
    if (!Info::has_package("INTEL")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"FixElectrode", "-log", "none", "-echo",  "screen", "-nocite", "-pk",
                         "intel",        "0",    "mode", "double", "omp",    "4",       "lrt",
                         "no",           "-sf",  "intel"};

    // cannot use more than 1 thread for dpd styles due to pRNG
    if (utils::strmatch(test_config.pair_style, "^dpd")) args[12] = "1";

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = init_lammps(args, test_config, true);

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
    double epsilon = 20 * test_config.epsilon;
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

    EXPECT_CHARGES("init_charges", lmp->atom, test_config.init_charges, epsilon);
    EXPECT_FORCES("init_forces", lmp->atom, test_config.init_forces, 10 * epsilon);
    EXPECT_STRESS("init_stress", pair->virial, test_config.init_stress, 10 * epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(pair->eng_coul, test_config.init_coul, epsilon);
    if (print_stats) std::cerr << "init_energy stats:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

TEST(FixElectrode, extract)
{
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"FixElectrode", "-log", "none", "-echo", "screen", "-nocite"};

    if (!verbose) ::testing::internal::CaptureStdout();
    LAMMPS *lmp = init_lammps(args, test_config, true);
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
