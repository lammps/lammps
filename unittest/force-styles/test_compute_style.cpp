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

// unit tests for compute styles intended for molecular systems

#include "error_stats.h"
#include "test_config.h"
#include "test_config_reader.h"
#include "test_main.h"
#include "yaml_reader.h"
#include "yaml_writer.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "atom.h"
#include "compute.h"
#include "fmt/format.h"
#include "force.h"
#include "info.h"
#include "input.h"
#include "kspace.h"
#include "lammps.h"
#include "modify.h"
#include "pair.h"
#include "platform.h"
#include "universe.h"
#include "update.h"
#include "utils.h"
#include "variable.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

using ::testing::HasSubstr;
using ::testing::StartsWith;

using namespace LAMMPS_NS;

void cleanup_lammps(LAMMPS *&lmp, const TestConfig &cfg)
{
    delete lmp;
    lmp = nullptr;
}

LAMMPS *init_lammps(LAMMPS::argv &args, const TestConfig &cfg, const bool use_respa)
{
    LAMMPS *lmp;

    lmp = new LAMMPS(args, MPI_COMM_WORLD);

    // check if prerequisite styles are available
    Info *info = new Info(lmp);
    int nfail  = 0;
    int reaxff_flag = 0;
    for (const auto &prerequisite : cfg.prerequisites) {
        std::string style = prerequisite.second;

        // this is a test for compute styles, so if the suffixed
        // version is not available, there is no reason to test.
        if (prerequisite.first == "compute") {
            if (lmp->suffix_enable) {
                style += "/";
                style += lmp->suffix;
            }
        }

        if (prerequisite.first == "pair" && prerequisite.second == "reaxff")
            reaxff_flag = 1;

        if (!info->has_style(prerequisite.first, style)) ++nfail;
    }
    delete info;
    if (nfail > 0) {
        cleanup_lammps(lmp, cfg);
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

    if (use_respa) command("run_style respa 2 1 bond 1 pair 2");

    if (!reaxff_flag) {
        // set up molecular system force field

        command("pair_style lj/cut 8.0");
        command("pair_coeff  1 1  0.02   2.5");
        command("pair_coeff  2 2  0.005  1.0");
        command("pair_coeff  2 4  0.005  0.5");
        command("pair_coeff  3 3  0.02   3.2");
        command("pair_coeff  4 4  0.015  3.1");
        command("pair_coeff  5 5  0.015  3.1");
        command("bond_style harmonic");
        command("bond_coeff  1 250.0 1.5");
        command("bond_coeff  2 300.0 1.1");
        command("bond_coeff  3 350.0 1.3");
        command("bond_coeff  4 650.0 1.2");
        command("bond_coeff  5 450.0 1.0");
        command("angle_style harmonic");
        command("angle_coeff  1  75.0 110.1");
        command("angle_coeff  2  45.0 111.0");
        command("angle_coeff  3  50.0 120.0");
        command("angle_coeff  4 100.0 108.5");
        command("group solute  molecule 1:2");
        command("group solvent molecule 3:5");
    }
    
    for (const auto &post_command : cfg.post_commands)
        command(post_command);

    command("timestep 0.25");
    command("run 0 post no");
    command("thermo 2");
    command("run 4 post no start 0 stop 8");
    return lmp;
}

// re-generate yaml file with current settings.

void generate_yaml_file(const char *outfile, const TestConfig &config)
{
    // initialize system geometry
    LAMMPS::argv args = {"ComputeStyle", "-log", "none", "-echo", "screen", "-nocite"};
    LAMMPS *lmp       = nullptr;
    try {
        lmp = init_lammps(args, config, false);
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

    // natoms
    writer.emit("natoms", natoms);

    auto *icompute = lmp->modify->get_compute_by_id("test");
    if (!icompute) {
        std::cerr << "ERROR: no compute defined with compute ID 'test'\n";
        exit(1);
    } else {

        // global scalar
        if (icompute->scalar_flag) {
            double value = icompute->compute_scalar();
            writer.emit("global_scalar", value);
        }

        // global vector
        if (icompute->vector_flag) {
            int num = icompute->size_vector;
            block   = std::to_string(num);
            icompute->compute_vector();
            for (int i = 0; i < num; ++i)
                block += fmt::format(" {}", icompute->vector[i]);
            writer.emit_block("global_vector", block);
        }

        // per-atom vector
        if (icompute->peratom_flag && icompute->size_peratom_cols == 0) {

            block = std::to_string(icompute->size_peratom_cols);
            icompute->compute_peratom();
         
            for (int i = 1; i <= natoms; ++i) {
                const int j = lmp->atom->map(i);
                if (j >= 0 && j < lmp->atom->nlocal)
                    block += fmt::format(" {} {}\n", i, icompute->vector_atom[j]);
            }
            writer.emit_block("peratom_vector", block);
        }

        // per-atom array
        if (icompute->peratom_flag && icompute->size_peratom_cols > 0) {

            block = std::to_string(icompute->size_peratom_cols) + "\n";
            icompute->compute_peratom();

            for (int i = 1; i <= natoms; ++i) {
                const int j = lmp->atom->map(i);
                if (j >= 0 && j < lmp->atom->nlocal) {
                    block += fmt::format(" {}", i);
                    for (int k = 0; k < icompute->size_peratom_cols; ++k)
                        block += fmt::format(" {}", icompute->array_atom[j][k]);
                    block += "\n";
                }
            }
            writer.emit_block("peratom_array", block);
        }

    }

    cleanup_lammps(lmp, config);
}

TEST(ComputeStyle, plain)
{
    if (!LAMMPS::is_installed_pkg("MOLECULE")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
#if defined(USING_STATIC_LIBS)
    if (test_config.skip_tests.count("static")) GTEST_SKIP();
#endif

    LAMMPS::argv args = {"ComputeStyle", "-log", "none", "-echo", "screen", "-nocite"};

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config, false);
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
    EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms, nlocal);

    double epsilon = test_config.epsilon;
    // relax test precision when using pppm and single precision FFTs
#if defined(FFT_SINGLE)
    if (lmp->force->kspace && utils::strmatch(lmp->force->kspace_style, "^pppm")) epsilon *= 2.0e8;
#endif

    ErrorStats stats;

    auto *icompute = lmp->modify->get_compute_by_id("test");
    if (!icompute) {
        FAIL() << "ERROR: no compute defined with compute ID 'test'\n";
    } else {

        stats.reset();
        // global scalar
        if (icompute->scalar_flag) {
            double value = icompute->compute_scalar();
            EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
        }

        // global vector
        if (icompute->vector_flag) {
            int num = icompute->size_vector;
            EXPECT_EQ(num, test_config.global_vector.size());

            icompute->compute_vector();

            for (int i = 0; i < num; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], icompute->vector[i],
                                      epsilon);
        }

        // per-atom vector
        if (icompute->peratom_flag && icompute->size_peratom_cols == 0) {
            icompute->compute_peratom();

            for (const auto &entry : test_config.peratom_vector) {
                int index = entry.first;
                double value = entry.second;
                
                int j = lmp->atom->map(index);
                if (j >= 0 && j < lmp->atom->nlocal) {
                    EXPECT_FP_LE_WITH_EPS(value, icompute->vector_atom[j], epsilon);
                }
            }
        }

        // per-atom array
        if (icompute->peratom_flag && icompute->size_peratom_cols > 0) {
            icompute->compute_peratom();
            int ncols = icompute->size_peratom_cols;
            
            for (const auto &entry : test_config.peratom_array) {
                int i = entry.first;
                const auto &values = entry.second;
                
                int j = lmp->atom->map(i);
                if (j >= 0 && j < lmp->atom->nlocal) {
                    ASSERT_EQ(values.size(), ncols);
                    for (int k = 0; k < ncols; ++k) {
                        EXPECT_FP_LE_WITH_EPS(values[k], icompute->array_atom[j][k], epsilon);
                    }
                }
            }
        }

        if (print_stats && stats.has_data())
            std::cerr << "global_data, normal run, verlet: " << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

TEST(ComputeStyle, kokkos_omp)
{
    if (!LAMMPS::is_installed_pkg("KOKKOS")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
    if (!Info::has_accelerator_feature("KOKKOS", "api", "openmp")) GTEST_SKIP();
    // if KOKKOS has GPU support enabled, it *must* be used. We cannot test OpenMP only.
    if (Info::has_accelerator_feature("KOKKOS", "api", "cuda") ||
        Info::has_accelerator_feature("KOKKOS", "api", "hip") ||
        Info::has_accelerator_feature("KOKKOS", "api", "sycl")) {
        GTEST_SKIP() << "Cannot test KOKKOS/OpenMP with GPU support enabled";
    }
    LAMMPS::argv args = {"ComputeStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-k",          "on",   "t",    "4",     "-sf",    "kk"};

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config, false);
    } catch (std::exception &e) {
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        FAIL() << e.what();
    }
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    if (!lmp) {
        std::cerr << "One or more prerequisite styles with /kk suffix\n"
                     "are not available in this LAMMPS configuration:\n";
        for (auto &prerequisite : test_config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "/kk\n";
        }
        GTEST_SKIP();
    }

    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms, nlocal);

    // relax error a bit for KOKKOS package
    double epsilon = 10.0 * test_config.epsilon;
    // relax test precision when using pppm and single precision FFTs
#if defined(FFT_SINGLE)
    if (lmp->force->kspace && utils::strmatch(lmp->force->kspace_style, "^pppm")) epsilon *= 2.0e8;
#endif

    ErrorStats stats;

    auto *icompute = lmp->modify->get_compute_by_id("test");

    if (!icompute) {
        FAIL() << "ERROR: no compute defined with compute ID 'test'\n";
    } else {

        stats.reset();
        // global scalar
        if (icompute->scalar_flag) {
            double value = icompute->compute_scalar();
            EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
        }

        // global vector
        if (icompute->vector_flag) {
            int num = icompute->size_vector;
            EXPECT_EQ(num, test_config.global_vector.size());

            icompute->compute_vector();

            for (int i = 0; i < num; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], icompute->vector[i],
                                      epsilon);
        }

        // per-atom vector
        if (icompute->peratom_flag && icompute->size_peratom_cols == 0) {
            icompute->compute_peratom();

            for (const auto &entry : test_config.peratom_vector) {
                int index = entry.first;
                double value = entry.second;
                
                int j = lmp->atom->map(index);
                if (j >= 0 && j < lmp->atom->nlocal) {
                    EXPECT_FP_LE_WITH_EPS(value, icompute->vector_atom[j], epsilon);
                }
            }
        }

        // per-atom array
        if (icompute->peratom_flag && icompute->size_peratom_cols > 0) {
            icompute->compute_peratom();
            int ncols = icompute->size_peratom_cols;
            
            for (const auto &entry : test_config.peratom_array) {
                int i = entry.first;
                const auto &values = entry.second;
                
                int j = lmp->atom->map(i);
                if (j >= 0 && j < lmp->atom->nlocal) {
                    ASSERT_EQ(values.size(), ncols);
                    for (int k = 0; k < ncols; ++k) {
                        EXPECT_FP_LE_WITH_EPS(values[k], icompute->array_atom[j][k], epsilon);
                    }
                }
            }
        }

        if (print_stats && stats.has_data())
            std::cerr << "global_data, normal run, verlet: " << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};
