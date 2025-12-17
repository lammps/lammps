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
#include "grid3d.h"
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
    int snap_flag = 0;
    int pace_flag = 0;
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

        if (prerequisite.first == "pair" && prerequisite.second == "snap")
            snap_flag = 1;

        if (prerequisite.first == "compute" && prerequisite.second == "pace")
            pace_flag = 1;

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

    if ( !cfg.pair_style.empty() && cfg.pair_style != "zero" ) {

        command("pair_style " + cfg.pair_style);

        for (const auto &pair_coeff : cfg.pair_coeff)
            command("pair_coeff " + pair_coeff);

    } else if (!pace_flag) {
        // set up molecular system force field as default
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
    if (!reaxff_flag && !snap_flag && !pace_flag)
        command("run 4 post no start 0 stop 8");
    return lmp;
}

// re-generate yaml file with current settings.

void generate_yaml_file(const char *outfile, const TestConfig &config)
{
    // initialize system geometry

    LAMMPS::argv args = {"ComputeStyle", "-log", "none", "-echo", "screen", "-nocite"};

    //LAMMPS::argv args = {"ComputeStyle", "-log", "none", "-echo", "screen", "-nocite", "-k",          "on",   "t",    "4",     "-sf",    "kk"};

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

    if ( !config.pair_style.empty() && config.pair_style != "zero" ) {
        // pair_style
        writer.emit("pair_style", config.pair_style);

        // pair_coeff
        block.clear();
        for (auto pair_coeff : config.pair_coeff) {
            block += pair_coeff + "\n";
        }
        writer.emit_block("pair_coeff", block);
    }

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

        // global array
        if (icompute->array_flag) {
            icompute->compute_array();
            int nrows = icompute->size_array_rows;
            int ncols = icompute->size_array_cols;
            block = fmt::format("{} {}\n", nrows, ncols);
            for (int i = 0; i < nrows; ++i) {
                for (int j = 0; j < ncols; ++j)
                    block += fmt::format(" {}", icompute->array[i][j]);
                block += "\n";
            }
            writer.emit_block("global_array", block);
        }

        // per-atom vector
        if (icompute->peratom_flag && icompute->size_peratom_cols == 0) {
            block.clear();
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
            icompute->compute_peratom();
            block = std::to_string(icompute->size_peratom_cols) + "\n";
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

        // local vector
        if (icompute->local_flag && icompute->size_local_cols == 0) {
            icompute->compute_local();
            int nrows = icompute->size_local_rows;
            block = fmt::format("{}\n", nrows);
            for (int i = 0; i < nrows; ++i)
                block += fmt::format(" {}", icompute->vector_local[i]);
            writer.emit_block("local_vector", block);
        }

        // local array
        if (icompute->local_flag && icompute->size_local_cols > 0) {
            icompute->compute_local();
            int nrows = icompute->size_local_rows;
            int ncols = icompute->size_local_cols;
            block = fmt::format("{} {}\n", nrows, ncols);
            for (int i = 0; i < nrows; ++i) {
                for (int j = 0; j < ncols; ++j)
                    block += fmt::format(" {}", icompute->array_local[i][j]);
                block += "\n";
            }
            writer.emit_block("local_array", block);
        }

        // pergrid vector/array
        if (icompute->pergrid_flag) {
            icompute->compute_pergrid();
            writer.emit("pergrid_name", config.pergrid_name);
            writer.emit("pergrid_data", config.pergrid_data);
            int dimension, nx, ny, nz, nxlo, nxhi, nylo, nyhi, nzlo, nzhi, ncol;
            int index = icompute->get_grid_by_name(config.pergrid_name, dimension);
            if( dimension == 2 ) {
                std::cerr << "ERROR: only Grid3D supported for unit testing.\n";
                exit(1);
            }
            Grid3d *grid3d = static_cast<Grid3d*>(icompute->get_grid_by_index(index));
            grid3d->get_size(nx, ny, nz);
            grid3d->get_bounds_owned(nxlo, nxhi, nylo, nyhi, nzlo, nzhi);
            icompute->get_griddata_by_name(index, config.pergrid_data, ncol);

            if (ncol == 0) {
                block = fmt::format("{} {} {} {}\n", dimension, nx, ny, nz);
                double ***vec3d = static_cast<double ***>(icompute->get_griddata_by_index(index));
                for (int iz = nzlo; iz <= nzhi; iz++)
                    for (int iy = nylo; iy <= nyhi; iy++)
                        for (int ix = nxlo; ix <= nxhi; ix++)
                            block += fmt::format(" {}\n", vec3d[iz][iy][ix]);
                writer.emit_block("pergrid_vector", block);
            } else {
                block = fmt::format("{} {} {} {} {}\n", dimension, nx, ny, nz, ncol);
                double ****array3d = static_cast<double ****>(icompute->get_griddata_by_index(index));
                for (int iz = nzlo; iz <= nzhi; iz++)
                    for (int iy = nylo; iy <= nyhi; iy++)
                        for (int ix = nxlo; ix <= nxhi; ix++) {
                            for (int n = 0; n < ncol; n++)
                                block += fmt::format(" {}", array3d[iz][iy][ix][n]);
                            block += "\n";
                        }
                writer.emit_block("pergrid_array", block);
            }
        }

    }

    cleanup_lammps(lmp, config);
}

TEST(ComputeStyle, plain)
{
    if (!Info::has_package("MOLECULE")) GTEST_SKIP();
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

        // global array
        if (icompute->array_flag) {
            icompute->compute_array();
            int nrows = icompute->size_array_rows;
            int ncols = icompute->size_array_cols;
            ASSERT_EQ(test_config.global_array.size(), nrows);
            for (int i = 0; i < nrows; ++i) {
            const auto values = test_config.global_array[i];
                ASSERT_EQ(values.size(), ncols);
                for (int j = 0; j < ncols; ++j) {
                
                    fprintf(stderr, "*** values[%i] %f icompute->array[%i][%i] %f abs\n",
                            j, values[j], i, j, icompute->array[i][j], std::abs(values[j] - icompute->array[i][j]));
                            
                    EXPECT_FP_LE_WITH_EPS(values[j], icompute->array[i][j], epsilon);
                    
                }
            }
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

        // local vector
        if (icompute->local_flag && icompute->size_local_cols == 0) {
            icompute->compute_local();
            int nrows = icompute->size_local_rows;
            ASSERT_EQ(test_config.local_vector.size(), nrows);
            for (int i = 0; i < nrows; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.local_vector[i], icompute->vector_local[i], epsilon);
        }

        // local array
        if (icompute->local_flag && icompute->size_local_cols > 0) {
            icompute->compute_local();
            int nrows = icompute->size_local_rows;
            int ncols = icompute->size_local_cols;
            ASSERT_EQ(test_config.local_array.size(), nrows);
            for (int i = 0; i < nrows; ++i) {
                const auto values = test_config.local_array[i];
                ASSERT_EQ(values.size(), ncols);
                for (int j = 0; j < ncols; ++j) {
                    EXPECT_FP_LE_WITH_EPS(values[j], icompute->array_local[i][j], epsilon);
                }
            }
        }

        // pergrid vector/array
        if (icompute->pergrid_flag) {
            icompute->compute_pergrid();
            int dimension, nx, ny, nz, nxlo, nxhi, nylo, nyhi, nzlo, nzhi, ncol;
            int index = icompute->get_grid_by_name(test_config.pergrid_name, dimension);
            if( dimension == 2 ) {
                std::cerr << "ERROR: only Grid3D supported for unit testing.\n";
                exit(1);
            }
            Grid3d *grid3d = static_cast<Grid3d*>(icompute->get_grid_by_index(index));
            grid3d->get_size(nx, ny, nz);
            grid3d->get_bounds_owned(nxlo, nxhi, nylo, nyhi, nzlo, nzhi);
            icompute->get_griddata_by_name(index, test_config.pergrid_data, ncol);
            if (ncol == 0) {
                double ***vec3d = static_cast<double ***>(icompute->get_griddata_by_index(index));
                for (int iz = nzlo; iz <= nzhi; iz++)
                    for (int iy = nylo; iy <= nyhi; iy++)
                        for (int ix = nxlo; ix <= nxhi; ix++)
                          EXPECT_FP_LE_WITH_EPS(
                            test_config.pergrid_vector[(iz*ny + iy)*nx + ix],
                            vec3d[iz][iy][ix], epsilon );
            } else {
                double ****array3d = static_cast<double ****>(icompute->get_griddata_by_index(index));
                for (int iz = nzlo; iz <= nzhi; iz++)
                    for (int iy = nylo; iy <= nyhi; iy++)
                        for (int ix = nxlo; ix <= nxhi; ix++) {
                            for (int n = 0; n < ncol; n++)
                                EXPECT_FP_LE_WITH_EPS(
                                    test_config.pergrid_array[(iz*ny + iy)*nx + ix][n],
                                    array3d[iz][iy][ix][n], epsilon );
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
    if (!Info::has_package("KOKKOS")) GTEST_SKIP();
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

        // global array
        if (icompute->array_flag) {
            icompute->compute_array();
            int nrows = icompute->size_array_rows;
            int ncols = icompute->size_array_cols;
            ASSERT_EQ(test_config.global_array.size(), nrows);
            for (int i = 0; i < nrows; ++i) {
            const auto values = test_config.global_array[i];
                ASSERT_EQ(values.size(), ncols);
                for (int j = 0; j < ncols; ++j) {
                    EXPECT_FP_LE_WITH_EPS(values[j], icompute->array[i][j], epsilon);
                }
            }
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

        // local vector
        if (icompute->local_flag && icompute->size_local_cols == 0) {
            icompute->compute_local();
            int nrows = icompute->size_local_rows;
            ASSERT_EQ(test_config.local_vector.size(), nrows);
            for (int i = 0; i < nrows; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.local_vector[i], icompute->vector_local[i], epsilon);
        }

        // local array (need to sort expected/actual local arrays to match in KOKKOS)
        if (icompute->local_flag && icompute->size_local_cols > 0) {
            icompute->compute_local();
            int nrows = icompute->size_local_rows;
            int ncols = icompute->size_local_cols;
            ASSERT_EQ(static_cast<int>(test_config.local_array.size()), nrows);

            // Create sorted copies of both arrays
            std::vector<std::vector<double>> expected_sorted = test_config.local_array;
            std::vector<std::vector<double>> actual_sorted;

            // Copy actual array to vector for sorting
            for (int i = 0; i < nrows; ++i) {
                std::vector<double> row(ncols);
                for (int j = 0; j < ncols; ++j) {
                    row[j] = icompute->array_local[i][j];
                }
                actual_sorted.push_back(row);
            }

            // Sort both arrays
            std::sort(expected_sorted.begin(), expected_sorted.end());
            std::sort(actual_sorted.begin(), actual_sorted.end());

            // Compare sorted arrays
            for (int i = 0; i < nrows; ++i) {
                ASSERT_EQ(expected_sorted[i].size(), static_cast<size_t>(ncols));
                for (int j = 0; j < ncols; ++j) {
                    EXPECT_FP_LE_WITH_EPS(expected_sorted[i][j], actual_sorted[i][j], epsilon);
                }
            }
        }

        // pergrid vector/array
        if (icompute->pergrid_flag) {
            icompute->compute_pergrid();
            int dimension, nx, ny, nz, nxlo, nxhi, nylo, nyhi, nzlo, nzhi, ncol;
            int index = icompute->get_grid_by_name(test_config.pergrid_name, dimension);
            if( dimension == 2 ) {
                std::cerr << "ERROR: only Grid3D supported for unit testing.\n";
                exit(1);
            }
            Grid3d *grid3d = static_cast<Grid3d*>(icompute->get_grid_by_index(index));
            grid3d->get_size(nx, ny, nz);
            grid3d->get_bounds_owned(nxlo, nxhi, nylo, nyhi, nzlo, nzhi);
            icompute->get_griddata_by_name(index, test_config.pergrid_data, ncol);
            if (ncol == 0) {
                double ***vec3d = static_cast<double ***>(icompute->get_griddata_by_index(index));
                for (int iz = nzlo; iz <= nzhi; iz++)
                    for (int iy = nylo; iy <= nyhi; iy++)
                        for (int ix = nxlo; ix <= nxhi; ix++)
                          EXPECT_FP_LE_WITH_EPS(
                            test_config.pergrid_vector[(iz*ny + iy)*nx + ix],
                            vec3d[iz][iy][ix], epsilon );
            } else {
                double ****array3d = static_cast<double ****>(icompute->get_griddata_by_index(index));
                for (int iz = nzlo; iz <= nzhi; iz++)
                    for (int iy = nylo; iy <= nyhi; iy++)
                        for (int ix = nxlo; ix <= nxhi; ix++) {
                            for (int n = 0; n < ncol; n++)
                                EXPECT_FP_LE_WITH_EPS(
                                    test_config.pergrid_array[(iz*ny + iy)*nx + ix][n],
                                    array3d[iz][iy][ix][n], epsilon );
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
