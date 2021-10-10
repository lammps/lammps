/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// unit tests for pair styles intended for molecular systems

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
#include "fix.h"
#include "fmt/format.h"
#include "force.h"
#include "info.h"
#include "input.h"
#include "kspace.h"
#include "lammps.h"
#include "modify.h"
#include "pair.h"
#include "universe.h"
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

static void delete_file(const std::string &filename)
{
    remove(filename.c_str());
};

void cleanup_lammps(LAMMPS *lmp, const TestConfig &cfg)
{
    delete_file(cfg.basename + ".restart");
    delete lmp;
    lmp = nullptr;
}

LAMMPS *init_lammps(int argc, char **argv, const TestConfig &cfg, const bool use_respa = false)
{
    LAMMPS *lmp;

    lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);

    // check if prerequisite styles are available
    Info *info = new Info(lmp);
    int nfail  = 0;
    for (auto &prerequisite : cfg.prerequisites) {
        std::string style = prerequisite.second;

        // this is a test for fix styles, so if the suffixed
        // version is not available, there is no reason to test.
        if (prerequisite.first == "fix") {
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

    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line.c_str());
    };

    command("variable input_dir index " + INPUT_FOLDER);
    for (auto &pre_command : cfg.pre_commands)
        command(pre_command);

    std::string input_file = INPUT_FOLDER + PATH_SEP + cfg.input_file;
    lmp->input->file(input_file.c_str());

    if (use_respa) command("run_style respa 2 1 bond 1 pair 2");

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

    for (auto &post_command : cfg.post_commands)
        command(post_command);

    command("timestep 0.25");
    command("run 0 post no");
    command("thermo 2");
    command("run 4 post no start 0 stop 8");
    command("write_restart " + cfg.basename + ".restart");
    command("run 4 post no start 0 stop 8");
    return lmp;
}

void restart_lammps(LAMMPS *lmp, const TestConfig &cfg, bool use_rmass, bool use_respa)
{
    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line.c_str());
    };

    command("clear");
    command("read_restart " + cfg.basename + ".restart");

    if (use_rmass) {
        command("fix rmass all property/atom rmass ghost yes");
        for (int i = 0; i < lmp->atom->ntypes; ++i)
            command(fmt::format("set type {} mass {}", i + 1, lmp->atom->mass[i + 1]));
    }

    if (use_respa) command("run_style respa 2 1 bond 1 pair 2");

    for (auto &post_command : cfg.post_commands)
        command(post_command);

    command("thermo 2");
    command("run 4 post no start 0 stop 8");
}

// re-generate yaml file with current settings.

void generate_yaml_file(const char *outfile, const TestConfig &config)
{
    // initialize system geometry
    const char *args[] = {"FixIntegrate", "-log", "none", "-echo", "screen", "-nocite"};

    char **argv = (char **)args;
    int argc    = sizeof(args) / sizeof(char *);
    LAMMPS *lmp = init_lammps(argc, argv, config);
    if (!lmp) {
        std::cerr << "One or more prerequisite styles are not available "
                     "in this LAMMPS configuration:\n";
        for (auto prerequisite : config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        return;
    }

    const int natoms = lmp->atom->natoms;
    std::string block("");
    YamlWriter writer(outfile);

    // write yaml header
    write_yaml_header(&writer, &test_config, lmp->version);

    // natoms
    writer.emit("natoms", natoms);

    int ifix = lmp->modify->find_fix("test");
    if (ifix < 0) {
        std::cerr << "ERROR: no fix defined with fix ID 'test'\n";
        exit(1);
    } else {
        Fix *fix = lmp->modify->fix[ifix];

        // run_stress, if enabled
        if (fix->thermo_virial) {
            auto stress = fix->virial;
            block       = fmt::format("{:23.16e} {:23.16e} {:23.16e} "
                                "{:23.16e} {:23.16e} {:23.16e}",
                                stress[0], stress[1], stress[2], stress[3], stress[4], stress[5]);
            writer.emit_block("run_stress", block);
        }

        // global scalar
        if (fix->scalar_flag) {
            double value = fix->compute_scalar();
            writer.emit("global_scalar", value);
        }

        // global vector
        if (fix->vector_flag) {
            int num = fix->size_vector;
            block   = std::to_string(num);
            for (int i = 0; i < num; ++i)
                block += fmt::format(" {}", fix->compute_vector(i));
            writer.emit_block("global_vector", block);
        }
    }

    // run_pos
    block.clear();
    auto x = lmp->atom->x;
    for (int i = 1; i <= natoms; ++i) {
        const int j = lmp->atom->map(i);
        block += fmt::format("{:3} {:23.16e} {:23.16e} {:23.16e}\n", i, x[j][0], x[j][1], x[j][2]);
    }
    writer.emit_block("run_pos", block);

    // run_vel
    block.clear();
    auto v = lmp->atom->v;
    for (int i = 1; i <= natoms; ++i) {
        const int j = lmp->atom->map(i);
        block += fmt::format("{:3} {:23.16e} {:23.16e} {:23.16e}\n", i, v[j][0], v[j][1], v[j][2]);
    }
    writer.emit_block("run_vel", block);
    cleanup_lammps(lmp, config);
    return;
}

TEST(FixTimestep, plain)
{
    if (!LAMMPS::is_installed_pkg("MOLECULE")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
#if defined(USING_STATIC_LIBS)
    if (test_config.skip_tests.count("static")) GTEST_SKIP();
#endif

    const char *args[] = {"FixTimestep", "-log", "none", "-echo", "screen", "-nocite"};

    char **argv = (char **)args;
    int argc    = sizeof(args) / sizeof(char *);

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp        = init_lammps(argc, argv, test_config);
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

    auto tag = lmp->atom->tag;
    auto x   = lmp->atom->x;
    auto v   = lmp->atom->v;
    ErrorStats stats;
    stats.reset();
    const std::vector<coord_t> &x_ref = test_config.run_pos;
    ASSERT_EQ(nlocal + 1, x_ref.size());
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(x[i][0], x_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(x[i][1], x_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(x[i][2], x_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << "run_pos, normal run, verlet: " << stats << std::endl;

    const std::vector<coord_t> &v_ref = test_config.run_vel;
    stats.reset();
    ASSERT_EQ(nlocal + 1, v_ref.size());
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(v[i][0], v_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(v[i][1], v_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(v[i][2], v_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << "run_vel, normal run, verlet: " << stats << std::endl;

    int ifix = lmp->modify->find_fix("test");
    if (ifix < 0) {
        FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
    } else {
        Fix *fix = lmp->modify->fix[ifix];
        if (fix->thermo_virial) {
            stats.reset();
            auto stress = fix->virial;
            EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, epsilon);
            if (print_stats) std::cerr << "run_stress normal run, verlet: " << stats << std::endl;
        }

        stats.reset();
        // global scalar
        if (fix->scalar_flag) {
            double value = fix->compute_scalar();
            EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
        }

        // global vector
        if (fix->vector_flag) {
            int num = fix->size_vector;
            EXPECT_EQ(num, test_config.global_vector.size());

            for (int i = 0; i < num; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], fix->compute_vector(i),
                                      epsilon);
        }

        // check t_target for thermostats

        int dim     = -1;
        double *ptr = (double *)fix->extract("t_target", dim);
        if ((ptr != nullptr) && (dim == 0)) {
            int ivar = lmp->input->variable->find("t_target");
            if (ivar >= 0) {
                double t_ref    = atof(lmp->input->variable->retrieve("t_target"));
                double t_target = *ptr;
                EXPECT_FP_LE_WITH_EPS(t_target, t_ref, epsilon);
            }
        }
        if (print_stats && stats.has_data())
            std::cerr << "global_data, normal run, verlet: " << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    restart_lammps(lmp, test_config, false, false);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    tag = lmp->atom->tag;
    x   = lmp->atom->x;
    v   = lmp->atom->v;
    stats.reset();
    ASSERT_EQ(nlocal + 1, x_ref.size());
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(x[i][0], x_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(x[i][1], x_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(x[i][2], x_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << "run_pos, restart, verlet: " << stats << std::endl;

    stats.reset();
    ASSERT_EQ(nlocal + 1, v_ref.size());
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(v[i][0], v_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(v[i][1], v_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(v[i][2], v_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << "run_vel, restart, verlet: " << stats << std::endl;

    ifix = lmp->modify->find_fix("test");
    if (ifix < 0) {
        FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
    } else {
        Fix *fix = lmp->modify->fix[ifix];
        if (fix->thermo_virial) {
            stats.reset();
            auto stress = fix->virial;
            EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, epsilon);
            if (print_stats) std::cerr << "run_stress restart, verlet: " << stats << std::endl;
        }

        stats.reset();

        // global scalar
        if (fix->scalar_flag) {
            double value = fix->compute_scalar();
            EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
        }

        // global vector
        if (fix->vector_flag) {
            int num = fix->size_vector;
            EXPECT_EQ(num, test_config.global_vector.size());

            for (int i = 0; i < num; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], fix->compute_vector(i),
                                      epsilon);
        }
        if (print_stats && stats.has_data())
            std::cerr << "global_data, restart, verlet: " << stats << std::endl;
    }

    if (lmp->atom->rmass == nullptr) {
        if (!verbose) ::testing::internal::CaptureStdout();
        restart_lammps(lmp, test_config, true, false);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        x   = lmp->atom->x;
        tag = lmp->atom->tag;
        stats.reset();
        ASSERT_EQ(nlocal + 1, x_ref.size());
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(x[i][0], x_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(x[i][1], x_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(x[i][2], x_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "run_pos, rmass, verlet: " << stats << std::endl;

        v = lmp->atom->v;
        stats.reset();
        ASSERT_EQ(nlocal + 1, v_ref.size());
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(v[i][0], v_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(v[i][1], v_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(v[i][2], v_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "run_vel, rmass, verlet: " << stats << std::endl;

        ifix = lmp->modify->find_fix("test");
        if (ifix < 0) {
            FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
        } else {
            Fix *fix = lmp->modify->fix[ifix];
            if (fix->thermo_virial) {
                stats.reset();
                auto stress = fix->virial;
                EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, epsilon);
                if (print_stats) std::cerr << "run_stress rmass, verlet: " << stats << std::endl;
            }

            stats.reset();

            // global scalar
            if (fix->scalar_flag) {
                double value = fix->compute_scalar();
                EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
            }

            // global vector
            if (fix->vector_flag) {
                int num = fix->size_vector;
                EXPECT_EQ(num, test_config.global_vector.size());

                for (int i = 0; i < num; ++i)
                    EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], fix->compute_vector(i),
                                          epsilon);
            }
            if (print_stats && stats.has_data())
                std::cerr << "global_data, rmass, verlet: " << stats << std::endl;
        }
    }

    // rigid fixes need work to test properly with r-RESPA.
    // fix nve/limit cannot work with r-RESPA
    ifix = lmp->modify->find_fix("test");
    if (!utils::strmatch(lmp->modify->fix[ifix]->style, "^rigid") &&
        !utils::strmatch(lmp->modify->fix[ifix]->style, "^nve/limit")) {

        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        ::testing::internal::CaptureStdout();
        lmp    = init_lammps(argc, argv, test_config, true);
        output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;

        // lower required precision by two orders of magnitude to accommodate respa
        epsilon *= 100.0;

        tag = lmp->atom->tag;
        x   = lmp->atom->x;
        v   = lmp->atom->v;
        stats.reset();
        ASSERT_EQ(nlocal + 1, x_ref.size());
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(x[i][0], x_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(x[i][1], x_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(x[i][2], x_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "run_pos, normal run, respa: " << stats << std::endl;

        stats.reset();
        ASSERT_EQ(nlocal + 1, v_ref.size());
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(v[i][0], v_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(v[i][1], v_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(v[i][2], v_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "run_vel, normal run, respa: " << stats << std::endl;

        ifix = lmp->modify->find_fix("test");
        if (ifix < 0) {
            FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
        } else {
            Fix *fix = lmp->modify->fix[ifix];
            if (fix->thermo_virial) {
                stats.reset();
                auto stress = fix->virial;
                EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, 1000 * epsilon);
                if (print_stats)
                    std::cerr << "run_stress normal run, respa: " << stats << std::endl;
            }

            stats.reset();

            // global scalar
            if (fix->scalar_flag) {
                double value = fix->compute_scalar();
                EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, 10 * epsilon);
            }

            // global vector
            if (fix->vector_flag) {
                int num = fix->size_vector;
                EXPECT_EQ(num, test_config.global_vector.size());

                for (int i = 0; i < num; ++i)
                    EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], fix->compute_vector(i),
                                          10 * epsilon);
            }
            if (print_stats && stats.has_data())
                std::cerr << "global_data, normal run, respa: " << stats << std::endl;
        }

        if (!verbose) ::testing::internal::CaptureStdout();
        restart_lammps(lmp, test_config, false, true);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        tag = lmp->atom->tag;
        x   = lmp->atom->x;
        v   = lmp->atom->v;
        stats.reset();
        ASSERT_EQ(nlocal + 1, x_ref.size());
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(x[i][0], x_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(x[i][1], x_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(x[i][2], x_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "run_pos, restart, respa: " << stats << std::endl;

        stats.reset();
        ASSERT_EQ(nlocal + 1, v_ref.size());
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(v[i][0], v_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(v[i][1], v_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(v[i][2], v_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "run_vel, restart, respa: " << stats << std::endl;

        ifix = lmp->modify->find_fix("test");
        if (ifix < 0) {
            FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
        } else {
            Fix *fix = lmp->modify->fix[ifix];
            if (fix->thermo_virial) {
                stats.reset();
                auto stress = fix->virial;
                EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, 1000 * epsilon);
                if (print_stats) std::cerr << "run_stress restart, respa: " << stats << std::endl;
            }

            stats.reset();

            // global scalar
            if (fix->scalar_flag) {
                double value = fix->compute_scalar();
                EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, 10 * epsilon);
            }

            // global vector
            if (fix->vector_flag) {
                int num = fix->size_vector;
                EXPECT_EQ(num, test_config.global_vector.size());

                for (int i = 0; i < num; ++i)
                    EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], fix->compute_vector(i),
                                          10 * epsilon);
            }
            if (print_stats && stats.has_data())
                std::cerr << "global_data, restart, respa: " << stats << std::endl;
        }

        if (lmp->atom->rmass == nullptr) {
            if (!verbose) ::testing::internal::CaptureStdout();
            restart_lammps(lmp, test_config, true, true);
            if (!verbose) ::testing::internal::GetCapturedStdout();

            x   = lmp->atom->x;
            tag = lmp->atom->tag;
            stats.reset();
            ASSERT_EQ(nlocal + 1, x_ref.size());
            for (int i = 0; i < nlocal; ++i) {
                EXPECT_FP_LE_WITH_EPS(x[i][0], x_ref[tag[i]].x, epsilon);
                EXPECT_FP_LE_WITH_EPS(x[i][1], x_ref[tag[i]].y, epsilon);
                EXPECT_FP_LE_WITH_EPS(x[i][2], x_ref[tag[i]].z, epsilon);
            }
            if (print_stats) std::cerr << "run_pos, rmass, respa: " << stats << std::endl;

            v = lmp->atom->v;
            stats.reset();
            ASSERT_EQ(nlocal + 1, v_ref.size());
            for (int i = 0; i < nlocal; ++i) {
                EXPECT_FP_LE_WITH_EPS(v[i][0], v_ref[tag[i]].x, epsilon);
                EXPECT_FP_LE_WITH_EPS(v[i][1], v_ref[tag[i]].y, epsilon);
                EXPECT_FP_LE_WITH_EPS(v[i][2], v_ref[tag[i]].z, epsilon);
            }
            if (print_stats) std::cerr << "run_vel, rmass, respa: " << stats << std::endl;

            ifix = lmp->modify->find_fix("test");
            if (ifix < 0) {
                FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
            } else {
                Fix *fix = lmp->modify->fix[ifix];
                if (fix->thermo_virial) {
                    stats.reset();
                    auto stress = fix->virial;
                    EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, 1000 * epsilon);
                    EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, 1000 * epsilon);
                    EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, 1000 * epsilon);
                    EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, 1000 * epsilon);
                    EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, 1000 * epsilon);
                    EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, 1000 * epsilon);
                    if (print_stats) std::cerr << "run_stress rmass, respa: " << stats << std::endl;
                }

                stats.reset();

                // global scalar
                if (fix->scalar_flag) {
                    double value = fix->compute_scalar();
                    EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, 10 * epsilon);
                }

                // global vector
                if (fix->vector_flag) {
                    int num = fix->size_vector;
                    EXPECT_EQ(num, test_config.global_vector.size());

                    for (int i = 0; i < num; ++i)
                        EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], fix->compute_vector(i),
                                              10 * epsilon);
                }
                if (print_stats && stats.has_data())
                    std::cerr << "global_data, rmass, respa: " << stats << std::endl;
            }
        }
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

TEST(FixTimestep, omp)
{
    if (!LAMMPS::is_installed_pkg("OPENMP")) GTEST_SKIP();
    if (!LAMMPS::is_installed_pkg("MOLECULE")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
#if defined(USING_STATIC_LIBS)
    if (test_config.skip_tests.count("static")) GTEST_SKIP();
#endif

    const char *args[] = {"FixTimestep", "-log", "none", "-echo", "screen", "-nocite",
                          "-pk",         "omp",  "4",    "-sf",   "omp"};

    char **argv = (char **)args;
    int argc    = sizeof(args) / sizeof(char *);

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp        = init_lammps(argc, argv, test_config);
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

    auto tag = lmp->atom->tag;
    auto x   = lmp->atom->x;
    auto v   = lmp->atom->v;
    ErrorStats stats;
    stats.reset();
    const std::vector<coord_t> &x_ref = test_config.run_pos;
    ASSERT_EQ(nlocal + 1, x_ref.size());
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(x[i][0], x_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(x[i][1], x_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(x[i][2], x_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << "run_pos, normal run, verlet: " << stats << std::endl;

    const std::vector<coord_t> &v_ref = test_config.run_vel;
    stats.reset();
    ASSERT_EQ(nlocal + 1, v_ref.size());
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(v[i][0], v_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(v[i][1], v_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(v[i][2], v_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << "run_vel, normal run, verlet: " << stats << std::endl;

    int ifix = lmp->modify->find_fix("test");
    if (ifix < 0) {
        FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
    } else {
        Fix *fix = lmp->modify->fix[ifix];
        if (fix->thermo_virial) {
            stats.reset();
            auto stress = fix->virial;
            EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, epsilon);
            if (print_stats) std::cerr << "run_stress normal run, verlet: " << stats << std::endl;
        }

        stats.reset();
        // global scalar
        if (fix->scalar_flag) {
            double value = fix->compute_scalar();
            EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
        }

        // global vector
        if (fix->vector_flag) {
            int num = fix->size_vector;
            EXPECT_EQ(num, test_config.global_vector.size());

            for (int i = 0; i < num; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], fix->compute_vector(i),
                                      epsilon);
        }

        // check t_target for thermostats

        int dim     = -1;
        double *ptr = (double *)fix->extract("t_target", dim);
        if ((ptr != nullptr) && (dim == 0)) {
            int ivar = lmp->input->variable->find("t_target");
            if (ivar >= 0) {
                double t_ref    = atof(lmp->input->variable->retrieve("t_target"));
                double t_target = *ptr;
                EXPECT_FP_LE_WITH_EPS(t_target, t_ref, epsilon);
            }
        }
        if (print_stats && stats.has_data())
            std::cerr << "global_data, normal run, verlet: " << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    restart_lammps(lmp, test_config, false, false);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    tag = lmp->atom->tag;
    x   = lmp->atom->x;
    v   = lmp->atom->v;
    stats.reset();
    ASSERT_EQ(nlocal + 1, x_ref.size());
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(x[i][0], x_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(x[i][1], x_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(x[i][2], x_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << "run_pos, restart, verlet: " << stats << std::endl;

    stats.reset();
    ASSERT_EQ(nlocal + 1, v_ref.size());
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(v[i][0], v_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(v[i][1], v_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(v[i][2], v_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << "run_vel, restart, verlet: " << stats << std::endl;

    ifix = lmp->modify->find_fix("test");
    if (ifix < 0) {
        FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
    } else {
        Fix *fix = lmp->modify->fix[ifix];
        if (fix->thermo_virial) {
            stats.reset();
            auto stress = fix->virial;
            EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, epsilon);
            EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, epsilon);
            if (print_stats) std::cerr << "run_stress restart, verlet: " << stats << std::endl;
        }

        stats.reset();

        // global scalar
        if (fix->scalar_flag) {
            double value = fix->compute_scalar();
            EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
        }

        // global vector
        if (fix->vector_flag) {
            int num = fix->size_vector;
            EXPECT_EQ(num, test_config.global_vector.size());

            for (int i = 0; i < num; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], fix->compute_vector(i),
                                      epsilon);
        }
        if (print_stats && stats.has_data())
            std::cerr << "global_data, restart, verlet: " << stats << std::endl;
    }

    if (lmp->atom->rmass == nullptr) {
        if (!verbose) ::testing::internal::CaptureStdout();
        restart_lammps(lmp, test_config, true, false);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        x   = lmp->atom->x;
        tag = lmp->atom->tag;
        stats.reset();
        ASSERT_EQ(nlocal + 1, x_ref.size());
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(x[i][0], x_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(x[i][1], x_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(x[i][2], x_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "run_pos, rmass, verlet: " << stats << std::endl;

        v = lmp->atom->v;
        stats.reset();
        ASSERT_EQ(nlocal + 1, v_ref.size());
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(v[i][0], v_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(v[i][1], v_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(v[i][2], v_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "run_vel, rmass, verlet: " << stats << std::endl;

        ifix = lmp->modify->find_fix("test");
        if (ifix < 0) {
            FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
        } else {
            Fix *fix = lmp->modify->fix[ifix];
            if (fix->thermo_virial) {
                stats.reset();
                auto stress = fix->virial;
                EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, epsilon);
                if (print_stats) std::cerr << "run_stress rmass, verlet: " << stats << std::endl;
            }

            stats.reset();

            // global scalar
            if (fix->scalar_flag) {
                double value = fix->compute_scalar();
                EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
            }

            // global vector
            if (fix->vector_flag) {
                int num = fix->size_vector;
                EXPECT_EQ(num, test_config.global_vector.size());

                for (int i = 0; i < num; ++i)
                    EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], fix->compute_vector(i),
                                          epsilon);
            }
            if (print_stats && stats.has_data())
                std::cerr << "global_data, rmass, verlet: " << stats << std::endl;
        }
    }

    // rigid fixes need work to test properly with r-RESPA,
    // also, torque is not supported by respa/omp
    ifix = lmp->modify->find_fix("test");
    if (!utils::strmatch(lmp->modify->fix[ifix]->style, "^rigid") && !lmp->atom->torque) {

        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        ::testing::internal::CaptureStdout();
        lmp    = init_lammps(argc, argv, test_config, true);
        output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;

        // lower required precision by two orders of magnitude to accommodate respa
        epsilon *= 100.0;

        tag = lmp->atom->tag;
        x   = lmp->atom->x;
        v   = lmp->atom->v;
        stats.reset();
        ASSERT_EQ(nlocal + 1, x_ref.size());
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(x[i][0], x_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(x[i][1], x_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(x[i][2], x_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "run_pos, normal run, respa: " << stats << std::endl;
        printf("x1\n");

        stats.reset();
        ASSERT_EQ(nlocal + 1, v_ref.size());
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(v[i][0], v_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(v[i][1], v_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(v[i][2], v_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "run_vel, normal run, respa: " << stats << std::endl;

        ifix = lmp->modify->find_fix("test");
        if (ifix < 0) {
            FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
        } else {
            Fix *fix = lmp->modify->fix[ifix];
            if (fix->thermo_virial) {
                stats.reset();
                auto stress = fix->virial;
                EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, 1000 * epsilon);
                if (print_stats)
                    std::cerr << "run_stress normal run, respa: " << stats << std::endl;
            }

            stats.reset();

            // global scalar
            if (fix->scalar_flag) {
                double value = fix->compute_scalar();
                EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, 10 * epsilon);
            }

            // global vector
            if (fix->vector_flag) {
                int num = fix->size_vector;
                EXPECT_EQ(num, test_config.global_vector.size());

                for (int i = 0; i < num; ++i)
                    EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], fix->compute_vector(i),
                                          10 * epsilon);
            }
            if (print_stats && stats.has_data())
                std::cerr << "global_data, normal run, respa: " << stats << std::endl;
        }

        if (!verbose) ::testing::internal::CaptureStdout();
        restart_lammps(lmp, test_config, false, true);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        tag = lmp->atom->tag;
        x   = lmp->atom->x;
        v   = lmp->atom->v;
        stats.reset();
        ASSERT_EQ(nlocal + 1, x_ref.size());
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(x[i][0], x_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(x[i][1], x_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(x[i][2], x_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "run_pos, restart, respa: " << stats << std::endl;

        stats.reset();
        ASSERT_EQ(nlocal + 1, v_ref.size());
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(v[i][0], v_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(v[i][1], v_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(v[i][2], v_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "run_vel, restart, respa: " << stats << std::endl;

        ifix = lmp->modify->find_fix("test");
        if (ifix < 0) {
            FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
        } else {
            Fix *fix = lmp->modify->fix[ifix];
            if (fix->thermo_virial) {
                stats.reset();
                auto stress = fix->virial;
                EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, 1000 * epsilon);
                EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, 1000 * epsilon);
                if (print_stats) std::cerr << "run_stress restart, respa: " << stats << std::endl;
            }

            stats.reset();

            // global scalar
            if (fix->scalar_flag) {
                double value = fix->compute_scalar();
                EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, 10 * epsilon);
            }

            // global vector
            if (fix->vector_flag) {
                int num = fix->size_vector;
                EXPECT_EQ(num, test_config.global_vector.size());

                for (int i = 0; i < num; ++i)
                    EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], fix->compute_vector(i),
                                          10 * epsilon);
            }
            if (print_stats && stats.has_data())
                std::cerr << "global_data, restart, respa: " << stats << std::endl;
        }

        if (lmp->atom->rmass == nullptr) {
            if (!verbose) ::testing::internal::CaptureStdout();
            restart_lammps(lmp, test_config, true, true);
            if (!verbose) ::testing::internal::GetCapturedStdout();

            x   = lmp->atom->x;
            tag = lmp->atom->tag;
            stats.reset();
            ASSERT_EQ(nlocal + 1, x_ref.size());
            for (int i = 0; i < nlocal; ++i) {
                EXPECT_FP_LE_WITH_EPS(x[i][0], x_ref[tag[i]].x, epsilon);
                EXPECT_FP_LE_WITH_EPS(x[i][1], x_ref[tag[i]].y, epsilon);
                EXPECT_FP_LE_WITH_EPS(x[i][2], x_ref[tag[i]].z, epsilon);
            }
            if (print_stats) std::cerr << "run_pos, rmass, respa: " << stats << std::endl;

            v = lmp->atom->v;
            stats.reset();
            ASSERT_EQ(nlocal + 1, v_ref.size());
            for (int i = 0; i < nlocal; ++i) {
                EXPECT_FP_LE_WITH_EPS(v[i][0], v_ref[tag[i]].x, epsilon);
                EXPECT_FP_LE_WITH_EPS(v[i][1], v_ref[tag[i]].y, epsilon);
                EXPECT_FP_LE_WITH_EPS(v[i][2], v_ref[tag[i]].z, epsilon);
            }
            if (print_stats) std::cerr << "run_vel, rmass, respa: " << stats << std::endl;

            ifix = lmp->modify->find_fix("test");
            if (ifix < 0) {
                FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
            } else {
                Fix *fix = lmp->modify->fix[ifix];
                if (fix->thermo_virial) {
                    stats.reset();
                    auto stress = fix->virial;
                    EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, 1000 * epsilon);
                    EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, 1000 * epsilon);
                    EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, 1000 * epsilon);
                    EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, 1000 * epsilon);
                    EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, 1000 * epsilon);
                    EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, 1000 * epsilon);
                    if (print_stats) std::cerr << "run_stress rmass, respa: " << stats << std::endl;
                }

                stats.reset();

                // global scalar
                if (fix->scalar_flag) {
                    double value = fix->compute_scalar();
                    EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, 10 * epsilon);
                }

                // global vector
                if (fix->vector_flag) {
                    int num = fix->size_vector;
                    EXPECT_EQ(num, test_config.global_vector.size());

                    for (int i = 0; i < num; ++i)
                        EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], fix->compute_vector(i),
                                              10 * epsilon);
                }
                if (print_stats && stats.has_data())
                    std::cerr << "global_data, rmass, respa: " << stats << std::endl;
            }
        }
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};
