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
#include "fix.h"
#include "force.h"
#include "info.h"
#include "input.h"
#include "kspace.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <exception>
#include <iostream>
#include <set>
#include <utility>

using ::testing::HasSubstr;
using ::testing::StartsWith;

using namespace LAMMPS_NS;

void cleanup_lammps(LAMMPS *&lmp, const TestConfig &cfg)
{
    platform::unlink(cfg.basename + ".restart");
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
    for (const auto &prerequisite : cfg.prerequisites) {
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
        lmp->input->one(line);
    };

    command("variable input_dir index " + INPUT_FOLDER);
    for (const auto &pre_command : cfg.pre_commands)
        command(pre_command);

    std::string input_file = platform::path_join(INPUT_FOLDER, cfg.input_file);
    lmp->input->file(input_file.c_str());

    if (use_respa) command("run_style respa 2 1 bond 1 pair 2");

    // set up molecular system force field and groups from the coeffs file
    // indicated by the YAML file (the input template only defines the geometry)

    if (cfg.input_coeffs.empty()) {
        std::cerr << "ERROR: no 'input_coeffs' file given in the YAML file\n";
        cleanup_lammps(lmp, cfg);
        return nullptr;
    }
    std::string coeffs_file = platform::path_join(INPUT_FOLDER, cfg.input_coeffs);
    lmp->input->file(coeffs_file.c_str());

    for (const auto &post_command : cfg.post_commands)
        command(post_command);

    // the default timestep of 0.25 assumes the (real units) molecular input templates;
    // systems needing a different timestep (e.g. spin dynamics in metal units) set the
    // "timestep" keyword in their yaml file.
    command(fmt::format("timestep {}", (cfg.timestep > 0.0) ? cfg.timestep : 0.25));
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
        lmp->input->one(line);
    };

    command("clear");
    command("read_restart " + cfg.basename + ".restart");

    if (use_rmass) {
        command("fix rmass all property/atom rmass ghost yes");
        for (int i = 0; i < lmp->atom->ntypes; ++i)
            command(fmt::format("set type {} mass {}", i + 1, lmp->atom->mass[i + 1]));
    }

    if (use_respa) command("run_style respa 2 1 bond 1 pair 2");

    for (const auto &post_command : cfg.post_commands)
        command(post_command);

    auto *ifix = lmp->modify->get_fix_by_id("test");
    // styles tagged "no_reset_dt" reject a timestep change (Fix::reset_dt() raises
    // an error, e.g. fix move), so do not exercise it for them
    if (ifix && !test_config.has_tag("no_reset_dt")) {
        // must be set to trigger calling Fix::reset_dt() with timestep
        lmp->update->first_update = 1;
        // test validity of Fix::reset_dt(). With run_style respa there may be segfaults
        command(fmt::format("timestep {}", (cfg.timestep > 0.0) ? cfg.timestep : 0.25));
    }
    command("thermo 2");
    command("run 4 post no start 0 stop 8");
}

// re-generate yaml file with current settings.

void generate_yaml_file(const char *outfile, const TestConfig &config)
{
    // initialize system geometry
    LAMMPS::argv args = {"FixIntegrate", "-log", "none", "-echo", "screen", "-nocite"};
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

    auto *ifix = lmp->modify->get_fix_by_id("test");
    if (!ifix) {
        std::cerr << "ERROR: no fix defined with fix ID 'test'\n";
        exit(1);
    } else {
        // run_stress, if enabled
        if (ifix->thermo_virial) {
            auto *stress = ifix->virial;
            // avoid false positives on tiny stresses. force to zero instead.
            for (int i = 0; i < 6; ++i)
                if (fabs(stress[i]) < 1.0e-13) stress[i] = 0.0;
            block = fmt::format("{:23.16e} {:23.16e} {:23.16e} {:23.16e} {:23.16e} {:23.16e}",
                                stress[0], stress[1], stress[2], stress[3], stress[4], stress[5]);
            writer.emit_block("run_stress", block);
        }

        // global scalar
        if (ifix->scalar_flag) {
            double value = ifix->compute_scalar();
            // avoid false positives on tiny values. force to zero instead.
            if (fabs(value) < 1.0e-13) value = 0.0;
            writer.emit("global_scalar", value);
        }

        // global vector
        if (ifix->vector_flag) {
            int num = ifix->size_vector;
            block   = std::to_string(num);
            double value;
            for (int i = 0; i < num; ++i) {
                // avoid false positives on tiny values. force to zero instead.
                value = ifix->compute_vector(i);
                if (fabs(value) < 1.0e-13) value = 0.0;
                block += fmt::format(" {:23.16e}", value);
            }
            writer.emit_block("global_vector", block);
        }
    }

    // run_pos
    block.clear();
    auto *x = lmp->atom->x;
    for (int i = 1; i <= natoms; ++i) {
        const int j = lmp->atom->map(i);
        block += fmt::format("{:3} {:23.16e} {:23.16e} {:23.16e}\n", i, x[j][0], x[j][1], x[j][2]);
    }
    writer.emit_block("run_pos", block);

    // run_vel
    block.clear();
    auto *v = lmp->atom->v;
    for (int i = 1; i <= natoms; ++i) {
        const int j = lmp->atom->map(i);
        block += fmt::format("{:3} {:23.16e} {:23.16e} {:23.16e}\n", i, v[j][0], v[j][1], v[j][2]);
    }
    writer.emit_block("run_vel", block);

    // run_torque

    if (lmp->atom->torque_flag) {
        block.clear();
        auto *t = lmp->atom->torque;
        for (int i = 1; i <= natoms; ++i) {
            const int j = lmp->atom->map(i);
            block +=
                fmt::format("{:3} {:23.16e} {:23.16e} {:23.16e}\n", i, t[j][0], t[j][1], t[j][2]);
        }
        writer.emit_block("run_torque", block);
    }

    // run_spin and run_mag_forces (only for atom_style spin)

    if (lmp->atom->sp_flag) {
        block.clear();
        auto *sp = lmp->atom->sp;
        for (int i = 1; i <= natoms; ++i) {
            const int j = lmp->atom->map(i);
            block += fmt::format("{:3} {:23.16e} {:23.16e} {:23.16e} {:23.16e}\n", i, sp[j][0],
                                 sp[j][1], sp[j][2], sp[j][3]);
        }
        writer.emit_block("run_spin", block);

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

TEST(FixTimestep, plain)
{
    // the "atom <style>" entry in the yaml prerequisites covers required packages
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
#if defined(USING_STATIC_LIBS)
    if (test_config.skip_tests.count("static")) GTEST_SKIP();
#endif

    LAMMPS::argv args = {"FixTimestep", "-log", "none", "-echo", "screen", "-nocite"};

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

    EXPECT_POSITIONS("run_pos (normal run, verlet)", lmp->atom, test_config.run_pos, epsilon);
    EXPECT_VELOCITIES("run_vel (normal run, verlet)", lmp->atom, test_config.run_vel, epsilon);
    if (lmp->atom->torque_flag)
        EXPECT_TORQUES("run_torques (normal run, verlet)", lmp->atom, test_config.run_torque,
                       epsilon);
    EXPECT_SPINS("run_spin (normal run, verlet)", lmp->atom, test_config.run_spin, epsilon);
    EXPECT_MAG_FORCES("run_mag_forces (normal run, verlet)", lmp->atom,
                      test_config.run_mag_forces, epsilon);

    auto *ifix = lmp->modify->get_fix_by_id("test");
    if (!ifix) {
        FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
    } else {
        if (ifix->thermo_virial) {
            EXPECT_STRESS("run_stress (normal run, verlet)", ifix->virial, test_config.run_stress,
                          epsilon);
        }

        stats.reset();
        // global scalar
        if (ifix->scalar_flag) {
            double value = ifix->compute_scalar();
            EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
        }

        // global vector
        if (ifix->vector_flag) {
            int num = ifix->size_vector;
            EXPECT_EQ(num, test_config.global_vector.size());

            for (int i = 0; i < num; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
                                      epsilon);
        }

        // check t_target for thermostats

        int dim   = -1;
        auto *ptr = (double *)ifix->extract("t_target", dim);
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

    EXPECT_POSITIONS("run_pos (restart, verlet)", lmp->atom, test_config.run_pos, epsilon);
    EXPECT_VELOCITIES("run_vel (restart, verlet)", lmp->atom, test_config.run_vel, epsilon);
    if (lmp->atom->torque_flag)
        EXPECT_TORQUES("run_torque (restart, verlet)", lmp->atom, test_config.run_torque, epsilon);
    EXPECT_SPINS("run_spin (restart, verlet)", lmp->atom, test_config.run_spin, epsilon);
    EXPECT_MAG_FORCES("run_mag_forces (restart, verlet)", lmp->atom, test_config.run_mag_forces,
                      epsilon);

    ifix = lmp->modify->get_fix_by_id("test");
    if (!ifix) {
        FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
    } else {
        if (ifix->thermo_virial) {
            EXPECT_STRESS("run_stress (restart, verlet)", ifix->virial, test_config.run_stress,
                          epsilon);
        }

        stats.reset();

        // global scalar
        if (ifix->scalar_flag) {
            double value = ifix->compute_scalar();
            EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
        }

        // global vector
        if (ifix->vector_flag) {
            int num = ifix->size_vector;
            EXPECT_EQ(num, test_config.global_vector.size());

            for (int i = 0; i < num; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
                                      epsilon);
        }
        if (print_stats && stats.has_data())
            std::cerr << "global_data, restart, verlet: " << stats << std::endl;
    }

    if (lmp->atom->rmass == nullptr) {
        if (!verbose) ::testing::internal::CaptureStdout();
        restart_lammps(lmp, test_config, true, false);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        EXPECT_POSITIONS("run_pos (rmass, verlet)", lmp->atom, test_config.run_pos, epsilon);
        EXPECT_VELOCITIES("run_vel (rmass, verlet)", lmp->atom, test_config.run_vel, epsilon);

        ifix = lmp->modify->get_fix_by_id("test");
        if (!ifix) {
            FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
        } else {
            if (ifix->thermo_virial) {
                EXPECT_STRESS("run_stress (rmass, verlet)", ifix->virial, test_config.run_stress,
                              epsilon);
            }

            stats.reset();

            // global scalar
            if (ifix->scalar_flag) {
                double value = ifix->compute_scalar();
                EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
            }

            // global vector
            if (ifix->vector_flag) {
                int num = ifix->size_vector;
                EXPECT_EQ(num, test_config.global_vector.size());

                for (int i = 0; i < num; ++i)
                    EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
                                          epsilon);
            }
            if (print_stats && stats.has_data())
                std::cerr << "global_data, rmass, verlet: " << stats << std::endl;
        }
    }

    // styles tagged "no_respa" are not exercised under r-RESPA: rigid fixes need
    // work to test properly with r-RESPA, fix nve/limit and fix recenter do not
    // support it, and stochastic integrators/barostats (brownian, gjf,
    // press/langevin) draw their random numbers differently under r-RESPA so the
    // trajectories cannot match.  Adding the tag to a YAML file is all that is
    // needed for a future case; no change to this driver is required.
    ifix = lmp->modify->get_fix_by_id("test");
    if (ifix && !test_config.has_tag("no_respa")) {
        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();

        ::testing::internal::CaptureStdout();
        try {
            lmp = init_lammps(args, test_config, true);
        } catch (std::exception &e) {
            output = ::testing::internal::GetCapturedStdout();
            if (verbose) std::cout << output;
            FAIL() << e.what();
        }
        output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;

        // lower required precision by two orders of magnitude to accommodate respa
        epsilon *= 100.0;

        EXPECT_POSITIONS("run_pos (normal run, respa)", lmp->atom, test_config.run_pos, epsilon);
        EXPECT_VELOCITIES("run_vel (normal run, respa)", lmp->atom, test_config.run_vel, epsilon);

        ifix = lmp->modify->get_fix_by_id("test");
        if (!ifix) {
            FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
        } else {
            if (ifix->thermo_virial) {
                EXPECT_STRESS("run_stress (normal run, respa)", ifix->virial,
                              test_config.run_stress, 1000 * epsilon);
            }

            stats.reset();

            // global scalar
            if (ifix->scalar_flag) {
                double value = ifix->compute_scalar();
                EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, 10 * epsilon);
            }

            // global vector
            if (ifix->vector_flag) {
                int num = ifix->size_vector;
                EXPECT_EQ(num, test_config.global_vector.size());

                for (int i = 0; i < num; ++i)
                    EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
                                          10 * epsilon);
            }
            if (print_stats && stats.has_data())
                std::cerr << "global_data, normal run, respa: " << stats << std::endl;
        }

        if (!verbose) ::testing::internal::CaptureStdout();
        restart_lammps(lmp, test_config, false, true);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        EXPECT_POSITIONS("run_pos (restart, respa)", lmp->atom, test_config.run_pos, epsilon);
        EXPECT_VELOCITIES("run_vel (restart, respa)", lmp->atom, test_config.run_vel, epsilon);

        ifix = lmp->modify->get_fix_by_id("test");
        if (!ifix) {
            FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
        } else {
            if (ifix->thermo_virial) {
                EXPECT_STRESS("run_stress (restart, respa)", ifix->virial, test_config.run_stress,
                              1000 * epsilon);
            }

            stats.reset();

            // global scalar
            if (ifix->scalar_flag) {
                double value = ifix->compute_scalar();
                EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, 10 * epsilon);
            }

            // global vector
            if (ifix->vector_flag) {
                int num = ifix->size_vector;
                EXPECT_EQ(num, test_config.global_vector.size());

                for (int i = 0; i < num; ++i)
                    EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
                                          10 * epsilon);
            }
            if (print_stats && stats.has_data())
                std::cerr << "global_data, restart, respa: " << stats << std::endl;
        }

        if (lmp->atom->rmass == nullptr) {
            if (!verbose) ::testing::internal::CaptureStdout();
            restart_lammps(lmp, test_config, true, true);
            if (!verbose) ::testing::internal::GetCapturedStdout();

            EXPECT_POSITIONS("run_pos (rmass, respa)", lmp->atom, test_config.run_pos, epsilon);
            EXPECT_VELOCITIES("run_vel (rmass, respa)", lmp->atom, test_config.run_vel, epsilon);

            ifix = lmp->modify->get_fix_by_id("test");
            if (!ifix) {
                FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
            } else {
                if (ifix->thermo_virial) {
                    EXPECT_STRESS("run_stress (rmass, respa)", ifix->virial, test_config.run_stress,
                                  1000 * epsilon);
                }

                stats.reset();

                // global scalar
                if (ifix->scalar_flag) {
                    double value = ifix->compute_scalar();
                    EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, 10 * epsilon);
                }

                // global vector
                if (ifix->vector_flag) {
                    int num = ifix->size_vector;
                    EXPECT_EQ(num, test_config.global_vector.size());

                    for (int i = 0; i < num; ++i)
                        EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
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
    if (!Info::has_package("OPENMP")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();
#if defined(USING_STATIC_LIBS)
    if (test_config.skip_tests.count("static")) GTEST_SKIP();
#endif

    LAMMPS::argv args = {"FixTimestep", "-log", "none", "-echo", "screen", "-nocite",
                         "-pk",         "omp",  "4",    "-sf",   "omp"};

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
        std::cerr << "One or more prerequisite styles with /omp suffix are not available "
                     "in this LAMMPS configuration:\n";
        for (auto &prerequisite : test_config.prerequisites) {
            if (prerequisite.first == "atom") {
                std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
            } else {
                std::cerr << prerequisite.first << "_style " << prerequisite.second << "/omp\n";
            }
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

    EXPECT_POSITIONS("run_pos (normal run, verlet)", lmp->atom, test_config.run_pos, epsilon);
    EXPECT_VELOCITIES("run_vel (normal run, verlet)", lmp->atom, test_config.run_vel, epsilon);

    auto *ifix = lmp->modify->get_fix_by_id("test");
    if (!ifix) {
        FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
    } else {
        if (ifix->thermo_virial) {
            EXPECT_STRESS("run_stress (normal run, verlet)", ifix->virial, test_config.run_stress,
                          epsilon);
        }

        stats.reset();
        // global scalar
        if (ifix->scalar_flag) {
            double value = ifix->compute_scalar();
            EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
        }

        // global vector
        if (ifix->vector_flag) {
            int num = ifix->size_vector;
            EXPECT_EQ(num, test_config.global_vector.size());

            for (int i = 0; i < num; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
                                      epsilon);
        }

        // check t_target for thermostats

        int dim   = -1;
        auto *ptr = (double *)ifix->extract("t_target", dim);
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

    EXPECT_POSITIONS("run_pos (restart, verlet)", lmp->atom, test_config.run_pos, epsilon);
    EXPECT_VELOCITIES("run_vel (restart, verlet)", lmp->atom, test_config.run_vel, epsilon);

    ifix = lmp->modify->get_fix_by_id("test");
    if (!ifix) {
        FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
    } else {
        if (ifix->thermo_virial) {
            EXPECT_STRESS("run_stress (restart, verlet)", ifix->virial, test_config.run_stress,
                          epsilon);
        }

        stats.reset();

        // global scalar
        if (ifix->scalar_flag) {
            double value = ifix->compute_scalar();
            EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
        }

        // global vector
        if (ifix->vector_flag) {
            int num = ifix->size_vector;
            EXPECT_EQ(num, test_config.global_vector.size());

            for (int i = 0; i < num; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
                                      epsilon);
        }
        if (print_stats && stats.has_data())
            std::cerr << "global_data, restart, verlet: " << stats << std::endl;
    }

    if (lmp->atom->rmass == nullptr) {
        if (!verbose) ::testing::internal::CaptureStdout();
        restart_lammps(lmp, test_config, true, false);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        EXPECT_POSITIONS("run_pos (rmass, verlet)", lmp->atom, test_config.run_pos, epsilon);
        EXPECT_VELOCITIES("run_vel (rmass, verlet)", lmp->atom, test_config.run_vel, epsilon);

        ifix = lmp->modify->get_fix_by_id("test");
        if (!ifix) {
            FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
        } else {
            if (ifix->thermo_virial) {
                EXPECT_STRESS("run_stress (rmass, verlet)", ifix->virial, test_config.run_stress,
                              epsilon);
            }

            stats.reset();

            // global scalar
            if (ifix->scalar_flag) {
                double value = ifix->compute_scalar();
                EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
            }

            // global vector
            if (ifix->vector_flag) {
                int num = ifix->size_vector;
                EXPECT_EQ(num, test_config.global_vector.size());

                for (int i = 0; i < num; ++i)
                    EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
                                          epsilon);
            }
            if (print_stats && stats.has_data())
                std::cerr << "global_data, rmass, verlet: " << stats << std::endl;
        }
    }

    // skip the r-RESPA leg for styles tagged "no_respa" (same set as the plain
    // fixture); also, torque is not supported by respa/omp
    ifix = lmp->modify->get_fix_by_id("test");
    if (ifix && !test_config.has_tag("no_respa") && !lmp->atom->torque) {

        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();

        ::testing::internal::CaptureStdout();
        try {
            lmp = init_lammps(args, test_config, true);
        } catch (std::exception &e) {
            output = ::testing::internal::GetCapturedStdout();
            if (verbose) std::cout << output;
            FAIL() << e.what();
        }
        output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;

        // lower required precision by two orders of magnitude to accommodate respa
        epsilon *= 100.0;

        EXPECT_POSITIONS("run_pos (normal run, respa)", lmp->atom, test_config.run_pos, epsilon);
        EXPECT_VELOCITIES("run_vel (normal run, respa)", lmp->atom, test_config.run_vel, epsilon);

        ifix = lmp->modify->get_fix_by_id("test");
        if (!ifix) {
            FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
        } else {
            if (ifix->thermo_virial) {
                EXPECT_STRESS("run_stress (normal run, respa)", ifix->virial,
                              test_config.run_stress, 1000 * epsilon);
            }

            stats.reset();

            // global scalar
            if (ifix->scalar_flag) {
                double value = ifix->compute_scalar();
                EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, 10 * epsilon);
            }

            // global vector
            if (ifix->vector_flag) {
                int num = ifix->size_vector;
                EXPECT_EQ(num, test_config.global_vector.size());

                for (int i = 0; i < num; ++i)
                    EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
                                          10 * epsilon);
            }
            if (print_stats && stats.has_data())
                std::cerr << "global_data, normal run, respa: " << stats << std::endl;
        }

        if (!verbose) ::testing::internal::CaptureStdout();
        restart_lammps(lmp, test_config, false, true);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        EXPECT_POSITIONS("run_pos (restart, respa)", lmp->atom, test_config.run_pos, epsilon);
        EXPECT_VELOCITIES("run_vel (restart, respa)", lmp->atom, test_config.run_vel, epsilon);

        ifix = lmp->modify->get_fix_by_id("test");
        if (!ifix) {
            FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
        } else {
            if (ifix->thermo_virial) {
                EXPECT_STRESS("run_stress (restart, respa)", ifix->virial, test_config.run_stress,
                              1000 * epsilon);
            }

            stats.reset();

            // global scalar
            if (ifix->scalar_flag) {
                double value = ifix->compute_scalar();
                EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, 10 * epsilon);
            }

            // global vector
            if (ifix->vector_flag) {
                int num = ifix->size_vector;
                EXPECT_EQ(num, test_config.global_vector.size());

                for (int i = 0; i < num; ++i)
                    EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
                                          10 * epsilon);
            }
            if (print_stats && stats.has_data())
                std::cerr << "global_data, restart, respa: " << stats << std::endl;
        }

        if (lmp->atom->rmass == nullptr) {
            if (!verbose) ::testing::internal::CaptureStdout();
            restart_lammps(lmp, test_config, true, true);
            if (!verbose) ::testing::internal::GetCapturedStdout();

            EXPECT_POSITIONS("run_pos (rmass, respa)", lmp->atom, test_config.run_pos, epsilon);
            EXPECT_VELOCITIES("run_vel (rmass, respa)", lmp->atom, test_config.run_vel, epsilon);

            ifix = lmp->modify->get_fix_by_id("test");
            if (!ifix) {
                FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
            } else {
                if (ifix->thermo_virial) {
                    EXPECT_STRESS("run_stress (rmass, respa)", ifix->virial, test_config.run_stress,
                                  1000 * epsilon);
                }

                stats.reset();

                // global scalar
                if (ifix->scalar_flag) {
                    double value = ifix->compute_scalar();
                    EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, 10 * epsilon);
                }

                // global vector
                if (ifix->vector_flag) {
                    int num = ifix->size_vector;
                    EXPECT_EQ(num, test_config.global_vector.size());

                    for (int i = 0; i < num; ++i)
                        EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
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

// precision of the KOKKOS package as selected with -D KOKKOS_PREC at compile time
static std::string kokkos_precision()
{
    if (Info::has_accelerator_feature("KOKKOS", "precision", "mixed")) return "mixed";
    if (Info::has_accelerator_feature("KOKKOS", "precision", "single")) return "single";
    return "double";
}

static void run_kokkos_test(LAMMPS::argv &args)
{
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
    // relax error a lot for reduced precision KOKKOS builds
    const std::string kk_precision = kokkos_precision();
    if (kk_precision == "mixed")
        epsilon *= 2.0e9;
    else if (kk_precision == "single")
        epsilon *= 1.0e10;
    // relax test precision when using pppm and single precision FFTs
#if defined(FFT_SINGLE)
    if (lmp->force->kspace && utils::strmatch(lmp->force->kspace_style, "^pppm")) epsilon *= 2.0e8;
#endif

    ErrorStats stats;

    EXPECT_POSITIONS("run_pos (normal run, verlet)", lmp->atom, test_config.run_pos, epsilon);
    EXPECT_VELOCITIES("run_vel (normal run, verlet)", lmp->atom, test_config.run_vel, epsilon);
    if (lmp->atom->torque_flag)
        EXPECT_TORQUES("run_torque (normal run, verlet)", lmp->atom, test_config.run_torque,
                       epsilon);
    EXPECT_SPINS("run_spin (normal run, verlet)", lmp->atom, test_config.run_spin, epsilon);
    EXPECT_MAG_FORCES("run_mag_forces (normal run, verlet)", lmp->atom,
                      test_config.run_mag_forces, epsilon);

    auto *ifix = lmp->modify->get_fix_by_id("test");

    if (!ifix) {
        FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
    } else {

        if (ifix->thermo_virial) {
            EXPECT_STRESS("run_stress (normal run, verlet)", ifix->virial, test_config.run_stress,
                          epsilon);
        }

        stats.reset();
        // global scalar
        if (ifix->scalar_flag) {
            double value = ifix->compute_scalar();
            EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
        }

        // global vector
        if (ifix->vector_flag) {
            int num = ifix->size_vector;
            EXPECT_EQ(num, test_config.global_vector.size());

            for (int i = 0; i < num; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
                                      epsilon);
        }

        // check t_target for thermostats

        int dim     = -1;
        double *ptr = (double *)ifix->extract("t_target", dim);
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

    EXPECT_POSITIONS("run_pos (restart, verlet)", lmp->atom, test_config.run_pos, epsilon);
    EXPECT_VELOCITIES("run_vel (restart, verlet)", lmp->atom, test_config.run_vel, epsilon);
    if (lmp->atom->torque_flag)
        EXPECT_TORQUES("run_torque (restart, verlet)", lmp->atom, test_config.run_torque, epsilon);
    EXPECT_SPINS("run_spin (restart, verlet)", lmp->atom, test_config.run_spin, epsilon);
    EXPECT_MAG_FORCES("run_mag_forces (restart, verlet)", lmp->atom, test_config.run_mag_forces,
                      epsilon);

    ifix = lmp->modify->get_fix_by_id("test");
    if (!ifix) {
        FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
    } else {
        if (ifix->thermo_virial) {
            EXPECT_STRESS("run_stress (restart, verlet)", ifix->virial, test_config.run_stress,
                          epsilon);
        }

        stats.reset();

        // global scalar
        if (ifix->scalar_flag) {
            double value = ifix->compute_scalar();
            EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
        }

        // global vector
        if (ifix->vector_flag) {
            int num = ifix->size_vector;
            EXPECT_EQ(num, test_config.global_vector.size());

            for (int i = 0; i < num; ++i)
                EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
                                      epsilon);
        }
        if (print_stats && stats.has_data())
            std::cerr << "global_data, restart, verlet: " << stats << std::endl;
    }

    if (lmp->atom->rmass == nullptr) {
        if (!verbose) ::testing::internal::CaptureStdout();
        restart_lammps(lmp, test_config, true, false);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        EXPECT_POSITIONS("run_pos (rmass, verlet)", lmp->atom, test_config.run_pos, epsilon);
        EXPECT_VELOCITIES("run_vel (rmass, verlet)", lmp->atom, test_config.run_vel, epsilon);

        ifix = lmp->modify->get_fix_by_id("test");
        if (!ifix) {
            FAIL() << "ERROR: no fix defined with fix ID 'test'\n";
        } else {
            if (ifix->thermo_virial) {
                EXPECT_STRESS("run_stress (rmass, verlet)", ifix->virial, test_config.run_stress,
                              epsilon);
            }

            stats.reset();

            // global scalar
            if (ifix->scalar_flag) {
                double value = ifix->compute_scalar();
                EXPECT_FP_LE_WITH_EPS(test_config.global_scalar, value, epsilon);
            }

            // global vector
            if (ifix->vector_flag) {
                int num = ifix->size_vector;
                EXPECT_EQ(num, test_config.global_vector.size());

                for (int i = 0; i < num; ++i)
                    EXPECT_FP_LE_WITH_EPS(test_config.global_vector[i], ifix->compute_vector(i),
                                          epsilon);
            }
            if (print_stats && stats.has_data())
                std::cerr << "global_data, rmass, verlet: " << stats << std::endl;
        }
    }

    // skip RESPA tests for KOKKOS

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
}

TEST(FixTimestep, kokkos_omp)
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

    LAMMPS::argv args = {"FixTimestep", "-log", "none", "-echo", "screen", "-nocite",
                         "-k",          "on",   "t",    "4",     "-sf",    "kk"};

    run_kokkos_test(args);
};

TEST(FixTimestep, kokkos_serial)
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

    LAMMPS::argv args = {"FixTimestep", "-log", "none", "-echo", "screen", "-nocite",
                         "-k",          "on",   "t",    "1",     "-sf",    "kk"};

    run_kokkos_test(args);
};

TEST(FixTimestep, kokkos_gpu)
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
    LAMMPS::argv args = {"FixTimestep", "-log", "none",   "-echo", "screen", "-nocite", "-k", "on",
                         "g",           "1",    "-sf",    "kk",    "-pk",     "kokkos",  "neigh",
                         "half", "newton", "on"};

    run_kokkos_test(args);
};
