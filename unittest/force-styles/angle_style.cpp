/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// unit tests for angle styles intended for molecular systems

#include "error_stats.h"
#include "test_config.h"
#include "test_config_reader.h"
#include "test_main.h"
#include "yaml_reader.h"
#include "yaml_writer.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "angle.h"
#include "atom.h"
#include "compute.h"
#include "force.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "modify.h"
#include "universe.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <mpi.h>

#include <iostream>
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
    delete_file(cfg.basename + ".data");
    delete_file(cfg.basename + "-coeffs.in");
    delete lmp;
}

LAMMPS *init_lammps(int argc, char **argv, const TestConfig &cfg, const bool newton = true)
{
    LAMMPS *lmp;

    lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);

    // check if prerequisite styles are available
    Info *info = new Info(lmp);
    int nfail  = 0;
    for (auto &prerequisite : cfg.prerequisites) {
        std::string style = prerequisite.second;

        // this is a test for angle styles, so if the suffixed
        // version is not available, there is no reason to test.
        if (prerequisite.first == "angle") {
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
        lmp->input->one(line.c_str());
    };
    auto parse_input_script = [&](const std::string &filename) {
        lmp->input->file(filename.c_str());
    };

    if (newton) {
        command("variable newton_bond index on");
    } else {
        command("variable newton_bond index off");
    }

    command("variable input_dir index " + INPUT_FOLDER);

    for (auto &pre_command : cfg.pre_commands) {
        command(pre_command);
    }

    std::string input_file = INPUT_FOLDER + PATH_SEP + cfg.input_file;
    parse_input_script(input_file);

    command("angle_style " + cfg.angle_style);

    for (auto &angle_coeff : cfg.angle_coeff) {
        command("angle_coeff " + angle_coeff);
    }

    for (auto &post_command : cfg.post_commands) {
        command(post_command);
    }

    command("run 0 post no");
    command("write_restart " + cfg.basename + ".restart");
    command("write_data " + cfg.basename + ".data");
    command("write_coeff " + cfg.basename + "-coeffs.in");

    return lmp;
}

void run_lammps(LAMMPS *lmp)
{
    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line.c_str());
    };

    command("fix 1 all nve");
    command("compute pe all pe/atom angle");
    command("compute sum all reduce sum c_pe");
    command("thermo_style custom step temp pe press c_sum");
    command("thermo 2");
    command("run 4 post no");
}

void restart_lammps(LAMMPS *lmp, const TestConfig &cfg)
{
    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line.c_str());
    };

    command("clear");
    command("read_restart " + cfg.basename + ".restart");

    if (!lmp->force->angle) {
        command("angle_style " + cfg.angle_style);
    }

    if ((cfg.angle_style.substr(0, 6) == "hybrid") || !lmp->force->angle->writedata) {
        for (auto &angle_coeff : cfg.angle_coeff) {
            command("angle_coeff " + angle_coeff);
        }
    }

    for (auto &post_command : cfg.post_commands) {
        command(post_command);
    }

    command("run 0 post no");
}

void data_lammps(LAMMPS *lmp, const TestConfig &cfg)
{
    // utility lambdas to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line.c_str());
    };
    auto parse_input_script = [&](const std::string &filename) {
        lmp->input->file(filename.c_str());
    };

    command("clear");
    command("variable angle_style delete");
    command("variable data_file  delete");
    command("variable newton_bond delete");
    command("variable newton_bond index on");

    for (auto &pre_command : cfg.pre_commands) {
        command(pre_command);
    }

    command("variable angle_style index '" + cfg.angle_style + "'");
    command("variable data_file index " + cfg.basename + ".data");

    std::string input_file = INPUT_FOLDER + PATH_SEP + cfg.input_file;
    parse_input_script(input_file);

    for (auto &angle_coeff : cfg.angle_coeff) {
        command("angle_coeff " + angle_coeff);
    }
    for (auto &post_command : cfg.post_commands) {
        command(post_command);
    }
    command("run 0 post no");
}

// re-generate yaml file with current settings.

void generate_yaml_file(const char *outfile, const TestConfig &config)
{
    // initialize system geometry
    const char *args[] = {"AngleStyle", "-log", "none", "-echo", "screen", "-nocite"};

    char **argv = (char **)args;
    int argc    = sizeof(args) / sizeof(char *);
    LAMMPS *lmp = init_lammps(argc, argv, config);
    if (!lmp) {
        std::cerr << "One or more prerequisite styles are not available "
                     "in this LAMMPS configuration:\n";
        for (auto &prerequisite : config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        return;
    }

    const int natoms  = lmp->atom->natoms;
    const int bufsize = 256;
    char buf[bufsize];
    std::string block("");

    YamlWriter writer(outfile);

    // lammps_version
    writer.emit("lammps_version", lmp->universe->version);

    // date_generated
    std::time_t now = time(NULL);
    block           = ctime(&now);
    block           = block.substr(0, block.find("\n") - 1);
    writer.emit("date_generated", block);

    // epsilon
    writer.emit("epsilon", config.epsilon);

    // prerequisites
    block.clear();
    for (auto &prerequisite : config.prerequisites) {
        block += prerequisite.first + " " + prerequisite.second + "\n";
    }
    writer.emit_block("prerequisites", block);

    // pre_commands
    block.clear();
    for (auto &command : config.pre_commands) {
        block += command + "\n";
    }
    writer.emit_block("pre_commands", block);

    // post_commands
    block.clear();
    for (auto &command : config.post_commands) {
        block += command + "\n";
    }
    writer.emit_block("post_commands", block);

    // input_file
    writer.emit("input_file", config.input_file);

    // angle_style
    writer.emit("angle_style", config.angle_style);

    // angle_coeff
    block.clear();
    for (auto &angle_coeff : config.angle_coeff) {
        block += angle_coeff + "\n";
    }
    writer.emit_block("angle_coeff", block);

    // equilibrium angle
    std::stringstream eqstr;
    eqstr << lmp->atom->nangletypes;
    for (int i = 0; i < lmp->atom->nangletypes; ++i) {
        eqstr << " " << lmp->force->angle->equilibrium_angle(i + 1);
    }
    writer.emit("equilibrium", eqstr.str());

    // extract
    block.clear();
    std::stringstream outstr;
    for (auto &data : config.extract) {
        outstr << data.first << " " << data.second << std::endl;
    }
    writer.emit_block("extract", outstr.str());

    // natoms
    writer.emit("natoms", natoms);

    // init_energy
    writer.emit("init_energy", lmp->force->angle->energy);

    // init_stress
    double *stress = lmp->force->angle->virial;
    snprintf(buf, bufsize, "% 23.16e % 23.16e % 23.16e % 23.16e % 23.16e % 23.16e", stress[0],
             stress[1], stress[2], stress[3], stress[4], stress[5]);
    writer.emit_block("init_stress", buf);

    // init_forces
    block.clear();
    double **f  = lmp->atom->f;
    tagint *tag = lmp->atom->tag;
    for (int i = 0; i < natoms; ++i) {
        snprintf(buf, bufsize, "% 3d % 23.16e % 23.16e % 23.16e\n", (int)tag[i], f[i][0], f[i][1],
                 f[i][2]);
        block += buf;
    }
    writer.emit_block("init_forces", block);

    // do a few steps of MD
    run_lammps(lmp);

    // run_energy
    writer.emit("run_energy", lmp->force->angle->energy);

    // run_stress
    stress = lmp->force->angle->virial;
    snprintf(buf, bufsize, "% 23.16e % 23.16e % 23.16e % 23.16e % 23.16e % 23.16e", stress[0],
             stress[1], stress[2], stress[3], stress[4], stress[5]);
    writer.emit_block("run_stress", buf);

    block.clear();
    f   = lmp->atom->f;
    tag = lmp->atom->tag;
    for (int i = 0; i < natoms; ++i) {
        snprintf(buf, bufsize, "% 3d % 23.16e % 23.16e % 23.16e\n", (int)tag[i], f[i][0], f[i][1],
                 f[i][2]);
        block += buf;
    }
    writer.emit_block("run_forces", block);

    cleanup_lammps(lmp, config);
    return;
}

TEST(AngleStyle, plain)
{
    const char *args[] = {"AngleStyle", "-log", "none", "-echo", "screen", "-nocite"};

    char **argv = (char **)args;
    int argc    = sizeof(args) / sizeof(char *);

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = init_lammps(argc, argv, test_config, true);

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
    double **f     = lmp->atom->f;
    tagint *tag    = lmp->atom->tag;
    ErrorStats stats;
    stats.reset();
    const std::vector<coord_t> &f_ref = test_config.init_forces;
    ASSERT_EQ(nlocal + 1, f_ref.size());
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(f[i][0], f_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][1], f_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][2], f_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << "init_forces stats, newton on: " << stats << std::endl;

    Angle *angle   = lmp->force->angle;
    double *stress = angle->virial;
    stats.reset();
    EXPECT_FP_LE_WITH_EPS(stress[0], test_config.init_stress.xx, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[1], test_config.init_stress.yy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[2], test_config.init_stress.zz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[3], test_config.init_stress.xy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[4], test_config.init_stress.xz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[5], test_config.init_stress.yz, epsilon);
    if (print_stats) std::cerr << "init_stress stats, newton on: " << stats << std::endl;

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(angle->energy, test_config.init_energy, epsilon);
    if (print_stats) std::cerr << "init_energy stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_lammps(lmp);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    f                                 = lmp->atom->f;
    stress                            = angle->virial;
    const std::vector<coord_t> &f_run = test_config.run_forces;
    ASSERT_EQ(nlocal + 1, f_run.size());
    stats.reset();
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(f[i][0], f_run[tag[i]].x, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][1], f_run[tag[i]].y, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][2], f_run[tag[i]].z, 10 * epsilon);
    }
    if (print_stats) std::cerr << "run_forces  stats, newton on: " << stats << std::endl;

    stress = angle->virial;
    stats.reset();
    EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, epsilon);
    if (print_stats) std::cerr << "run_stress  stats, newton on: " << stats << std::endl;

    stats.reset();
    int id        = lmp->modify->find_compute("sum");
    double energy = lmp->modify->compute[id]->compute_scalar();
    EXPECT_FP_LE_WITH_EPS(angle->energy, test_config.run_energy, epsilon);
    EXPECT_FP_LE_WITH_EPS(angle->energy, energy, epsilon);
    if (print_stats) std::cerr << "run_energy  stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    lmp = init_lammps(argc, argv, test_config, false);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // skip over these tests if newton bond is forced to be on
    if (lmp->force->newton_bond == 0) {

        f   = lmp->atom->f;
        tag = lmp->atom->tag;
        stats.reset();
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(f[i][0], f_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(f[i][1], f_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(f[i][2], f_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "init_forces stats, newton off:" << stats << std::endl;

        angle  = lmp->force->angle;
        stress = angle->virial;
        stats.reset();
        EXPECT_FP_LE_WITH_EPS(stress[0], test_config.init_stress.xx, 2 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[1], test_config.init_stress.yy, 2 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[2], test_config.init_stress.zz, 2 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[3], test_config.init_stress.xy, 2 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[4], test_config.init_stress.xz, 2 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[5], test_config.init_stress.yz, 2 * epsilon);
        if (print_stats) std::cerr << "init_stress stats, newton off:" << stats << std::endl;

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(angle->energy, test_config.init_energy, epsilon);
        if (print_stats) std::cerr << "init_energy stats, newton off:" << stats << std::endl;

        if (!verbose) ::testing::internal::CaptureStdout();
        run_lammps(lmp);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        f      = lmp->atom->f;
        stress = angle->virial;
        stats.reset();
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(f[i][0], f_run[tag[i]].x, 10 * epsilon);
            EXPECT_FP_LE_WITH_EPS(f[i][1], f_run[tag[i]].y, 10 * epsilon);
            EXPECT_FP_LE_WITH_EPS(f[i][2], f_run[tag[i]].z, 10 * epsilon);
        }
        if (print_stats) std::cerr << "run_forces  stats, newton off:" << stats << std::endl;

        stress = angle->virial;
        stats.reset();
        EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, epsilon);
        if (print_stats) std::cerr << "run_stress  stats, newton off:" << stats << std::endl;

        stats.reset();
        id     = lmp->modify->find_compute("sum");
        energy = lmp->modify->compute[id]->compute_scalar();
        EXPECT_FP_LE_WITH_EPS(angle->energy, test_config.run_energy, epsilon);
        EXPECT_FP_LE_WITH_EPS(angle->energy, energy, epsilon);
        if (print_stats) std::cerr << "run_energy  stats, newton off:" << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    restart_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    f   = lmp->atom->f;
    tag = lmp->atom->tag;
    stats.reset();
    ASSERT_EQ(nlocal + 1, f_ref.size());
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(f[i][0], f_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][1], f_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][2], f_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << "restart_forces stats:" << stats << std::endl;

    angle  = lmp->force->angle;
    stress = angle->virial;
    stats.reset();
    EXPECT_FP_LE_WITH_EPS(stress[0], test_config.init_stress.xx, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[1], test_config.init_stress.yy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[2], test_config.init_stress.zz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[3], test_config.init_stress.xy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[4], test_config.init_stress.xz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[5], test_config.init_stress.yz, epsilon);
    if (print_stats) std::cerr << "restart_stress stats:" << stats << std::endl;

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(angle->energy, test_config.init_energy, epsilon);
    if (print_stats) std::cerr << "restart_energy stats:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    data_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    f   = lmp->atom->f;
    tag = lmp->atom->tag;
    stats.reset();
    ASSERT_EQ(nlocal + 1, f_ref.size());
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(f[i][0], f_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][1], f_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][2], f_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << "data_forces stats:" << stats << std::endl;

    angle  = lmp->force->angle;
    stress = angle->virial;
    stats.reset();
    EXPECT_FP_LE_WITH_EPS(stress[0], test_config.init_stress.xx, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[1], test_config.init_stress.yy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[2], test_config.init_stress.zz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[3], test_config.init_stress.xy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[4], test_config.init_stress.xz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[5], test_config.init_stress.yz, epsilon);
    if (print_stats) std::cerr << "data_stress stats:" << stats << std::endl;

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(angle->energy, test_config.init_energy, epsilon);
    if (print_stats) std::cerr << "data_energy stats:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

TEST(AngleStyle, omp)
{
    if (!LAMMPS::is_installed_pkg("USER-OMP")) GTEST_SKIP();
    const char *args[] = {"AngleStyle", "-log", "none", "-echo", "screen", "-nocite",
                          "-pk",        "omp",  "4",    "-sf",   "omp"};

    char **argv = (char **)args;
    int argc    = sizeof(args) / sizeof(char *);

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = init_lammps(argc, argv, test_config, true);

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

    // relax error a bit for USER-OMP package
    double epsilon                    = 5.0 * test_config.epsilon;
    double **f                        = lmp->atom->f;
    tagint *tag                       = lmp->atom->tag;
    const std::vector<coord_t> &f_ref = test_config.init_forces;
    ErrorStats stats;
    stats.reset();
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(f[i][0], f_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][1], f_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][2], f_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << "init_forces stats, newton on: " << stats << std::endl;

    Angle *angle   = lmp->force->angle;
    double *stress = angle->virial;

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(stress[0], test_config.init_stress.xx, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[1], test_config.init_stress.yy, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[2], test_config.init_stress.zz, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[3], test_config.init_stress.xy, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[4], test_config.init_stress.xz, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[5], test_config.init_stress.yz, 10 * epsilon);
    if (print_stats) std::cerr << "init_stress stats, newton on: " << stats << std::endl;

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(angle->energy, test_config.init_energy, epsilon);
    if (print_stats) std::cerr << "init_energy stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_lammps(lmp);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    f                                 = lmp->atom->f;
    stress                            = angle->virial;
    const std::vector<coord_t> &f_run = test_config.run_forces;
    ASSERT_EQ(nlocal + 1, f_run.size());
    stats.reset();
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(f[i][0], f_run[tag[i]].x, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][1], f_run[tag[i]].y, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][2], f_run[tag[i]].z, 10 * epsilon);
    }
    if (print_stats) std::cerr << "run_forces  stats, newton on: " << stats << std::endl;

    stress = angle->virial;
    stats.reset();
    EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, 10 * epsilon);
    if (print_stats) std::cerr << "run_stress  stats, newton on: " << stats << std::endl;

    stats.reset();
    int id        = lmp->modify->find_compute("sum");
    double energy = lmp->modify->compute[id]->compute_scalar();
    EXPECT_FP_LE_WITH_EPS(angle->energy, test_config.run_energy, epsilon);
    // TODO: this is currently broken for USER-OMP with angle style hybrid
    // needs to be fixed in the main code somewhere. Not sure where, though.
    if (test_config.angle_style.substr(0, 6) != "hybrid")
        EXPECT_FP_LE_WITH_EPS(angle->energy, energy, epsilon);
    if (print_stats) std::cerr << "run_energy  stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    lmp = init_lammps(argc, argv, test_config, false);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // skip over these tests if newton bond is forced to be on
    if (lmp->force->newton_bond == 0) {

        f   = lmp->atom->f;
        tag = lmp->atom->tag;
        stats.reset();
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(f[i][0], f_ref[tag[i]].x, epsilon);
            EXPECT_FP_LE_WITH_EPS(f[i][1], f_ref[tag[i]].y, epsilon);
            EXPECT_FP_LE_WITH_EPS(f[i][2], f_ref[tag[i]].z, epsilon);
        }
        if (print_stats) std::cerr << "init_forces stats, newton off:" << stats << std::endl;

        angle  = lmp->force->angle;
        stress = angle->virial;
        stats.reset();
        EXPECT_FP_LE_WITH_EPS(stress[0], test_config.init_stress.xx, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[1], test_config.init_stress.yy, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[2], test_config.init_stress.zz, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[3], test_config.init_stress.xy, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[4], test_config.init_stress.xz, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[5], test_config.init_stress.yz, 10 * epsilon);
        if (print_stats) std::cerr << "init_stress stats, newton off:" << stats << std::endl;

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(angle->energy, test_config.init_energy, epsilon);
        if (print_stats) std::cerr << "init_energy stats, newton off:" << stats << std::endl;

        if (!verbose) ::testing::internal::CaptureStdout();
        run_lammps(lmp);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        f = lmp->atom->f;
        stats.reset();
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(f[i][0], f_run[tag[i]].x, 10 * epsilon);
            EXPECT_FP_LE_WITH_EPS(f[i][1], f_run[tag[i]].y, 10 * epsilon);
            EXPECT_FP_LE_WITH_EPS(f[i][2], f_run[tag[i]].z, 10 * epsilon);
        }
        if (print_stats) std::cerr << "run_forces  stats, newton off:" << stats << std::endl;

        stress = angle->virial;
        stats.reset();
        EXPECT_FP_LE_WITH_EPS(stress[0], test_config.run_stress.xx, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[1], test_config.run_stress.yy, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[2], test_config.run_stress.zz, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[3], test_config.run_stress.xy, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[4], test_config.run_stress.xz, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[5], test_config.run_stress.yz, 10 * epsilon);
        if (print_stats) std::cerr << "run_stress  stats, newton off:" << stats << std::endl;

        stats.reset();
        id     = lmp->modify->find_compute("sum");
        energy = lmp->modify->compute[id]->compute_scalar();
        EXPECT_FP_LE_WITH_EPS(angle->energy, test_config.run_energy, epsilon);
        // TODO: this is currently broken for USER-OMP with angle style hybrid
        // needs to be fixed in the main code somewhere. Not sure where, though.
        if (test_config.angle_style.substr(0, 6) != "hybrid")
            EXPECT_FP_LE_WITH_EPS(angle->energy, energy, epsilon);
        if (print_stats) std::cerr << "run_energy  stats, newton off:" << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

TEST(AngleStyle, single)
{
    const char *args[] = {"AngleStyle", "-log", "none", "-echo", "screen", "-nocite"};

    char **argv = (char **)args;
    int argc    = sizeof(args) / sizeof(char *);

    // create a LAMMPS instance with standard settings to detect the number of atom types
    if (!verbose) ::testing::internal::CaptureStdout();
    LAMMPS *lmp = init_lammps(argc, argv, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    if (!lmp) {
        std::cerr << "One or more prerequisite styles are not available "
                     "in this LAMMPS configuration:\n";
        for (auto &prerequisite : test_config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        GTEST_SKIP();
    }

    // gather some information and skip if unsupported
    int nangletypes = lmp->atom->nangletypes;
    int molecular   = lmp->atom->molecular;
    if (molecular != 1) {
        std::cerr << "Only simple molecular atom styles are supported\n";
        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        GTEST_SKIP();
    }

    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line.c_str());
    };

    // now start over
    if (!verbose) ::testing::internal::CaptureStdout();
    command("clear");
    command("variable newton_bond delete");
    command("variable newton_bond index on");

    command("variable input_dir index " + INPUT_FOLDER);

    for (auto &pre_command : test_config.pre_commands) {
        command(pre_command);
    }

    command("atom_style molecular");
    command("units ${units}");
    command("boundary p p p");
    command("newton ${newton_pair} ${newton_bond}");
    command("special_bonds lj/coul "
            "${bond_factor} ${angle_factor} ${dihedral_factor}");

    command("atom_modify map array");
    command("region box block -10.0 10.0 -10.0 10.0 -10.0 10.0 units box");

    char buf[10];
    std::string cmd("create_box 1 box");
    cmd += " angle/types ";
    snprintf(buf, 10, "%d", nangletypes);
    cmd += buf;
    cmd += " extra/angle/per/atom 2";
    cmd += " extra/special/per/atom 2";
    command(cmd);

    command("pair_style zero 8.0");
    command("pair_coeff * *");

    command("angle_style " + test_config.angle_style);
    Angle *angle = lmp->force->angle;

    for (auto &angle_coeff : test_config.angle_coeff) {
        command("angle_coeff " + angle_coeff);
    }

    // create (only) three atoms and one angle
    command("mass * 1.0");
    command("create_atoms 1 single  5.0 -0.75  0.4 units box");
    command("create_atoms 1 single  5.5  0.25 -0.1 units box");
    command("create_atoms 1 single  5.0  0.75  0.4 units box");
    command("create_bonds single/angle 1 1 2 3");

    for (auto &post_command : test_config.post_commands) {
        command(post_command);
    }

    command("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    int idx1       = lmp->atom->map(1);
    int idx2       = lmp->atom->map(2);
    int idx3       = lmp->atom->map(3);
    double epsilon = test_config.epsilon;
    double eangle[4], esingle[4];

    eangle[0]  = angle->energy;
    esingle[0] = angle->single(1, idx1, idx2, idx3);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("displace_atoms all random 0.5 0.5 0.5 23456");
    command("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    idx1       = lmp->atom->map(1);
    idx2       = lmp->atom->map(2);
    idx3       = lmp->atom->map(3);
    eangle[1]  = angle->energy;
    esingle[1] = angle->single(1, idx1, idx2, idx3);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("displace_atoms all random 0.5 0.5 0.5 456963");
    command("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    idx1       = lmp->atom->map(1);
    idx2       = lmp->atom->map(2);
    idx3       = lmp->atom->map(3);
    eangle[2]  = angle->energy;
    esingle[2] = angle->single(1, idx1, idx2, idx3);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("displace_atoms all random 0.5 0.5 0.5 9726532");
    command("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    idx1       = lmp->atom->map(1);
    idx2       = lmp->atom->map(2);
    idx3       = lmp->atom->map(3);
    eangle[3]  = angle->energy;
    esingle[3] = angle->single(1, idx1, idx2, idx3);

    ErrorStats stats;
    EXPECT_FP_LE_WITH_EPS(eangle[0], esingle[0], epsilon);
    EXPECT_FP_LE_WITH_EPS(eangle[1], esingle[1], epsilon);
    EXPECT_FP_LE_WITH_EPS(eangle[2], esingle[2], epsilon);
    EXPECT_FP_LE_WITH_EPS(eangle[3], esingle[3], epsilon);
    if (print_stats) std::cerr << "single_energy  stats:" << stats << std::endl;

    int i = 0;
    for (auto &dist : test_config.equilibrium)
        EXPECT_NEAR(dist, angle->equilibrium_angle(++i), 0.00001);

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
}
