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

// unit tests for dihedral styles intended for molecular systems

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
#include "dihedral.h"
#include "fmt/format.h"
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

        // this is a test for dihedral styles, so if the suffixed
        // version is not available, there is no reason to test.
        if (prerequisite.first == "dihedral") {
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

    command("dihedral_style " + cfg.dihedral_style);

    for (auto &dihedral_coeff : cfg.dihedral_coeff) {
        command("dihedral_coeff " + dihedral_coeff);
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
    command("compute pe all pe/atom dihedral");
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

    if (!lmp->force->dihedral) {
        command("dihedral_style " + cfg.dihedral_style);
    }

    if ((cfg.dihedral_style.substr(0, 6) == "hybrid") || !lmp->force->dihedral->writedata) {
        for (auto &dihedral_coeff : cfg.dihedral_coeff) {
            command("dihedral_coeff " + dihedral_coeff);
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
    command("variable dihedral_style delete");
    command("variable data_file  delete");
    command("variable newton_bond delete");
    command("variable newton_bond index on");

    for (auto &pre_command : cfg.pre_commands) {
        command(pre_command);
    }

    command("variable dihedral_style index '" + cfg.dihedral_style + "'");
    command("variable data_file index " + cfg.basename + ".data");

    // special treatment for dihedral styles charmm and charmmfsw
    if (cfg.dihedral_style == "charmm") {
        command("variable pair_style delete");
        command("variable pair_style index 'lj/charmm/coul/charmm 7.0 8.0'");
    } else if (cfg.dihedral_style == "charmmfsw") {
        command("variable pair_style delete");
        command("variable pair_style index 'lj/charmmfsw/coul/charmmfsh 7.0 8.0'");
    }

    std::string input_file = INPUT_FOLDER + PATH_SEP + cfg.input_file;
    parse_input_script(input_file);

    for (auto &dihedral_coeff : cfg.dihedral_coeff) {
        command("dihedral_coeff " + dihedral_coeff);
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
    const char *args[] = {"DihedralStyle", "-log", "none", "-echo", "screen", "-nocite"};

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

    const int natoms = lmp->atom->natoms;
    std::string block("");

    YamlWriter writer(outfile);

    // write yaml header
    write_yaml_header(&writer, &test_config, lmp->version);

    // dihedral_style
    writer.emit("dihedral_style", config.dihedral_style);

    // dihedral_coeff
    block.clear();
    for (auto &dihedral_coeff : config.dihedral_coeff) {
        block += dihedral_coeff + "\n";
    }
    writer.emit_block("dihedral_coeff", block);

    // extract
    block.clear();
    for (auto data : config.extract)
        block += fmt::format("{} {}\n", data.first, data.second);
    writer.emit_block("extract", block);

    // natoms
    writer.emit("natoms", natoms);

    // init_energy
    writer.emit("init_energy", lmp->force->dihedral->energy);

    // init_stress
    auto stress = lmp->force->dihedral->virial;
    block = fmt::format("{:23.16e} {:23.16e} {:23.16e} {:23.16e} {:23.16e} {:23.16e}", stress[0],
                        stress[1], stress[2], stress[3], stress[4], stress[5]);
    writer.emit_block("init_stress", block);

    // init_forces
    block.clear();
    auto f = lmp->atom->f;
    for (int i = 1; i <= natoms; ++i) {
        const int j = lmp->atom->map(i);
        block += fmt::format("{:3} {:23.16e} {:23.16e} {:23.16e}\n", i, f[j][0], f[j][1], f[j][2]);
    }
    writer.emit_block("init_forces", block);

    // do a few steps of MD
    run_lammps(lmp);

    // run_energy
    writer.emit("run_energy", lmp->force->dihedral->energy);

    // run_stress
    stress = lmp->force->dihedral->virial;
    block  = fmt::format("{:23.16e} {:23.16e} {:23.16e} {:23.16e} {:23.16e} {:23.16e}", stress[0],
                        stress[1], stress[2], stress[3], stress[4], stress[5]);
    writer.emit_block("run_stress", block);

    block.clear();
    f = lmp->atom->f;
    for (int i = 1; i <= natoms; ++i) {
        const int j = lmp->atom->map(i);
        block += fmt::format("{:3} {:23.16e} {:23.16e} {:23.16e}\n", i, f[j][0], f[j][1], f[j][2]);
    }
    writer.emit_block("run_forces", block);

    cleanup_lammps(lmp, config);
    return;
}

TEST(DihedralStyle, plain)
{
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    const char *args[] = {"DihedralStyle", "-log", "none", "-echo", "screen", "-nocite"};

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

    auto f   = lmp->atom->f;
    auto tag = lmp->atom->tag;
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

    auto dihedral = lmp->force->dihedral;
    auto stress   = dihedral->virial;
    stats.reset();
    EXPECT_FP_LE_WITH_EPS(stress[0], test_config.init_stress.xx, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[1], test_config.init_stress.yy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[2], test_config.init_stress.zz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[3], test_config.init_stress.xy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[4], test_config.init_stress.xz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[5], test_config.init_stress.yz, epsilon);
    if (print_stats) std::cerr << "init_stress stats, newton on: " << stats << std::endl;

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(dihedral->energy, test_config.init_energy, epsilon);
    if (print_stats) std::cerr << "init_energy stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_lammps(lmp);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    f      = lmp->atom->f;
    tag    = lmp->atom->tag;
    stress = dihedral->virial;

    const std::vector<coord_t> &f_run = test_config.run_forces;
    ASSERT_EQ(nlocal + 1, f_run.size());
    stats.reset();
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(f[i][0], f_run[tag[i]].x, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][1], f_run[tag[i]].y, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][2], f_run[tag[i]].z, 10 * epsilon);
    }
    if (print_stats) std::cerr << "run_forces  stats, newton on: " << stats << std::endl;

    stress = dihedral->virial;
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
    EXPECT_FP_LE_WITH_EPS(dihedral->energy, test_config.run_energy, epsilon);
    EXPECT_FP_LE_WITH_EPS(dihedral->energy, energy, epsilon);
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

        dihedral = lmp->force->dihedral;
        stress   = dihedral->virial;
        stats.reset();
        EXPECT_FP_LE_WITH_EPS(stress[0], test_config.init_stress.xx, 2 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[1], test_config.init_stress.yy, 2 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[2], test_config.init_stress.zz, 2 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[3], test_config.init_stress.xy, 2 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[4], test_config.init_stress.xz, 2 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[5], test_config.init_stress.yz, 2 * epsilon);
        if (print_stats) std::cerr << "init_stress stats, newton off:" << stats << std::endl;

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(dihedral->energy, test_config.init_energy, epsilon);
        if (print_stats) std::cerr << "init_energy stats, newton off:" << stats << std::endl;

        if (!verbose) ::testing::internal::CaptureStdout();
        run_lammps(lmp);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        f      = lmp->atom->f;
        tag    = lmp->atom->tag;
        stress = dihedral->virial;
        stats.reset();
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(f[i][0], f_run[tag[i]].x, 10 * epsilon);
            EXPECT_FP_LE_WITH_EPS(f[i][1], f_run[tag[i]].y, 10 * epsilon);
            EXPECT_FP_LE_WITH_EPS(f[i][2], f_run[tag[i]].z, 10 * epsilon);
        }
        if (print_stats) std::cerr << "run_forces  stats, newton off:" << stats << std::endl;

        stress = dihedral->virial;
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
        EXPECT_FP_LE_WITH_EPS(dihedral->energy, test_config.run_energy, epsilon);
        EXPECT_FP_LE_WITH_EPS(dihedral->energy, energy, epsilon);
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

    dihedral = lmp->force->dihedral;
    stress   = dihedral->virial;
    stats.reset();
    EXPECT_FP_LE_WITH_EPS(stress[0], test_config.init_stress.xx, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[1], test_config.init_stress.yy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[2], test_config.init_stress.zz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[3], test_config.init_stress.xy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[4], test_config.init_stress.xz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[5], test_config.init_stress.yz, epsilon);
    if (print_stats) std::cerr << "restart_stress stats:" << stats << std::endl;

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(dihedral->energy, test_config.init_energy, epsilon);
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

    dihedral = lmp->force->dihedral;
    stress   = dihedral->virial;
    stats.reset();
    EXPECT_FP_LE_WITH_EPS(stress[0], test_config.init_stress.xx, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[1], test_config.init_stress.yy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[2], test_config.init_stress.zz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[3], test_config.init_stress.xy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[4], test_config.init_stress.xz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[5], test_config.init_stress.yz, epsilon);
    if (print_stats) std::cerr << "data_stress stats:" << stats << std::endl;

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(dihedral->energy, test_config.init_energy, epsilon);
    if (print_stats) std::cerr << "data_energy stats:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

TEST(DihedralStyle, omp)
{
    if (!LAMMPS::is_installed_pkg("OPENMP")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    const char *args[] = {"DihedralStyle", "-log", "none", "-echo", "screen", "-nocite",
                          "-pk",           "omp",  "4",    "-sf",   "omp"};

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

    // relax error a bit for OPENMP package
    double epsilon = 5.0 * test_config.epsilon;

    auto f   = lmp->atom->f;
    auto tag = lmp->atom->tag;

    const std::vector<coord_t> &f_ref = test_config.init_forces;
    ErrorStats stats;
    stats.reset();
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(f[i][0], f_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][1], f_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][2], f_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << "init_forces stats, newton on: " << stats << std::endl;

    auto dihedral = lmp->force->dihedral;
    auto stress   = dihedral->virial;

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(stress[0], test_config.init_stress.xx, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[1], test_config.init_stress.yy, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[2], test_config.init_stress.zz, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[3], test_config.init_stress.xy, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[4], test_config.init_stress.xz, 10 * epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[5], test_config.init_stress.yz, 10 * epsilon);
    if (print_stats) std::cerr << "init_stress stats, newton on: " << stats << std::endl;

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(dihedral->energy, test_config.init_energy, epsilon);
    if (print_stats) std::cerr << "init_energy stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_lammps(lmp);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    f      = lmp->atom->f;
    tag    = lmp->atom->tag;
    stress = dihedral->virial;

    const std::vector<coord_t> &f_run = test_config.run_forces;
    ASSERT_EQ(nlocal + 1, f_run.size());
    stats.reset();
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(f[i][0], f_run[tag[i]].x, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][1], f_run[tag[i]].y, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][2], f_run[tag[i]].z, 10 * epsilon);
    }
    if (print_stats) std::cerr << "run_forces  stats, newton on: " << stats << std::endl;

    stress = dihedral->virial;
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
    EXPECT_FP_LE_WITH_EPS(dihedral->energy, test_config.run_energy, epsilon);
    // TODO: this is currently broken for OPENMP with dihedral style hybrid
    // needs to be fixed in the main code somewhere. Not sure where, though.
    if (test_config.dihedral_style.substr(0, 6) != "hybrid")
        EXPECT_FP_LE_WITH_EPS(dihedral->energy, energy, epsilon);
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

        dihedral = lmp->force->dihedral;
        stress   = dihedral->virial;
        stats.reset();
        EXPECT_FP_LE_WITH_EPS(stress[0], test_config.init_stress.xx, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[1], test_config.init_stress.yy, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[2], test_config.init_stress.zz, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[3], test_config.init_stress.xy, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[4], test_config.init_stress.xz, 10 * epsilon);
        EXPECT_FP_LE_WITH_EPS(stress[5], test_config.init_stress.yz, 10 * epsilon);
        if (print_stats) std::cerr << "init_stress stats, newton off:" << stats << std::endl;

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(dihedral->energy, test_config.init_energy, epsilon);
        if (print_stats) std::cerr << "init_energy stats, newton off:" << stats << std::endl;

        if (!verbose) ::testing::internal::CaptureStdout();
        run_lammps(lmp);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        f   = lmp->atom->f;
        tag = lmp->atom->tag;
        stats.reset();
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(f[i][0], f_run[tag[i]].x, 10 * epsilon);
            EXPECT_FP_LE_WITH_EPS(f[i][1], f_run[tag[i]].y, 10 * epsilon);
            EXPECT_FP_LE_WITH_EPS(f[i][2], f_run[tag[i]].z, 10 * epsilon);
        }
        if (print_stats) std::cerr << "run_forces  stats, newton off:" << stats << std::endl;

        stress = dihedral->virial;
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
        EXPECT_FP_LE_WITH_EPS(dihedral->energy, test_config.run_energy, epsilon);
        // TODO: this is currently broken for OPENMP with dihedral style hybrid
        // needs to be fixed in the main code somewhere. Not sure where, though.
        if (test_config.dihedral_style.substr(0, 6) != "hybrid")
            EXPECT_FP_LE_WITH_EPS(dihedral->energy, energy, epsilon);
        if (print_stats) std::cerr << "run_energy  stats, newton off:" << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};
