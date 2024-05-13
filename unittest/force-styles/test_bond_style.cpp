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

// unit tests for bond styles intended for molecular systems

#include "error_stats.h"
#include "test_config.h"
#include "test_config_reader.h"
#include "test_main.h"
#include "yaml_reader.h"
#include "yaml_writer.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "atom.h"
#include "bond.h"
#include "compute.h"
#include "exceptions.h"
#include "fix.h"
#include "fmt/format.h"
#include "force.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "modify.h"
#include "platform.h"
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

void cleanup_lammps(LAMMPS *lmp, const TestConfig &cfg)
{
    platform::unlink(cfg.basename + ".restart");
    platform::unlink(cfg.basename + ".data");
    platform::unlink(cfg.basename + "-coeffs.in");
    delete lmp;
}

LAMMPS *init_lammps(LAMMPS::argv &args, const TestConfig &cfg, const bool newton = true)
{
    LAMMPS *lmp;

    lmp = new LAMMPS(args, MPI_COMM_WORLD);

    // check if prerequisite styles are available
    Info *info = new Info(lmp);
    int nfail  = 0;
    for (auto &prerequisite : cfg.prerequisites) {
        std::string style = prerequisite.second;

        // this is a test for bond styles, so if the suffixed
        // version is not available, there is no reason to test.
        if (prerequisite.first == "bond") {
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
        command("variable newton_bond index on");
    } else {
        command("variable newton_bond index off");
    }

    command("variable input_dir index " + INPUT_FOLDER);

    for (auto &pre_command : cfg.pre_commands) {
        command(pre_command);
    }

    std::string input_file = platform::path_join(INPUT_FOLDER, cfg.input_file);
    parse_input_script(input_file);

    command("bond_style " + cfg.bond_style);

    for (auto &bond_coeff : cfg.bond_coeff) {
        command("bond_coeff " + bond_coeff);
    }

    for (auto &post_command : cfg.post_commands) {
        command(post_command);
    }

    command("run 0 post no");
    command("write_restart " + cfg.basename + ".restart");
    command("write_data " + cfg.basename + ".data nofix");
    command("write_coeff " + cfg.basename + "-coeffs.in");

    return lmp;
}

void run_lammps(LAMMPS *lmp)
{
    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    command("fix 1 all nve");
    command("compute pe all pe/atom bond");
    command("compute sum all reduce sum c_pe");
    command("thermo_style custom step temp pe press c_sum");
    command("thermo 2");
    command("run 4 post no");
}

void restart_lammps(LAMMPS *lmp, const TestConfig &cfg)
{
    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    command("clear");
    command("read_restart " + cfg.basename + ".restart");

    if (!lmp->force->bond) {
        command("bond_style " + cfg.bond_style);
    }

    if ((cfg.bond_style.substr(0, 6) == "hybrid") || !lmp->force->bond->writedata) {
        for (auto &bond_coeff : cfg.bond_coeff) {
            command("bond_coeff " + bond_coeff);
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
        lmp->input->one(line);
    };
    auto parse_input_script = [&](const std::string &filename) {
        lmp->input->file(filename.c_str());
    };

    command("clear");
    command("variable bond_style delete");
    command("variable data_file  delete");
    command("variable newton_bond delete");
    command("variable newton_bond index on");

    for (auto &pre_command : cfg.pre_commands) {
        command(pre_command);
    }

    command("variable bond_style index '" + cfg.bond_style + "'");
    command("variable data_file index " + cfg.basename + ".data");

    std::string input_file = platform::path_join(INPUT_FOLDER, cfg.input_file);
    parse_input_script(input_file);

    for (auto &bond_coeff : cfg.bond_coeff) {
        command("bond_coeff " + bond_coeff);
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
    LAMMPS::argv args = {"BondStyle", "-log", "none", "-echo", "screen", "-nocite"};

    LAMMPS *lmp = init_lammps(args, config);
    if (!lmp) {
        std::cerr << "One or more prerequisite styles are not available "
                     "in this LAMMPS configuration:\n";
        for (auto &prerequisite : config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        return;
    }

    const int natoms = lmp->atom->natoms;
    std::string block;

    YamlWriter writer(outfile);

    // write yaml header
    write_yaml_header(&writer, &test_config, lmp->version);

    // bond_style
    writer.emit("bond_style", config.bond_style);

    // bond_coeff
    block.clear();
    for (auto &bond_coeff : config.bond_coeff) {
        block += bond_coeff + "\n";
    }
    writer.emit_block("bond_coeff", block);

    // equilibrium distance
    block = fmt::format("{}", lmp->atom->nbondtypes);
    for (int i = 0; i < lmp->atom->nbondtypes; ++i)
        block += fmt::format(" {}", lmp->force->bond->equilibrium_distance(i + 1));
    writer.emit("equilibrium", block);

    // extract
    block.clear();
    for (auto data : config.extract)
        block += fmt::format("{} {}\n", data.first, data.second);
    writer.emit_block("extract", block);

    // natoms
    writer.emit("natoms", natoms);

    // init_energy
    writer.emit("init_energy", lmp->force->bond->energy);

    // init_stress
    auto stress = lmp->force->bond->virial;
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
    writer.emit("run_energy", lmp->force->bond->energy);

    // run_stress
    stress = lmp->force->bond->virial;
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
}

TEST(BondStyle, plain)
{
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"BondStyle", "-log", "none", "-echo", "screen", "-nocite"};

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

    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms, nlocal);

    double epsilon = test_config.epsilon;

    ErrorStats stats;
    auto bond = lmp->force->bond;

    EXPECT_FORCES("init_forces (newton on)", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_STRESS("init_stress (newton on)", bond->virial, test_config.init_stress, epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(bond->energy, test_config.init_energy, epsilon);
    if (print_stats) std::cerr << "init_energy stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_lammps(lmp);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    EXPECT_FORCES("run_forces (newton on)", lmp->atom, test_config.run_forces, 10 * epsilon);
    EXPECT_STRESS("run_stress (newton on)", bond->virial, test_config.run_stress, epsilon);

    stats.reset();
    int id        = lmp->modify->find_compute("sum");
    double energy = lmp->modify->compute[id]->compute_scalar();
    EXPECT_FP_LE_WITH_EPS(bond->energy, test_config.run_energy, epsilon);
    EXPECT_FP_LE_WITH_EPS(bond->energy, energy, epsilon);
    if (print_stats) std::cerr << "run_energy  stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    lmp = init_lammps(args, test_config, false);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // skip over these tests if newton bond is forced to be on
    if (lmp->force->newton_bond == 0) {
        bond = lmp->force->bond;

        EXPECT_FORCES("init_forces (newton off)", lmp->atom, test_config.init_forces, epsilon);
        EXPECT_STRESS("init_stress (newton off)", bond->virial, test_config.init_stress,
                      2 * epsilon);

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(bond->energy, test_config.init_energy, epsilon);
        if (print_stats) std::cerr << "init_energy stats, newton off:" << stats << std::endl;

        if (!verbose) ::testing::internal::CaptureStdout();
        run_lammps(lmp);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        EXPECT_FORCES("run_forces (newton off)", lmp->atom, test_config.run_forces, 10 * epsilon);
        EXPECT_STRESS("run_stress (newton off)", bond->virial, test_config.run_stress, epsilon);

        stats.reset();
        id     = lmp->modify->find_compute("sum");
        energy = lmp->modify->compute[id]->compute_scalar();
        EXPECT_FP_LE_WITH_EPS(bond->energy, test_config.run_energy, epsilon);
        EXPECT_FP_LE_WITH_EPS(bond->energy, energy, epsilon);
        if (print_stats) std::cerr << "run_energy  stats, newton off:" << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    restart_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    bond = lmp->force->bond;

    EXPECT_FORCES("restart_forces", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_STRESS("restart_stress", bond->virial, test_config.init_stress, epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(bond->energy, test_config.init_energy, epsilon);
    if (print_stats) std::cerr << "restart_energy stats:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    data_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    bond = lmp->force->bond;

    EXPECT_FORCES("data_forces", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_STRESS("data_stress", bond->virial, test_config.init_stress, epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(bond->energy, test_config.init_energy, epsilon);
    if (print_stats) std::cerr << "data_energy stats:" << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

TEST(BondStyle, omp)
{
    if (!LAMMPS::is_installed_pkg("OPENMP")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"BondStyle", "-log", "none", "-echo", "screen", "-nocite",
                         "-pk",       "omp",  "4",    "-sf",   "omp"};

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = init_lammps(args, test_config, true);

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

    ErrorStats stats;
    auto bond = lmp->force->bond;

    EXPECT_FORCES("init_forces (newton on)", lmp->atom, test_config.init_forces, epsilon);
    EXPECT_STRESS("init_stress (newton on)", bond->virial, test_config.init_stress, 10 * epsilon);

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(bond->energy, test_config.init_energy, epsilon);
    if (print_stats) std::cerr << "init_energy stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_lammps(lmp);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    EXPECT_FORCES("run_forces (newton on)", lmp->atom, test_config.run_forces, 10 * epsilon);
    EXPECT_STRESS("run_stress (newton on)", bond->virial, test_config.run_stress, 10 * epsilon);

    stats.reset();
    int id        = lmp->modify->find_compute("sum");
    double energy = lmp->modify->compute[id]->compute_scalar();
    EXPECT_FP_LE_WITH_EPS(bond->energy, test_config.run_energy, epsilon);
    // TODO: this is currently broken for OPENMP with bond style hybrid
    // needs to be fixed in the main code somewhere. Not sure where, though.
    if (test_config.bond_style.substr(0, 6) != "hybrid")
        EXPECT_FP_LE_WITH_EPS(bond->energy, energy, epsilon);
    if (print_stats) std::cerr << "run_energy  stats, newton on: " << stats << std::endl;

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    lmp = init_lammps(args, test_config, false);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    // skip over these tests if newton bond is forced to be on
    if (lmp->force->newton_bond == 0) {
        bond = lmp->force->bond;

        EXPECT_FORCES("init_forces (newton off)", lmp->atom, test_config.init_forces, epsilon);
        EXPECT_STRESS("init_stress (newton off)", bond->virial, test_config.init_stress,
                      10 * epsilon);

        stats.reset();
        EXPECT_FP_LE_WITH_EPS(bond->energy, test_config.init_energy, epsilon);
        if (print_stats) std::cerr << "init_energy stats, newton off:" << stats << std::endl;

        if (!verbose) ::testing::internal::CaptureStdout();
        run_lammps(lmp);
        if (!verbose) ::testing::internal::GetCapturedStdout();

        EXPECT_FORCES("run_forces (newton off)", lmp->atom, test_config.run_forces, 10 * epsilon);
        EXPECT_STRESS("run_stress (newton off)", bond->virial, test_config.run_stress,
                      10 * epsilon);

        stats.reset();
        id     = lmp->modify->find_compute("sum");
        energy = lmp->modify->compute[id]->compute_scalar();
        EXPECT_FP_LE_WITH_EPS(bond->energy, test_config.run_energy, epsilon);
        // TODO: this is currently broken for OPENMP with bond style hybrid
        // needs to be fixed in the main code somewhere. Not sure where, though.
        if (test_config.bond_style.substr(0, 6) != "hybrid")
            EXPECT_FP_LE_WITH_EPS(bond->energy, energy, epsilon);
        if (print_stats) std::cerr << "run_energy  stats, newton off:" << stats << std::endl;
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
};


TEST(BondStyle, numdiff)
{
    if (!LAMMPS::is_installed_pkg("EXTRA-FIX")) GTEST_SKIP();
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"BondStyle", "-log", "none", "-echo", "screen", "-nocite"};

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

    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms, nlocal);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("fix diff all numdiff 2 6.05504e-6");
    lmp->input->one("run 2 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    Fix *ifix = lmp->modify->get_fix_by_id("diff");
    if (ifix) {
        double epsilon = test_config.epsilon * 5.0e8;
        ErrorStats stats;
        double **f1 = lmp->atom->f;
        double **f2 = ifix->array_atom;
        SCOPED_TRACE("EXPECT FORCES: numdiff");
        for (int i = 0; i < nlocal; ++i) {
            EXPECT_FP_LE_WITH_EPS(f1[i][0], f2[i][0], epsilon);
            EXPECT_FP_LE_WITH_EPS(f1[i][1], f2[i][1], epsilon);
            EXPECT_FP_LE_WITH_EPS(f1[i][2], f2[i][2], epsilon);
        }
        if (print_stats)
            std::cerr << "numdiff  stats: " << stats << " epsilon: " << epsilon << std::endl;
    }
    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
}

TEST(BondStyle, single)
{
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"BondStyle", "-log", "none", "-echo", "screen", "-nocite"};

    // create a LAMMPS instance with standard settings to detect the number of atom types
    if (!verbose) ::testing::internal::CaptureStdout();
    LAMMPS *lmp = init_lammps(args, test_config);
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
    int nbondtypes = lmp->atom->nbondtypes;
    int molecular  = lmp->atom->molecular;
    if (molecular != Atom::MOLECULAR) {
        std::cerr << "Only simple molecular atom styles are supported\n";
        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        GTEST_SKIP();
    }

    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
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

    command(fmt::format("create_box 1 box bond/types {} "
                        "extra/bond/per/atom 2 extra/special/per/atom 2",
                        nbondtypes));

    command("pair_style zero 8.0");
    command("pair_coeff * *");

    command("bond_style " + test_config.bond_style);
    auto bond = lmp->force->bond;

    for (auto &bond_coeff : test_config.bond_coeff) {
        command("bond_coeff " + bond_coeff);
    }

    // create (only) four atoms and two bonds
    command("mass * 1.0");
    command("create_atoms 1 single  5.0 -0.75  0.4 units box");
    command("create_atoms 1 single  5.5  0.25 -0.1 units box");
    command("create_atoms 1 single -5.0  0.75  0.4 units box");
    command("create_atoms 1 single -5.5 -0.25 -0.1 units box");
    command("create_bonds single/bond 1 1 2");
    command("create_bonds single/bond 2 3 4");

    for (auto &post_command : test_config.post_commands) {
        command(post_command);
    }

    command("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    int idx1       = lmp->atom->map(1);
    int idx2       = lmp->atom->map(2);
    int idx3       = lmp->atom->map(3);
    int idx4       = lmp->atom->map(4);
    double epsilon = test_config.epsilon;
    double **f     = lmp->atom->f;
    double **x     = lmp->atom->x;
    double delx1   = x[idx2][0] - x[idx1][0];
    double dely1   = x[idx2][1] - x[idx1][1];
    double delz1   = x[idx2][2] - x[idx1][2];
    double rsq1    = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
    double delx2   = x[idx4][0] - x[idx3][0];
    double dely2   = x[idx4][1] - x[idx3][1];
    double delz2   = x[idx4][2] - x[idx3][2];
    double rsq2    = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
    double fsingle = 0.0;
    double ebond[4], esngl[4];
    ErrorStats stats;

    ebond[0] = bond->energy;
    esngl[0] = bond->single(1, rsq1, idx1, idx2, fsingle);
    EXPECT_FP_LE_WITH_EPS(f[idx1][0], -fsingle * delx1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx1][1], -fsingle * dely1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx1][2], -fsingle * delz1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx2][0], fsingle * delx1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx2][1], fsingle * dely1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx2][2], fsingle * delz1, epsilon);

    esngl[0] += bond->single(2, rsq2, idx3, idx4, fsingle);
    EXPECT_FP_LE_WITH_EPS(f[idx3][0], -fsingle * delx2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx3][1], -fsingle * dely2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx3][2], -fsingle * delz2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx4][0], fsingle * delx2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx4][1], fsingle * dely2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx4][2], fsingle * delz2, epsilon);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("displace_atoms all random 0.5 0.5 0.5 23456");
    command("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    f       = lmp->atom->f;
    x       = lmp->atom->x;
    idx1    = lmp->atom->map(1);
    idx2    = lmp->atom->map(2);
    idx3    = lmp->atom->map(3);
    idx4    = lmp->atom->map(4);
    delx1   = x[idx2][0] - x[idx1][0];
    dely1   = x[idx2][1] - x[idx1][1];
    delz1   = x[idx2][2] - x[idx1][2];
    rsq1    = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
    delx2   = x[idx4][0] - x[idx3][0];
    dely2   = x[idx4][1] - x[idx3][1];
    delz2   = x[idx4][2] - x[idx3][2];
    rsq2    = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
    fsingle = 0.0;

    ebond[1] = bond->energy;
    esngl[1] = bond->single(1, rsq1, idx1, idx2, fsingle);
    EXPECT_FP_LE_WITH_EPS(f[idx1][0], -fsingle * delx1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx1][1], -fsingle * dely1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx1][2], -fsingle * delz1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx2][0], fsingle * delx1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx2][1], fsingle * dely1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx2][2], fsingle * delz1, epsilon);

    esngl[1] += bond->single(2, rsq2, idx3, idx4, fsingle);
    EXPECT_FP_LE_WITH_EPS(f[idx3][0], -fsingle * delx2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx3][1], -fsingle * dely2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx3][2], -fsingle * delz2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx4][0], fsingle * delx2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx4][1], fsingle * dely2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx4][2], fsingle * delz2, epsilon);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("displace_atoms all random 0.5 0.5 0.5 456963");
    command("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    f       = lmp->atom->f;
    x       = lmp->atom->x;
    idx1    = lmp->atom->map(1);
    idx2    = lmp->atom->map(2);
    idx3    = lmp->atom->map(3);
    idx4    = lmp->atom->map(4);
    delx1   = x[idx2][0] - x[idx1][0];
    dely1   = x[idx2][1] - x[idx1][1];
    delz1   = x[idx2][2] - x[idx1][2];
    rsq1    = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
    delx2   = x[idx4][0] - x[idx3][0];
    dely2   = x[idx4][1] - x[idx3][1];
    delz2   = x[idx4][2] - x[idx3][2];
    rsq2    = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
    fsingle = 0.0;

    ebond[2] = bond->energy;
    esngl[2] = bond->single(1, rsq1, idx1, idx2, fsingle);
    EXPECT_FP_LE_WITH_EPS(f[idx1][0], -fsingle * delx1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx1][1], -fsingle * dely1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx1][2], -fsingle * delz1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx2][0], fsingle * delx1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx2][1], fsingle * dely1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx2][2], fsingle * delz1, epsilon);

    esngl[2] += bond->single(2, rsq2, idx3, idx4, fsingle);
    EXPECT_FP_LE_WITH_EPS(f[idx3][0], -fsingle * delx2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx3][1], -fsingle * dely2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx3][2], -fsingle * delz2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx4][0], fsingle * delx2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx4][1], fsingle * dely2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx4][2], fsingle * delz2, epsilon);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("displace_atoms all random 0.5 0.5 0.5 9726532");
    command("run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    f       = lmp->atom->f;
    x       = lmp->atom->x;
    idx1    = lmp->atom->map(1);
    idx2    = lmp->atom->map(2);
    idx3    = lmp->atom->map(3);
    idx4    = lmp->atom->map(4);
    delx1   = x[idx2][0] - x[idx1][0];
    dely1   = x[idx2][1] - x[idx1][1];
    delz1   = x[idx2][2] - x[idx1][2];
    rsq1    = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
    delx2   = x[idx4][0] - x[idx3][0];
    dely2   = x[idx4][1] - x[idx3][1];
    delz2   = x[idx4][2] - x[idx3][2];
    rsq2    = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
    fsingle = 0.0;

    ebond[3] = bond->energy;
    esngl[3] = bond->single(1, rsq1, idx1, idx2, fsingle);
    EXPECT_FP_LE_WITH_EPS(f[idx1][0], -fsingle * delx1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx1][1], -fsingle * dely1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx1][2], -fsingle * delz1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx2][0], fsingle * delx1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx2][1], fsingle * dely1, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx2][2], fsingle * delz1, epsilon);

    esngl[3] += bond->single(2, rsq2, idx3, idx4, fsingle);
    EXPECT_FP_LE_WITH_EPS(f[idx3][0], -fsingle * delx2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx3][1], -fsingle * dely2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx3][2], -fsingle * delz2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx4][0], fsingle * delx2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx4][1], fsingle * dely2, epsilon);
    EXPECT_FP_LE_WITH_EPS(f[idx4][2], fsingle * delz2, epsilon);
    if (print_stats) std::cerr << "single_force  stats:" << stats << std::endl;

    stats.reset();
    EXPECT_FP_LE_WITH_EPS(ebond[0], esngl[0], epsilon);
    EXPECT_FP_LE_WITH_EPS(ebond[1], esngl[1], epsilon);
    EXPECT_FP_LE_WITH_EPS(ebond[2], esngl[2], epsilon);
    EXPECT_FP_LE_WITH_EPS(ebond[3], esngl[3], epsilon);
    if (print_stats) std::cerr << "single_energy  stats:" << stats << std::endl;

    int i = 0;
    for (auto &dist : test_config.equilibrium)
        EXPECT_NEAR(dist, bond->equilibrium_distance(++i), 0.00001);

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
}

TEST(BondStyle, extract)
{
    if (test_config.skip_tests.count(test_info_->name())) GTEST_SKIP();

    LAMMPS::argv args = {"BondStyle", "-log", "none", "-echo", "screen", "-nocite"};

    if (!verbose) ::testing::internal::CaptureStdout();
    LAMMPS *lmp = init_lammps(args, test_config, true);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    if (!lmp) {
        std::cerr << "One or more prerequisite styles are not available "
                     "in this LAMMPS configuration:\n";
        for (auto prerequisite : test_config.prerequisites) {
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        }
        GTEST_SKIP();
    }

    auto bond = lmp->force->bond;
    void *ptr = nullptr;
    int dim   = 0;
    for (auto extract : test_config.extract) {
        ptr = bond->extract(extract.first.c_str(), dim);
        EXPECT_NE(ptr, nullptr);
        EXPECT_EQ(dim, extract.second);
    }
    ptr = bond->extract("does_not_exist", dim);
    EXPECT_EQ(ptr, nullptr);

    for (int i = 1; i <= lmp->atom->nbondtypes; ++i)
        EXPECT_GE(bond->equilibrium_distance(i), 0.0);

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
}
