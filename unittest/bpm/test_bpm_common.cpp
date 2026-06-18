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

// Shared implementation for the BPM pair+bond verification test drivers.  Every
// driver builds its system entirely from the YAML file: a 'variables' block
// provides ${var} substitution, 'pre_commands' create the geometry and the
// per-atom storage, 'pair_style'/'pair_coeff' set the contact model,
// 'bond_style'/'bond_coeff' set the bond model (the system under test), and
// 'post_commands' create the bonds (after the pair style and bond style exist)
// and add the integrator.  The trajectory is run in segments ('run_segments')
// and per-atom positions, velocities and forces are compared against the
// recorded reference after each segment, together with any analytic model.

#include "test_bpm_common.h"

#include "test_analytic_models.h"
#include "test_config.h"
#include "test_main.h"
#include "../yaml_writer.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "atom.h"
#include "force.h"
#include "info.h"
#include "input.h"

#include "fmt/format.h"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

using ::testing::HasSubstr;
using ::testing::StartsWith;

using namespace LAMMPS_NS;

static void cleanup_lammps(LAMMPS *&lmp, const TestConfig &cfg)
{
    (void) cfg;
    delete lmp;
    lmp = nullptr;
}

static LAMMPS *init_lammps(LAMMPS::argv &args, const TestConfig &cfg, const bool newton)
{
    LAMMPS *lmp = new LAMMPS(args, MPI_COMM_WORLD);

    // check if prerequisite styles are available
    Info *info = new Info(lmp);
    int nfail  = 0;
    for (const auto &prerequisite : cfg.prerequisites) {
        if (!info->has_style(prerequisite.first, prerequisite.second)) ++nfail;
    }
    delete info;
    if (nfail > 0) {
        cleanup_lammps(lmp, cfg);
        return nullptr;
    }

    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    // BPM bond styles mandate newton bond off; the pair newton flag follows the
    // requested setting.  newton must precede box creation in pre_commands.
    command(newton ? "newton on off" : "newton off off");
    command("variable input_dir index " + INPUT_FOLDER);

    for (const auto &var : cfg.variables)
        command("variable " + var.first + " index " + var.second);

    // geometry, atom storage (fix property/atom), and the initial special_bonds
    for (const auto &pre_command : cfg.pre_commands)
        command(pre_command);

    // contact pair style
    command("pair_style " + cfg.pair_style);
    for (const auto &pair_coeff : cfg.pair_coeff)
        command("pair_coeff " + pair_coeff);

    // create the bonds (under the build-time special weights) and restore the
    // run-time special weights; the integrator and any loading are added here too
    for (const auto &post_command : cfg.post_commands)
        command(post_command);

    // bond style (the system under test) -- defined after create_bonds so the BPM
    // run-time special-bond weights (coul 1,1,1) are already in effect
    command("bond_style " + cfg.bond_style);
    for (const auto &bond_coeff : cfg.bond_coeff)
        command("bond_coeff " + bond_coeff);

    command("run 0 post no");
    return lmp;
}

// run all segments and compare per-atom state against the reference after each.
// lmp is taken by reference so the optional restart round trip can swap in a
// fresh instance rebuilt from a restart file; the continued trajectory must
// still match the (uninterrupted) reference, which exercises the restart of the
// per-type coefficients, the property/atom fields and the bond history.
static void run_and_check(LAMMPS *&lmp, LAMMPS::argv &args, const TestConfig &cfg, double epsilon,
                          const std::string &label)
{
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    for (std::size_t i = 0; i < cfg.run_segments.size(); ++i) {
        command("run " + std::to_string(cfg.run_segments[i]) + " post no");
        const std::string tag = label + ", seg " + std::to_string(i);

        if (i < cfg.seg_pos.size())
            EXPECT_POSITIONS("run_pos (" + tag + ")", lmp->atom, cfg.seg_pos[i], epsilon);
        if (i < cfg.seg_vel.size())
            EXPECT_VELOCITIES("run_vel (" + tag + ")", lmp->atom, cfg.seg_vel[i], epsilon);
        if (i < cfg.seg_force.size())
            EXPECT_FORCES("run_force (" + tag + ")", lmp->atom, cfg.seg_force[i], epsilon);

        check_analytic_model(cfg, lmp, (int) i);

        // optional restart round trip: write a restart, rebuild from it, and
        // continue.  restart_commands restore the non-restart setup (lattice,
        // integrator) that read_restart does not.  command() uses lmp by
        // reference, so it targets the rebuilt instance after the swap.
        if ((cfg.restart_segment >= 0) && ((int) i == cfg.restart_segment)) {
            const std::string rfile = "bpm_restart_" + cfg.basename + ".tmp";
            command("write_restart " + rfile);
            delete lmp;
            lmp = new LAMMPS(args, MPI_COMM_WORLD);
            command("read_restart " + rfile);
            for (const auto &rc : cfg.restart_commands) command(rc);
            command("run 0 post no");
        }
    }
}

// re-generate yaml file with current settings.

void generate_yaml_file(const char *outfile, const TestConfig &config)
{
    LAMMPS::argv args = {"BPM", "-log", "none", "-echo", "screen", "-nocite"};
    LAMMPS *lmp       = nullptr;
    try {
        lmp = init_lammps(args, config, false);
    } catch (std::exception &e) {
        FAIL() << e.what();
    }
    if (!lmp) {
        std::cerr << "One or more prerequisite styles are not available "
                     "in this LAMMPS configuration:\n";
        for (auto prerequisite : config.prerequisites)
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        return;
    }

    const int natoms = lmp->atom->natoms;
    YamlWriter writer(outfile);

    write_yaml_header(&writer, &test_config, lmp->version);
    writer.emit("natoms", natoms);

    std::string block;
    for (const auto &var : config.variables)
        block += var.first + " " + var.second + "\n";
    writer.emit_block("variables", block);

    writer.emit("pair_style", config.pair_style);
    block.clear();
    for (const auto &pair_coeff : config.pair_coeff)
        block += pair_coeff + "\n";
    writer.emit_block("pair_coeff", block);

    writer.emit("bond_style", config.bond_style);
    block.clear();
    for (const auto &bond_coeff : config.bond_coeff)
        block += bond_coeff + "\n";
    writer.emit_block("bond_coeff", block);

    block.clear();
    for (std::size_t i = 0; i < config.run_segments.size(); ++i) {
        if (i) block += " ";
        block += std::to_string(config.run_segments[i]);
    }
    writer.emit_block("run_segments", block);

    if (config.analytic_enable) {
        writer.emit("analytic_enable", std::string("yes"));
        writer.emit("analytic_model", config.analytic_model);
        writer.emit("analytic_tol", config.analytic_tol);
        writer.emit("analytic_segment", (long) config.analytic_segment);
        if (config.analytic_only) writer.emit("analytic_only", std::string("yes"));
    }

    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    if (config.analytic_only) {
        for (std::size_t i = 0; i < config.run_segments.size(); ++i)
            command("run " + std::to_string(config.run_segments[i]) + " post no");
        cleanup_lammps(lmp, config);
        return;
    }

    std::string pos_block, vel_block, force_block;

    // iterate over local atoms by tag; rows are keyed by (segment, tag), so the
    // emission order is irrelevant when read back.
    for (std::size_t i = 0; i < config.run_segments.size(); ++i) {
        command("run " + std::to_string(config.run_segments[i]) + " post no");
        auto *x         = lmp->atom->x;
        auto *v         = lmp->atom->v;
        auto *f         = lmp->atom->f;
        tagint *tag     = lmp->atom->tag;
        const int local = lmp->atom->nlocal;
        for (int j = 0; j < local; ++j) {
            const tagint id = tag[j];
            pos_block += fmt::format("{:3} {:3} {:23.16e} {:23.16e} {:23.16e}\n", i, id, x[j][0],
                                     x[j][1], x[j][2]);
            vel_block += fmt::format("{:3} {:3} {:23.16e} {:23.16e} {:23.16e}\n", i, id, v[j][0],
                                     v[j][1], v[j][2]);
            force_block += fmt::format("{:3} {:3} {:23.16e} {:23.16e} {:23.16e}\n", i, id, f[j][0],
                                       f[j][1], f[j][2]);
        }
    }
    writer.emit_block("run_pos", pos_block);
    writer.emit_block("run_vel", vel_block);
    writer.emit_block("run_force", force_block);

    cleanup_lammps(lmp, config);
}

void run_bpm_trajectory_test(bool newton, const std::string &label)
{
    LAMMPS::argv args = {"BPM", "-log", "none", "-echo", "screen", "-nocite"};

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp = nullptr;
    try {
        lmp = init_lammps(args, test_config, newton);
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
        for (auto &prerequisite : test_config.prerequisites)
            std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
        GTEST_SKIP();
    }

    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    EXPECT_THAT(output, HasSubstr("Loop time"));

    // abort if running in parallel and not all atoms are local
    const int nlocal = lmp->atom->nlocal;
    ASSERT_EQ(lmp->atom->natoms, nlocal);

    double epsilon = test_config.epsilon;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_and_check(lmp, args, test_config, epsilon, label);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
}
