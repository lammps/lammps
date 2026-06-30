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

// Shared implementation for the DEM verification test drivers
// (test_dem_01 ... test_dem_11).  Every driver builds its system entirely
// from the YAML file: a 'variables' block provides ${var} substitution,
// 'pre_commands' create the geometry, 'pair_style'/'pair_coeff' set the
// contact model, and 'post_commands' add the integrator, gravity and walls.
// The trajectory is run in segments ('run_segments') and per-atom positions,
// velocities, torques and angular velocities are compared against the
// recorded reference after each segment, together with any analytic model
// selected in the YAML.  Driver-specific behaviour lives in the YAML files
// and in test_analytic_models.cpp, so the drivers themselves only differ in
// the GoogleTest suite name.

#include "test_dem_common.h"

#include "test_analytic_models.h"
#include "test_config.h"
#include "test_main.h"
#include "yaml_writer.h"

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

    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    // newton must be set before the simulation box is created in pre_commands
    command(newton ? "newton on" : "newton off");
    command("variable input_dir index " + INPUT_FOLDER);

    // user-defined variables enable ${var} substitution in all command strings
    for (const auto &var : cfg.variables)
        command("variable " + var.first + " index " + var.second);

    // build the system geometry (units, atom_style, region, create_atoms, set, ...)
    for (const auto &pre_command : cfg.pre_commands)
        command(pre_command);

    // contact model
    command("pair_style " + cfg.pair_style);
    for (const auto &pair_coeff : cfg.pair_coeff)
        command("pair_coeff " + pair_coeff);

    // integrator, gravity, walls, ...
    for (const auto &post_command : cfg.post_commands)
        command(post_command);

    command("run 0 post no");
    return lmp;
}

// run all segments and compare per-atom state against the reference after each
static void run_and_check(LAMMPS *lmp, const TestConfig &cfg, double epsilon,
                          const std::string &label)
{
    auto command = [&](const std::string &line) {
        lmp->input->one(line);
    };

    const bool has_torque = lmp->atom->torque_flag;
    const bool has_omega  = lmp->atom->omega_flag;
    const bool has_angmom = lmp->atom->angmom_flag;

    for (std::size_t i = 0; i < cfg.run_segments.size(); ++i) {
        command("run " + std::to_string(cfg.run_segments[i]) + " post no");
        const std::string tag = label + ", seg " + std::to_string(i);

        if (i < cfg.seg_pos.size())
            EXPECT_POSITIONS("run_pos (" + tag + ")", lmp->atom, cfg.seg_pos[i], epsilon);
        if (i < cfg.seg_vel.size())
            EXPECT_VELOCITIES("run_vel (" + tag + ")", lmp->atom, cfg.seg_vel[i], epsilon);
        if (has_torque && (i < cfg.seg_torque.size()))
            EXPECT_TORQUES("run_torque (" + tag + ")", lmp->atom, cfg.seg_torque[i], epsilon);
        if (has_omega && (i < cfg.seg_omega.size()))
            EXPECT_OMEGA("run_omega (" + tag + ")", lmp->atom, cfg.seg_omega[i], epsilon);
        if (has_angmom && (i < cfg.seg_angmom.size()))
            EXPECT_ANGMOM("run_angmom (" + tag + ")", lmp->atom, cfg.seg_angmom[i], epsilon);

        check_analytic_model(cfg, lmp, (int) i);
    }
}

// re-generate yaml file with current settings.

void generate_yaml_file(const char *outfile, const TestConfig &config)
{
    LAMMPS::argv args = {"DEM", "-log", "none", "-echo", "screen", "-nocite"};
    LAMMPS *lmp       = nullptr;
    try {
        lmp = init_lammps(args, config, true);
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

    // write yaml header
    write_yaml_header(&writer, &test_config, lmp->version);

    // natoms
    writer.emit("natoms", natoms);

    // variables block (echo back verbatim)
    std::string block;
    for (const auto &var : config.variables)
        block += var.first + " " + var.second + "\n";
    writer.emit_block("variables", block);

    // pair style and coefficients
    writer.emit("pair_style", config.pair_style);
    block.clear();
    for (const auto &pair_coeff : config.pair_coeff)
        block += pair_coeff + "\n";
    writer.emit_block("pair_coeff", block);

    // run segments
    block.clear();
    for (std::size_t i = 0; i < config.run_segments.size(); ++i) {
        if (i) block += " ";
        block += std::to_string(config.run_segments[i]);
    }
    writer.emit_block("run_segments", block);

    // optional analytic check flags (echo back if enabled)
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

    // For analytic-only (chaotic bulk) tests -- e.g. angle of repose -- the
    // trajectory is not bit-reproducible across newton on/off or platforms, so
    // we record no per-atom reference.  Only the bulk analytic observable is
    // robust; it is checked live during the test (check_analytic_model).  We
    // still execute the run segments so the reference reflects a completed run.
    if (config.analytic_only) {
        for (std::size_t i = 0; i < config.run_segments.size(); ++i)
            command("run " + std::to_string(config.run_segments[i]) + " post no");
        cleanup_lammps(lmp, config);
        return;
    }

    const bool has_torque = lmp->atom->torque_flag;
    const bool has_omega  = lmp->atom->omega_flag;
    const bool has_angmom = lmp->atom->angmom_flag;
    std::string pos_block, vel_block, torque_block, omega_block, angmom_block;

    // iterate over local atoms by tag; granular/atomic systems have no atom map.
    // rows are keyed by (segment, tag), so the emission order is irrelevant.
    for (std::size_t i = 0; i < config.run_segments.size(); ++i) {
        command("run " + std::to_string(config.run_segments[i]) + " post no");
        auto *x         = lmp->atom->x;
        auto *v         = lmp->atom->v;
        auto *t         = lmp->atom->torque;
        auto *w         = lmp->atom->omega;
        auto *angmom    = lmp->atom->angmom;
        tagint *tag     = lmp->atom->tag;
        const int local = lmp->atom->nlocal;
        for (int j = 0; j < local; ++j) {
            const tagint id = tag[j];
            pos_block += fmt::format("{:3} {:3} {:23.16e} {:23.16e} {:23.16e}\n", i, id, x[j][0],
                                     x[j][1], x[j][2]);
            vel_block += fmt::format("{:3} {:3} {:23.16e} {:23.16e} {:23.16e}\n", i, id, v[j][0],
                                     v[j][1], v[j][2]);
            if (has_torque)
                torque_block += fmt::format("{:3} {:3} {:23.16e} {:23.16e} {:23.16e}\n", i, id,
                                            t[j][0], t[j][1], t[j][2]);
            if (has_omega)
                omega_block += fmt::format("{:3} {:3} {:23.16e} {:23.16e} {:23.16e}\n", i, id,
                                           w[j][0], w[j][1], w[j][2]);
            if (has_angmom)
                angmom_block += fmt::format("{:3} {:3} {:23.16e} {:23.16e} {:23.16e}\n", i, id,
                                            angmom[j][0], angmom[j][1], angmom[j][2]);
        }
    }
    writer.emit_block("run_pos", pos_block);
    writer.emit_block("run_vel", vel_block);
    if (has_torque) writer.emit_block("run_torque", torque_block);
    if (has_omega) writer.emit_block("run_omega", omega_block);
    if (has_angmom) writer.emit_block("run_angmom", angmom_block);

    cleanup_lammps(lmp, config);
}

void run_dem_trajectory_test(bool newton, const std::string &label)
{
    LAMMPS::argv args = {"DEM", "-log", "none", "-echo", "screen", "-nocite"};

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

    // when requesting newton off, only compare if the pair style honored it
    if (!newton && (lmp->force->newton_pair != 0)) {
        if (!verbose) ::testing::internal::CaptureStdout();
        cleanup_lammps(lmp, test_config);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        GTEST_SKIP() << "pair style forces newton pair on";
    }

    double epsilon = test_config.epsilon;

    if (!verbose) ::testing::internal::CaptureStdout();
    run_and_check(lmp, test_config, epsilon, label);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    if (!verbose) ::testing::internal::CaptureStdout();
    cleanup_lammps(lmp, test_config);
    if (!verbose) ::testing::internal::GetCapturedStdout();
}
