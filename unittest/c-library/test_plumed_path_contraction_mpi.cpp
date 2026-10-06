// Coordinate-map and conservative-force tests for partially contracted paths.
#define LAMMPS_LIB_MPI 1
#include "atom.h"
#include "fix.h"
#include "lammps.h"
#include "library.h"
#include "modify.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

#include "../testing/test_mpi_main.h"
#include "gtest/gtest.h"

namespace {
struct State {
    double energy         = 0.0;
    int ownership_changes = 0;
    std::array<double, 6> force{}, virial{}, position{};
};

double basis(int mode, int bead, int beads)
{
    const double norm = 1.0 / std::sqrt(double(beads));
    if (mode == 0) return norm;
    if (beads % 2 == 0 && mode == beads / 2) return (bead % 2 ? -1.0 : 1.0) * norm;
    const double phase = 2.0 * std::acos(-1.0) * mode * bead / beads;
    return std::sqrt(2.0) * norm * (mode <= beads / 2 ? std::cos(phase) : std::sin(phase));
}

State evaluate(int beads, double fraction, bool normal_modes, bool biased = true,
               int displaced = -1, int component = 0, double delta = 0.0, double strain = 0.0,
               bool crossing = false, bool centroid = false, int steps = 0, bool restart = false,
               double shear = 0.0, bool radial = false, bool opes = false, bool learn = false)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const int ranks_per_bead = size / beads;
    const int bead           = rank / ranks_per_bead;
    const std::string stem   = "test_contracted_" + std::to_string(size);
    if (rank == 0) {
        std::ofstream out(stem + ".dat");
        if (radial) {
            out << "d: DISTANCE ATOMS=1,2\n"
                << "s: CUSTOM ARG=d FUNC=x+x*x*x PERIODIC=NO\n"
                << "t: CUSTOM ARG=d FUNC=x*x PERIODIC=NO\n";
        } else {
            out << "d: DISTANCE ATOMS=1,2 COMPONENTS\n"
                << "s: CUSTOM ARG=d.x FUNC=x+x*x*x PERIODIC=NO\n"
                << "t: CUSTOM ARG=d.y,d.z FUNC=x*x+0.3*y PERIODIC=NO\n";
        }
        if (!centroid) out << "m: ENSEMBLE ARG=s,t\n";
        out << "v: CUSTOM ARG=" << (centroid ? "s,t" : "m.s,m.t")
            << " FUNC=0.5*(x-0.7)^2+0.3*x*y+0.4*y*y PERIODIC=NO\n";
        if (biased && opes) {
            out << "b: OPES_METAD ARG=m.s,m.t PACE=2 BARRIER=4 TEMP=1 SIGMA=0.5,0.5 FIXED_SIGMA "
                   "FMT=%24.17g FILE="
                << stem << ".kernels";
            if (learn)
                out << " STATE_WFILE=" << stem << ".state STATE_WSTRIDE=2 RESTART=NO\n";
            else
                out << " STATE_RFILE=" << stem << ".state UPDATE_UNTIL=0 RESTART=YES\n";
        } else if (biased)
            out << "b: BIASVALUE ARG=v\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);
    LAMMPS_NS::LAMMPS::argv args = {"test_contracted",
                                    "-screen",
                                    "none",
                                    "-log",
                                    "none",
                                    "-partition",
                                    std::to_string(beads) + "x" + std::to_string(ranks_per_bead),
                                    "-in",
                                    "none",
                                    "-nocite"};
    auto *lmp                    = new LAMMPS_NS::LAMMPS(args, MPI_COMM_WORLD);
    auto command                 = [&](const std::string &line) {
        lammps_command(lmp, line.c_str());
    };
    command("units lj");
    command("atom_style atomic");
    command("atom_modify map array");
    if (ranks_per_bead == 2) command("processors 2 1 1");
    const double scale = 1.0 + strain;
    std::ostringstream box;
    box << std::setprecision(17) << "region box " << (shear == 0.0 ? "block" : "prism") << " 0 "
        << 8 * scale << " 0 " << 8 * scale << " 0 " << 8 * scale;
    if (shear != 0.0) box << ' ' << 8 * scale * shear << " 0 0";
    command(box.str());
    command("create_box 1 box");
    const std::array<double, 3> anchor = {crossing ? 7.3 : 3.5, 4.0, 4.0};
    for (int atom = 0; atom < 2; ++atom) {
        std::array<double, 3> x = anchor;
        if (atom) {
            x[0] += 1.0 + 0.05 * bead;
            x[1] += 0.1 - 0.03 * bead;
            x[2] += 0.15 + 0.05 * bead;
            if (bead == displaced) x[component] += delta;
        }
        std::array<int, 3> images{};
        for (int d = 0; d < 3; ++d) {
            images[d] = int(std::floor(x[d] / 8.0));
            x[d]      = scale * (x[d] - 8 * images[d]);
        }
        x[0] += shear * x[1];
        std::ostringstream create;
        create << std::setprecision(17) << "create_atoms 1 single " << x[0] << ' ' << x[1] << ' '
               << x[2] << " units box";
        command(create.str());
        command("set atom " + std::to_string(atom + 1) + " image " + std::to_string(images[0]) +
                " " + std::to_string(images[1]) + " " + std::to_string(images[2]));
    }
    command("mass 1 1.0");
    command("pair_style lj/cut 2.5");
    command("pair_coeff * * 0.05 1.0");
    command("timestep 0.001");
    command("neigh_modify every 1 delay 0 check no");
    command(steps ? "velocity all set 200 0 0" : "velocity all set 0 0 0");
    const std::string pimd =
        std::string("fix fpimd all pimd/langevin method ") + (normal_modes ? "nmpimd" : "pimd") +
        " ensemble nvt integrator obabo thermostat PILE_L 2468 tau 1.0 temp 1.0 fixcom no";
    command(pimd);
    std::string fix = "fix bias all plumed plumedfile " + stem + ".dat outfile " + stem +
                      ".log path_integral " + (centroid ? "centroid" : "bead_mean") +
                      " pimd_fix fpimd";
    if (fraction >= 0.0) fix += " path_contraction " + std::to_string(fraction);
    command(fix);
    auto owners = [&]() {
        std::array<int, 2> local{}, global{};
        for (int i = 0; i < lmp->atom->nlocal; ++i)
            local[lmp->atom->tag[i] - 1] = rank + 1;
        MPI_Allreduce(local.data(), global.data(), 2, MPI_INT, MPI_MAX, lmp->world);
        return global;
    };
    const auto initial_owners = owners();
    command("run " + std::to_string(steps) + " post no");
    const auto final_owners = owners();
    if (restart) {
        const std::string checkpoint = stem + ".restart." + std::to_string(bead);
        command("write_restart " + checkpoint);
        command("clear");
        command("read_restart " + checkpoint);
        command(pimd);
        command(fix);
        command("run 0 post no");
    }
    State result;
    for (int i = 0; i < 2; ++i)
        result.ownership_changes += initial_owners[i] != final_owners[i];
    // The legacy C gather_atoms API rejects BIGBIG even for this two-atom
    // fixture. Gather the fixed-size output by tag with the native Atom API.
    auto gather = [&](double **values, std::array<double, 6> &output) {
        std::array<double, 6> local{};
        for (int i = 0; i < lmp->atom->nlocal; ++i)
            for (int d = 0; d < 3; ++d)
                local[3 * (lmp->atom->tag[i] - 1) + d] = values[i][d];
        MPI_Allreduce(local.data(), output.data(), 6, MPI_DOUBLE, MPI_SUM, lmp->world);
    };
    gather(lmp->atom->f, result.force);
    gather(lmp->atom->x, result.position);
    auto *plumed              = lmp->modify->get_fix_by_id("bias");
    const double contribution = rank % ranks_per_bead == 0 ? plumed->compute_scalar() : 0.0;
    MPI_Allreduce(&contribution, &result.energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(plumed->virial, result.virial.data(), 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(lammps_has_error(lmp), 0);
    delete lmp;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::remove((stem + ".dat").c_str());
        std::remove((stem + ".log").c_str());
        for (int b = 0; b < beads; ++b) {
            std::remove((stem + ".log." + std::to_string(b)).c_str());
            std::remove((stem + ".restart." + std::to_string(b)).c_str());
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return result;
}
} // namespace

TEST(MPI, contracted_path_endpoints_and_physical_force_preservation)
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int beads : {2, size}) {
        for (bool normal : {false, true}) {
            for (bool crossing : {false, true}) {
                const auto native_mean = evaluate(beads, -1, normal, true, -1, 0, 0, 0, crossing);
                const auto native_centroid =
                    evaluate(beads, -1, normal, true, -1, 0, 0, 0, crossing, true);
                for (double fraction : {0.0, 1.0}) {
                    const auto mapped =
                        evaluate(beads, fraction, normal, true, -1, 0, 0, 0, crossing);
                    const auto &reference = fraction == 0 ? native_centroid : native_mean;
                    EXPECT_NEAR(mapped.energy, reference.energy, 1e-11);
                    for (int d = 0; d < 6; ++d) {
                        EXPECT_NEAR(mapped.force[d], reference.force[d], 2e-10);
                        EXPECT_NEAR(mapped.position[d], reference.position[d], 1e-12);
                        EXPECT_NEAR(mapped.virial[d], reference.virial[d], 2e-10);
                    }
                }
                const auto zero        = evaluate(beads, -1, normal, false, -1, 0, 0, 0, crossing);
                const auto mapped_zero = evaluate(beads, 0.5, normal, false, -1, 0, 0, 0, crossing);
                for (int d = 0; d < 6; ++d)
                    EXPECT_NEAR(zero.force[d], mapped_zero.force[d], 1e-12);
            }
        }
    }
}

TEST(MPI, contracted_path_cartesian_derivatives_and_strain)
{
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Odd P with multiple spatial ranks when launched with six ranks.
    const int beads = size / 2;
    const int bead  = rank / 2;
    for (bool normal : {false, true}) {
        const auto state = evaluate(beads, 0.5, normal);
        const auto zero  = evaluate(beads, 0.5, normal, false);
        // Closed-form moments for the arithmetic-progression path. This
        // checks the force independently of subtracting nearby bias energies.
        const double cx  = 1.0 + 0.025 * (beads - 1);
        const double cy  = 0.1 - 0.015 * (beads - 1);
        const double cz  = 0.15 + 0.025 * (beads - 1);
        const double mx2 = 0.0025 * (beads * beads - 1) / 12.0;
        const double my2 = 0.0009 * (beads * beads - 1) / 12.0;
        const double cv0 = cx + cx * cx * cx + 0.75 * cx * mx2;
        const double cv1 = cy * cy + 0.25 * my2 + 0.3 * cz;
        const double g0 = cv0 - 0.7 + 0.3 * cv1, g1 = 0.3 * cv0 + 0.8 * cv1;
        EXPECT_NEAR(state.energy,
                    0.5 * (cv0 - 0.7) * (cv0 - 0.7) + 0.3 * cv0 * cv1 + 0.4 * cv1 * cv1, 1e-11);
        std::array<double, 3> exact{};
        for (int b = 0; b < beads; ++b) {
            const double u                    = 1.0 + 0.05 * b - cx;
            const std::array<double, 3> force = {
                -g0 * (1 + 3 * cx * cx + 0.75 * mx2 + 1.5 * cx * u + 0.375 * (u * u - mx2)),
                -g1 * (2 * cy + 0.5 * (0.1 - 0.03 * b - cy)), -0.3 * g1};
            const double coefficient = normal ? basis(bead, b, beads) : (bead == b);
            for (int d = 0; d < 3; ++d)
                exact[d] += coefficient * force[d];
        }
        for (int d = 0; d < 3; ++d) {
            EXPECT_NEAR(state.force[3 + d] - zero.force[3 + d], exact[d], 2e-11);
            EXPECT_NEAR(state.force[d] - zero.force[d], -exact[d], 2e-11);
        }
        for (double step : {1e-4, 1e-5, 1e-6}) {
            SCOPED_TRACE("normal=" + std::to_string(normal) + " step=" + std::to_string(step));
            for (int component = 0; component < 3; ++component) {
                SCOPED_TRACE("component=" + std::to_string(component));
                double expected = 0.0;
                for (int displaced = 0; displaced < beads; ++displaced) {
                    const auto plus =
                        evaluate(beads, 0.5, normal, true, displaced, component, step);
                    const auto minus =
                        evaluate(beads, 0.5, normal, true, displaced, component, -step);
                    const double force = -beads * (plus.energy - minus.energy) / (2 * step);
                    expected +=
                        force * (normal ? basis(bead, displaced, beads) : (bead == displaced));
                }
                if (rank % 2 == 0)
                    std::cout << std::setprecision(17) << "CONTRACTION_FD normal=" << normal
                              << " bead=" << bead << " component=" << component << " step=" << step
                              << " force=" << state.force[3 + component] - zero.force[3 + component]
                              << " fd=" << expected << std::endl;
                // Require agreement on two consecutive predeclared stencils.
                // The smallest stencil is retained as a roundoff diagnostic:
                // energy subtraction can amplify sub-picounit evaluation error.
                if (step >= 1e-5) {
                    EXPECT_NEAR(state.force[3 + component] - zero.force[3 + component], expected,
                                2e-7);
                    EXPECT_NEAR(state.force[component] - zero.force[component], -expected, 2e-7);
                }
                EXPECT_TRUE(std::isfinite(expected));
            }
            const auto plus       = evaluate(beads, 0.5, normal, true, -1, 0, 0, step);
            const auto minus      = evaluate(beads, 0.5, normal, true, -1, 0, 0, -step);
            const double expected = -beads * (plus.energy - minus.energy) / (2 * step);
            EXPECT_NEAR(state.virial[0] + state.virial[1] + state.virial[2], expected, 2e-6);
        }
    }
}

TEST(MPI, contracted_path_migration_and_restart)
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const int beads = size / 2;
    for (bool normal : {false, true}) {
        for (double fraction : {0.0, 0.5, 1.0}) {
            const auto before =
                evaluate(beads, fraction, normal, true, -1, 0, 0, 0, false, false, 4);
            const auto resumed =
                evaluate(beads, fraction, normal, true, -1, 0, 0, 0, false, false, 4, true);
            int local_changes = before.ownership_changes, total_changes = 0;
            MPI_Allreduce(&local_changes, &total_changes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            // A moving atom may cross a spatial boundary in only some beads.
            EXPECT_GT(total_changes, 0);
            EXPECT_NEAR(before.energy, resumed.energy, 1e-11);
            for (int d = 0; d < 6; ++d) {
                EXPECT_NEAR(before.position[d], resumed.position[d], 1e-12);
                EXPECT_NEAR(before.force[d], resumed.force[d], 2e-9);
            }
            if (fraction != 0.5) {
                const auto control =
                    evaluate(beads, -1, normal, true, -1, 0, 0, 0, false, fraction == 0, 4);
                EXPECT_NEAR(before.energy, control.energy, 1e-10);
                for (int d = 0; d < 6; ++d) {
                    EXPECT_NEAR(before.position[d], control.position[d], 1e-11);
                    EXPECT_NEAR(before.force[d], control.force[d], 2e-9);
                }
            }
        }
    }
}

TEST(MPI, contracted_path_shear_derivative)
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const int beads = size / 2;
    // LAMMPS reports six tensor components. A rotationally invariant field
    // has symmetric virial, making xy equal to the xy-strain derivative.
    // Cartesian-component fields can have a nonsymmetric nine-component
    // tensor and must not be tested with the transposed six-component entry.
    const double shear = 0.07;
    const auto state =
        evaluate(beads, 0.5, false, true, -1, 0, 0, 0, false, false, 0, false, shear, true);
    for (double step : {1e-4, 1e-5}) {
        const auto plus  = evaluate(beads, 0.5, false, true, -1, 0, 0, 0, false, false, 0, false,
                                    shear + step, true);
        const auto minus = evaluate(beads, 0.5, false, true, -1, 0, 0, 0, false, false, 0, false,
                                    shear - step, true);
        EXPECT_NEAR(state.virial[3], -beads * (plus.energy - minus.energy) / (2 * step), 2e-7);
    }
}

TEST(MPI, contracted_path_input_contract)
{
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const std::string stem = "test_contraction_guard_" + std::to_string(size);
    if (rank == 0) {
        std::ofstream out(stem + ".dat");
        out << "d: DISTANCE ATOMS=1,2\nm: ENSEMBLE ARG=d\n"
            << "b: RESTRAINT ARG=m.d AT=1 KAPPA=1\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int scenario = 0; scenario < 6; ++scenario) {
        SCOPED_TRACE("scenario=" + std::to_string(scenario));
        LAMMPS_NS::LAMMPS::argv args = {"test_guard",
                                        "-screen",
                                        "none",
                                        "-log",
                                        "none",
                                        "-partition",
                                        std::to_string(size) + "x1",
                                        "-in",
                                        "none",
                                        "-nocite"};
        auto *lmp                    = new LAMMPS_NS::LAMMPS(args, MPI_COMM_WORLD);
        auto command                 = [&](const std::string &line) {
            lammps_command(lmp, line.c_str());
        };
        lammps_set_show_error(lmp, 0);
        command("units lj");
        command("atom_style atomic");
        command("atom_modify map array");
        command("region box block 0 8 0 8 0 8");
        command("create_box 1 box");
        command("create_atoms 1 single 3 4 4");
        command("create_atoms 1 single 4 4 4");
        command("mass 1 1");
        command("pair_style zero 2.5");
        command("pair_coeff * *");
        std::string pimd = "fix fpimd all pimd/langevin method nmpimd ensemble ";
        pimd += scenario == 4 ? "npt" : "nvt";
        pimd += " integrator obabo thermostat PILE_L 2468 tau 1 temp 1 fixcom no";
        if (scenario == 4) pimd += " iso 1 barostat BZP taup 1";
        command(pimd);
        if (scenario == 5) command("fix moving all deform 1 x erate 0.01 remap x");
        std::string fix = "fix bias all plumed plumedfile " + stem + ".dat outfile " + stem +
                          ".log path_integral ";
        fix += scenario == 2 ? "centroid" : "bead_mean";
        fix += " pimd_fix fpimd path_contraction ";
        fix += scenario == 0                ? "-0.1"
               : scenario == 1              ? "1.1"
               : scenario == 3 && rank == 0 ? "0.25"
                                            : "0.5";
        command(fix);
        if (scenario >= 3) command("run 0 post no");
        const int failed = lammps_has_error(lmp);
        EXPECT_EQ(failed, 1) << "scenario=" << scenario;
        if (failed) {
            char message[1024];
            lammps_get_last_error_message(lmp, message, sizeof(message));
            EXPECT_NE(std::string(message).find("path_contraction"), std::string::npos);
        }
        delete lmp;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (rank == 0) {
        std::remove((stem + ".dat").c_str());
        std::remove((stem + ".log").c_str());
        for (int b = 0; b < size; ++b)
            std::remove((stem + ".log." + std::to_string(b)).c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

TEST(MPI, contracted_path_frozen_opes_derivatives_and_state)
{
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int beads = size / 2, bead = rank / 2;
    const std::string stem = "test_contracted_" + std::to_string(size);
    auto contents          = [](const std::string &file) {
        std::ifstream in(file);
        return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    };
    evaluate(beads, 0.5, false, true, -1, 0, 0, 0, false, false, 6, false, 0, false, true, true);
    const std::string original = contents(stem + ".state." + std::to_string(bead));
    EXPECT_FALSE(original.empty());
    EXPECT_EQ(original, contents(stem + ".state.0"));
    int counter = -1;
    std::istringstream rows(original);
    std::string line;
    while (std::getline(rows, line)) {
        if (line.rfind("#! SET counter", 0) != 0) continue;
        std::istringstream row(line);
        std::string marker, set, key;
        row >> marker >> set >> key >> counter;
    }
    EXPECT_EQ(counter, 4) << original;
    const auto state =
        evaluate(beads, 0.5, false, true, -1, 0, 0, 0, false, false, 0, false, 0, false, true);
    const auto zero = evaluate(beads, 0.5, false, false);
    EXPECT_GT(std::abs(state.energy), 1e-6);
    for (double step : {1e-4, 1e-5}) {
        for (int b = 0; b < beads; ++b) {
            for (int d = 0; d < 3; ++d) {
                const auto plus  = evaluate(beads, 0.5, false, true, b, d, step, 0, false, false, 0,
                                            false, 0, false, true);
                const auto minus = evaluate(beads, 0.5, false, true, b, d, -step, 0, false, false,
                                            0, false, 0, false, true);
                if (b == bead)
                    EXPECT_NEAR(state.force[3 + d] - zero.force[3 + d],
                                -beads * (plus.energy - minus.energy) / (2 * step), 2e-6);
            }
        }
    }
    EXPECT_EQ(original, contents(stem + ".state." + std::to_string(bead)));
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0 && !::testing::Test::HasFailure()) {
        for (int b = 0; b < beads; ++b) {
            for (const std::string &ext : {".state.", ".kernels."}) {
                const std::string file = stem + ext + std::to_string(b);
                std::remove(file.c_str());
                std::remove(("bck.last." + file).c_str());
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
