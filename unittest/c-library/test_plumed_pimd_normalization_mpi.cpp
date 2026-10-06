// Physical-bias normalization in the beta/P ring-polymer Hamiltonian.

#define LAMMPS_LIB_MPI 1
#include "lammps.h"
#include "library.h"
#include "utils.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "../testing/test_mpi_main.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {
struct BiasState {
    double physical_bias   = 0.0;
    double energy_residual = 0.0;
    std::array<double, 6> force{};
};

BiasState evaluate(const char *mode, int beads, int displaced_bead, double displacement,
                   bool biased, bool normal_modes, bool thermo_energy = true)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const bool probability      = std::string(mode).find("probability_") == 0;
    const bool conditional      = std::string(mode).find("conditional_") == 0;
    const double coupling       = std::string(mode) == "conditional_path" ? 0.5 : 0.0;
    const int ranks_per_bead    = size / beads;
    const int bead              = rank / ranks_per_bead;
    const std::string partition = std::to_string(beads) + "x" + std::to_string(ranks_per_bead);
    const std::string input     = "test_physical_bias_" + std::string(mode) + ".dat";
    if (rank == 0) {
        std::ofstream out(input);
        out << "d: DISTANCE ATOMS=1,2 NOPBC\n"
            << "s: CUSTOM ARG=d FUNC=x*x PERIODIC=NO\n";
        if (probability) {
            // Core-function oracle: actual optional actions are covered by PLUMED
            // native tests and installed-kernel integration qualification.
            out << "v: CUSTOM ARG=d FUNC=0.4*x*x+0.2*x*x*x*x PERIODIC=NO\n"
                << "ratio: CUSTOM ARG=v FUNC=exp(-x) PERIODIC=NO\n"
                << "meanratio: ENSEMBLE ARG=ratio\n"
                << "va: CUSTOM ARG=meanratio.ratio FUNC=-log(x) PERIODIC=NO\n"
                << "dx: DISTANCE ATOMS=1,2 COMPONENTS NOPBC\n"
                << "mean: ENSEMBLE ARG=dx.x\n"
                << "vc: CUSTOM ARG=mean.dx.x FUNC=0.7*x*x PERIODIC=NO\n";
            const std::string name(mode);
            if (name == "probability_mixed")
                out << "energy: CUSTOM ARG=vc,va FUNC=-log(0.5*exp(-x)+0.5*exp(-y-0.2)) "
                       "PERIODIC=NO\n";
            if (biased)
                out << "b: BIASVALUE ARG="
                    << (name == "probability_pure"   ? "va"
                        : name == "probability_zero" ? "vc"
                                                     : "energy")
                    << "\n";
        } else if (conditional) {
            // Build the nonlinear CV after averaging Cartesian components.
            // The existing adapter owns all bead factors and normal-mode transforms.
            out << "dx: DISTANCE ATOMS=1,2 COMPONENTS NOPBC\n"
                << "mean: ENSEMBLE ARG=dx.x\n"
                << "c: CUSTOM ARG=mean.dx.x FUNC=x*x PERIODIC=NO\n"
                << "h: CUSTOM ARG=d FUNC=exp(-x*x/2) PERIODIC=NO\n"
                << "fraction: ENSEMBLE ARG=h\n"
                << "a: CUSTOM ARG=fraction.h FUNC=0.1+x PERIODIC=NO\n"
                << "m: CUSTOM ARG=c FUNC=exp(-1+0.2*x+0.03*x*x) PERIODIC=NO\n"
                << "ratio: CUSTOM ARG=a,m FUNC=log(" << 1.0 - coupling << "+" << coupling
                << "*x/y) PERIODIC=NO\n"
                << "energy: CUSTOM ARG=c,ratio FUNC=0.5*(x-0.7)^2-y PERIODIC=NO\n";
            if (biased) out << "b: BIASVALUE ARG=energy\n";
        } else {
            if (std::string(mode) == "bead_mean") out << "m: ENSEMBLE ARG=s\n";
            if (biased)
                out << "b: RESTRAINT ARG=" << (std::string(mode) == "bead_mean" ? "m.s" : "s")
                    << " AT=0.7 KAPPA=1.0\n";
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    LAMMPS_NS::LAMMPS::argv args = {"test_physical_bias", "-screen", "none", "-log", "none",
                                    "-partition",         partition, "-in",  "none", "-nocite"};
    auto *lmp                    = new LAMMPS_NS::LAMMPS(args, MPI_COMM_WORLD);
    auto command                 = [&](const std::string &line) {
        lammps_command(lmp, line.c_str());
    };
    command("units lj");
    command("atom_style atomic");
    command("atom_modify map array");
    command("region box block 0 8 0 8 0 8");
    command("create_box 1 box");
    command("create_atoms 1 single 3.5 4 4 units box");
    std::ostringstream position;
    position << std::setprecision(17) << "create_atoms 1 single "
             << 4.5 + 0.2 * bead + (bead == displaced_bead ? displacement : 0.0)
             << " 4 4 units box";
    command(position.str());
    command("mass 1 1.0");
    command("pair_style lj/cut 2.5");
    command("pair_coeff * * 0.05 1.0");
    command("velocity all set 0 0 0");
    command(std::string("fix fpimd all pimd/langevin method ") +
            (normal_modes ? "nmpimd" : "pimd") +
            " ensemble nvt integrator obabo "
            "thermostat PILE_L 2468 tau 1.0 temp 1.0 fixcom no");
    command("fix bias all plumed plumedfile " + input +
            " outfile test_physical_bias.log "
            "path_integral " +
            ((conditional || probability) ? "bead_mean" : mode) + " pimd_fix fpimd");
    command(std::string("fix_modify bias energy ") + (thermo_energy ? "yes" : "no"));
    command("run 0 post no");

    BiasState result;
    lammps_gather_atoms(lmp, "f", 1, 3, result.force.data());
    auto *energy = static_cast<double *>(
        lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, 0, 0));
    EXPECT_NE(energy, nullptr);
    double contribution = energy && rank % ranks_per_bead == 0 ? *energy : 0.0;
    if (energy) lammps_free(energy);
    MPI_Allreduce(&contribution, &result.physical_bias, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (displaced_bead == -1) {
        // Cartesian PIMD updates energy estimators at end_of_step, not run 0.
        // Compare all terms at the same completed step, after coordinates move.
        command("compute physical_energy all pe pair");
        command("run 1 post no");
        auto *potential = static_cast<double *>(
            lammps_extract_fix(lmp, "fpimd", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, 2, 0));
        auto *current_bias = static_cast<double *>(
            lammps_extract_fix(lmp, "bias", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, 0, 0));
        auto *physical = static_cast<double *>(
            lammps_extract_compute(lmp, "physical_energy", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR));
        EXPECT_NE(potential, nullptr);
        EXPECT_NE(current_bias, nullptr);
        EXPECT_NE(physical, nullptr);
        contribution = potential && current_bias && physical && rank % ranks_per_bead == 0
                           ? *potential - *physical - beads * (*current_bias)
                           : 0.0;
        MPI_Allreduce(&contribution, &result.energy_residual, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        if (potential) lammps_free(potential);
        if (current_bias) lammps_free(current_bias);
    }
    EXPECT_EQ(lammps_has_error(lmp), 0);
    delete lmp;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::remove(input.c_str());
        std::remove("test_physical_bias.log");
        for (int b = 0; b < beads; ++b)
            std::remove(("test_physical_bias.log." + std::to_string(b)).c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return result;
}
// Orthonormal real Fourier basis for the ring polymer. This is an analytic
// coordinate definition, independent of the integrator's transformation arrays.
double normal_mode_coefficient(int mode, int bead, int beads)
{
    const double normalization = 1.0 / std::sqrt(double(beads));
    if (mode == 0) return normalization;
    if (beads % 2 == 0 && mode == beads / 2) return (bead % 2 ? -1.0 : 1.0) * normalization;
    const double phase = 2.0 * std::acos(-1.0) * mode * bead / beads;
    return std::sqrt(2.0) * normalization * (mode <= beads / 2 ? std::cos(phase) : std::sin(phase));
}

void check_physical_bias(bool normal_modes)
{
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ASSERT_EQ(size, 4);
    constexpr double delta = 1.0e-5;
    for (int beads : {1, 2, 4}) {
        const int bead = rank / (size / beads);
        for (const char *mode :
             {"centroid", "bead_mean", "bead_density", "conditional_zero", "conditional_path",
              "probability_pure", "probability_zero", "probability_mixed"}) {
            if (beads == 1 && std::string(mode) != "centroid") continue;
            SCOPED_TRACE(std::string(mode) + " P=" + std::to_string(beads));
            const auto zero      = evaluate(mode, beads, -1, 0.0, false, normal_modes);
            const auto state     = evaluate(mode, beads, -1, 0.0, true, normal_modes);
            double mean_distance = 0.0, mean_squared = 0.0, mean_bias = 0.0;
            for (int b = 0; b < beads; ++b) {
                const double distance = 1.0 + 0.2 * b;
                mean_distance += distance / beads;
                mean_squared += distance * distance / beads;
                mean_bias += 0.5 * std::pow(distance * distance - 0.7, 2) / beads;
            }
            double expected_bias = std::string(mode) == "centroid"
                                       ? 0.5 * std::pow(mean_distance * mean_distance - 0.7, 2)
                                   : std::string(mode) == "bead_mean"
                                       ? 0.5 * std::pow(mean_squared - 0.7, 2)
                                       : mean_bias;
            if (std::string(mode).find("conditional_") == 0) {
                double score = 0.1;
                for (int b = 0; b < beads; ++b) {
                    const double x = 1.0 + 0.2 * b;
                    score += std::exp(-x * x / 2) / beads;
                }
                const double c          = mean_distance * mean_distance;
                const double normalizer = std::exp(-1 + 0.2 * c + 0.03 * c * c);
                const double lambda     = std::string(mode) == "conditional_path" ? 0.5 : 0.0;
                expected_bias =
                    0.5 * std::pow(c - 0.7, 2) - std::log(1 - lambda + lambda * score / normalizer);
            }
            if (std::string(mode).find("probability_") == 0) {
                double ratio = 0;
                for (int b = 0; b < beads; ++b) {
                    const double x = 1.0 + 0.2 * b;
                    ratio += std::exp(-0.4 * x * x - 0.2 * std::pow(x, 4)) / beads;
                }
                const double vc = 0.7 * mean_distance * mean_distance;
                expected_bias = std::string(mode) == "probability_pure" ? -std::log(ratio)
                                : std::string(mode) == "probability_zero"
                                    ? vc
                                    : -std::log(0.5 * std::exp(-vc) + 0.5 * ratio * std::exp(-0.2));
            }
            EXPECT_NEAR(state.physical_bias, expected_bias, 1.0e-12);

            EXPECT_NEAR(zero.energy_residual, 0.0, 1.0e-11);
            EXPECT_NEAR(state.energy_residual, 0.0, 1.0e-11);
            const auto without_thermo = evaluate(mode, beads, -1, 0.0, true, normal_modes, false);
            EXPECT_NEAR(without_thermo.physical_bias, state.physical_bias, 1.0e-12);
            EXPECT_NEAR(without_thermo.energy_residual, 0.0, 1.0e-11);
            for (int d = 0; d < 6; ++d)
                EXPECT_NEAR(without_thermo.force[d], state.force[d], 1.0e-12);

            double expected_force = 0.0;
            for (int displaced = 0; displaced < beads; ++displaced) {
                const auto plus  = evaluate(mode, beads, displaced, delta, true, normal_modes);
                const auto minus = evaluate(mode, beads, displaced, -delta, true, normal_modes);
                // Differentiate measured physical B, not a replica of the force code.
                // At beta/P the dynamical potential is H_ring + P*B.
                const double force =
                    -beads * (plus.physical_bias - minus.physical_bias) / (2 * delta);
                const double coefficient = normal_modes
                                               ? normal_mode_coefficient(bead, displaced, beads)
                                               : (bead == displaced ? 1.0 : 0.0);
                expected_force += coefficient * force;
            }
            EXPECT_NEAR(state.force[3] - zero.force[3], expected_force, 2.0e-8);
            EXPECT_NEAR(state.force[0] - zero.force[0], -expected_force, 2.0e-8);
            for (int d : {1, 2, 4, 5})
                EXPECT_NEAR(state.force[d] - zero.force[d], 0.0, 1.0e-12);
        }
    }
}

} // namespace

TEST(MPI, plumed_pimd_physical_bias_finite_difference)
{
    check_physical_bias(false);
}

TEST(MPI, plumed_nmpimd_physical_bias_finite_difference)
{
    check_physical_bias(true);
}
