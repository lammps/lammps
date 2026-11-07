// SPDX-License-Identifier: GPL-2.0-or-later OR LGPL-2.1-or-later

#include "info.h"
#include "input.h"
#include "lammps.h"
#include "library.h"
#include "lmptype.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

using namespace LAMMPS_NS;

namespace {

enum class Backend { SERIAL, OPENMP, KOKKOS_HALF, KOKKOS_FULL };
enum class Potential { D3, LJ, D3_LJ, LJ_D3 };

bool is_kokkos(Backend backend)
{
    return (backend == Backend::KOKKOS_HALF) || (backend == Backend::KOKKOS_FULL);
}

bool uses_gpu()
{
    return Info::has_accelerator_feature("KOKKOS", "api", "cuda") ||
        Info::has_accelerator_feature("KOKKOS", "api", "hip") ||
        Info::has_accelerator_feature("KOKKOS", "api", "sycl");
}

bool available(Backend backend)
{
    if (backend == Backend::OPENMP) return Info::has_package("OPENMP");
    if (!is_kokkos(backend)) return true;
    return Info::has_package("KOKKOS") && (!uses_gpu() || Info::has_kokkos_gpu_device());
}

bool reduced_precision(Backend backend)
{
    return is_kokkos(backend) &&
        !Info::has_accelerator_feature("KOKKOS", "precision", "double");
}

std::unique_ptr<LAMMPS> make_lammps(Backend backend)
{
    LAMMPS::argv args = {"DispersionD3Regression", "-log", "none", "-screen", "none", "-nocite"};
    if (backend == Backend::OPENMP) args.insert(args.end(), {"-sf", "omp", "-pk", "omp", "4"});
    if (is_kokkos(backend)) {
        const bool full = backend == Backend::KOKKOS_FULL;
        args.insert(args.end(), {"-k", "on", uses_gpu() ? "g" : "t", "1", "-sf", "kk",
                                 "-pk", "kokkos", "neigh", full ? "full" : "half",
                                 "newton", full ? "off" : "on"});
    }
    return std::make_unique<LAMMPS>(args, MPI_COMM_WORLD);
}

void commands(LAMMPS &lmp, const std::vector<std::string> &lines)
{
    for (const auto &line : lines) lmp.input->one(line);
}

struct Result {
    double energy;
    std::array<double, 6> virial;
    std::vector<double> forces;
};

Result measure(LAMMPS &lmp)
{
    lmp.input->one("run 0 post no");
    Result result;
    result.energy = lammps_get_thermo(&lmp, "pe");
    const char *fields[] = {"pxx", "pyy", "pzz", "pxy", "pxz", "pyz"};
    const auto nktv2p = *static_cast<double *>(lammps_extract_global(&lmp, "nktv2p"));
    const double scale = lammps_get_thermo(&lmp, "vol") / nktv2p;
    // Compare virials in energy units, without magnifying reduction roundoff
    // through the conversion from eV/Angstrom^3 to bar.
    for (int i = 0; i < 6; ++i) result.virial[i] = scale * lammps_get_thermo(&lmp, fields[i]);

    const int natoms = static_cast<int>(lammps_get_natoms(&lmp));
    const int nlocal = lammps_extract_setting(&lmp, "nlocal");
    const auto *ids = static_cast<tagint *>(lammps_extract_atom(&lmp, "id"));
    const auto *forces = static_cast<double **>(lammps_extract_atom(&lmp, "f"));
    std::vector<double> local(3 * natoms, 0.0);
    for (int i = 0; i < nlocal; ++i)
        for (int d = 0; d < 3; ++d) local[3 * (ids[i] - 1) + d] = forces[i][d];
    result.forces.resize(3 * natoms);
    MPI_Allreduce(local.data(), result.forces.data(), 3 * natoms, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return result;
}

void expect_match(const Result &reference, const Result &candidate, Backend backend)
{
    const double tolerance = reduced_precision(backend) ? 3.0e-5 : 1.0e-10;
    ASSERT_TRUE(std::isfinite(reference.energy));
    ASSERT_TRUE(std::isfinite(candidate.energy));
    EXPECT_NEAR(reference.energy, candidate.energy,
                tolerance * std::max(1.0, std::abs(reference.energy)));
    for (std::size_t i = 0; i < reference.virial.size(); ++i) {
        ASSERT_TRUE(std::isfinite(reference.virial[i]));
        ASSERT_TRUE(std::isfinite(candidate.virial[i]));
        EXPECT_NEAR(reference.virial[i], candidate.virial[i],
                    tolerance * std::max(1.0, std::abs(reference.virial[i])));
    }
    ASSERT_EQ(reference.forces.size(), candidate.forces.size());
    for (std::size_t i = 0; i < reference.forces.size(); ++i) {
        ASSERT_TRUE(std::isfinite(reference.forces[i]));
        ASSERT_TRUE(std::isfinite(candidate.forces[i]));
        EXPECT_NEAR(reference.forces[i], candidate.forces[i],
                    tolerance * std::max(1.0, std::abs(reference.forces[i]))) << "force component " << i;
    }
}

Result crystal(Backend backend, const std::string &damping, double lattice, bool mixed,
               Potential potential = Potential::D3)
{
    auto lmp = make_lammps(backend);
    commands(*lmp, {"units metal", "atom_style atomic", "atom_modify map array",
                    "lattice fcc " + std::to_string(lattice), "region box block 0 3 0 3 0 3",
                    std::string("create_box ") + (mixed ? "8" : "1") + " box",
                    "create_atoms 1 box", "mass * 28.0855", "group probe id 1",
                    "displace_atoms probe move 0.03 -0.02 0.01 units box"});
    if (mixed)
        for (int i = 1; i <= 108; ++i)
            lmp->input->one("set atom " + std::to_string(i) + " type " + std::to_string(1 + i % 8));
    const std::string elements = mixed ? "C H O N Si C O Si" : "Si";
    const std::string d3 = "dispersion/d3 " + damping + " pbe 12.0 8.0";
    if (potential == Potential::D3) {
        commands(*lmp, {"pair_style " + d3, "pair_coeff * * " + elements});
    } else if (potential == Potential::LJ) {
        commands(*lmp, {"pair_style lj/cut 6.0", "pair_coeff * * 0.1 3.0"});
    } else {
        const std::string styles = potential == Potential::D3_LJ ?
            d3 + " lj/cut 6.0" : "lj/cut 6.0 " + d3;
        commands(*lmp, {"pair_style hybrid/overlay " + styles,
                        "pair_coeff * * dispersion/d3 " + elements,
                        "pair_coeff * * lj/cut 0.1 3.0"});
    }
    commands(*lmp, {"neighbor 0.5 bin", "neigh_modify every 1 delay 0 check no",
                    "thermo_style custom step pe pxx pyy pzz pxy pxz pyz"});
    return measure(*lmp);
}

class DispersionD3Test : public ::testing::TestWithParam<Backend> {};

TEST_P(DispersionD3Test, HighCoordinationMatchesCPU)
{
    if (!available(GetParam())) GTEST_SKIP() << "Backend unavailable";
    for (const double lattice : {3.2, 3.4, 3.84, 3.90}) {
        SCOPED_TRACE(lattice);
        const auto reference = crystal(Backend::SERIAL, "bj", lattice, false);
        const auto candidate = crystal(GetParam(), "bj", lattice, false);
        expect_match(reference, candidate, GetParam());
    }
}

TEST_P(DispersionD3Test, DampingAndElementMappingMatchCPU)
{
    if (!available(GetParam())) GTEST_SKIP() << "Backend unavailable";
    for (const std::string damping : {"original", "zerom", "bj", "bjm"}) {
        SCOPED_TRACE(damping);
        const auto reference = crystal(Backend::SERIAL, damping, 5.43, true);
        const auto candidate = crystal(GetParam(), damping, 5.43, true);
        expect_match(reference, candidate, GetParam());
    }
}

TEST_P(DispersionD3Test, HybridOverlayIsAdditive)
{
    if (!available(GetParam())) GTEST_SKIP() << "Backend unavailable";
    auto expected = crystal(GetParam(), "bj", 5.43, false, Potential::D3);
    const auto lj = crystal(GetParam(), "bj", 5.43, false, Potential::LJ);
    expected.energy += lj.energy;
    for (std::size_t i = 0; i < expected.virial.size(); ++i) expected.virial[i] += lj.virial[i];
    for (std::size_t i = 0; i < expected.forces.size(); ++i) expected.forces[i] += lj.forces[i];
    for (const auto potential : {Potential::D3_LJ, Potential::LJ_D3}) {
        const auto actual = crystal(GetParam(), "bj", 5.43, false, potential);
        expect_match(expected, actual, GetParam());
    }
}

TEST_P(DispersionD3Test, ScaledPairEnergyHasConsistentForces)
{
    if (!available(GetParam()) || !Info::has_package("MOLECULE"))
        GTEST_SKIP() << "Backend or MOLECULE package unavailable";
    const double delta = reduced_precision(GetParam()) ? 1.0e-3 : 1.0e-5;
    const double tolerance = reduced_precision(GetParam()) ? 5.0e-5 : 1.0e-7;
    for (const std::string weights : {"1.0 1.0 1.0", "0.5 0.8 1.0", "0.0 0.0 1.0"}) {
        SCOPED_TRACE(weights);
        auto lmp = make_lammps(GetParam());
        commands(*lmp, {"units metal", "atom_style bond", "atom_modify map array", "boundary f f f",
                        "region box block -5 15 -5 15 -5 15",
                        "create_box 1 box bond/types 1 extra/bond/per/atom 2 extra/special/per/atom 3",
                        "mass 1 28.0855", "create_atoms 1 single 0 0 0",
                        "create_atoms 1 single 2.3 0.1 0", "create_atoms 1 single 4.6 0.2 0.4",
                        "create_atoms 1 single 2.0 2.4 0.3", "create_bonds single/bond 1 1 2",
                        "create_bonds single/bond 1 2 3", "bond_style zero", "bond_coeff 1",
                        "special_bonds lj " + weights, "pair_style dispersion/d3 bj pbe 8.0 8.0",
                        "pair_coeff * * Si", "neighbor 0.5 bin",
                        "neigh_modify every 1 delay 0 check no", "group probe id 1",
                        "thermo_style custom step pe pxx pyy pzz pxy pxz pyz"});
        for (int d = 0; d < 3; ++d) {
            const auto reference = measure(*lmp);
            auto displace = [&](double step) {
                std::array<double, 3> move = {0.0, 0.0, 0.0};
                move[d] = step;
                lmp->input->one("displace_atoms probe move " + std::to_string(move[0]) + " " +
                               std::to_string(move[1]) + " " + std::to_string(move[2]) + " units box");
            };
            displace(delta);
            const double plus = measure(*lmp).energy;
            displace(-2.0 * delta);
            const double minus = measure(*lmp).energy;
            displace(delta);
            ASSERT_TRUE(std::isfinite(reference.forces[d]));
            ASSERT_TRUE(std::isfinite(plus));
            ASSERT_TRUE(std::isfinite(minus));
            EXPECT_NEAR(reference.forces[d], -(plus - minus) / (2.0 * delta), tolerance)
                << "force component " << d;
        }
    }
}

INSTANTIATE_TEST_SUITE_P(Backends, DispersionD3Test,
                        ::testing::Values(Backend::SERIAL, Backend::OPENMP,
                                          Backend::KOKKOS_HALF, Backend::KOKKOS_FULL));

} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    lammps_kokkos_finalize();
    MPI_Finalize();
    return result;
}
