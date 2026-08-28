/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// regression tests for the inertial spin dynamics (tspin) styles of the SPIN
// package.  these cover the invariants that are easy to break by accident and
// that the force-style tests cannot express: ownership of the internal
// per-atom state, ownership of the spin mass, the keyword guards, and the
// treatment of non-magnetic atoms.

#include "../testing/core.h"

#include "atom.h"
#include "fix.h"
#include "force.h"
#include "info.h"
#include "math_const.h"
#include "modify.h"
#include "pair.h"
#include "platform.h"
#include "utils.h"

#include "fmt/format.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <mpi.h>
#include <string>
#include <vector>

using ::testing::ContainsRegex;

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

class TSpinTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "TSpinTest";
        LAMMPSTest::SetUp();
        if (!info->has_style("fix", "nve/tspin")) GTEST_SKIP() << "SPIN package not available";
    }

    // 4x4x4 bcc iron; the atoms of group "cold" are left non-magnetic
    void build(bool nonmagnetic = false)
    {
        BEGIN_HIDE_OUTPUT();
        command("units metal");
        command("atom_style spin");
        command("lattice bcc 2.8665");
        command("region box block 0 2 0 2 0 2");
        command("create_box 1 box");
        command("create_atoms 1 box");
        command("mass 1 55.845");
        if (nonmagnetic) {
            command("group cold id 1 2 3 4");
            command("group warm subtract all cold");
            command("set group warm spin/atom/random 31 2.2");
        } else {
            command("set group all spin/atom/random 31 2.2");
        }
        command("pair_style zero 4.0");
        command("pair_coeff * *");
        command("fix pin all spring/tspin 1.0 2.2");
        END_HIDE_OUTPUT();
    }

    double *smass()
    {
        int flag, cols;
        const int index = lmp->atom->find_custom("tspin_smass", flag, cols);
        return (index < 0) ? nullptr : lmp->atom->dvector[index];
    }

    double **vspin()
    {
        int flag, cols;
        const int index = lmp->atom->find_custom("tspin_vs", flag, cols);
        return (index < 0) ? nullptr : lmp->atom->darray[index];
    }
};

// the integrator creates the internal state fix and the two custom properties

TEST_F(TSpinTest, CreatesInternalState)
{
    build();
    BEGIN_HIDE_OUTPUT();
    command("fix 1 all nve/tspin spinmass 0.0075");
    END_HIDE_OUTPUT();

    Fix *state = lmp->modify->get_fix_by_id("TSPIN_STATE");
    ASSERT_NE(state, nullptr);
    ASSERT_THAT(std::string(state->style), ContainsRegex("^property/atom"));

    int flag, cols;
    ASSERT_GE(lmp->atom->find_custom("tspin_vs", flag, cols), 0);
    ASSERT_EQ(flag, 1);
    ASSERT_EQ(cols, 3);
    ASSERT_GE(lmp->atom->find_custom("tspin_smass", flag, cols), 0);
    ASSERT_EQ(flag, 1);
    ASSERT_EQ(cols, 0);
}

// removing the internal state must give a controlled error, not a crash

TEST_F(TSpinTest, RemovedStateIsDiagnosed)
{
    build();
    BEGIN_HIDE_OUTPUT();
    command("fix 1 all nve/tspin spinmass 0.0075");
    command("run 0 post no");
    command("unfix TSPIN_STATE");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*needs the internal per-atom state of the tspin styles.*",
                 command("run 0 post no"););
}

// the reserved fix ID may not be taken by something else

TEST_F(TSpinTest, ReservedFixIdIsProtected)
{
    build();
    BEGIN_HIDE_OUTPUT();
    command("fix TSPIN_STATE all nve");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*is reserved for the internal state of the tspin styles.*",
                 command("fix 1 all nve/tspin spinmass 0.0075"););
}

// the spin degrees of freedom and the chain masses are fixed at setup, so a
// dynamic group would silently give the wrong temperature

TEST_F(TSpinTest, DynamicGroupIsRejected)
{
    build();
    BEGIN_HIDE_OUTPUT();
    command("group dyn dynamic all region box");
    command("fix 1 dyn nvt/tspin temp 300 300 0.1 spinmass 0.0075");
    END_HIDE_OUTPUT();

    // the check happens when the fix list is initialized, not when it is created
    TEST_FAILURE(".*does not allow use with a dynamic group.*", command("run 0 post no"););
}

// the tspin keywords must not be accepted by the plain Nose-Hoover styles

TEST_F(TSpinTest, KeywordsDoNotLeakIntoPlainNH)
{
    build();
    TEST_FAILURE(".*Unknown fix nvt keyword: spinmass.*",
                 command("fix x all nvt temp 300 300 0.1 spinmass 0.01"););
    TEST_FAILURE(".*Unknown fix npt keyword: lattice.*",
                 command("fix y all npt temp 300 300 0.1 iso 0 0 1 lattice frozen"););
}

// dilate takes a group ID, which may be spelled like one of our keywords

TEST_F(TSpinTest, DilateGroupNameIsNotEaten)
{
    build();
    BEGIN_HIDE_OUTPUT();
    command("group spin region box");
    command("fix 1 all npt/tspin temp 300 300 0.1 iso 0 0 1 dilate spin spinmass 0.0075");
    command("run 0 post no");
    END_HIDE_OUTPUT();
    ASSERT_NE(lmp->modify->get_fix_by_id("1"), nullptr);
}

// the fix only owns the spin mass when it was given the spinmass keyword

TEST_F(TSpinTest, SpinMassOwnership)
{
    build();
    BEGIN_HIDE_OUTPUT();
    command("fix 1 all nve/tspin");
    command("velocity/tspin all create 300.0 12345 spinmass 0.02");
    command("run 0 post no");
    END_HIDE_OUTPUT();
    ASSERT_NEAR(smass()[0], 0.02 * 55.845, 1.0e-12);

    BEGIN_HIDE_OUTPUT();
    command("unfix 1");
    command("fix 2 all nve/tspin spinmass 0.0075");
    command("run 0 post no");
    END_HIDE_OUTPUT();
    ASSERT_NEAR(smass()[0], 0.0075 * 55.845, 1.0e-12);
}

// mom yes zeroes the total spin momentum, not the plain sum of velocities

TEST_F(TSpinTest, MomentumIsMassWeighted)
{
    build();
    BEGIN_HIDE_OUTPUT();
    command("fix 1 all nve/tspin spinmass 0.0075");
    command("velocity/tspin all create 300.0 12345 mom yes");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    double **v_s = vspin();
    double *m_s  = smass();
    double p[3]  = {0.0, 0.0, 0.0};
    for (int i = 0; i < lmp->atom->nlocal; ++i)
        for (int d = 0; d < 3; ++d) p[d] += m_s[i] * v_s[i][d];
    double pall[3];
    MPI_Allreduce(p, pall, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (double pd : pall) ASSERT_NEAR(pd, 0.0, 1.0e-10);
}

// non-magnetic atoms get no spin mass, no spin velocity and no error

TEST_F(TSpinTest, NonMagneticAtomsAreSkipped)
{
    build(true);
    BEGIN_HIDE_OUTPUT();
    command("fix 1 all nve/tspin");
    command("velocity/tspin all create 300.0 12345 spinmass 0.0075");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    double **sp   = lmp->atom->sp;
    double **v_s  = vspin();
    double *m_s   = smass();
    for (int i = 0; i < lmp->atom->nlocal; ++i) {
        if (sp[i][3] > 1.0e-8) {
            ASSERT_GT(m_s[i], 0.0);
        } else {
            ASSERT_DOUBLE_EQ(m_s[i], 0.0);
            for (int d = 0; d < 3; ++d) ASSERT_DOUBLE_EQ(v_s[i][d], 0.0);
        }
    }
}

// a dynamical spin remains active when its modulus becomes small or passes
// exactly through zero; a positive spin mass, not the instantaneous modulus,
// identifies a dynamical spin

TEST_F(TSpinTest, SpinCrossesZeroModulus)
{
    build();
    BEGIN_HIDE_OUTPUT();
    command("unfix pin");
    command("fix 1 all nve/tspin spinmass 0.0075");
    command("set atom * d2_tspin_vs[1] 0.0");
    command("set atom * d2_tspin_vs[2] 0.0");
    command("set atom * d2_tspin_vs[3] 0.0");
    command("set atom 1 spin/atom 3.0e-8 1.0 0.0 0.0");
    command("set atom 1 d2_tspin_vs[1] -1.5e-8");
    command("timestep 1.0");
    command("run 3 post no");
    END_HIDE_OUTPUT();

    int iatom = -1;
    for (int i = 0; i < lmp->atom->nlocal; ++i)
        if (lmp->atom->tag[i] == 1) iatom = i;
    ASSERT_GE(iatom, 0);

    double **sp = lmp->atom->sp;
    ASSERT_NEAR(sp[iatom][3] * sp[iatom][0], -1.5e-8, 1.0e-18);
    ASSERT_NEAR(sp[iatom][3], 1.5e-8, 1.0e-18);
    ASSERT_NEAR(vspin()[iatom][0], -1.5e-8, 1.0e-18);
}

// the same zero crossing remains active in the FixNH-derived integrator; the
// nearly uncoupled chain isolates the spin drift while exercising that path

TEST_F(TSpinTest, SpinCrossesZeroModulusNVT)
{
    build();
    BEGIN_HIDE_OUTPUT();
    command("unfix pin");
    command("fix 1 all nvt/tspin temp 1.0e-20 1.0e-20 1.0e20 lattice frozen spinmass 0.0075");
    command("set atom * d2_tspin_vs[1] 0.0");
    command("set atom * d2_tspin_vs[2] 0.0");
    command("set atom * d2_tspin_vs[3] 0.0");
    command("set atom 1 spin/atom 3.0e-8 1.0 0.0 0.0");
    command("set atom 1 d2_tspin_vs[1] -1.5e-8");
    command("timestep 1.0");
    command("run 3 post no");
    END_HIDE_OUTPUT();

    int iatom = -1;
    for (int i = 0; i < lmp->atom->nlocal; ++i)
        if (lmp->atom->tag[i] == 1) iatom = i;
    ASSERT_GE(iatom, 0);

    double **sp = lmp->atom->sp;
    ASSERT_NEAR(sp[iatom][3] * sp[iatom][0], -1.5e-8, 1.0e-18);
    ASSERT_NEAR(sp[iatom][3], 1.5e-8, 1.0e-18);
    ASSERT_NEAR(vspin()[iatom][0], -1.5e-8, 1.0e-18);
}

// the internal state is restart data, not data-file data

TEST_F(TSpinTest, InternalStateIsNotWrittenToDataFile)
{
    build();
    BEGIN_HIDE_OUTPUT();
    command("fix 1 all nve/tspin spinmass 0.0075");
    command("velocity/tspin all create 300.0 12345");
    command("run 0 post no");
    command("write_data tspin_test.data");
    END_HIDE_OUTPUT();

    FILE *fp = fopen("tspin_test.data", "r");
    ASSERT_NE(fp, nullptr);
    char line[512];
    bool found = false;
    while (fgets(line, sizeof(line), fp))
        if (strstr(line, "TSPIN_STATE")) found = true;
    fclose(fp);
    platform::unlink("tspin_test.data");
    ASSERT_FALSE(found);
}

// a restart round trip must bring the spin velocities and masses back

TEST_F(TSpinTest, RestartRestoresSpinState)
{
    build();
    BEGIN_HIDE_OUTPUT();
    command("fix 1 all nve/tspin spinmass 0.0075");
    command("velocity/tspin all create 300.0 12345");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    const int nlocal = lmp->atom->nlocal;
    std::vector<double> vref(3 * nlocal), mref(nlocal);
    for (int i = 0; i < nlocal; ++i) {
        mref[i] = smass()[i];
        for (int d = 0; d < 3; ++d) vref[3 * i + d] = vspin()[i][d];
    }

    BEGIN_HIDE_OUTPUT();
    command("write_restart tspin_test.restart");
    command("clear");
    command("read_restart tspin_test.restart");
    command("pair_style zero 4.0");
    command("pair_coeff * *");
    command("fix pin all spring/tspin 1.0 2.2");
    command("fix 1 all nve/tspin");
    command("run 0 post no");
    END_HIDE_OUTPUT();
    platform::unlink("tspin_test.restart");

    ASSERT_EQ(lmp->atom->nlocal, nlocal);
    for (int i = 0; i < nlocal; ++i) {
        ASSERT_DOUBLE_EQ(smass()[i], mref[i]);
        for (int d = 0; d < 3; ++d) ASSERT_DOUBLE_EQ(vspin()[i][d], vref[3 * i + d]);
    }
}

// the magnetic pair styles express an energy on the unit sphere and cannot
// drive inertial dynamics

TEST_F(TSpinTest, SpinPairStylesAreRejected)
{
    if (!info->has_style("pair", "spin/exchange")) GTEST_SKIP() << "pair spin/exchange missing";
    build();
    BEGIN_HIDE_OUTPUT();
    command("pair_style spin/exchange 3.5");
    command("pair_coeff * * exchange 3.4 0.02726 0.2171 1.841");
    command("fix 1 all nve/tspin spinmass 0.0075");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*does not provide the derivative of the energy with respect to the full spin "
                 "vector.*",
                 command("run 0 post no"););
}

// spin/dipole/cut defines an energy on the complete spin vectors.  Its fm
// output must equal |S|/hbar times the negative Cartesian spin gradient, and
// its mechanical force must remain active when tspin moves the lattice.

TEST_F(TSpinTest, DipoleCutProvidesFullGradient)
{
    if (!info->has_style("pair", "spin/dipole/cut"))
        GTEST_SKIP() << "pair spin/dipole/cut missing";

    constexpr double m1 = 2.0;
    constexpr double m2 = 3.0;
    constexpr double r  = 1.5;
    const double isqrt2 = 1.0 / std::sqrt(2.0);
    const std::array<double, 3> s1 = {m1 * isqrt2, m1 * isqrt2, 0.0};

    BEGIN_HIDE_OUTPUT();
    command("units metal");
    command("atom_style spin");
    command("boundary f f f");
    command("region box block -5 5 -5 5 -5 5 units box");
    command("create_box 1 box");
    command("create_atoms 1 single -0.75 0 0 units box");
    command("create_atoms 1 single 0.75 0 0 units box");
    command("mass 1 1.0");
    command(fmt::format("set atom 1 spin/atom {} {} {} 0.0", m1, isqrt2, isqrt2));
    command("set atom 2 spin/atom 3.0 0.6 0.0 0.8");
    command("pair_style spin/dipole/cut 4.0");
    command("pair_coeff * * 4.0");
    command("fix 1 all nve/tspin lattice moving spinmass 1.0");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    int i1 = -1;
    for (int i = 0; i < lmp->atom->nlocal; ++i)
        if (lmp->atom->tag[i] == 1) i1 = i;
    ASSERT_GE(i1, 0);

    const double hbar = lmp->force->hplanck / MathConst::MY_2PI;
    std::array<double, 3> force_from_fm;
    std::array<double, 3> lattice_force;
    for (int d = 0; d < 3; ++d) {
        force_from_fm[d] = hbar * lmp->atom->fm[i1][d] / m1;
        lattice_force[d] = lmp->atom->f[i1][d];
    }

    const double prefactor = 9.274e-4 * 9.274e-4 * 784.15 / (4.0 * MathConst::MY_PI);
    const double expected_energy =
        prefactor * m1 * m2 / (r * r * r) * (0.6 * isqrt2 - 3.0 * 0.6 * isqrt2);
    ASSERT_NEAR(lmp->force->pair->eng_vdwl, expected_energy, 1.0e-15);

    const std::array<double, 3> expected_spin_force = {
        prefactor * m2 / (r * r * r) * 1.2, 0.0,
        -prefactor * m2 / (r * r * r) * 0.8};
    const std::array<double, 3> expected_lattice_force = {
        3.0 * prefactor * m1 * m2 / (r * r * r * r) * 0.6 * std::sqrt(2.0),
        -3.0 * prefactor * m1 * m2 / (r * r * r * r) * 0.3 * std::sqrt(2.0),
        -3.0 * prefactor * m1 * m2 / (r * r * r * r) * 0.4 * std::sqrt(2.0)};
    for (int d = 0; d < 3; ++d) {
        ASSERT_NEAR(force_from_fm[d], expected_spin_force[d], 1.0e-15);
        ASSERT_NEAR(lattice_force[d], expected_lattice_force[d], 1.0e-15);
    }

    auto spin_energy = [&](const std::array<double, 3> &spin) {
        const double smag =
            std::sqrt(spin[0] * spin[0] + spin[1] * spin[1] + spin[2] * spin[2]);
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("set atom 1 spin/atom {:.17g} {:.17g} {:.17g} {:.17g}", smag,
                            spin[0] / smag, spin[1] / smag, spin[2] / smag));
        command("run 0 post no");
        END_HIDE_OUTPUT();
        return lmp->force->pair->eng_vdwl;
    };

    constexpr double delta = 1.0e-6;
    for (int d = 0; d < 3; ++d) {
        auto plus = s1;
        auto minus = s1;
        plus[d] += delta;
        minus[d] -= delta;
        const double numerical_force = -(spin_energy(plus) - spin_energy(minus)) / (2.0 * delta);
        ASSERT_NEAR(force_from_fm[d], numerical_force, 1.0e-12);
    }

    const double radial_force =
        isqrt2 * (force_from_fm[0] + force_from_fm[1]);
    ASSERT_GT(std::abs(radial_force), 1.0e-6);

    spin_energy(s1);

    auto position_energy = [&](int dim, double value) {
        const char axis = "xyz"[dim];
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("set atom 1 {} {:.17g}", axis, value));
        command("run 0 post no");
        END_HIDE_OUTPUT();
        return lmp->force->pair->eng_vdwl;
    };

    const std::array<double, 3> x1 = {-0.75, 0.0, 0.0};
    for (int d = 0; d < 3; ++d) {
        const double numerical_force =
            -(position_energy(d, x1[d] + delta) - position_energy(d, x1[d] - delta)) /
            (2.0 * delta);
        ASSERT_NEAR(lattice_force[d], numerical_force, 1.0e-12);
        position_energy(d, x1[d]);
    }
}

// With a fixed separation along x and spring/tspin centered at zero, the
// two-spin dipole Hamiltonian is quadratic in the Cartesian spin vectors.  For
// equal spins and masses, the symmetric x and y modes are independent harmonic
// oscillators with known frequencies.

TEST_F(TSpinTest, DipoleCoupledModesConvergeToAnalyticSolution)
{
    if (!info->has_style("pair", "spin/dipole/cut"))
        GTEST_SKIP() << "pair spin/dipole/cut missing";

    constexpr double r = 1.5;
    constexpr double spring_k = 1.0e-4;
    constexpr double final_time = 1.0;
    const double isqrt2 = 1.0 / std::sqrt(2.0);

    auto run_mode = [&](double dt, int nsteps) {
        BEGIN_HIDE_OUTPUT();
        command("clear");
        command("units metal");
        command("atom_style spin");
        command("boundary f f f");
        command("region box block -5 5 -5 5 -5 5 units box");
        command("create_box 1 box");
        command("create_atoms 1 single -0.75 0 0 units box");
        command("create_atoms 1 single 0.75 0 0 units box");
        command("mass 1 1.0");
        command(fmt::format("set group all spin/atom {} {} {} 0.0", std::sqrt(2.0),
                            isqrt2, isqrt2));
        command("pair_style spin/dipole/cut 4.0");
        command("pair_coeff * * 4.0");
        command(fmt::format("fix pin all spring/tspin {:.17g} 0.0", spring_k));
        command("fix 1 all nve/tspin lattice frozen spinmass 1.0");
        command("set group all d2_tspin_vs[1] 0.0");
        command("set group all d2_tspin_vs[2] 0.0");
        command("set group all d2_tspin_vs[3] 0.0");
        command(fmt::format("timestep {:.17g}", dt));
        command(fmt::format("run {} post no", nsteps));
        END_HIDE_OUTPUT();

        const double prefactor = 9.274e-4 * 9.274e-4 * 784.15 / (4.0 * MathConst::MY_PI);
        const double dipole_k = prefactor / (r * r * r);
        const double omega_x = std::sqrt(lmp->force->ftm2v * (spring_k - 2.0 * dipole_k));
        const double omega_y = std::sqrt(lmp->force->ftm2v * (spring_k + dipole_k));
        const std::array<double, 3> spin_exact = {
            std::cos(omega_x * final_time), std::cos(omega_y * final_time), 0.0};
        const std::array<double, 3> velocity_exact = {
            -omega_x * std::sin(omega_x * final_time),
            -omega_y * std::sin(omega_y * final_time), 0.0};

        double error = 0.0;
        int count = 0;
        double **sp = lmp->atom->sp;
        double **vs = vspin();
        for (int i = 0; i < lmp->atom->nlocal; ++i) {
            if ((lmp->atom->tag[i] != 1) && (lmp->atom->tag[i] != 2)) continue;
            count++;
            for (int d = 0; d < 3; ++d) {
                const double spin_value = sp[i][3] * sp[i][d];
                error = std::max(error, std::abs(spin_value - spin_exact[d]));
                error = std::max(error, std::abs(vs[i][d] - velocity_exact[d]));
            }
        }
        EXPECT_EQ(count, 2);
        return error;
    };

    const double error_dt = run_mode(0.002, 500);
    const double error_half = run_mode(0.001, 1000);
    const double error_quarter = run_mode(0.0005, 2000);

    EXPECT_LT(error_dt, 1.0e-6);
    EXPECT_NEAR(error_dt / error_half, 4.0, 0.05);
    EXPECT_NEAR(error_half / error_quarter, 4.0, 0.05);
}

// lattice and spin select the sectors independently; both frozen is an error

TEST_F(TSpinTest, LatticeAndSpinSwitches)
{
    build();
    BEGIN_HIDE_OUTPUT();
    command("fix 1 all nvt/tspin temp 300 300 0.1 spinmass 0.0075 lattice frozen");
    command("velocity/tspin all create 300.0 12345");
    command("run 0 post no");
    command("unfix 1");
    command("fix 2 all nvt/tspin temp 300 300 0.1 spinmass 0.0075 spin frozen");
    command("run 0 post no");
    command("unfix 2");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*nothing to integrate with both lattice and spin frozen.*",
                 command("fix 3 all nvt/tspin temp 300 300 0.1 lattice frozen spin frozen"););
}

// rRESPA has no per-level magnetic force and must be refused

TEST_F(TSpinTest, RespaIsRejected)
{
    build();
    BEGIN_HIDE_OUTPUT();
    command("fix 1 all nve/tspin spinmass 0.0075");
    command("run_style respa 2 2");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*is not compatible with run_style respa.*", command("run 0 post no"););
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (LAMMPS_NS::platform::mpi_vendor() == "Open MPI" && !LAMMPS_NS::Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. "
                     "Death tests will be skipped\n";

    for (int i = 1; i < argc; ++i)
        if (strcmp(argv[i], "-v") == 0) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
