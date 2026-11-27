// unit tests for the parallel (sectoring) algorithm of fix nve/spin
//
// With a frozen lattice and only exchange interactions, every single-spin
// update of the sequential algorithm is a rotation about the local exchange
// field, which conserves the total energy to round-off. If two interacting
// spins were updated at the same time (e.g. because of a wrong sector
// assignment or stale ghost spins), the energy would no longer be conserved.

#include "exceptions.h"
#include "fmt/format.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "library.h"

#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../testing/test_mpi_main.h"

using ::testing::HasSubstr;

namespace LAMMPS_NS {

// max. allowed deviation of the total energy (eV); round-off is ~1e-13,
// broken sectoring gives errors of order 0.1 eV
static constexpr double EPSILON_ENERGY = 1.0e-9;

class FixNVESpinMPITest : public ::testing::Test {
public:
    void command(const std::string &line) { lmp->input->one(line); }

protected:
    const char *testbinary = "LAMMPSTest";
    LAMMPS *lmp            = nullptr;

    void SetUp() override
    {
        int nprocs;
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        if (nprocs != 4) GTEST_SKIP() << "Test requires 4 MPI processes";
        if (!Info::has_package("SPIN")) GTEST_SKIP() << "Test requires the SPIN package";

        LAMMPS::argv args = {testbinary, "-log", "none", "-echo", "screen", "-nocite"};
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(args, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        lmp = nullptr;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    // bcc iron with random spins, exchange interaction only, frozen lattice

    void init_bcc(int nx, int ny, int nz, const std::string &procs)
    {
        command("units           metal");
        command("atom_style      spin");
        command("atom_modify     map array");
        command("boundary        p p p");
        command("processors      " + procs);
        command("lattice         bcc 2.8665");
        command(fmt::format("region          box block 0 {} 0 {} 0 {}", nx, ny, nz));
        command("create_box      1 box");
        command("create_atoms    1 box");
        command("mass            1 55.845");
        command("set             group all spin/random 31 2.2");
        command("pair_style      spin/exchange 3.5");
        command("pair_coeff      * * exchange 3.4 0.02726 0.2171 1.841");
        command("neighbor        0.1 bin");
        command("fix             1 all nve/spin lattice frozen");
        command("timestep        0.0005");
        command("thermo_style    custom step pe etotal");
    }

    // 2d hexagonal lattice with random spins, exchange interaction only, frozen lattice

    void init_hex2d(int nx, int ny, const std::string &procs)
    {
        command("units           metal");
        command("dimension       2");
        command("atom_style      spin");
        command("atom_modify     map array");
        command("boundary        p p p");
        command("processors      " + procs);
        command("lattice         hex 2.5");
        command(fmt::format("region          box block 0 {} 0 {} -0.5 0.5", nx, ny));
        command("create_box      1 box");
        command("create_atoms    1 box");
        command("mass            1 55.845");
        command("set             group all spin/random 31 2.2");
        command("pair_style      spin/exchange 3.0");
        command("pair_coeff      * * exchange 3.0 0.02726 0.2171 1.841");
        command("neighbor        0.1 bin");
        command("fix             1 all nve/spin lattice frozen");
        command("timestep        0.0005");
        command("thermo_style    custom step pe etotal");
    }

    // run nchunk x nstep steps and return the largest deviation of the total energy

    double max_energy_error(int nchunk, int nstep)
    {
        command("run 0 post no");
        const double e0 = lammps_get_thermo(lmp, "etotal");
        double maxerr   = 0.0;
        for (int i = 0; i < nchunk; ++i) {
            command(fmt::format("run {} post no", nstep));
            maxerr = std::max(maxerr, std::fabs(lammps_get_thermo(lmp, "etotal") - e0));
        }
        return maxerr;
    }
};

// sectors in x and y; z is not split and its periodic images are local ghosts

TEST_F(FixNVESpinMPITest, film_3d)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    init_bcc(8, 8, 3, "2 2 1");
    const double err = max_energy_error(4, 50);
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_LT(err, EPSILON_ENERGY);
}

// sectors in z only

TEST_F(FixNVESpinMPITest, bar_3d)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    init_bcc(4, 4, 12, "1 1 4");
    const double err = max_energy_error(4, 50);
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_LT(err, EPSILON_ENERGY);
}

// 2d system: sectors in x and y, z is never split

TEST_F(FixNVESpinMPITest, hex_2d)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    init_hex2d(12, 12, "2 2 1");
    const double err = max_energy_error(4, 50);
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_LT(err, EPSILON_ENERGY);
}

// sub-domains too thin for two sectors in x: must stop with an error

TEST_F(FixNVESpinMPITest, too_thin)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    init_bcc(6, 6, 6, "4 1 1");
    std::string mesg;
    try {
        command("run 0 post no");
    } catch (LAMMPSException &e) {
        mesg = e.what();
    }
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(mesg, HasSubstr("Processor sub-domain of fix nve/spin is"));
    EXPECT_THAT(mesg, HasSubstr("wide in x direction"));
}
} // namespace LAMMPS_NS
