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

#include "domain.h"
#include "fmt/format.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "lattice.h"
#include "region.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

#if defined(OMPI_MAJOR_VERSION)
const bool have_openmpi = true;
#else
const bool have_openmpi = false;
#endif

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::ExitedWithCode;
using ::testing::MatchesRegex;
using ::testing::StrEq;

#define TEST_FAILURE(errmsg, ...)                                 \
    if (Info::has_exceptions()) {                                 \
        ::testing::internal::CaptureStdout();                     \
        ASSERT_ANY_THROW({__VA_ARGS__});                          \
        auto mesg = ::testing::internal::GetCapturedStdout();     \
        ASSERT_THAT(mesg, MatchesRegex(errmsg));                  \
    } else {                                                      \
        if (!have_openmpi) {                                      \
            ::testing::internal::CaptureStdout();                 \
            ASSERT_DEATH({__VA_ARGS__}, "");                      \
            auto mesg = ::testing::internal::GetCapturedStdout(); \
            ASSERT_THAT(mesg, MatchesRegex(errmsg));              \
        }                                                         \
    }

class LatticeRegionTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {"LatticeRegionTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        lmp->input->one("units metal");
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(LatticeRegionTest, lattice_none)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("lattice none 2.0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    auto lattice = lmp->domain->lattice;
    ASSERT_EQ(lattice->style, Lattice::NONE);
    ASSERT_EQ(lattice->xlattice, 2.0);
    ASSERT_EQ(lattice->ylattice, 2.0);
    ASSERT_EQ(lattice->zlattice, 2.0);
    ASSERT_EQ(lattice->nbasis, 0);
    ASSERT_EQ(lattice->basis, nullptr);

    TEST_FAILURE(".*ERROR: Illegal lattice command.*", lmp->input->one("lattice"););
    TEST_FAILURE(".*ERROR: Illegal lattice command.*", lmp->input->one("lattice xxx"););
    TEST_FAILURE(".*ERROR: Illegal lattice command.*", lmp->input->one("lattice none 1.0 origin"););
    TEST_FAILURE(".*ERROR: Expected floating point.*", lmp->input->one("lattice none xxx"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units lj");
    lmp->input->one("lattice none 1.0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lattice->xlattice, 1.0);
    ASSERT_EQ(lattice->ylattice, 1.0);
    ASSERT_EQ(lattice->zlattice, 1.0);
}

TEST_F(LatticeRegionTest, lattice_sc)
{
    ::testing::internal::CaptureStdout();
    lmp->input->one("lattice sc 2.0");
    auto output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    ASSERT_THAT(output, MatchesRegex(".*Lattice spacing in x,y,z = 2.0* 2.0* 2.0*.*"));

    auto lattice = lmp->domain->lattice;
    ASSERT_EQ(lattice->style, Lattice::SC);
    ASSERT_EQ(lattice->xlattice, 2.0);
    ASSERT_EQ(lattice->ylattice, 2.0);
    ASSERT_EQ(lattice->zlattice, 2.0);
    ASSERT_EQ(lattice->nbasis, 1);
    ASSERT_NE(lattice->basis, nullptr);
    ASSERT_EQ(lattice->a1[0], 1.0);
    ASSERT_EQ(lattice->a1[1], 0.0);
    ASSERT_EQ(lattice->a1[2], 0.0);
    ASSERT_EQ(lattice->a2[0], 0.0);
    ASSERT_EQ(lattice->a2[1], 1.0);
    ASSERT_EQ(lattice->a2[2], 0.0);
    ASSERT_EQ(lattice->a3[0], 0.0);
    ASSERT_EQ(lattice->a3[1], 0.0);
    ASSERT_EQ(lattice->a3[2], 1.0);
    ASSERT_EQ(lattice->basis[0][0], 0.0);
    ASSERT_EQ(lattice->basis[0][1], 0.0);
    ASSERT_EQ(lattice->basis[0][2], 0.0);

    TEST_FAILURE(".*ERROR: Illegal lattice command.*",
                 lmp->input->one("lattice sc 1.0 origin 1.0 1.0 1.0"););
    TEST_FAILURE(".*ERROR: Illegal lattice command.*",
                 lmp->input->one("lattice sc 1.0 origin 1.0"););
    TEST_FAILURE(".*ERROR: Expected floating point.*",
                 lmp->input->one("lattice sc 1.0 origin xxx 1.0 1.0"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("units lj");
    lmp->input->one("lattice sc 2.0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_DOUBLE_EQ(lattice->xlattice, pow(0.5, 1.0 / 3.0));
    ASSERT_DOUBLE_EQ(lattice->ylattice, pow(0.5, 1.0 / 3.0));
    ASSERT_DOUBLE_EQ(lattice->zlattice, pow(0.5, 1.0 / 3.0));

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("dimension 2");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Lattice style incompatible with simulation dimension.*",
                 lmp->input->one("lattice sc 1.0"););
}

TEST_F(LatticeRegionTest, lattice_bcc)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("lattice bcc 4.2 orient x 1 1 0 orient y -1 1 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    auto lattice = lmp->domain->lattice;
    ASSERT_EQ(lattice->style, Lattice::BCC);
    ASSERT_DOUBLE_EQ(lattice->xlattice, sqrt(2.0) * 4.2);
    ASSERT_DOUBLE_EQ(lattice->ylattice, sqrt(2.0) * 4.2);
    ASSERT_DOUBLE_EQ(lattice->zlattice, 4.2);
    ASSERT_EQ(lattice->nbasis, 2);
    ASSERT_EQ(lattice->basis[0][0], 0.0);
    ASSERT_EQ(lattice->basis[0][1], 0.0);
    ASSERT_EQ(lattice->basis[0][2], 0.0);
    ASSERT_EQ(lattice->basis[1][0], 0.5);
    ASSERT_EQ(lattice->basis[1][1], 0.5);
    ASSERT_EQ(lattice->basis[1][2], 0.5);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("dimension 2");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Lattice style incompatible with simulation dimension.*",
                 lmp->input->one("lattice bcc 1.0"););
}

TEST_F(LatticeRegionTest, lattice_fcc)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("lattice fcc 3.5 origin 0.5 0.5 0.5");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    auto lattice = lmp->domain->lattice;
    ASSERT_EQ(lattice->style, Lattice::FCC);
    ASSERT_DOUBLE_EQ(lattice->xlattice, 3.5);
    ASSERT_DOUBLE_EQ(lattice->ylattice, 3.5);
    ASSERT_DOUBLE_EQ(lattice->zlattice, 3.5);
    ASSERT_EQ(lattice->nbasis, 4);
    ASSERT_EQ(lattice->basis[0][0], 0.0);
    ASSERT_EQ(lattice->basis[0][1], 0.0);
    ASSERT_EQ(lattice->basis[0][2], 0.0);
    ASSERT_EQ(lattice->basis[1][0], 0.5);
    ASSERT_EQ(lattice->basis[1][1], 0.5);
    ASSERT_EQ(lattice->basis[1][2], 0.0);
    ASSERT_EQ(lattice->basis[2][0], 0.5);
    ASSERT_EQ(lattice->basis[2][1], 0.0);
    ASSERT_EQ(lattice->basis[2][2], 0.5);
    ASSERT_EQ(lattice->basis[3][0], 0.0);
    ASSERT_EQ(lattice->basis[3][1], 0.5);
    ASSERT_EQ(lattice->basis[3][2], 0.5);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("dimension 2");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Lattice style incompatible with simulation dimension.*",
                 lmp->input->one("lattice fcc 1.0"););
}

TEST_F(LatticeRegionTest, lattice_sq)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("dimension 2");
    lmp->input->one("lattice sq 3.0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    auto lattice = lmp->domain->lattice;
    ASSERT_EQ(lattice->style, Lattice::SQ);
    ASSERT_DOUBLE_EQ(lattice->xlattice, 3.0);
    ASSERT_DOUBLE_EQ(lattice->ylattice, 3.0);
    ASSERT_DOUBLE_EQ(lattice->zlattice, 3.0);
    ASSERT_EQ(lattice->nbasis, 1);
    ASSERT_EQ(lattice->basis[0][0], 0.0);
    ASSERT_EQ(lattice->basis[0][1], 0.0);
    ASSERT_EQ(lattice->basis[0][2], 0.0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("dimension 3");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Lattice style incompatible with simulation dimension.*",
                 lmp->input->one("lattice sq 1.0"););
}

TEST_F(LatticeRegionTest, lattice_sq2)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("dimension 2");
    lmp->input->one("lattice sq2 2.0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    auto lattice = lmp->domain->lattice;
    ASSERT_EQ(lattice->style, Lattice::SQ2);
    ASSERT_DOUBLE_EQ(lattice->xlattice, 2.0);
    ASSERT_DOUBLE_EQ(lattice->ylattice, 2.0);
    ASSERT_DOUBLE_EQ(lattice->zlattice, 2.0);
    ASSERT_EQ(lattice->nbasis, 2);
    ASSERT_EQ(lattice->basis[0][0], 0.0);
    ASSERT_EQ(lattice->basis[0][1], 0.0);
    ASSERT_EQ(lattice->basis[0][2], 0.0);
    ASSERT_EQ(lattice->basis[1][0], 0.5);
    ASSERT_EQ(lattice->basis[1][1], 0.5);
    ASSERT_EQ(lattice->basis[1][2], 0.0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("dimension 3");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Lattice style incompatible with simulation dimension.*",
                 lmp->input->one("lattice sq2 1.0"););
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (have_openmpi && !LAMMPS_NS::Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. "
                     "Death tests will be skipped\n";

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
