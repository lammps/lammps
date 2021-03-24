/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom.h"
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

    void command(const std::string &cmd) { lmp->input->one(cmd); }
};

TEST_F(LatticeRegionTest, lattice_none)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice none 2.0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    auto lattice = lmp->domain->lattice;
    ASSERT_EQ(lattice->style, Lattice::NONE);
    ASSERT_EQ(lattice->xlattice, 2.0);
    ASSERT_EQ(lattice->ylattice, 2.0);
    ASSERT_EQ(lattice->zlattice, 2.0);
    ASSERT_EQ(lattice->nbasis, 0);
    ASSERT_EQ(lattice->basis, nullptr);

    TEST_FAILURE(".*ERROR: Illegal lattice command.*", command("lattice"););
    TEST_FAILURE(".*ERROR: Illegal lattice command.*", command("lattice xxx"););
    TEST_FAILURE(".*ERROR: Illegal lattice command.*", command("lattice none 1.0 origin"););
    TEST_FAILURE(".*ERROR: Expected floating point.*", command("lattice none xxx"););

    if (!verbose) ::testing::internal::CaptureStdout();
    command("units lj");
    command("lattice none 1.0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    lattice = lmp->domain->lattice;
    ASSERT_EQ(lattice->xlattice, 1.0);
    ASSERT_EQ(lattice->ylattice, 1.0);
    ASSERT_EQ(lattice->zlattice, 1.0);
}

TEST_F(LatticeRegionTest, lattice_sc)
{
    ::testing::internal::CaptureStdout();
    command("lattice sc 1.0 spacing 1.5 2.0 3.0");
    auto output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    ASSERT_THAT(output, MatchesRegex(".*Lattice spacing in x,y,z = 1.50* 2.0* 3.0*.*"));

    auto lattice = lmp->domain->lattice;
    ASSERT_EQ(lattice->xlattice, 1.5);
    ASSERT_EQ(lattice->ylattice, 2.0);
    ASSERT_EQ(lattice->zlattice, 3.0);

    ::testing::internal::CaptureStdout();
    command("lattice sc 2.0");
    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    ASSERT_THAT(output, MatchesRegex(".*Lattice spacing in x,y,z = 2.0* 2.0* 2.0*.*"));

    lattice = lmp->domain->lattice;
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
                 command("lattice sc 1.0 origin 1.0 1.0 1.0"););
    TEST_FAILURE(".*ERROR: Illegal lattice command.*", command("lattice sc 1.0 origin 1.0"););
    TEST_FAILURE(".*ERROR: Illegal lattice command.*",
                 command("lattice sc 1.0 origin 0.0 0.0 0.0 xxx"););
    TEST_FAILURE(".*ERROR: Expected floating point.*",
                 command("lattice sc 1.0 origin xxx 1.0 1.0"););
    TEST_FAILURE(".*ERROR: Lattice orient vectors are not orthogonal.*",
                 command("lattice sc 1.0 orient x 2 2 0"););
    TEST_FAILURE(".*ERROR: Lattice orient vectors are not right-handed.*",
                 command("lattice sc 1.0 orient y 0 -1 0"););
    TEST_FAILURE(".*ERROR: Lattice spacings are invalid.*",
                 command("lattice sc 1.0 spacing 0.0 1.0 1.0"););
    TEST_FAILURE(".*ERROR: Lattice spacings are invalid.*",
                 command("lattice sc 1.0 spacing 1.0 -0.1 1.0"););

    if (!verbose) ::testing::internal::CaptureStdout();
    command("units lj");
    command("lattice sc 2.0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    lattice = lmp->domain->lattice;
    ASSERT_DOUBLE_EQ(lattice->xlattice, pow(0.5, 1.0 / 3.0));
    ASSERT_DOUBLE_EQ(lattice->ylattice, pow(0.5, 1.0 / 3.0));
    ASSERT_DOUBLE_EQ(lattice->zlattice, pow(0.5, 1.0 / 3.0));

    if (!verbose) ::testing::internal::CaptureStdout();
    command("dimension 2");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Lattice style incompatible with simulation dimension.*",
                 command("lattice sc 1.0"););
}

TEST_F(LatticeRegionTest, lattice_bcc)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice bcc 4.2 orient x 1 1 0 orient y -1 1 0");
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
    command("dimension 2");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Lattice style incompatible with simulation dimension.*",
                 command("lattice bcc 1.0"););
}

TEST_F(LatticeRegionTest, lattice_fcc)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice fcc 3.5 origin 0.5 0.5 0.5");
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

    TEST_FAILURE(".*ERROR: Invalid option in lattice command for non-custom style.*",
                 command("lattice fcc 1.0 basis 0.0 0.0 0.0"););
    TEST_FAILURE(".*ERROR: Invalid option in lattice command for non-custom style.*",
                 command("lattice fcc 1.0 a1 0.0 1.0 0.0"););
    TEST_FAILURE(".*ERROR: Illegal lattice command.*", command("lattice fcc 1.0 orient w 1 0 0"););

    if (!verbose) ::testing::internal::CaptureStdout();
    command("dimension 2");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Lattice style incompatible with simulation dimension.*",
                 command("lattice fcc 1.0"););
}

TEST_F(LatticeRegionTest, lattice_hcp)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice hcp 3.0 orient z 0 0 1");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    auto lattice = lmp->domain->lattice;
    ASSERT_EQ(lattice->style, Lattice::HCP);
    ASSERT_DOUBLE_EQ(lattice->xlattice, 3.0);
    ASSERT_DOUBLE_EQ(lattice->ylattice, 3.0 * sqrt(3.0));
    ASSERT_DOUBLE_EQ(lattice->zlattice, 2.0 * sqrt(6.0));
    ASSERT_EQ(lattice->nbasis, 4);
    ASSERT_EQ(lattice->basis[0][0], 0.0);
    ASSERT_EQ(lattice->basis[0][1], 0.0);
    ASSERT_EQ(lattice->basis[0][2], 0.0);
    ASSERT_EQ(lattice->basis[1][0], 0.5);
    ASSERT_EQ(lattice->basis[1][1], 0.5);
    ASSERT_EQ(lattice->basis[1][2], 0.0);
    ASSERT_EQ(lattice->basis[2][0], 0.5);
    ASSERT_DOUBLE_EQ(lattice->basis[2][1], 5.0 / 6.0);
    ASSERT_EQ(lattice->basis[2][2], 0.5);
    ASSERT_EQ(lattice->basis[3][0], 0.0);
    ASSERT_DOUBLE_EQ(lattice->basis[3][1], 1.0 / 3.0);
    ASSERT_EQ(lattice->basis[3][2], 0.5);
    ASSERT_EQ(lattice->a1[0], 1.0);
    ASSERT_EQ(lattice->a1[1], 0.0);
    ASSERT_EQ(lattice->a1[2], 0.0);
    ASSERT_EQ(lattice->a2[0], 0.0);
    ASSERT_DOUBLE_EQ(lattice->a2[1], sqrt(3.0));
    ASSERT_EQ(lattice->a2[2], 0.0);
    ASSERT_EQ(lattice->a3[0], 0.0);
    ASSERT_EQ(lattice->a3[1], 0.0);
    ASSERT_DOUBLE_EQ(lattice->a3[2], sqrt(8.0 / 3.0));

    TEST_FAILURE(".*ERROR: Invalid option in lattice command for non-custom style.*",
                 command("lattice hcp 1.0 a2 0.0 1.0 0.0"););
    TEST_FAILURE(".*ERROR: Invalid option in lattice command for non-custom style.*",
                 command("lattice hcp 1.0 a3 0.0 1.0 0.0"););
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dimension 2");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Lattice style incompatible with simulation dimension.*",
                 command("lattice hcp 1.0"););
}

TEST_F(LatticeRegionTest, lattice_diamond)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice diamond 4.1 orient x 1 1 2 orient y -1 1 0 orient z -1 -1 1");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    auto lattice = lmp->domain->lattice;
    ASSERT_EQ(lattice->style, Lattice::DIAMOND);
    ASSERT_DOUBLE_EQ(lattice->xlattice, 6.6952719636073539);
    ASSERT_DOUBLE_EQ(lattice->ylattice, 5.7982756057296889);
    ASSERT_DOUBLE_EQ(lattice->zlattice, 7.1014083110323973);
    ASSERT_EQ(lattice->nbasis, 8);
    ASSERT_EQ(lattice->basis[0][0], 0.0);
    ASSERT_EQ(lattice->basis[0][1], 0.0);
    ASSERT_EQ(lattice->basis[0][2], 0.0);
    ASSERT_EQ(lattice->basis[1][0], 0.0);
    ASSERT_EQ(lattice->basis[1][1], 0.5);
    ASSERT_EQ(lattice->basis[1][2], 0.5);
    ASSERT_EQ(lattice->basis[2][0], 0.5);
    ASSERT_EQ(lattice->basis[2][1], 0.0);
    ASSERT_EQ(lattice->basis[2][2], 0.5);
    ASSERT_EQ(lattice->basis[3][0], 0.5);
    ASSERT_EQ(lattice->basis[3][1], 0.5);
    ASSERT_EQ(lattice->basis[3][2], 0.0);
    ASSERT_EQ(lattice->basis[4][0], 0.25);
    ASSERT_EQ(lattice->basis[4][1], 0.25);
    ASSERT_EQ(lattice->basis[4][2], 0.25);
    ASSERT_EQ(lattice->basis[5][0], 0.25);
    ASSERT_EQ(lattice->basis[5][1], 0.75);
    ASSERT_EQ(lattice->basis[5][2], 0.75);
    ASSERT_EQ(lattice->basis[6][0], 0.75);
    ASSERT_EQ(lattice->basis[6][1], 0.25);
    ASSERT_EQ(lattice->basis[6][2], 0.75);
    ASSERT_EQ(lattice->basis[7][0], 0.75);
    ASSERT_EQ(lattice->basis[7][1], 0.75);
    ASSERT_EQ(lattice->basis[7][2], 0.25);
    ASSERT_EQ(lattice->a1[0], 1.0);
    ASSERT_EQ(lattice->a1[1], 0.0);
    ASSERT_EQ(lattice->a1[2], 0.0);
    ASSERT_EQ(lattice->a2[0], 0.0);
    ASSERT_EQ(lattice->a2[1], 1.0);
    ASSERT_EQ(lattice->a2[2], 0.0);
    ASSERT_EQ(lattice->a3[0], 0.0);
    ASSERT_EQ(lattice->a3[1], 0.0);
    ASSERT_EQ(lattice->a3[2], 1.0);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("dimension 2");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Lattice style incompatible with simulation dimension.*",
                 command("lattice diamond 1.0"););
}

TEST_F(LatticeRegionTest, lattice_sq)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dimension 2");
    command("lattice sq 3.0");
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

    TEST_FAILURE(".*ERROR: Lattice settings are not compatible with 2d simulation.*",
                 command("lattice sq 1.0 orient x 1 1 2 orient y -1 1 0 orient z -1 -1 1"););

    if (!verbose) ::testing::internal::CaptureStdout();
    command("dimension 3");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Lattice style incompatible with simulation dimension.*",
                 command("lattice sq 1.0"););
}

TEST_F(LatticeRegionTest, lattice_sq2)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dimension 2");
    command("lattice sq2 2.0");
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
    command("dimension 3");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Lattice style incompatible with simulation dimension.*",
                 command("lattice sq2 1.0"););
}

TEST_F(LatticeRegionTest, lattice_hex)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("dimension 2");
    command("lattice hex 2.0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    auto lattice = lmp->domain->lattice;
    ASSERT_EQ(lattice->style, Lattice::HEX);
    ASSERT_DOUBLE_EQ(lattice->xlattice, 2.0);
    ASSERT_DOUBLE_EQ(lattice->ylattice, 3.4641016151377544);
    ASSERT_DOUBLE_EQ(lattice->zlattice, 2.0);
    ASSERT_EQ(lattice->nbasis, 2);
    ASSERT_EQ(lattice->basis[0][0], 0.0);
    ASSERT_EQ(lattice->basis[0][1], 0.0);
    ASSERT_EQ(lattice->basis[0][2], 0.0);
    ASSERT_EQ(lattice->basis[1][0], 0.5);
    ASSERT_EQ(lattice->basis[1][1], 0.5);
    ASSERT_EQ(lattice->basis[1][2], 0.0);
    ASSERT_EQ(lattice->a1[0], 1.0);
    ASSERT_EQ(lattice->a1[1], 0.0);
    ASSERT_EQ(lattice->a1[2], 0.0);
    ASSERT_EQ(lattice->a2[0], 0.0);
    ASSERT_DOUBLE_EQ(lattice->a2[1], sqrt(3.0));
    ASSERT_EQ(lattice->a2[2], 0.0);
    ASSERT_EQ(lattice->a3[0], 0.0);
    ASSERT_EQ(lattice->a3[1], 0.0);
    ASSERT_EQ(lattice->a3[2], 1.0);

    if (!verbose) ::testing::internal::CaptureStdout();
    command("dimension 3");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Lattice style incompatible with simulation dimension.*",
                 command("lattice hex 1.0"););
}

TEST_F(LatticeRegionTest, lattice_custom)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("variable a equal 4.34");
    command("variable b equal $a*sqrt(3.0)");
    command("variable c equal $a*sqrt(8.0/3.0)");
    command("variable t equal 1.0/3.0");
    command("variable f equal 5.0/6.0");
    command("lattice custom  1.0     "
            "a1      $a   0.0  0.0   "
            "a2      0.0  $b   0.0   "
            "a3      0.0  0.0  $c    "
            "basis   0.0  0.0  0.0   "
            "basis   0.5  0.5  0.0   "
            "basis   $t   0.0  0.5   "
            "basis   $f   0.5  0.5   "
            "basis   0.0  0.0  0.625 "
            "basis   0.5  0.5  0.625 "
            "basis   $t   0.0  0.125 "
            "basis   $f   0.5  0.125 ");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    auto lattice = lmp->domain->lattice;
    ASSERT_EQ(lattice->style, Lattice::CUSTOM);
    ASSERT_DOUBLE_EQ(lattice->xlattice, 4.34);
    ASSERT_DOUBLE_EQ(lattice->ylattice, 4.34 * sqrt(3.0));
    ASSERT_DOUBLE_EQ(lattice->zlattice, 4.34 * sqrt(8.0 / 3.0));
    ASSERT_EQ(lattice->nbasis, 8);
    ASSERT_DOUBLE_EQ(lattice->basis[0][0], 0.0);
    ASSERT_DOUBLE_EQ(lattice->basis[0][1], 0.0);
    ASSERT_DOUBLE_EQ(lattice->basis[0][2], 0.0);
    ASSERT_DOUBLE_EQ(lattice->basis[1][0], 0.5);
    ASSERT_DOUBLE_EQ(lattice->basis[1][1], 0.5);
    ASSERT_DOUBLE_EQ(lattice->basis[1][2], 0.0);
    ASSERT_NEAR(lattice->basis[2][0], 1.0 / 3.0, 1.0e-14);
    ASSERT_DOUBLE_EQ(lattice->basis[2][1], 0.0);
    ASSERT_DOUBLE_EQ(lattice->basis[2][2], 0.5);
    ASSERT_DOUBLE_EQ(lattice->basis[3][0], 5.0 / 6.0);
    ASSERT_DOUBLE_EQ(lattice->basis[3][1], 0.5);
    ASSERT_DOUBLE_EQ(lattice->basis[3][2], 0.5);
    ASSERT_DOUBLE_EQ(lattice->basis[4][0], 0.0);
    ASSERT_DOUBLE_EQ(lattice->basis[4][1], 0.0);
    ASSERT_DOUBLE_EQ(lattice->basis[4][2], 0.625);
    ASSERT_DOUBLE_EQ(lattice->basis[5][0], 0.5);
    ASSERT_DOUBLE_EQ(lattice->basis[5][1], 0.5);
    ASSERT_DOUBLE_EQ(lattice->basis[5][2], 0.625);
    ASSERT_NEAR(lattice->basis[6][0], 1.0 / 3.0, 1.0e-14);
    ASSERT_DOUBLE_EQ(lattice->basis[6][1], 0.0);
    ASSERT_DOUBLE_EQ(lattice->basis[6][2], 0.125);
    ASSERT_DOUBLE_EQ(lattice->basis[7][0], 5.0 / 6.0);
    ASSERT_DOUBLE_EQ(lattice->basis[7][1], 0.5);
    ASSERT_DOUBLE_EQ(lattice->basis[7][2], 0.125);
    ASSERT_DOUBLE_EQ(lattice->a1[0], 4.34);
    ASSERT_DOUBLE_EQ(lattice->a1[1], 0.0);
    ASSERT_DOUBLE_EQ(lattice->a1[2], 0.0);
    ASSERT_DOUBLE_EQ(lattice->a2[0], 0.0);
    ASSERT_DOUBLE_EQ(lattice->a2[1], 4.34 * sqrt(3.0));
    ASSERT_DOUBLE_EQ(lattice->a2[2], 0.0);
    ASSERT_DOUBLE_EQ(lattice->a3[0], 0.0);
    ASSERT_DOUBLE_EQ(lattice->a3[1], 0.0);
    ASSERT_DOUBLE_EQ(lattice->a3[2], 4.34 * sqrt(8.0 / 3.0));

    TEST_FAILURE(".*ERROR: Illegal lattice command.*",
                 command("lattice custom 1.0 basis -0.1 0 0"););
    TEST_FAILURE(".*ERROR: Illegal lattice command.*",
                 command("lattice custom 1.0 basis 0.0 1.0 0"););

    if (!verbose) ::testing::internal::CaptureStdout();
    command("dimension 2");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: No basis atoms in lattice.*", command("lattice custom 1.0"););
    TEST_FAILURE(".*ERROR: Lattice settings are not compatible with 2d simulation.*",
                 command("lattice custom 1.0 origin 0.5 0.5 0.5 basis 0.0 0.0 0.0"););
    TEST_FAILURE(".*ERROR: Lattice settings are not compatible with 2d simulation.*",
                 command("lattice custom 1.0 a1 1.0 1.0 1.0 basis 0.0 0.0 0.0"););
    TEST_FAILURE(".*ERROR: Lattice settings are not compatible with 2d simulation.*",
                 command("lattice custom 1.0 a2 1.0 1.0 1.0 basis 0.0 0.0 0.0"););
    TEST_FAILURE(".*ERROR: Lattice settings are not compatible with 2d simulation.*",
                 command("lattice custom 1.0 a3 1.0 1.0 1.0 basis 0.0 0.0 0.0"););
}

TEST_F(LatticeRegionTest, region_fail)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice none 2.0");
    command("region box block 0 1 0 1 0 1");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*ERROR: Create_atoms command before simulation box is defined.*",
                 command("create_atoms 1 box"););
    if (!verbose) ::testing::internal::CaptureStdout();
    command("create_box 1 box");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Cannot create atoms with undefined lattice.*",
                 command("create_atoms 1 box"););
}

TEST_F(LatticeRegionTest, region_block_lattice)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice sc 1.5");
    command("region box block 0 2 0 2 0 2 units lattice");
    command("create_box 1 box");
    command("create_atoms 1 box");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_EQ(lmp->domain->triclinic, 0);
    auto x = lmp->atom->x;
    ASSERT_EQ(lmp->atom->natoms, 8);
    ASSERT_DOUBLE_EQ(x[0][0], 0.0);
    ASSERT_DOUBLE_EQ(x[0][1], 0.0);
    ASSERT_DOUBLE_EQ(x[0][2], 0.0);
    ASSERT_DOUBLE_EQ(x[1][0], 1.5);
    ASSERT_DOUBLE_EQ(x[1][1], 0.0);
    ASSERT_DOUBLE_EQ(x[1][2], 0.0);
    ASSERT_DOUBLE_EQ(x[2][0], 0.0);
    ASSERT_DOUBLE_EQ(x[2][1], 1.5);
    ASSERT_DOUBLE_EQ(x[2][2], 0.0);
    ASSERT_DOUBLE_EQ(x[3][0], 1.5);
    ASSERT_DOUBLE_EQ(x[3][1], 1.5);
    ASSERT_DOUBLE_EQ(x[3][2], 0.0);
}

TEST_F(LatticeRegionTest, region_block_box)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice sc 1.5 origin 0.75 0.75 0.75");
    command("region box block 0 2 0 2 0 2 units box");
    command("create_box 1 box");
    command("create_atoms 1 box");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->domain->triclinic, 0);

    auto x = lmp->atom->x;
    ASSERT_EQ(lmp->atom->natoms, 1);
    ASSERT_DOUBLE_EQ(x[0][0], 1.125);
    ASSERT_DOUBLE_EQ(x[0][1], 1.125);
    ASSERT_DOUBLE_EQ(x[0][2], 1.125);
}

TEST_F(LatticeRegionTest, region_cone)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice fcc 2.5 origin 0.5 0.5 0.5");
    command("region box cone x 1.0 1.0 0.5 2.1 0.0 2.0");
    command("create_box 1 box");
    command("create_atoms 1 region box");
    command("write_dump all atom init.lammpstrj");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->domain->triclinic, 0);
    ASSERT_EQ(lmp->atom->natoms, 42);
}

TEST_F(LatticeRegionTest, region_cylinder)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice fcc 2.5 origin 0.5 0.5 0.5");
    command("region box cylinder z 1.0 1.0 2.1 0.0 2.0 ");
    command("create_box 1 box");
    command("create_atoms 1 region box");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->domain->triclinic, 0);
    ASSERT_EQ(lmp->atom->natoms, 114);
}

TEST_F(LatticeRegionTest, region_prism)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice bcc 2.5 origin 0.75 0.75 0.75");
    command("region box prism 0 2 0 2 0 2 0.5 0.0 0.0");
    command("create_box 1 box");
    command("create_atoms 1 box");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->domain->triclinic, 1);
    ASSERT_EQ(lmp->atom->natoms, 16);
}

TEST_F(LatticeRegionTest, region_sphere)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice fcc 2.5 origin 0.5 0.5 0.5");
    command("region box sphere 1.0 1.0 1.0 1.1");
    command("create_box 1 box");
    command("create_atoms 1 region box");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->domain->triclinic, 0);
    ASSERT_EQ(lmp->atom->natoms, 14);
}

TEST_F(LatticeRegionTest, region_union)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice fcc 2.5 origin 0.5 0.5 0.5");
    command("region part1 sphere 2.0 1.0 1.0 1.1");
    command("region part2 block 0.0 2.0 0.0 2.0 0.0 2.0");
    command("region box union 2 part1 part2");
    command("create_box 1 box");
    command("create_atoms 1 region box");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->domain->triclinic, 0);
    ASSERT_EQ(lmp->atom->natoms, 67);
}

TEST_F(LatticeRegionTest, region_intersect)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice fcc 2.5 origin 0.5 0.5 0.5");
    command("region part1 sphere 2.0 1.0 1.0 1.8");
    command("region part2 block 0.0 2.0 0.0 2.0 0.0 2.0");
    command("region box intersect 2 part1 part2");
    command("create_box 1 box");
    command("create_atoms 1 region box");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->domain->triclinic, 0);
    ASSERT_EQ(lmp->atom->natoms, 21);
}

TEST_F(LatticeRegionTest, region_plane)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    command("lattice fcc 2.5 origin 0.5 0.5 0.5");
    command("region box block 0.0 2.0 0.0 2.0 0.0 2.0");
    command("region part1 plane 0.5 1.0 0.0  0.75 0.0 0.0");
    command("region part2 plane 1.5 1.0 0.0  0.75 0.0 0.0 side out");
    command("region atoms intersect 2 part1 part2");
    command("create_box 1 box");
    command("create_atoms 1 region atoms");
    command("write_dump all atom init.lammpstrj");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->domain->triclinic, 0);
    ASSERT_EQ(lmp->atom->natoms, 16);
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
