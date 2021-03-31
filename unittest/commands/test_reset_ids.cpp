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
#include "fmt/format.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "output.h"
#include "update.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
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
using ::testing::MatchesRegex;

#define GETIDX(i) lmp->atom->map(i)

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

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

class ResetIDsTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {"ResetIDsTest", "-log", "none", "-nocite", "-echo", "screen"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp        = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        Info *info = new Info(lmp);
        if (info->has_style("atom", "full")) {
            lmp->input->one("variable input_dir index " STRINGIFY(TEST_INPUT_FOLDER));
            lmp->input->one("include ${input_dir}/in.fourmol");
        }
        delete info;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(ResetIDsTest, MolIDAll)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();

    auto molid = lmp->atom->molecule;
    ASSERT_EQ(molid[GETIDX(1)], 1);
    ASSERT_EQ(molid[GETIDX(2)], 1);
    ASSERT_EQ(molid[GETIDX(3)], 1);
    ASSERT_EQ(molid[GETIDX(4)], 1);
    ASSERT_EQ(molid[GETIDX(5)], 1);
    ASSERT_EQ(molid[GETIDX(6)], 1);
    ASSERT_EQ(molid[GETIDX(7)], 1);
    ASSERT_EQ(molid[GETIDX(8)], 2);
    ASSERT_EQ(molid[GETIDX(9)], 2);
    ASSERT_EQ(molid[GETIDX(10)], 2);
    ASSERT_EQ(molid[GETIDX(11)], 2);
    ASSERT_EQ(molid[GETIDX(12)], 2);
    ASSERT_EQ(molid[GETIDX(13)], 2);
    ASSERT_EQ(molid[GETIDX(14)], 2);
    ASSERT_EQ(molid[GETIDX(15)], 2);
    ASSERT_EQ(molid[GETIDX(16)], 2);
    ASSERT_EQ(molid[GETIDX(17)], 2);
    ASSERT_EQ(molid[GETIDX(18)], 3);
    ASSERT_EQ(molid[GETIDX(19)], 3);
    ASSERT_EQ(molid[GETIDX(20)], 3);
    ASSERT_EQ(molid[GETIDX(21)], 4);
    ASSERT_EQ(molid[GETIDX(22)], 4);
    ASSERT_EQ(molid[GETIDX(23)], 4);
    ASSERT_EQ(molid[GETIDX(24)], 5);
    ASSERT_EQ(molid[GETIDX(25)], 5);
    ASSERT_EQ(molid[GETIDX(26)], 5);
    ASSERT_EQ(molid[GETIDX(27)], 6);
    ASSERT_EQ(molid[GETIDX(28)], 6);
    ASSERT_EQ(molid[GETIDX(29)], 6);

    // the original data file has two different molecule IDs
    // for two residues of the same molecule/fragment.
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("reset_mol_ids all");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_EQ(molid[GETIDX(1)], 1);
    ASSERT_EQ(molid[GETIDX(2)], 1);
    ASSERT_EQ(molid[GETIDX(3)], 1);
    ASSERT_EQ(molid[GETIDX(4)], 1);
    ASSERT_EQ(molid[GETIDX(5)], 1);
    ASSERT_EQ(molid[GETIDX(6)], 1);
    ASSERT_EQ(molid[GETIDX(7)], 1);
    ASSERT_EQ(molid[GETIDX(8)], 1);
    ASSERT_EQ(molid[GETIDX(9)], 1);
    ASSERT_EQ(molid[GETIDX(10)], 1);
    ASSERT_EQ(molid[GETIDX(11)], 1);
    ASSERT_EQ(molid[GETIDX(12)], 1);
    ASSERT_EQ(molid[GETIDX(13)], 1);
    ASSERT_EQ(molid[GETIDX(14)], 1);
    ASSERT_EQ(molid[GETIDX(15)], 1);
    ASSERT_EQ(molid[GETIDX(16)], 1);
    ASSERT_EQ(molid[GETIDX(17)], 1);
    ASSERT_EQ(molid[GETIDX(18)], 2);
    ASSERT_EQ(molid[GETIDX(19)], 2);
    ASSERT_EQ(molid[GETIDX(20)], 2);
    ASSERT_EQ(molid[GETIDX(21)], 3);
    ASSERT_EQ(molid[GETIDX(22)], 3);
    ASSERT_EQ(molid[GETIDX(23)], 3);
    ASSERT_EQ(molid[GETIDX(24)], 4);
    ASSERT_EQ(molid[GETIDX(25)], 4);
    ASSERT_EQ(molid[GETIDX(26)], 4);
    ASSERT_EQ(molid[GETIDX(27)], 5);
    ASSERT_EQ(molid[GETIDX(28)], 5);
    ASSERT_EQ(molid[GETIDX(29)], 5);
}

TEST_F(ResetIDsTest, DeletePlusAtomID)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();

    auto molid = lmp->atom->molecule;

    // delete two water molecules
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("group allwater molecule 3:6");
    lmp->input->one("group twowater molecule 4:6:2");
    lmp->input->one("delete_atoms group twowater compress no bond yes");
    lmp->input->one("reset_mol_ids all");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->atom->natoms, 23);
    ASSERT_EQ(lmp->atom->map_tag_max, 26);

    ASSERT_EQ(molid[GETIDX(1)], 1);
    ASSERT_EQ(molid[GETIDX(2)], 1);
    ASSERT_EQ(molid[GETIDX(3)], 1);
    ASSERT_EQ(molid[GETIDX(4)], 1);
    ASSERT_EQ(molid[GETIDX(5)], 1);
    ASSERT_EQ(molid[GETIDX(6)], 1);
    ASSERT_EQ(molid[GETIDX(7)], 1);
    ASSERT_EQ(molid[GETIDX(8)], 1);
    ASSERT_EQ(molid[GETIDX(9)], 1);
    ASSERT_EQ(molid[GETIDX(10)], 1);
    ASSERT_EQ(molid[GETIDX(11)], 1);
    ASSERT_EQ(molid[GETIDX(12)], 1);
    ASSERT_EQ(molid[GETIDX(13)], 1);
    ASSERT_EQ(molid[GETIDX(14)], 1);
    ASSERT_EQ(molid[GETIDX(15)], 1);
    ASSERT_EQ(molid[GETIDX(16)], 1);
    ASSERT_EQ(molid[GETIDX(17)], 1);
    ASSERT_EQ(molid[GETIDX(18)], 2);
    ASSERT_EQ(molid[GETIDX(19)], 2);
    ASSERT_EQ(molid[GETIDX(20)], 2);
    ASSERT_EQ(molid[GETIDX(24)], 3);
    ASSERT_EQ(molid[GETIDX(25)], 3);
    ASSERT_EQ(molid[GETIDX(26)], 3);

    // now also check and reset the atom ids

    ASSERT_GE(GETIDX(1), 0);
    ASSERT_GE(GETIDX(2), 0);
    ASSERT_GE(GETIDX(3), 0);
    ASSERT_GE(GETIDX(4), 0);
    ASSERT_GE(GETIDX(5), 0);
    ASSERT_GE(GETIDX(6), 0);
    ASSERT_GE(GETIDX(7), 0);
    ASSERT_GE(GETIDX(8), 0);
    ASSERT_GE(GETIDX(9), 0);
    ASSERT_GE(GETIDX(10), 0);
    ASSERT_GE(GETIDX(11), 0);
    ASSERT_GE(GETIDX(12), 0);
    ASSERT_GE(GETIDX(13), 0);
    ASSERT_GE(GETIDX(14), 0);
    ASSERT_GE(GETIDX(15), 0);
    ASSERT_GE(GETIDX(16), 0);
    ASSERT_GE(GETIDX(17), 0);
    ASSERT_GE(GETIDX(18), 0);
    ASSERT_GE(GETIDX(19), 0);
    ASSERT_GE(GETIDX(20), 0);
    ASSERT_EQ(GETIDX(21), -1);
    ASSERT_EQ(GETIDX(22), -1);
    ASSERT_EQ(GETIDX(23), -1);
    ASSERT_GE(GETIDX(24), 0);
    ASSERT_GE(GETIDX(25), 0);
    ASSERT_GE(GETIDX(26), 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("reset_atom_ids");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_EQ(lmp->atom->map_tag_max, 23);
    for (int i = 1; i <= 23; ++i)
        ASSERT_GE(GETIDX(i), 0);
}

TEST_F(ResetIDsTest, PartialOffset)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();

    auto molid = lmp->atom->molecule;

    // delete two water molecules
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("group allwater molecule 3:6");
    lmp->input->one("group nowater subtract all allwater");
    lmp->input->one("reset_mol_ids allwater offset 4");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->atom->natoms, 29);
    ASSERT_EQ(lmp->atom->map_tag_max, 29);

    ASSERT_EQ(molid[GETIDX(1)], 1);
    ASSERT_EQ(molid[GETIDX(2)], 1);
    ASSERT_EQ(molid[GETIDX(3)], 1);
    ASSERT_EQ(molid[GETIDX(4)], 1);
    ASSERT_EQ(molid[GETIDX(5)], 1);
    ASSERT_EQ(molid[GETIDX(6)], 1);
    ASSERT_EQ(molid[GETIDX(7)], 1);
    ASSERT_EQ(molid[GETIDX(8)], 2);
    ASSERT_EQ(molid[GETIDX(9)], 2);
    ASSERT_EQ(molid[GETIDX(10)], 2);
    ASSERT_EQ(molid[GETIDX(11)], 2);
    ASSERT_EQ(molid[GETIDX(12)], 2);
    ASSERT_EQ(molid[GETIDX(13)], 2);
    ASSERT_EQ(molid[GETIDX(14)], 2);
    ASSERT_EQ(molid[GETIDX(15)], 2);
    ASSERT_EQ(molid[GETIDX(16)], 2);
    ASSERT_EQ(molid[GETIDX(17)], 2);
    ASSERT_EQ(molid[GETIDX(18)], 5);
    ASSERT_EQ(molid[GETIDX(19)], 5);
    ASSERT_EQ(molid[GETIDX(20)], 5);
    ASSERT_EQ(molid[GETIDX(21)], 6);
    ASSERT_EQ(molid[GETIDX(22)], 6);
    ASSERT_EQ(molid[GETIDX(23)], 6);
    ASSERT_EQ(molid[GETIDX(24)], 7);
    ASSERT_EQ(molid[GETIDX(25)], 7);
    ASSERT_EQ(molid[GETIDX(26)], 7);
    ASSERT_EQ(molid[GETIDX(27)], 8);
    ASSERT_EQ(molid[GETIDX(28)], 8);
    ASSERT_EQ(molid[GETIDX(29)], 8);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("reset_mol_ids nowater offset 0");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_EQ(molid[GETIDX(1)], 1);
    ASSERT_EQ(molid[GETIDX(2)], 1);
    ASSERT_EQ(molid[GETIDX(3)], 1);
    ASSERT_EQ(molid[GETIDX(4)], 1);
    ASSERT_EQ(molid[GETIDX(5)], 1);
    ASSERT_EQ(molid[GETIDX(6)], 1);
    ASSERT_EQ(molid[GETIDX(7)], 1);
    ASSERT_EQ(molid[GETIDX(8)], 1);
    ASSERT_EQ(molid[GETIDX(9)], 1);
    ASSERT_EQ(molid[GETIDX(10)], 1);
    ASSERT_EQ(molid[GETIDX(11)], 1);
    ASSERT_EQ(molid[GETIDX(12)], 1);
    ASSERT_EQ(molid[GETIDX(13)], 1);
    ASSERT_EQ(molid[GETIDX(14)], 1);
    ASSERT_EQ(molid[GETIDX(15)], 1);
    ASSERT_EQ(molid[GETIDX(16)], 1);
    ASSERT_EQ(molid[GETIDX(17)], 1);
    ASSERT_EQ(molid[GETIDX(18)], 5);
    ASSERT_EQ(molid[GETIDX(19)], 5);
    ASSERT_EQ(molid[GETIDX(20)], 5);
    ASSERT_EQ(molid[GETIDX(21)], 6);
    ASSERT_EQ(molid[GETIDX(22)], 6);
    ASSERT_EQ(molid[GETIDX(23)], 6);
    ASSERT_EQ(molid[GETIDX(24)], 7);
    ASSERT_EQ(molid[GETIDX(25)], 7);
    ASSERT_EQ(molid[GETIDX(26)], 7);
    ASSERT_EQ(molid[GETIDX(27)], 8);
    ASSERT_EQ(molid[GETIDX(28)], 8);
    ASSERT_EQ(molid[GETIDX(29)], 8);
}

TEST_F(ResetIDsTest, DeleteAdd)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();

    auto molid = lmp->atom->molecule;

    // delete two water molecules
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("group allwater molecule 3:6");
    lmp->input->one("group twowater molecule 4:6:2");
    lmp->input->one("group nowater subtract all allwater");
    lmp->input->one("delete_atoms group twowater compress no bond yes mol yes");
    lmp->input->one("reset_mol_ids allwater");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->atom->natoms, 23);
    ASSERT_EQ(lmp->atom->map_tag_max, 26);

    ASSERT_EQ(molid[GETIDX(1)], 1);
    ASSERT_EQ(molid[GETIDX(2)], 1);
    ASSERT_EQ(molid[GETIDX(3)], 1);
    ASSERT_EQ(molid[GETIDX(4)], 1);
    ASSERT_EQ(molid[GETIDX(5)], 1);
    ASSERT_EQ(molid[GETIDX(6)], 1);
    ASSERT_EQ(molid[GETIDX(7)], 1);
    ASSERT_EQ(molid[GETIDX(8)], 2);
    ASSERT_EQ(molid[GETIDX(9)], 2);
    ASSERT_EQ(molid[GETIDX(10)], 2);
    ASSERT_EQ(molid[GETIDX(11)], 2);
    ASSERT_EQ(molid[GETIDX(12)], 2);
    ASSERT_EQ(molid[GETIDX(13)], 2);
    ASSERT_EQ(molid[GETIDX(14)], 2);
    ASSERT_EQ(molid[GETIDX(15)], 2);
    ASSERT_EQ(molid[GETIDX(16)], 2);
    ASSERT_EQ(molid[GETIDX(17)], 2);
    ASSERT_EQ(molid[GETIDX(18)], 3);
    ASSERT_EQ(molid[GETIDX(19)], 3);
    ASSERT_EQ(molid[GETIDX(20)], 3);
    ASSERT_EQ(molid[GETIDX(24)], 4);
    ASSERT_EQ(molid[GETIDX(25)], 4);
    ASSERT_EQ(molid[GETIDX(26)], 4);

    // now also check and reset the atom ids

    ASSERT_GE(GETIDX(1), 0);
    ASSERT_GE(GETIDX(2), 0);
    ASSERT_GE(GETIDX(3), 0);
    ASSERT_GE(GETIDX(4), 0);
    ASSERT_GE(GETIDX(5), 0);
    ASSERT_GE(GETIDX(6), 0);
    ASSERT_GE(GETIDX(7), 0);
    ASSERT_GE(GETIDX(8), 0);
    ASSERT_GE(GETIDX(9), 0);
    ASSERT_GE(GETIDX(10), 0);
    ASSERT_GE(GETIDX(11), 0);
    ASSERT_GE(GETIDX(12), 0);
    ASSERT_GE(GETIDX(13), 0);
    ASSERT_GE(GETIDX(14), 0);
    ASSERT_GE(GETIDX(15), 0);
    ASSERT_GE(GETIDX(16), 0);
    ASSERT_GE(GETIDX(17), 0);
    ASSERT_GE(GETIDX(18), 0);
    ASSERT_GE(GETIDX(19), 0);
    ASSERT_GE(GETIDX(20), 0);
    ASSERT_EQ(GETIDX(21), -1);
    ASSERT_EQ(GETIDX(22), -1);
    ASSERT_EQ(GETIDX(23), -1);
    ASSERT_GE(GETIDX(24), 0);
    ASSERT_GE(GETIDX(25), 0);
    ASSERT_GE(GETIDX(26), 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("reset_atom_ids sort yes");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_EQ(lmp->atom->map_tag_max, 23);
    for (int i = 1; i <= 23; ++i)
        ASSERT_GE(GETIDX(i), 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("reset_mol_ids nowater offset 1");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_EQ(molid[GETIDX(1)], 2);
    ASSERT_EQ(molid[GETIDX(2)], 2);
    ASSERT_EQ(molid[GETIDX(3)], 2);
    ASSERT_EQ(molid[GETIDX(4)], 2);
    ASSERT_EQ(molid[GETIDX(5)], 2);
    ASSERT_EQ(molid[GETIDX(6)], 2);
    ASSERT_EQ(molid[GETIDX(7)], 2);
    ASSERT_EQ(molid[GETIDX(8)], 2);
    ASSERT_EQ(molid[GETIDX(9)], 2);
    ASSERT_EQ(molid[GETIDX(10)], 2);
    ASSERT_EQ(molid[GETIDX(11)], 2);
    ASSERT_EQ(molid[GETIDX(12)], 2);
    ASSERT_EQ(molid[GETIDX(13)], 3);
    ASSERT_EQ(molid[GETIDX(14)], 3);
    ASSERT_EQ(molid[GETIDX(15)], 3);
    ASSERT_EQ(molid[GETIDX(16)], 2);
    ASSERT_EQ(molid[GETIDX(17)], 2);
    ASSERT_EQ(molid[GETIDX(18)], 2);
    ASSERT_EQ(molid[GETIDX(19)], 2);
    ASSERT_EQ(molid[GETIDX(20)], 2);
    ASSERT_EQ(molid[GETIDX(21)], 4);
    ASSERT_EQ(molid[GETIDX(22)], 4);
    ASSERT_EQ(molid[GETIDX(23)], 4);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("create_atoms 1 single 0.0 0.0 0.0");
    lmp->input->one("create_atoms 2 single 1.0 0.0 0.0");
    lmp->input->one("create_atoms 3 single 2.0 0.0 0.0");
    lmp->input->one("create_atoms 4 single 3.0 0.0 0.0");
    lmp->input->one("reset_mol_ids all single yes");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->atom->natoms, 27);
    ASSERT_EQ(lmp->atom->map_tag_max, 27);

    ASSERT_EQ(molid[GETIDX(21)], 3);
    ASSERT_EQ(molid[GETIDX(22)], 3);
    ASSERT_EQ(molid[GETIDX(23)], 3);
    ASSERT_EQ(molid[GETIDX(24)], 4);
    ASSERT_EQ(molid[GETIDX(25)], 5);
    ASSERT_EQ(molid[GETIDX(26)], 6);
    ASSERT_EQ(molid[GETIDX(27)], 7);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("reset_mol_ids all single no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_EQ(molid[GETIDX(21)], 3);
    ASSERT_EQ(molid[GETIDX(22)], 3);
    ASSERT_EQ(molid[GETIDX(23)], 3);
    ASSERT_EQ(molid[GETIDX(24)], 0);
    ASSERT_EQ(molid[GETIDX(25)], 0);
    ASSERT_EQ(molid[GETIDX(26)], 0);
    ASSERT_EQ(molid[GETIDX(27)], 0);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("reset_mol_ids all compress no single yes");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_EQ(molid[GETIDX(21)], 21);
    ASSERT_EQ(molid[GETIDX(22)], 21);
    ASSERT_EQ(molid[GETIDX(23)], 21);
    ASSERT_EQ(molid[GETIDX(24)], 24);
    ASSERT_EQ(molid[GETIDX(25)], 25);
    ASSERT_EQ(molid[GETIDX(26)], 26);
    ASSERT_EQ(molid[GETIDX(27)], 27);
}

TEST_F(ResetIDsTest, TopologyData)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();

    // delete two water molecules
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("group allwater molecule 3:6");
    lmp->input->one("group twowater molecule 4:6:2");
    lmp->input->one("group nowater subtract all allwater");
    lmp->input->one("delete_atoms group twowater compress no bond yes mol yes");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->atom->natoms, 23);
    ASSERT_EQ(lmp->atom->map_tag_max, 26);

    auto num_bond     = lmp->atom->num_bond;
    auto num_angle    = lmp->atom->num_angle;
    auto bond_atom    = lmp->atom->bond_atom;
    auto angle_atom1  = lmp->atom->angle_atom1;
    auto angle_atom2  = lmp->atom->angle_atom2;
    auto angle_atom3  = lmp->atom->angle_atom3;
    ASSERT_EQ(num_bond[GETIDX(1)], 2);
    ASSERT_EQ(bond_atom[GETIDX(1)][0], 2);
    ASSERT_EQ(bond_atom[GETIDX(1)][1], 3);
    ASSERT_EQ(num_bond[GETIDX(2)], 0);
    ASSERT_EQ(num_bond[GETIDX(3)], 3);
    ASSERT_EQ(bond_atom[GETIDX(3)][0], 4);
    ASSERT_EQ(bond_atom[GETIDX(3)][1], 5);
    ASSERT_EQ(bond_atom[GETIDX(3)][2], 6);
    ASSERT_EQ(num_bond[GETIDX(4)], 0);
    ASSERT_EQ(num_bond[GETIDX(5)], 0);
    ASSERT_EQ(num_bond[GETIDX(6)], 2);
    ASSERT_EQ(bond_atom[GETIDX(6)][0], 8);
    ASSERT_EQ(bond_atom[GETIDX(6)][1], 7);
    ASSERT_EQ(num_bond[GETIDX(7)], 0);
    ASSERT_EQ(num_bond[GETIDX(8)], 2);
    ASSERT_EQ(bond_atom[GETIDX(8)][0], 9);
    ASSERT_EQ(bond_atom[GETIDX(8)][1], 10);
    ASSERT_EQ(num_bond[GETIDX(9)], 0);
    ASSERT_EQ(num_bond[GETIDX(10)], 3);
    ASSERT_EQ(bond_atom[GETIDX(10)][0], 11);
    ASSERT_EQ(bond_atom[GETIDX(10)][1], 12);
    ASSERT_EQ(bond_atom[GETIDX(10)][2], 16);
    ASSERT_EQ(num_bond[GETIDX(11)], 0);
    ASSERT_EQ(num_bond[GETIDX(12)], 3);
    ASSERT_EQ(bond_atom[GETIDX(12)][0], 13);
    ASSERT_EQ(bond_atom[GETIDX(12)][1], 14);
    ASSERT_EQ(bond_atom[GETIDX(12)][2], 15);
    ASSERT_EQ(num_bond[GETIDX(13)], 0);
    ASSERT_EQ(num_bond[GETIDX(14)], 0);
    ASSERT_EQ(num_bond[GETIDX(15)], 0);
    ASSERT_EQ(num_bond[GETIDX(16)], 1);
    ASSERT_EQ(bond_atom[GETIDX(16)][0], 17);
    ASSERT_EQ(num_bond[GETIDX(17)], 0);
    ASSERT_EQ(num_bond[GETIDX(18)], 2);
    ASSERT_EQ(bond_atom[GETIDX(18)][0], 19);
    ASSERT_EQ(bond_atom[GETIDX(18)][1], 20);
    ASSERT_EQ(num_bond[GETIDX(19)], 0);
    ASSERT_EQ(num_bond[GETIDX(20)], 0);
    ASSERT_EQ(num_bond[GETIDX(24)], 2);
    ASSERT_EQ(bond_atom[GETIDX(24)][0], 25);
    ASSERT_EQ(bond_atom[GETIDX(24)][1], 26);
    ASSERT_EQ(num_bond[GETIDX(25)], 0);
    ASSERT_EQ(num_bond[GETIDX(26)], 0);

    ASSERT_EQ(num_angle[GETIDX(1)], 1);
    ASSERT_EQ(angle_atom1[GETIDX(1)][0], 2);
    ASSERT_EQ(angle_atom2[GETIDX(1)][0], 1);
    ASSERT_EQ(angle_atom3[GETIDX(1)][0], 3);
    ASSERT_EQ(num_angle[GETIDX(2)], 0);
    ASSERT_EQ(num_angle[GETIDX(3)], 6);
    ASSERT_EQ(angle_atom1[GETIDX(3)][0], 1);
    ASSERT_EQ(angle_atom2[GETIDX(3)][0], 3);
    ASSERT_EQ(angle_atom3[GETIDX(3)][0], 5);
    ASSERT_EQ(angle_atom1[GETIDX(3)][1], 1);
    ASSERT_EQ(angle_atom2[GETIDX(3)][1], 3);
    ASSERT_EQ(angle_atom3[GETIDX(3)][1], 4);
    ASSERT_EQ(angle_atom1[GETIDX(3)][2], 1);
    ASSERT_EQ(angle_atom2[GETIDX(3)][2], 3);
    ASSERT_EQ(angle_atom3[GETIDX(3)][2], 6);
    ASSERT_EQ(angle_atom1[GETIDX(3)][3], 4);
    ASSERT_EQ(angle_atom2[GETIDX(3)][3], 3);
    ASSERT_EQ(angle_atom3[GETIDX(3)][3], 5);
    ASSERT_EQ(angle_atom1[GETIDX(3)][4], 5);
    ASSERT_EQ(angle_atom2[GETIDX(3)][4], 3);
    ASSERT_EQ(angle_atom3[GETIDX(3)][4], 6);
    ASSERT_EQ(num_angle[GETIDX(18)], 1);
    ASSERT_EQ(angle_atom1[GETIDX(18)][0], 19);
    ASSERT_EQ(angle_atom2[GETIDX(18)][0], 18);
    ASSERT_EQ(angle_atom3[GETIDX(18)][0], 20);
    ASSERT_EQ(num_angle[GETIDX(24)], 1);
    ASSERT_EQ(angle_atom1[GETIDX(24)][0], 25);
    ASSERT_EQ(angle_atom2[GETIDX(24)][0], 24);
    ASSERT_EQ(angle_atom3[GETIDX(24)][0], 26);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("reset_atom_ids sort yes");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    num_bond     = lmp->atom->num_bond;
    num_angle    = lmp->atom->num_angle;
    bond_atom    = lmp->atom->bond_atom;
    angle_atom1  = lmp->atom->angle_atom1;
    angle_atom2  = lmp->atom->angle_atom2;
    angle_atom3  = lmp->atom->angle_atom3;
    ASSERT_EQ(num_bond[GETIDX(1)], 2);
    ASSERT_EQ(bond_atom[GETIDX(1)][0], 3);
    ASSERT_EQ(bond_atom[GETIDX(1)][1], 2);
    ASSERT_EQ(num_bond[GETIDX(2)], 0);
    ASSERT_EQ(num_bond[GETIDX(3)], 2);
    ASSERT_EQ(bond_atom[GETIDX(3)][0], 16);
    ASSERT_EQ(bond_atom[GETIDX(3)][1], 5);
    ASSERT_EQ(num_bond[GETIDX(4)], 0);
    ASSERT_EQ(num_bond[GETIDX(5)], 3);
    ASSERT_EQ(bond_atom[GETIDX(5)][0], 4);
    ASSERT_EQ(bond_atom[GETIDX(5)][1], 8);
    ASSERT_EQ(bond_atom[GETIDX(5)][2], 18);
    ASSERT_EQ(num_bond[GETIDX(6)], 0);
    ASSERT_EQ(num_bond[GETIDX(7)], 0);
    ASSERT_EQ(num_bond[GETIDX(8)], 3);
    ASSERT_EQ(bond_atom[GETIDX(8)][0], 9);
    ASSERT_EQ(bond_atom[GETIDX(8)][1], 6);
    ASSERT_EQ(bond_atom[GETIDX(8)][2], 7);
    ASSERT_EQ(num_bond[GETIDX(9)], 0);
    ASSERT_EQ(num_bond[GETIDX(10)], 0);
    ASSERT_EQ(num_bond[GETIDX(11)], 3);
    ASSERT_EQ(bond_atom[GETIDX(11)][0], 10);
    ASSERT_EQ(bond_atom[GETIDX(11)][1], 19);
    ASSERT_EQ(bond_atom[GETIDX(11)][2], 1);
    ASSERT_EQ(num_bond[GETIDX(12)], 0);
    ASSERT_EQ(num_bond[GETIDX(13)], 0);
    ASSERT_EQ(num_bond[GETIDX(14)], 2);
    ASSERT_EQ(bond_atom[GETIDX(14)][0], 13);
    ASSERT_EQ(bond_atom[GETIDX(14)][1], 15);
    ASSERT_EQ(num_bond[GETIDX(15)], 0);
    ASSERT_EQ(num_bond[GETIDX(16)], 0);
    ASSERT_EQ(num_bond[GETIDX(17)], 0);
    ASSERT_EQ(num_bond[GETIDX(18)], 1);
    ASSERT_EQ(bond_atom[GETIDX(18)][0], 17);
    ASSERT_EQ(num_bond[GETIDX(19)], 0);
    ASSERT_EQ(num_bond[GETIDX(20)], 2);
    ASSERT_EQ(bond_atom[GETIDX(20)][0], 12);
    ASSERT_EQ(bond_atom[GETIDX(20)][1], 11);
    ASSERT_EQ(num_bond[GETIDX(21)], 0);
    ASSERT_EQ(num_bond[GETIDX(22)], 2);
    ASSERT_EQ(bond_atom[GETIDX(22)][0], 21);
    ASSERT_EQ(bond_atom[GETIDX(22)][1], 23);
    ASSERT_EQ(num_bond[GETIDX(23)], 0);

    ASSERT_EQ(num_angle[GETIDX(1)], 3);
    ASSERT_EQ(angle_atom1[GETIDX(1)][0], 11);
    ASSERT_EQ(angle_atom2[GETIDX(1)][0], 1);
    ASSERT_EQ(angle_atom3[GETIDX(1)][0], 2);
    ASSERT_EQ(angle_atom1[GETIDX(1)][1], 11);
    ASSERT_EQ(angle_atom2[GETIDX(1)][1], 1);
    ASSERT_EQ(angle_atom3[GETIDX(1)][1], 3);
    ASSERT_EQ(angle_atom1[GETIDX(1)][2], 2);
    ASSERT_EQ(angle_atom2[GETIDX(1)][2], 1);
    ASSERT_EQ(angle_atom3[GETIDX(1)][2], 3);
    ASSERT_EQ(num_angle[GETIDX(2)], 0);
    ASSERT_EQ(num_angle[GETIDX(5)], 6);
    ASSERT_EQ(angle_atom1[GETIDX(5)][0], 3);
    ASSERT_EQ(angle_atom2[GETIDX(5)][0], 5);
    ASSERT_EQ(angle_atom3[GETIDX(5)][0], 4);
    ASSERT_EQ(angle_atom1[GETIDX(5)][1], 3);
    ASSERT_EQ(angle_atom2[GETIDX(5)][1], 5);
    ASSERT_EQ(angle_atom3[GETIDX(5)][1], 18);
    ASSERT_EQ(angle_atom1[GETIDX(5)][2], 4);
    ASSERT_EQ(angle_atom2[GETIDX(5)][2], 5);
    ASSERT_EQ(angle_atom3[GETIDX(5)][2], 8);
    ASSERT_EQ(angle_atom1[GETIDX(5)][3], 8);
    ASSERT_EQ(angle_atom2[GETIDX(5)][3], 5);
    ASSERT_EQ(angle_atom3[GETIDX(5)][3], 18);
    ASSERT_EQ(angle_atom1[GETIDX(5)][4], 3);
    ASSERT_EQ(angle_atom2[GETIDX(5)][4], 5);
    ASSERT_EQ(angle_atom3[GETIDX(5)][4], 8);
    ASSERT_EQ(angle_atom1[GETIDX(5)][5], 4);
    ASSERT_EQ(angle_atom2[GETIDX(5)][5], 5);
    ASSERT_EQ(angle_atom3[GETIDX(5)][5], 18);
    ASSERT_EQ(num_angle[GETIDX(20)], 1);
    ASSERT_EQ(angle_atom1[GETIDX(20)][0], 12);
    ASSERT_EQ(angle_atom2[GETIDX(20)][0], 20);
    ASSERT_EQ(angle_atom3[GETIDX(20)][0], 11);
    ASSERT_EQ(num_angle[GETIDX(22)], 1);
    ASSERT_EQ(angle_atom1[GETIDX(22)][0], 21);
    ASSERT_EQ(angle_atom2[GETIDX(22)][0], 22);
    ASSERT_EQ(angle_atom3[GETIDX(22)][0], 23);
}

TEST_F(ResetIDsTest, DeathTests)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();

    TEST_FAILURE(".*ERROR: Illegal reset_mol_ids command.*", lmp->input->one("reset_mol_ids"););
    TEST_FAILURE(".*ERROR: Illegal reset_mol_ids command.*",
                 lmp->input->one("reset_mol_ids all offset 1 1"););
    TEST_FAILURE(".*ERROR: Illegal reset_mol_ids command.*",
                 lmp->input->one("reset_mol_ids all offset -2"););
    TEST_FAILURE(".*ERROR on proc 0: Expected integer.*",
                 lmp->input->one("reset_mol_ids all offset xxx"););
    TEST_FAILURE(".*ERROR on proc 0: Expected integer.*",
                 lmp->input->one("reset_mol_ids all compress yes single no offset xxx"););
    TEST_FAILURE(".*ERROR: Illegal reset_mol_ids command.*",
                 lmp->input->one("reset_mol_ids all offset"););
    TEST_FAILURE(".*ERROR: Illegal reset_mol_ids command.*",
                 lmp->input->one("reset_mol_ids all compress"););

    TEST_FAILURE(".*ERROR: Illegal reset_mol_ids command.*",
                 lmp->input->one("reset_mol_ids all compress xxx"););
    TEST_FAILURE(".*ERROR: Illegal reset_mol_ids command.*",
                 lmp->input->one("reset_mol_ids all single"););
    TEST_FAILURE(".*ERROR: Illegal reset_mol_ids command.*",
                 lmp->input->one("reset_mol_ids all single xxx"););
}

TEST(ResetMolIds, CMDFail)
{
    LAMMPS *lmp;
    const char *args[] = {"ResetIDsTest", "-log", "none", "-nocite", "-echo", "screen"};
    char **argv        = (char **)args;
    int argc           = sizeof(args) / sizeof(char *);
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
    if (!verbose) ::testing::internal::GetCapturedStdout();

    TEST_FAILURE(".*ERROR: Reset_mol_ids command before simulation box is.*",
                 lmp->input->one("reset_mol_ids all"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("atom_modify id no");
    lmp->input->one("region box block 0 1 0 1 0 1");
    lmp->input->one("create_box 1 box");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Cannot use reset_mol_ids unless.*",
                 lmp->input->one("reset_mol_ids all"););

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("clear");
    lmp->input->one("region box block 0 1 0 1 0 1");
    lmp->input->one("create_box 1 box");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    TEST_FAILURE(".*ERROR: Can only use reset_mol_ids.*", lmp->input->one("reset_mol_ids all"););

    if (!verbose) ::testing::internal::CaptureStdout();
    delete lmp;
    if (!verbose) ::testing::internal::GetCapturedStdout();
};

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
