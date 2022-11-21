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

#include "../testing/core.h"
#include "atom.h"
#include "fmt/format.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "library.h"
#include "output.h"
#include "update.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

#define GETIDX(i) lmp->atom->map(i)

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

class ResetAtomsIDTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "ResetAtomsIDTest";
        LAMMPSTest::SetUp();
        if (info->has_style("atom", "full")) {
            BEGIN_HIDE_OUTPUT();
            command("variable input_dir index \"" STRINGIFY(TEST_INPUT_FOLDER) "\"");
            command("include \"${input_dir}/in.fourmol\"");
            END_HIDE_OUTPUT();
        }
    }
};

TEST_F(ResetAtomsIDTest, MolIDAll)
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
    BEGIN_HIDE_OUTPUT();
    command("reset_atoms mol all");
    END_HIDE_OUTPUT();

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

TEST_F(ResetAtomsIDTest, DeletePlusAtomID)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();

    auto molid = lmp->atom->molecule;

    // delete two water molecules
    BEGIN_HIDE_OUTPUT();
    command("group allwater molecule 3:6");
    command("group twowater molecule 4:6:2");
    command("delete_atoms group twowater compress no bond yes");
    command("reset_atoms mol all");
    END_HIDE_OUTPUT();
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

    BEGIN_HIDE_OUTPUT();
    command("reset_atoms id");
    END_HIDE_OUTPUT();

    ASSERT_EQ(lmp->atom->map_tag_max, 23);
    for (int i = 1; i <= 23; ++i)
        ASSERT_GE(GETIDX(i), 0);
}

TEST_F(ResetAtomsIDTest, PartialOffset)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();

    auto molid = lmp->atom->molecule;

    // delete two water molecules
    BEGIN_HIDE_OUTPUT();
    command("group allwater molecule 3:6");
    command("group nowater subtract all allwater");
    command("reset_atoms mol allwater offset 4");
    END_HIDE_OUTPUT();
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

    BEGIN_HIDE_OUTPUT();
    command("reset_atoms mol nowater offset 0");
    END_HIDE_OUTPUT();

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

TEST_F(ResetAtomsIDTest, DeleteAdd)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();

    auto molid = lmp->atom->molecule;

    // delete two water molecules
    BEGIN_HIDE_OUTPUT();
    command("group allwater molecule 3:6");
    command("group twowater molecule 4:6:2");
    command("group nowater subtract all allwater");
    command("delete_atoms group twowater compress no bond yes mol yes");
    command("reset_atoms mol allwater");
    END_HIDE_OUTPUT();
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

    BEGIN_HIDE_OUTPUT();
    command("reset_atoms id sort yes");
    END_HIDE_OUTPUT();

    ASSERT_EQ(lmp->atom->map_tag_max, 23);
    for (int i = 1; i <= 23; ++i)
        ASSERT_GE(GETIDX(i), 0);

    BEGIN_HIDE_OUTPUT();
    command("reset_atoms mol nowater offset 1");
    END_HIDE_OUTPUT();

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

    BEGIN_HIDE_OUTPUT();
    command("create_atoms 1 single 0.0 0.0 0.0");
    command("create_atoms 2 single 1.0 0.0 0.0");
    command("create_atoms 3 single 2.0 0.0 0.0");
    command("create_atoms 4 single 3.0 0.0 0.0");
    command("reset_atoms mol all single yes");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->atom->natoms, 27);
    ASSERT_EQ(lmp->atom->map_tag_max, 27);

    ASSERT_EQ(molid[GETIDX(21)], 3);
    ASSERT_EQ(molid[GETIDX(22)], 3);
    ASSERT_EQ(molid[GETIDX(23)], 3);
    ASSERT_EQ(molid[GETIDX(24)], 4);
    ASSERT_EQ(molid[GETIDX(25)], 5);
    ASSERT_EQ(molid[GETIDX(26)], 6);
    ASSERT_EQ(molid[GETIDX(27)], 7);

    BEGIN_HIDE_OUTPUT();
    command("reset_atoms mol all single no");
    END_HIDE_OUTPUT();

    ASSERT_EQ(molid[GETIDX(21)], 3);
    ASSERT_EQ(molid[GETIDX(22)], 3);
    ASSERT_EQ(molid[GETIDX(23)], 3);
    ASSERT_EQ(molid[GETIDX(24)], 0);
    ASSERT_EQ(molid[GETIDX(25)], 0);
    ASSERT_EQ(molid[GETIDX(26)], 0);
    ASSERT_EQ(molid[GETIDX(27)], 0);

    BEGIN_HIDE_OUTPUT();
    command("reset_atoms mol all compress no single yes");
    END_HIDE_OUTPUT();

    ASSERT_EQ(molid[GETIDX(21)], 21);
    ASSERT_EQ(molid[GETIDX(22)], 21);
    ASSERT_EQ(molid[GETIDX(23)], 21);
    ASSERT_EQ(molid[GETIDX(24)], 24);
    ASSERT_EQ(molid[GETIDX(25)], 25);
    ASSERT_EQ(molid[GETIDX(26)], 26);
    ASSERT_EQ(molid[GETIDX(27)], 27);
}

TEST_F(ResetAtomsIDTest, TopologyData)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();

    // delete two water molecules
    BEGIN_HIDE_OUTPUT();
    command("group allwater molecule 3:6");
    command("group twowater molecule 4:6:2");
    command("group nowater subtract all allwater");
    command("delete_atoms group twowater compress no bond yes mol yes");
    END_HIDE_OUTPUT();
    ASSERT_EQ(lmp->atom->natoms, 23);
    ASSERT_EQ(lmp->atom->map_tag_max, 26);

    auto num_bond    = lmp->atom->num_bond;
    auto num_angle   = lmp->atom->num_angle;
    auto bond_atom   = lmp->atom->bond_atom;
    auto angle_atom1 = lmp->atom->angle_atom1;
    auto angle_atom2 = lmp->atom->angle_atom2;
    auto angle_atom3 = lmp->atom->angle_atom3;
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

    BEGIN_HIDE_OUTPUT();
    command("reset_atoms id sort yes");
    END_HIDE_OUTPUT();

    num_bond    = lmp->atom->num_bond;
    num_angle   = lmp->atom->num_angle;
    bond_atom   = lmp->atom->bond_atom;
    angle_atom1 = lmp->atom->angle_atom1;
    angle_atom2 = lmp->atom->angle_atom2;
    angle_atom3 = lmp->atom->angle_atom3;
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

TEST_F(ResetAtomsIDTest, DeathTests)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();

    TEST_FAILURE(".*ERROR: Illegal reset_atoms mol command.*", command("reset_atoms mol"););
    TEST_FAILURE(".*ERROR: Unknown reset_atoms mol keyword: 1.*",
                 command("reset_atoms mol all offset 1 1"););
    TEST_FAILURE(".*ERROR: Illegal reset_atoms mol offset: -2.*",
                 command("reset_atoms mol all offset -2"););
    TEST_FAILURE(".*ERROR on proc 0: Expected integer.*",
                 command("reset_atoms mol all offset xxx"););
    TEST_FAILURE(".*ERROR on proc 0: Expected integer.*",
                 command("reset_atoms mol all compress yes single no offset xxx"););
    TEST_FAILURE(".*ERROR: Illegal reset_atoms mol offset command: missing argument.*",
                 command("reset_atoms mol all offset"););
    TEST_FAILURE(".*ERROR: Illegal reset_atoms mol compress command: missing argument.*",
                 command("reset_atoms mol all compress"););

    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of 'xxx'.*",
                 command("reset_atoms mol all compress xxx"););
    TEST_FAILURE(".*ERROR: Illegal reset_atoms mol single command: missing argument.*",
                 command("reset_atoms mol all single"););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of 'xxx'.*",
                 command("reset_atoms mol all single xxx"););

    TEST_FAILURE(".*ERROR: Illegal reset_atoms image command: missing argument.*",
                 command("reset_atoms image"););
    TEST_FAILURE(".*ERROR: Unknown reset_atoms image keyword: xxx.*",
                 command("reset_atoms image all xxx"););
    TEST_FAILURE(".*ERROR: Could not find reset_atoms image group xxx.*",
                 command("reset_atoms image xxx"););
}

class ResetAtomsMolTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "ResetAtomsMolTest";
        LAMMPSTest::SetUp();
    }
};

TEST_F(ResetAtomsMolTest, FailBeforeBox)
{
    TEST_FAILURE(".*ERROR: Reset_atoms id command before simulation box is.*",
                 command("reset_atoms id"););
    TEST_FAILURE(".*ERROR: Reset_atoms mol command before simulation box is.*",
                 command("reset_atoms mol all"););
    TEST_FAILURE(".*ERROR: Reset_atoms image command before simulation box is.*",
                 command("reset_atoms image all"););
}

TEST_F(ResetAtomsMolTest, FailMissingId)
{
    BEGIN_HIDE_OUTPUT();
    command("atom_modify id no");
    command("region box block 0 1 0 1 0 1");
    command("create_box 1 box");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*ERROR: Cannot use reset_atoms mol unless.*", command("reset_atoms mol all"););
    TEST_FAILURE(".*ERROR: Cannot use reset_atoms image unless.*",
                 command("reset_atoms image all"););
}

TEST_F(ResetAtomsMolTest, FailOnlyMolecular)
{
    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("region box block 0 1 0 1 0 1");
    command("create_box 1 box");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*ERROR: Can only use reset_atoms mol.*", command("reset_atoms mol all"););
}

class ResetAtomsImageTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "ResetAtomsImageTest";
        LAMMPSTest::SetUp();
        if (info->has_style("atom", "full")) {
            BEGIN_HIDE_OUTPUT();
            command("variable input_dir index \"" STRINGIFY(TEST_INPUT_FOLDER) "\"");
            command("include \"${input_dir}/in.fourmol\"");
            command("create_atoms 1 single 1.0 1.0 1.0");
            command("create_atoms 1 single 2.0 1.0 1.0");
            command("create_atoms 1 single 1.0 2.0 1.0");
            command("set atom 1*7 image 1 1 -1");
            command("set atom 8*9 image 2 -1 0");
            command("set atom 8*9 image 2 -1 0");
            command("set atom 10 image 2 0 0");
            command("set atom 11*12 image 2 0 -1");
            command("set atom 13*15 image 1 0 0");
            command("set atom 16*17 image 0 1 1");
            command("set atom 18*19 image 0 1 0");
            command("set atom 20 image 0 2 0");
            command("set atom 21*23 image 20 -1 0");
            command("set atom 27*28 image 1 0 0");
            command("set atom 29 image 1 2 0");
            command("set atom 31 image -2 0 1");
            command("set atom 32 image 0 20 0");
            END_HIDE_OUTPUT();
        }
    }
};

TEST_F(ResetAtomsImageTest, ResetAtomsImage)
{
    if (lmp->atom->natoms == 0) GTEST_SKIP();
    EXPECT_EQ(lmp->atom->image[GETIDX(1)], lammps_encode_image_flags(1, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(2)], lammps_encode_image_flags(1, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(3)], lammps_encode_image_flags(1, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(4)], lammps_encode_image_flags(1, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(5)], lammps_encode_image_flags(1, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(6)], lammps_encode_image_flags(1, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(7)], lammps_encode_image_flags(1, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(8)], lammps_encode_image_flags(2, -1, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(9)], lammps_encode_image_flags(2, -1, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(10)], lammps_encode_image_flags(2, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(11)], lammps_encode_image_flags(2, 0, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(12)], lammps_encode_image_flags(2, 0, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(13)], lammps_encode_image_flags(1, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(14)], lammps_encode_image_flags(1, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(15)], lammps_encode_image_flags(1, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(16)], lammps_encode_image_flags(0, 1, 1));
    EXPECT_EQ(lmp->atom->image[GETIDX(17)], lammps_encode_image_flags(0, 1, 1));
    EXPECT_EQ(lmp->atom->image[GETIDX(18)], lammps_encode_image_flags(0, 1, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(19)], lammps_encode_image_flags(0, 1, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(20)], lammps_encode_image_flags(0, 2, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(21)], lammps_encode_image_flags(20, -1, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(22)], lammps_encode_image_flags(20, -1, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(23)], lammps_encode_image_flags(20, -1, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(24)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(25)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(26)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(27)], lammps_encode_image_flags(1, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(28)], lammps_encode_image_flags(1, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(29)], lammps_encode_image_flags(1, 2, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(30)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(31)], lammps_encode_image_flags(-2, 0, 1));
    EXPECT_EQ(lmp->atom->image[GETIDX(32)], lammps_encode_image_flags(0, 20, 0));
    BEGIN_HIDE_OUTPUT();
    command("group subset id 5:32");
    command("reset_atoms image subset");
    END_HIDE_OUTPUT();
    EXPECT_EQ(lmp->atom->image[GETIDX(1)], lammps_encode_image_flags(1, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(2)], lammps_encode_image_flags(1, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(3)], lammps_encode_image_flags(1, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(4)], lammps_encode_image_flags(1, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(5)], lammps_encode_image_flags(0, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(6)], lammps_encode_image_flags(0, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(7)], lammps_encode_image_flags(0, 1, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(8)], lammps_encode_image_flags(1, -1, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(9)], lammps_encode_image_flags(1, -1, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(10)], lammps_encode_image_flags(1, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(11)], lammps_encode_image_flags(1, 0, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(12)], lammps_encode_image_flags(1, 0, -1));
    EXPECT_EQ(lmp->atom->image[GETIDX(13)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(14)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(15)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(16)], lammps_encode_image_flags(-1, 1, 1));
    EXPECT_EQ(lmp->atom->image[GETIDX(17)], lammps_encode_image_flags(-1, 1, 1));
    EXPECT_EQ(lmp->atom->image[GETIDX(18)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(19)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(20)], lammps_encode_image_flags(0, 1, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(21)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(22)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(23)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(24)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(25)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(26)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(27)], lammps_encode_image_flags(0, -1, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(28)], lammps_encode_image_flags(0, -1, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(29)], lammps_encode_image_flags(0, 1, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(30)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(30)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(31)], lammps_encode_image_flags(0, 0, 0));
    EXPECT_EQ(lmp->atom->image[GETIDX(32)], lammps_encode_image_flags(0, 0, 0));
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (LAMMPS_NS::platform::mpi_vendor() == "Open MPI" && !Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. Death tests will be skipped\n";

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
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
