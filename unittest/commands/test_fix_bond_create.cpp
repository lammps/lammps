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
#include "fix.h"
#include "info.h"
#include "lammps.h"
#include "library.h"
#include "modify.h"
#include "utils.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include <string>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

using ::testing::StrEq;

// fix bond/create and fix bond/create/angle are tested with small systems that
// are created on-the-fly.  there is no time integration fix, so the atoms do
// not move during a run and the bonds that are created are fully determined by
// the geometry and the fix settings.

class FixBondCreateTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "FixBondCreateTest";
        LAMMPSTest::SetUp();
    }

    bool have_styles() const
    {
        return info->has_style("atom", "molecular") && info->has_style("fix", "bond/create");
    }

    // periodic 10x10x10 box with 3 atom types and pair style cutoff 2.0.
    // newton_bond must be set before the box is created.

    void create_system(const std::string &newton_bond = "on")
    {
        BEGIN_HIDE_OUTPUT();
        command("newton on " + newton_bond);
        command("units lj");
        command("atom_style molecular");
        command("region box block 0 10 0 10 0 10");
        command("create_box 3 box bond/types 2 angle/types 1 dihedral/types 1 improper/types 1 "
                "extra/bond/per/atom 4 extra/angle/per/atom 6 extra/dihedral/per/atom 9 "
                "extra/improper/per/atom 9 extra/special/per/atom 12");
        command("mass * 1.0");
        command("pair_style zero 2.0");
        command("pair_coeff * *");
        command("bond_style zero");
        command("bond_coeff *");
        command("angle_style zero");
        command("angle_coeff *");
        command("dihedral_style zero");
        command("dihedral_coeff *");
        command("improper_style zero");
        command("improper_coeff *");
        END_HIDE_OUTPUT();
    }

    // atom IDs are assigned consecutively starting from 1

    void add_atom(int type, double x, double y, double z)
    {
        command(fmt::format("create_atoms {} single {} {} {} units box", type, x, y, z));
    }

    // sorted list of all bonds, angles, dihedrals, or impropers as strings
    // "type:atom1-...-atomN".  except for impropers, the order of the atoms
    // is reversed when needed so that the first atom ID is lower than the last.

    std::string topology_list(const std::string &kind)
    {
        int natom  = 4;
        bigint num = 0;
        if (kind == "bonds") {
            natom = 2;
            num   = lmp->atom->nbonds;
        } else if (kind == "angles") {
            natom = 3;
            num   = lmp->atom->nangles;
        } else if (kind == "dihedrals") {
            num = lmp->atom->ndihedrals;
        } else if (kind == "impropers") {
            num = lmp->atom->nimpropers;
        }

        // with newton_bond off and inconsistently stored copies, the gather
        // functions may return up to natom entries for each of the num items.
        // reserve space for that, so that such a bug fails a test instead of
        // corrupting memory.

        std::vector<tagint> data((natom + 1) * natom * num);
        if (kind == "bonds") lammps_gather_bonds(lmp, data.data());
        if (kind == "angles") lammps_gather_angles(lmp, data.data());
        if (kind == "dihedrals") lammps_gather_dihedrals(lmp, data.data());
        if (kind == "impropers") lammps_gather_impropers(lmp, data.data());

        std::vector<std::string> entries;
        for (bigint i = 0; i < num; ++i) {
            auto *entry = &data[i * (natom + 1)];
            if ((kind != "impropers") && (entry[1] > entry[natom]))
                std::reverse(entry + 1, entry + natom + 1);
            auto text = fmt::format("{}:{}", entry[0], entry[1]);
            for (int j = 2; j <= natom; ++j)
                text += fmt::format("-{}", entry[j]);
            entries.push_back(text);
        }
        std::sort(entries.begin(), entries.end());
        return utils::join_words(entries, " ");
    }

    std::string bond_list() { return topology_list("bonds"); }

    std::vector<int> atom_types()
    {
        std::vector<int> types(lmp->atom->natoms);
        lammps_gather_atoms(lmp, "type", 0, 1, types.data());
        return types;
    }

    // number of bonds created on the last step (n = 0) and in total (n = 1)

    double created(int n) { return lmp->modify->get_fix_by_id("bc")->compute_vector(n); }
};

TEST_F(FixBondCreateTest, CreateBonds)
{
    if (!have_styles()) GTEST_SKIP();

    create_system();
    BEGIN_HIDE_OUTPUT();
    // eligible pair within the cutoff
    add_atom(1, 2.0, 2.0, 2.0);
    add_atom(2, 3.5, 2.0, 2.0);
    // eligible pair within the cutoff across a periodic boundary
    add_atom(1, 9.6, 5.0, 5.0);
    add_atom(2, 0.4, 5.0, 5.0);
    // eligible pair outside the cutoff
    add_atom(1, 2.0, 7.0, 7.0);
    add_atom(2, 2.0, 7.0, 8.7);
    // pair within the cutoff with the wrong atom types
    add_atom(1, 7.0, 2.0, 2.0);
    add_atom(1, 7.0, 3.0, 2.0);
    command("fix bc all bond/create 1 1 2 1.6 1");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_EQ(lmp->atom->nbonds, 2);
    EXPECT_THAT(bond_list(), StrEq("1:1-2 1:3-4"));
    EXPECT_EQ(created(0), 2.0);
    EXPECT_EQ(created(1), 2.0);

    // the same pair of atoms must not be bonded twice

    BEGIN_HIDE_OUTPUT();
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_EQ(lmp->atom->nbonds, 2);
    EXPECT_EQ(created(0), 0.0);
    EXPECT_EQ(created(1), 2.0);
}

TEST_F(FixBondCreateTest, TypeLabels)
{
    if (!have_styles()) GTEST_SKIP();

    create_system();
    BEGIN_HIDE_OUTPUT();
    command("labelmap atom 1 C 2 N 3 X");
    command("labelmap bond 1 C-C 2 C-N");
    add_atom(1, 2.0, 2.0, 2.0);
    add_atom(2, 3.0, 2.0, 2.0);
    command("fix bc all bond/create 1 C N 1.5 C-N");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("2:1-2"));
}

// atom 1 has two eligible partners: on each step an atom can form only one
// bond, and it is formed with the closest eligible partner

TEST_F(FixBondCreateTest, ClosestPartnerFirst)
{
    if (!have_styles()) GTEST_SKIP();

    create_system();
    BEGIN_HIDE_OUTPUT();
    add_atom(1, 5.0, 5.0, 5.0);
    add_atom(2, 6.0, 5.0, 5.0);
    add_atom(2, 5.0, 6.2, 5.0);
    command("fix bc all bond/create 1 1 2 1.5 1 atype 1");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-2"));
    EXPECT_EQ(lmp->atom->nangles, 0);

    BEGIN_HIDE_OUTPUT();
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-2 1:1-3"));
    EXPECT_THAT(topology_list("angles"), StrEq("1:2-1-3"));
}

// same geometry, but each atom may have only one bond and changes its type
// when it has reached that limit.  the atom types are collected with
// lammps_gather_atoms(), which is not compatible with -DLAMMPS_BIGBIG.

TEST_F(FixBondCreateTest, MaxBondsAndNewType)
{
#if defined(LAMMPS_BIGBIG)
    GTEST_SKIP() << "lammps_gather_atoms() is not compatible with -DLAMMPS_BIGBIG";
#endif
    if (!have_styles()) GTEST_SKIP();

    create_system();
    BEGIN_HIDE_OUTPUT();
    add_atom(1, 5.0, 5.0, 5.0);
    add_atom(2, 6.0, 5.0, 5.0);
    add_atom(2, 5.0, 6.2, 5.0);
    command("fix bc all bond/create 1 1 2 1.5 1 iparam 1 3 jparam 1 3");
    command("run 2 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-2"));
    EXPECT_EQ(created(1), 1.0);
    EXPECT_THAT(atom_types(), ::testing::ElementsAre(3, 3, 2));
}

// bonds of the bond type to be created that exist when the run starts count
// toward the limit, bonds of other types do not

TEST_F(FixBondCreateTest, ExistingBondsCount)
{
    if (!have_styles()) GTEST_SKIP();

    create_system();
    BEGIN_HIDE_OUTPUT();
    add_atom(1, 5.0, 5.0, 5.0);
    add_atom(2, 6.0, 5.0, 5.0);
    add_atom(3, 4.0, 5.0, 5.0);
    command("create_bonds single/bond 1 1 3");
    command("fix bc all bond/create 1 1 2 1.5 1 iparam 1 1");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-3"));

    BEGIN_HIDE_OUTPUT();
    command("unfix bc");
    command("fix bc all bond/create 1 1 2 1.5 2 iparam 1 1");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-3 2:1-2"));
}

TEST_F(FixBondCreateTest, Nevery)
{
    if (!have_styles()) GTEST_SKIP();

    create_system();
    BEGIN_HIDE_OUTPUT();
    add_atom(1, 5.0, 5.0, 5.0);
    add_atom(2, 6.0, 5.0, 5.0);
    command("fix bc all bond/create 5 1 2 1.5 1");
    command("run 4 post no");
    END_HIDE_OUTPUT();

    EXPECT_EQ(lmp->atom->nbonds, 0);

    BEGIN_HIDE_OUTPUT();
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-2"));
}

TEST_F(FixBondCreateTest, Probability)
{
    if (!have_styles()) GTEST_SKIP();

    create_system();
    BEGIN_HIDE_OUTPUT();
    add_atom(1, 5.0, 5.0, 5.0);
    add_atom(2, 6.0, 5.0, 5.0);
    command("fix bc all bond/create 1 1 2 1.5 1 prob 0.0 8273");
    command("run 5 post no");
    END_HIDE_OUTPUT();

    EXPECT_EQ(lmp->atom->nbonds, 0);

    BEGIN_HIDE_OUTPUT();
    command("unfix bc");
    command("fix bc all bond/create 1 1 2 1.5 1 prob 1.0 8273");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-2"));
}

// atoms 1 and 2 are in the same molecule, atoms 3 and 4 are not

TEST_F(FixBondCreateTest, MoleculeKeyword)
{
    if (!have_styles()) GTEST_SKIP();

    const std::vector<std::pair<std::string, std::string>> settings = {
        {"off", "1:1-2 1:3-4"}, {"inter", "1:3-4"}, {"intra", "1:1-2"}};

    for (const auto &[keyword, bonds] : settings) {
        BEGIN_HIDE_OUTPUT();
        command("clear");
        END_HIDE_OUTPUT();
        create_system();
        BEGIN_HIDE_OUTPUT();
        add_atom(1, 2.0, 5.0, 5.0);
        add_atom(2, 3.0, 5.0, 5.0);
        add_atom(1, 7.0, 5.0, 5.0);
        add_atom(2, 8.0, 5.0, 5.0);
        command("set atom 1*2 mol 1");
        command("set atom 3 mol 2");
        command("set atom 4 mol 3");
        command("fix bc all bond/create 1 1 2 1.5 1 molecule " + keyword);
        command("run 1 post no");
        END_HIDE_OUTPUT();

        EXPECT_THAT(bond_list(), StrEq(bonds)) << "molecule " << keyword;
    }
}

// the new bond 2-3 connects the existing bonds 1-2, 2-5, and 3-4 and thus
// creates the angles 1-2-3, 3-2-5, and 2-3-4, the dihedrals 1-2-3-4 and
// 5-2-3-4, and one improper with atom 2 as central atom.  with newton_bond
// off, each of them is stored with all of its atoms.  the topology and the
// energies must be the same for both settings.  the reference energies are
// from the same topology defined with create_bonds.

TEST_F(FixBondCreateTest, AnglesDihedralsImpropers)
{
    if (!have_styles() || !info->has_style("dihedral", "harmonic") ||
        !info->has_style("improper", "harmonic"))
        GTEST_SKIP();

    for (const std::string newton_bond : {"on", "off"}) {
        BEGIN_HIDE_OUTPUT();
        command("clear");
        END_HIDE_OUTPUT();
        create_system(newton_bond);
        BEGIN_HIDE_OUTPUT();
        command("dihedral_style harmonic");
        command("dihedral_coeff 1 1.0 1 2");
        command("improper_style harmonic");
        command("improper_coeff 1 1.0 30.0");
        add_atom(3, 2.4, 5.8, 5.0);
        add_atom(1, 3.0, 5.0, 5.0);
        add_atom(2, 4.2, 5.0, 5.0);
        add_atom(3, 4.6, 5.0, 5.8);
        add_atom(3, 3.0, 5.9, 4.4);
        command("create_bonds single/bond 1 1 2");
        command("create_bonds single/bond 1 2 5");
        command("create_bonds single/bond 1 3 4");
        command("fix bc all bond/create 1 1 2 1.5 1 atype 1 dtype 1 itype 1");
        command("run 1 post no");
        command("run 0 post no");
        END_HIDE_OUTPUT();

        EXPECT_THAT(bond_list(), StrEq("1:1-2 1:2-3 1:2-5 1:3-4")) << "newton " << newton_bond;
        EXPECT_THAT(topology_list("angles"), StrEq("1:1-2-3 1:2-3-4 1:3-2-5"))
            << "newton " << newton_bond;
        EXPECT_THAT(topology_list("dihedrals"), StrEq("1:1-2-3-4 1:4-3-2-5"))
            << "newton " << newton_bond;
        EXPECT_THAT(topology_list("impropers"), StrEq("1:2-1-5-3")) << "newton " << newton_bond;
        EXPECT_NEAR(lammps_get_thermo(lmp, "edihed"), 0.123076923076923, 1.0e-12)
            << "newton " << newton_bond;
        EXPECT_NEAR(lammps_get_thermo(lmp, "eimp"), 0.000636261188646685, 1.0e-12)
            << "newton " << newton_bond;
    }
}

TEST_F(FixBondCreateTest, InvalidArguments)
{
    if (!have_styles()) GTEST_SKIP();

    create_system();

    TEST_FAILURE(".*ERROR: Illegal fix bond/create command: missing argument.*",
                 command("fix bc all bond/create 1 1 2 1.5"););
    TEST_FAILURE(".*ERROR: Invalid fix bond/create nevery value 0.*",
                 command("fix bc all bond/create 0 1 2 1.5 1"););
    TEST_FAILURE(".*ERROR: Invalid atom type 0 in fix bond/create command.*",
                 command("fix bc all bond/create 1 0 2 1.5 1"););
    TEST_FAILURE(".*ERROR: Invalid atom type 4 in fix bond/create command.*",
                 command("fix bc all bond/create 1 1 4 1.5 1"););
    TEST_FAILURE(".*ERROR: Invalid fix bond/create cutoff value -1.*",
                 command("fix bc all bond/create 1 1 2 -1.0 1"););
    TEST_FAILURE(".*ERROR: Invalid bond type 3 in fix bond/create command.*",
                 command("fix bc all bond/create 1 1 2 1.5 3"););
    TEST_FAILURE(".*ERROR: Illegal fix bond/create iparam command: missing argument.*",
                 command("fix bc all bond/create 1 1 2 1.5 1 iparam 1"););
    TEST_FAILURE(".*ERROR: Invalid fix bond/create iparam maxbond value -1.*",
                 command("fix bc all bond/create 1 1 2 1.5 1 iparam -1 1"););
    TEST_FAILURE(".*ERROR: Invalid atom type 4 in fix bond/create jparam command.*",
                 command("fix bc all bond/create 1 1 2 1.5 1 jparam 1 4"););
    TEST_FAILURE(".*ERROR: Invalid fix bond/create prob fraction value 1.5.*",
                 command("fix bc all bond/create 1 1 2 1.5 1 prob 1.5 8273"););
    TEST_FAILURE(".*ERROR: Invalid fix bond/create prob seed value 0.*",
                 command("fix bc all bond/create 1 1 2 1.5 1 prob 0.5 0"););
    TEST_FAILURE(".*ERROR: Invalid fix bond/create atype value -1.*",
                 command("fix bc all bond/create 1 1 2 1.5 1 atype -1"););
    TEST_FAILURE(".*ERROR: Unknown option both for fix bond/create molecule.*",
                 command("fix bc all bond/create 1 1 2 1.5 1 molecule both"););
    TEST_FAILURE(".*ERROR: Unknown fix bond/create keyword xyz.*",
                 command("fix bc all bond/create 1 1 2 1.5 1 xyz 1"););
    TEST_FAILURE(".*ERROR: Fix bond/create iparam and jparam values must be the same when both "
                 "atom types are the same.*",
                 command("fix bc all bond/create 1 1 1 1.5 1 iparam 2 1 jparam 1 1"););
    TEST_FAILURE(".*ERROR: Keyword aconstrain requires fix bond/create/angle.*",
                 command("fix bc all bond/create 1 1 2 1.5 1 aconstrain 90 180"););

    TEST_FAILURE(".*ERROR: Illegal fix bond/create/angle aconstrain command: missing argument.*",
                 command("fix bc all bond/create/angle 1 1 2 1.5 1 aconstrain 90"););
    TEST_FAILURE(".*ERROR: Fix bond/create/angle aconstrain amin value 120 must be smaller than "
                 "amax value 90.*",
                 command("fix bc all bond/create/angle 1 1 2 1.5 1 aconstrain 120 90"););
    TEST_FAILURE(".*ERROR: Fix bond/create/angle aconstrain amin value -10 must be between 0 and "
                 "180 degrees.*",
                 command("fix bc all bond/create/angle 1 1 2 1.5 1 aconstrain -10 90"););
    TEST_FAILURE(".*ERROR: Fix bond/create/angle aconstrain amax value 200 must be between 0 and "
                 "180 degrees.*",
                 command("fix bc all bond/create/angle 1 1 2 1.5 1 aconstrain 90 200"););

    // errors detected when a run is set up

    BEGIN_HIDE_OUTPUT();
    command("fix bc all bond/create 1 1 2 2.5 1");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*ERROR: Fix bond/create cutoff 2.5 is longer than the pair style cutoff 2 for "
                 "atom types 1 and 2.*",
                 command("run 0 post no"););

    BEGIN_HIDE_OUTPUT();
    command("unfix bc");
    command("fix bc all bond/create 1 1 2 1.5 1 atype 2");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*ERROR: Fix bond/create atype value 2 is larger than the number of angle "
                 "types 1.*",
                 command("run 0 post no"););

    BEGIN_HIDE_OUTPUT();
    command("unfix bc");
    command("fix bc all bond/create 1 1 2 1.5 1");
    command("pair_style none");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*ERROR: Fix bond/create requires a pair style.*", command("run 0 post no"););
}

TEST_F(FixBondCreateTest, FixBondBreakInvalidArguments)
{
    if (!have_styles() || !info->has_style("fix", "bond/break")) GTEST_SKIP();

    create_system();

    TEST_FAILURE(".*ERROR: Illegal fix bond/break command: missing argument.*",
                 command("fix bb all bond/break 1 1"););
    TEST_FAILURE(".*ERROR: Invalid fix bond/break nevery value 0.*",
                 command("fix bb all bond/break 0 1 1.5"););
    TEST_FAILURE(".*ERROR: Invalid bond type 3 in fix bond/break command.*",
                 command("fix bb all bond/break 1 3 1.5"););
    TEST_FAILURE(".*ERROR: Invalid fix bond/break cutoff value -1.*",
                 command("fix bb all bond/break 1 1 -1.0"););
    TEST_FAILURE(".*ERROR: Illegal fix bond/break prob command: missing argument.*",
                 command("fix bb all bond/break 1 1 1.5 prob 0.5"););
    TEST_FAILURE(".*ERROR: Invalid fix bond/break prob fraction value -0.5.*",
                 command("fix bb all bond/break 1 1 1.5 prob -0.5 8273"););
    TEST_FAILURE(".*ERROR: Invalid fix bond/break prob seed value -5.*",
                 command("fix bb all bond/break 1 1 1.5 prob 0.5 -5"););
    TEST_FAILURE(".*ERROR: Unknown fix bond/break keyword xyz.*",
                 command("fix bb all bond/break 1 1 1.5 xyz 1"););
}

TEST_F(FixBondCreateTest, NonMolecularSystem)
{
    if (!have_styles()) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("atom_style atomic");
    command("region box block 0 10 0 10 0 10");
    command("create_box 2 box");
    END_HIDE_OUTPUT();

    TEST_FAILURE(".*ERROR: Cannot use fix bond/create with non-molecular systems.*",
                 command("fix bc all bond/create 1 1 2 1.5 1"););
    TEST_FAILURE(".*ERROR: Cannot use fix bond/create/angle with non-molecular systems.*",
                 command("fix bc all bond/create/angle 1 1 2 1.5 1 aconstrain 90 180"););
    if (info->has_style("fix", "bond/break")) {
        TEST_FAILURE(".*ERROR: Cannot use fix bond/break with non-molecular systems.*",
                     command("fix bb all bond/break 1 1 1.5"););
    }
}

// the angle constraint is applied at the atom of a prospective pair that
// already has a bond.  the angle 3-1-2 is about 116.6 degrees and the angle
// 6-4-5 is about 66.0 degrees.

TEST_F(FixBondCreateTest, AngleConstraintOneSide)
{
    if (!have_styles()) GTEST_SKIP();

    create_system();
    BEGIN_HIDE_OUTPUT();
    add_atom(1, 2.0, 2.0, 2.0);
    add_atom(2, 2.5, 3.0, 2.0);
    add_atom(3, 1.0, 2.0, 2.0);
    add_atom(1, 7.0, 7.0, 7.0);
    add_atom(2, 6.6, 7.9, 7.0);
    add_atom(3, 6.0, 7.0, 7.0);
    command("create_bonds single/bond 1 1 3");
    command("create_bonds single/bond 1 4 6");
    command("fix bc all bond/create/angle 1 1 2 1.5 1 aconstrain 90 180");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-2 1:1-3 1:4-6"));

    BEGIN_HIDE_OUTPUT();
    command("unfix bc");
    command("fix bc all bond/create/angle 1 1 2 1.5 1 aconstrain 0 90");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-2 1:1-3 1:4-5 1:4-6"));
}

// when both atoms of a prospective pair already have a bond, the angle
// constraint must be met at both of them.  atoms 3-1-2-4 are on a line, while
// the angle 5-6-8 is about 68.2 degrees.

TEST_F(FixBondCreateTest, AngleConstraintBothSides)
{
    if (!have_styles()) GTEST_SKIP();

    create_system();
    BEGIN_HIDE_OUTPUT();
    add_atom(1, 2.0, 2.0, 2.0);
    add_atom(2, 3.0, 2.0, 2.0);
    add_atom(3, 1.0, 2.0, 2.0);
    add_atom(3, 4.0, 2.0, 2.0);
    add_atom(1, 2.0, 7.0, 7.0);
    add_atom(2, 3.0, 7.0, 7.0);
    add_atom(3, 1.0, 7.0, 7.0);
    add_atom(3, 2.6, 8.0, 7.0);
    command("create_bonds single/bond 1 1 3");
    command("create_bonds single/bond 1 2 4");
    command("create_bonds single/bond 1 5 7");
    command("create_bonds single/bond 1 6 8");
    command("fix bc all bond/create/angle 1 1 2 1.5 1 aconstrain 90 180");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-2 1:1-3 1:2-4 1:5-7 1:6-8"));
}

// the existing bonds 1-3 and 4-6 cross a periodic boundary and the angles
// 2-1-3 and 5-4-6 are 180 degrees when the closest images of atoms 3 and 6
// are used.  with any other image, the angles are 0 degrees.

TEST_F(FixBondCreateTest, AngleConstraintPeriodic)
{
    if (!have_styles()) GTEST_SKIP();

    create_system();
    BEGIN_HIDE_OUTPUT();
    add_atom(1, 9.5, 5.0, 5.0);
    add_atom(2, 8.5, 5.0, 5.0);
    add_atom(3, 0.5, 5.0, 5.0);
    add_atom(2, 5.0, 0.4, 5.0);
    add_atom(1, 5.0, 1.4, 5.0);
    add_atom(3, 5.0, 9.6, 5.0);
    command("set atom 3 image 1 0 0");
    command("set atom 6 image 0 -1 0");
    command("create_bonds single/bond 1 1 3");
    command("create_bonds single/bond 1 4 6");
    command("fix bc all bond/create/angle 1 1 2 1.5 1 aconstrain 0 90");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-3 1:4-6"));

    BEGIN_HIDE_OUTPUT();
    command("unfix bc");
    command("fix bc all bond/create/angle 1 1 2 1.5 1 aconstrain 90 180");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-2 1:1-3 1:4-5 1:4-6"));
}

// atoms 3-1-2 are on a line.  for these positions, the computed cosine of the
// angle is slightly below -1 with IEEE 754 double precision math, so it must be
// clamped to get an angle of 180 degrees instead of NaN.  the amax setting
// is not an integer.

TEST_F(FixBondCreateTest, AngleConstraintLinear)
{
    if (!have_styles()) GTEST_SKIP();

    create_system();
    BEGIN_HIDE_OUTPUT();
    add_atom(1, 5.0, 5.0, 5.0);
    add_atom(2, 5.1, 5.1, 5.1);
    add_atom(3, 4.8, 4.8, 4.8);
    command("create_bonds single/bond 1 1 3");
    command("fix bc all bond/create/angle 1 1 2 1.5 1 aconstrain 90 179.5");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-3"));

    BEGIN_HIDE_OUTPUT();
    command("unfix bc");
    command("fix bc all bond/create/angle 1 1 2 1.5 1 aconstrain 90 180");
    command("run 1 post no");
    END_HIDE_OUTPUT();

    EXPECT_THAT(bond_list(), StrEq("1:1-2 1:1-3"));
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") verbose = true;
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
