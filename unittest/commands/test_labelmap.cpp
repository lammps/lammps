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

#include "lammps.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "info.h"
#include "label_map.h"
#include "library.h"
#include "math_const.h"
#include "modify.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using ::testing::ContainsRegex;
using ::testing::ExitedWithCode;
using ::testing::StrEq;

namespace LAMMPS_NS {
class LabelMapTest : public LAMMPSTest {
protected:
    Atom *atom;
    Domain *domain;
    void SetUp() override
    {
        testbinary = "LabelMapTest";
        args       = {"-log", "none", "-echo", "screen", "-nocite", "-v", "num", "1"};
        LAMMPSTest::SetUp();
        atom   = lmp->atom;
        domain = lmp->domain;
    }

    void TearDown() override { LAMMPSTest::TearDown(); }
};

TEST_F(LabelMapTest, Atoms)
{
    EXPECT_EQ(atom->natoms, 0);
    EXPECT_EQ(domain->box_exist, 0);
    EXPECT_EQ(atom->labelmapflag, 0);
    EXPECT_EQ(atom->types_style, Atom::NUMERIC);
    ASSERT_EQ(atom->lmap, nullptr);
    TEST_FAILURE(".*ERROR: Labelmap command before simulation box is.*",
                 command("labelmap atom 3 C1"););

    BEGIN_HIDE_OUTPUT();
    command("region box block 0 2 0 2 0 2");
    command("create_box 4 box");
    END_HIDE_OUTPUT();
    EXPECT_EQ(domain->box_exist, 1);
    EXPECT_EQ(atom->lmap, nullptr);
    EXPECT_EQ(atom->labelmapflag, 0);
    EXPECT_EQ(utils::expand_type(FLERR, "C1", Atom::ATOM, lmp), nullptr);

    BEGIN_HIDE_OUTPUT();
    command("labelmap atom 2 N1");
    command("labelmap atom 3 O1 4 H1");
    command("mass * 1.0");
    command("mass O1 3.0");
    command("mass N1 2.0");
    command("mass H1 4.0");
    END_HIDE_OUTPUT();
    EXPECT_EQ(atom->labelmapflag, 1);
    ASSERT_NE(atom->lmap, nullptr);
    EXPECT_FALSE(atom->lmap->is_complete(Atom::ATOM));
    EXPECT_DOUBLE_EQ(atom->mass[1], 1.0);
    EXPECT_DOUBLE_EQ(atom->mass[2], 2.0);
    EXPECT_DOUBLE_EQ(atom->mass[3], 3.0);
    EXPECT_DOUBLE_EQ(atom->mass[4], 4.0);

    BEGIN_HIDE_OUTPUT();
    command("labelmap atom 1 C1 2 N2 3 ' O#' 1 C1 4 H# 2 N3"); // second '#' starts comment
    command("mass \"O#\" 10.0");
    END_HIDE_OUTPUT();
    EXPECT_TRUE(atom->lmap->is_complete(Atom::ATOM));
    EXPECT_EQ(atom->lmap->find_type("C1", Atom::ATOM), 1);
    EXPECT_EQ(atom->lmap->find_type("N2", Atom::ATOM), 2);
    EXPECT_EQ(atom->lmap->find_type("O#", Atom::ATOM), 3);
    EXPECT_EQ(atom->lmap->find_type("H", Atom::ATOM), 4);
    EXPECT_EQ(atom->lmap->find_type("X", Atom::ATOM), -1);
    EXPECT_EQ(atom->lmap->find_type("", Atom::ATOM), -1);
    EXPECT_THAT(atom->lmap->find_label(1, Atom::ATOM), StrEq("C1"));
    EXPECT_THAT(atom->lmap->find_label(2, Atom::ATOM), StrEq("N2"));
    EXPECT_THAT(atom->lmap->find_label(3, Atom::ATOM), StrEq("O#"));
    EXPECT_THAT(atom->lmap->find_label(4, Atom::ATOM), StrEq("H"));
    EXPECT_THAT(atom->lmap->find_label(-1, Atom::ATOM), StrEq(""));
    EXPECT_THAT(atom->lmap->find_label(5, Atom::ATOM), StrEq(""));
    EXPECT_DOUBLE_EQ(atom->mass[3], 10.0);

    EXPECT_EQ(utils::expand_type(FLERR, "1", Atom::ATOM, lmp), nullptr);
    EXPECT_EQ(utils::expand_type(FLERR, "*3", Atom::ATOM, lmp), nullptr);
    EXPECT_EQ(utils::expand_type(FLERR, "1*2", Atom::ATOM, lmp), nullptr);
    EXPECT_EQ(utils::expand_type(FLERR, "*", Atom::ATOM, lmp), nullptr);
    EXPECT_EQ(utils::expand_type(FLERR, "**", Atom::ATOM, lmp), nullptr);
    EXPECT_EQ(utils::expand_type(FLERR, "1*2*", Atom::ATOM, lmp), nullptr);

    auto *expanded = utils::expand_type(FLERR, "C1", Atom::ATOM, lmp);
    EXPECT_THAT(expanded, StrEq("1"));
    delete[] expanded;
    expanded = utils::expand_type(FLERR, "O#", Atom::ATOM, lmp);
    EXPECT_THAT(expanded, StrEq("3"));
    delete[] expanded;
    TEST_FAILURE(".*ERROR: Atom type string XX not found in labelmap.*",
                 utils::expand_type(FLERR, "XX", Atom::ATOM, lmp););

    TEST_FAILURE(".*ERROR: Labelmap atom type 0 must be within 1-4.*",
                 command("labelmap atom 0 C1"););
    TEST_FAILURE(".*ERROR: Labelmap atom type 5 must be within 1-4.*",
                 command("labelmap atom 5 C1"););
    TEST_FAILURE(".*ERROR: Type label string 1C for atom type 1 is invalid.*",
                 command("labelmap atom 1 1C"););
    TEST_FAILURE(".*ERROR: Type label string #C for atom type 1 is invalid.*",
                 command("labelmap atom 1 '#C'"););
    TEST_FAILURE(".*ERROR: Type label string CA CB for atom type 1 is invalid.*",
                 command("labelmap atom 1 ' CA CB '"););
    TEST_FAILURE(".*ERROR: Type label string \\*C for atom type 1 is invalid.*",
                 command("labelmap atom 1 *C"););
    TEST_FAILURE(".*ERROR: The atom type label N2 is already in use for type 2.*",
                 command("labelmap atom 1 N2"););

    TEST_FAILURE(".*ERROR: No bond types allowed with current box settings.*",
                 command("labelmap bond 1 C1-C1"););
    TEST_FAILURE(".*ERROR: No angle types allowed with current box settings.*",
                 command("labelmap angle 1 C1-C1-C1"););
    TEST_FAILURE(".*ERROR: No dihedral types allowed with current box settings.*",
                 command("labelmap dihedral 1 C1-C1-C1-C1"););
    TEST_FAILURE(".*ERROR: No improper types allowed with current box settings.*",
                 command("labelmap improper 1 C1-C1-C1-C1"););

    TEST_FAILURE(".*ERROR: Incorrect number of arguments for labelmap command.*",
                 command("labelmap atom 1 C1 2"););
    TEST_FAILURE(".*ERROR: Incorrect number of arguments for labelmap command.*",
                 command("labelmap atom 1 C1 atom 2 C2"););
    TEST_FAILURE(".*ERROR: Incorrect number of arguments for labelmap clear command.*",
                 command("labelmap clear atom"););
    TEST_FAILURE(".*ERROR: Incorrect number of arguments for labelmap clear command.*",
                 command("labelmap clear atom bond"););
    TEST_FAILURE(".*ERROR: Incorrect number of arguments for labelmap write command.*",
                 command("labelmap write"););
    TEST_FAILURE(".*ERROR: Incorrect number of arguments for labelmap write command.*",
                 command("labelmap write filename xxx"););
    TEST_FAILURE(".*ERROR: Illegal labelmap atom command: missing argument.*",
                 command("labelmap atom 1"););
    TEST_FAILURE(".*ERROR: Illegal labelmap atom command: missing argument.*",
                 command("labelmap atom"););

    BEGIN_HIDE_OUTPUT();
    command("labelmap clear");
    command("labelmap atom 3 C1 2 N2");
    END_HIDE_OUTPUT();
    EXPECT_FALSE(atom->lmap->is_complete(Atom::ATOM));
    EXPECT_EQ(atom->lmap->find_type("C1", Atom::ATOM), 3);
    EXPECT_EQ(atom->lmap->find_type("N2", Atom::ATOM), 2);

    BEGIN_HIDE_OUTPUT();
    command("labelmap clear");
    command(R"(labelmap atom 1 "C1'" 2 'C2"' 3 """C1'-C2" """ 4 """ C2"-C1'""")");
    END_HIDE_OUTPUT();
    EXPECT_TRUE(atom->lmap->is_complete(Atom::ATOM));
    EXPECT_EQ(atom->lmap->find_type("C1'", Atom::ATOM), 1);
    EXPECT_EQ(atom->lmap->find_type(R"(C2")", Atom::ATOM), 2);
    EXPECT_EQ(atom->lmap->find_type(R"(C1'-C2")", Atom::ATOM), 3);
    EXPECT_EQ(atom->lmap->find_type(R"(C2"-C1')", Atom::ATOM), 4);
    EXPECT_THAT(atom->lmap->find_label(1, Atom::ATOM), StrEq("C1'"));
    EXPECT_THAT(atom->lmap->find_label(2, Atom::ATOM), StrEq(R"(C2")"));
    EXPECT_THAT(atom->lmap->find_label(3, Atom::ATOM), StrEq(R"(C1'-C2")"));
    EXPECT_THAT(atom->lmap->find_label(4, Atom::ATOM), StrEq(R"(C2"-C1')"));
}

TEST_F(LabelMapTest, Topology)
{
    if (!info->has_style("atom", "full")) GTEST_SKIP();

    EXPECT_EQ(atom->natoms, 0);
    EXPECT_EQ(atom->nbonds, 0);
    EXPECT_EQ(atom->nangles, 0);
    EXPECT_EQ(atom->ndihedrals, 0);
    EXPECT_EQ(atom->nimpropers, 0);
    EXPECT_EQ(domain->box_exist, 0);
    EXPECT_EQ(atom->labelmapflag, 0);
    ASSERT_EQ(atom->lmap, nullptr);
    TEST_FAILURE(".*ERROR: Labelmap command before simulation box is.*",
                 command("labelmap atom 3 C1"););

    BEGIN_HIDE_OUTPUT();
    command("atom_style full");
    command("region box block 0 2 0 2 0 2");
    command("create_box 2 box bond/types 3 angle/types 2 dihedral/types 1 improper/types 1");
    command("labelmap atom 1 C1");
    END_HIDE_OUTPUT();
    EXPECT_EQ(atom->labelmapflag, 1);
    ASSERT_NE(atom->lmap, nullptr);
    EXPECT_FALSE(atom->lmap->is_complete(Atom::ATOM));
    EXPECT_FALSE(atom->lmap->is_complete(Atom::BOND));
    EXPECT_FALSE(atom->lmap->is_complete(Atom::ANGLE));
    EXPECT_FALSE(atom->lmap->is_complete(Atom::DIHEDRAL));
    EXPECT_FALSE(atom->lmap->is_complete(Atom::IMPROPER));

    BEGIN_HIDE_OUTPUT();
    command("labelmap atom 2 \"N2'\"");
    command("labelmap bond 1 C1-N2 2 [C1][C1] 3 N2=N2");
    command(R"(labelmap angle 1 C1-N2-C1 2 """ N2'-C1"-N2' """)");
    command("labelmap dihedral 1 'C1-N2-C1-N2'");
    command("labelmap improper 1 \"C1-N2-C1-N2\"");
    command("mass C1 12.0");
    command("mass \"N2'\" 14.0");
    command("labelmap write labelmap_topology.inc");
    END_HIDE_OUTPUT();

    EXPECT_TRUE(atom->lmap->is_complete(Atom::ATOM));
    EXPECT_TRUE(atom->lmap->is_complete(Atom::BOND));
    EXPECT_TRUE(atom->lmap->is_complete(Atom::ANGLE));
    EXPECT_TRUE(atom->lmap->is_complete(Atom::DIHEDRAL));
    EXPECT_TRUE(atom->lmap->is_complete(Atom::IMPROPER));
    EXPECT_EQ(atom->lmap->find_type("C1", Atom::ATOM), 1);
    EXPECT_EQ(atom->lmap->find_type("N2'", Atom::ATOM), 2);
    EXPECT_EQ(atom->lmap->find_type("C1-N2", Atom::BOND), 1);
    EXPECT_EQ(atom->lmap->find_type("[C1][C1]", Atom::BOND), 2);
    EXPECT_EQ(atom->lmap->find_type("N2=N2", Atom::BOND), 3);
    EXPECT_EQ(atom->lmap->find_type("C1-N2-C1", Atom::ANGLE), 1);
    EXPECT_EQ(atom->lmap->find_type("N2'-C1\"-N2'", Atom::ANGLE), 2);
    EXPECT_EQ(atom->lmap->find_type("C1-N2-C1-N2", Atom::DIHEDRAL), 1);
    EXPECT_EQ(atom->lmap->find_type("C1-N2-C1-N2", Atom::IMPROPER), 1);
    EXPECT_EQ(atom->lmap->find_type("X", Atom::ATOM), -1);
    EXPECT_EQ(atom->lmap->find_type("N2'-C1\"-N2'", Atom::BOND), -1);

    EXPECT_THAT(atom->lmap->find_label(1, Atom::BOND), StrEq("C1-N2"));
    EXPECT_THAT(atom->lmap->find_label(2, Atom::BOND), StrEq("[C1][C1]"));
    EXPECT_THAT(atom->lmap->find_label(3, Atom::BOND), StrEq("N2=N2"));
    EXPECT_THAT(atom->lmap->find_label(1, Atom::ANGLE), StrEq("C1-N2-C1"));
    EXPECT_THAT(atom->lmap->find_label(2, Atom::ANGLE), StrEq(R"(N2'-C1"-N2')"));
    EXPECT_THAT(atom->lmap->find_label(1, Atom::DIHEDRAL), StrEq("C1-N2-C1-N2"));
    EXPECT_THAT(atom->lmap->find_label(1, Atom::IMPROPER), StrEq("C1-N2-C1-N2"));
    EXPECT_THAT(atom->lmap->find_label(0, Atom::BOND), StrEq(""));
    EXPECT_THAT(atom->lmap->find_label(4, Atom::BOND), StrEq(""));
    EXPECT_THAT(atom->lmap->find_label(-1, Atom::ANGLE), StrEq(""));
    EXPECT_THAT(atom->lmap->find_label(3, Atom::ANGLE), StrEq(""));
    EXPECT_THAT(atom->lmap->find_label(0, Atom::DIHEDRAL), StrEq(""));
    EXPECT_THAT(atom->lmap->find_label(2, Atom::DIHEDRAL), StrEq(""));
    EXPECT_THAT(atom->lmap->find_label(-1, Atom::IMPROPER), StrEq(""));
    EXPECT_THAT(atom->lmap->find_label(0, Atom::IMPROPER), StrEq(""));
    EXPECT_DOUBLE_EQ(atom->mass[1], 12.0);
    EXPECT_DOUBLE_EQ(atom->mass[2], 14.0);

    BEGIN_HIDE_OUTPUT();
    command("labelmap clear");
    command("labelmap atom 1 C1");
    END_HIDE_OUTPUT();
    EXPECT_EQ(atom->labelmapflag, 1);
    ASSERT_NE(atom->lmap, nullptr);
    EXPECT_FALSE(atom->lmap->is_complete(Atom::ATOM));
    EXPECT_FALSE(atom->lmap->is_complete(Atom::BOND));
    EXPECT_FALSE(atom->lmap->is_complete(Atom::ANGLE));
    EXPECT_FALSE(atom->lmap->is_complete(Atom::DIHEDRAL));
    EXPECT_FALSE(atom->lmap->is_complete(Atom::IMPROPER));

    BEGIN_HIDE_OUTPUT();
    command("include labelmap_topology.inc");
    END_HIDE_OUTPUT();

    EXPECT_TRUE(atom->lmap->is_complete(Atom::ATOM));
    EXPECT_TRUE(atom->lmap->is_complete(Atom::BOND));
    EXPECT_TRUE(atom->lmap->is_complete(Atom::ANGLE));
    EXPECT_TRUE(atom->lmap->is_complete(Atom::DIHEDRAL));
    EXPECT_TRUE(atom->lmap->is_complete(Atom::IMPROPER));
    EXPECT_EQ(atom->lmap->find_type("C1", Atom::ATOM), 1);
    EXPECT_EQ(atom->lmap->find_type("N2'", Atom::ATOM), 2);
    EXPECT_EQ(atom->lmap->find_type("C1-N2", Atom::BOND), 1);
    EXPECT_EQ(atom->lmap->find_type("[C1][C1]", Atom::BOND), 2);
    EXPECT_EQ(atom->lmap->find_type("N2=N2", Atom::BOND), 3);
    EXPECT_EQ(atom->lmap->find_type("C1-N2-C1", Atom::ANGLE), 1);
    EXPECT_EQ(atom->lmap->find_type("N2'-C1\"-N2'", Atom::ANGLE), 2);
    EXPECT_EQ(atom->lmap->find_type("C1-N2-C1-N2", Atom::DIHEDRAL), 1);
    EXPECT_EQ(atom->lmap->find_type("C1-N2-C1-N2", Atom::IMPROPER), 1);
    EXPECT_EQ(atom->lmap->find_type("X", Atom::ATOM), -1);
    EXPECT_EQ(atom->lmap->find_type("N2'-C1\"-N2'", Atom::BOND), -1);
    platform::unlink("labelmap_topology.inc");

    auto *expanded = utils::expand_type(FLERR, "N2'", Atom::ATOM, lmp);
    EXPECT_THAT(expanded, StrEq("2"));
    delete[] expanded;
    expanded = utils::expand_type(FLERR, "[C1][C1]", Atom::BOND, lmp);
    EXPECT_THAT(expanded, StrEq("2"));
    delete[] expanded;
    expanded = utils::expand_type(FLERR, "C1-N2-C1", Atom::ANGLE, lmp);
    EXPECT_THAT(expanded, StrEq("1"));
    delete[] expanded;
    expanded = utils::expand_type(FLERR, "C1-N2-C1-N2", Atom::DIHEDRAL, lmp);
    EXPECT_THAT(expanded, StrEq("1"));
    delete[] expanded;
    expanded = utils::expand_type(FLERR, "C1-N2-C1-N2", Atom::IMPROPER, lmp);
    EXPECT_THAT(expanded, StrEq("1"));
    delete[] expanded;
    TEST_FAILURE(".*ERROR: Bond type string XX not found in labelmap.*",
                 utils::expand_type(FLERR, "XX", Atom::BOND, lmp););
    TEST_FAILURE(".*ERROR: Angle type string XX not found in labelmap.*",
                 utils::expand_type(FLERR, "XX", Atom::ANGLE, lmp););
    TEST_FAILURE(".*ERROR: Dihedral type string XX not found in labelmap.*",
                 utils::expand_type(FLERR, "XX", Atom::DIHEDRAL, lmp););
    TEST_FAILURE(".*ERROR: Improper type string XX not found in labelmap.*",
                 utils::expand_type(FLERR, "XX", Atom::IMPROPER, lmp););

    // check for inference

    BEGIN_HIDE_OUTPUT();
    command(R"(labelmap atom 2 "N2")");
    command(R"(labelmap bond 2 "C1-C1")");
    command(R"(labelmap bond 3 "N2-N2")");
    command(R"(labelmap angle 2 "N2-C1-N2")");
    END_HIDE_OUTPUT();

    EXPECT_EQ(atom->lmap->infer_bondtype(1, 1), 2);
    EXPECT_EQ(atom->lmap->infer_bondtype(1, 2), 1);
    EXPECT_EQ(atom->lmap->infer_bondtype(2, 2), 3);
    EXPECT_EQ(atom->lmap->infer_bondtype(2, 3), 0);
    EXPECT_EQ(atom->lmap->infer_bondtype({"N2", "N2"}), 3);
    EXPECT_EQ(atom->lmap->infer_bondtype({"C1", "N2"}), 1);
    EXPECT_EQ(atom->lmap->infer_bondtype({"C1", "C1"}), 2);
    EXPECT_EQ(atom->lmap->infer_bondtype({"N2", "C1"}), -1);
    EXPECT_EQ(atom->lmap->infer_bondtype({"N3", "C1"}), 0);
    EXPECT_EQ(atom->lmap->infer_angletype(1, 2, 1), 1);
    EXPECT_EQ(atom->lmap->infer_angletype(2, 1, 2), 2);
    EXPECT_EQ(atom->lmap->infer_angletype(1, 1, 1), 0);
    EXPECT_EQ(atom->lmap->infer_angletype(3, 1, 1), 0);
    EXPECT_EQ(atom->lmap->infer_angletype(1, 3, 1), 0);
    EXPECT_EQ(atom->lmap->infer_angletype(1, 1, 3), 0);
    EXPECT_EQ(atom->lmap->infer_angletype({"C1", "N2", "C1"}), 1);
    EXPECT_EQ(atom->lmap->infer_angletype({"N2", "C1", "N2"}), 2);
    EXPECT_EQ(atom->lmap->infer_angletype({"C1", "N2", "N2"}), 0);
    EXPECT_EQ(atom->lmap->infer_angletype({"C1", "N3", "C1"}), 0);
    EXPECT_EQ(atom->lmap->infer_dihedraltype(1, 2, 1, 2), 1);
    EXPECT_EQ(atom->lmap->infer_dihedraltype(2, 1, 2, 2), 0);
    EXPECT_EQ(atom->lmap->infer_dihedraltype(2, 1, 2, 1), -1);
    EXPECT_EQ(atom->lmap->infer_dihedraltype(3, 1, 2, 1), 0);
    EXPECT_EQ(atom->lmap->infer_dihedraltype(2, 3, 2, 1), 0);
    EXPECT_EQ(atom->lmap->infer_dihedraltype(2, 1, 3, 1), 0);
    EXPECT_EQ(atom->lmap->infer_dihedraltype(2, 1, 2, 3), 0);
    EXPECT_EQ(atom->lmap->infer_dihedraltype({"C1", "N2", "C1", "N2"}), 1);
    EXPECT_EQ(atom->lmap->infer_dihedraltype({"N2", "C1", "N2", "C1"}), -1);
    EXPECT_EQ(atom->lmap->infer_dihedraltype({"C1", "N2", "N2", "C1"}), 0);
    EXPECT_EQ(atom->lmap->infer_dihedraltype({"N3", "C1", "N2", "C1"}), 0);
    EXPECT_EQ(atom->lmap->infer_impropertype(1, 2, 1, 2), 1);
    EXPECT_EQ(atom->lmap->infer_impropertype(2, 1, 2, 2), 0);
    EXPECT_EQ(atom->lmap->infer_impropertype(2, 1, 2, 1), -1);
    EXPECT_EQ(atom->lmap->infer_impropertype(3, 1, 2, 1), 0);
    EXPECT_EQ(atom->lmap->infer_impropertype(2, 3, 2, 1), 0);
    EXPECT_EQ(atom->lmap->infer_impropertype(2, 1, 3, 1), 0);
    EXPECT_EQ(atom->lmap->infer_impropertype(2, 1, 2, 3), 0);
    EXPECT_EQ(atom->lmap->infer_impropertype({"C1", "N2", "C1", "N2"}), 1);
    EXPECT_EQ(atom->lmap->infer_impropertype({"N2", "C1", "N2", "C1"}), -1);
    EXPECT_EQ(atom->lmap->infer_impropertype({"C1", "N2", "N2", "C1"}), -1);
    EXPECT_EQ(atom->lmap->infer_impropertype({"N3", "C1", "N2", "C1"}), 0);
}

TEST_F(LabelMapTest, CheckLabelsWriteData)
{
    if (!info->has_style("atom", "full")) GTEST_SKIP();
    if (!info->has_style("bond", "harmonic")) GTEST_SKIP();
    if (!info->has_style("angle", "harmonic")) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("atom_style full");
    command("region box block -5 5 -5 5 -5 5");
    command("create_box 2 box bond/types 1 angle/types 1 extra/bond/per/atom 2 "
            "extra/angle/per/atom 2 extra/special/per/atom 4");
    command("labelmap atom 1 C 2 N");
    command("labelmap bond 1 C-N");
    command("labelmap angle 1 C-N-C");
    command("mass * 1.0");
    command("pair_style lj/cut 2.5");
    command("pair_coeff * * 0.0 1.0");
    command("bond_style harmonic");
    command("bond_coeff * 0.0 1.0");
    command("angle_style harmonic");
    command("angle_coeff * 0.0 180.0");
    command("create_atoms C single 0.0 0.0 0.0");
    command("create_atoms N single 1.0 0.0 0.0");
    command("create_atoms C single 2.0 0.0 0.0");
    command("create_atoms N single 3.0 0.0 0.0");
    command("create_bonds single/bond 1 1 2");
    command("create_bonds single/bond 1 2 3");
    command("create_bonds single/bond 1 3 4");
    command("create_bonds single/angle 1 1 2 3");
    command("velocity all set 1.0 0.0 0.0");
    command("fix 1 all nve");
    command("labelmap check_labels a");
    END_HIDE_OUTPUT();

    BEGIN_CAPTURE_OUTPUT();
    command("run 10 post no");
    auto text = END_CAPTURE_OUTPUT();
    EXPECT_THAT(text, ContainsRegex(".*All angles in the simulation have self-consistent type "
                                    "labels.*"));

    // atoms 2-3-4 are N-C-N, which does not match the C-N-C label of angle type 1

    BEGIN_HIDE_OUTPUT();
    command("create_bonds single/angle 1 2 3 4");
    command("labelmap check_labels a");
    END_HIDE_OUTPUT();

    BEGIN_CAPTURE_OUTPUT();
    command("run 0 post no");
    text = END_CAPTURE_OUTPUT();
    EXPECT_THAT(text, ContainsRegex(".*WARNING: Angle between atoms 2, 3, 4 has constituent atom "
                                    "types \\(N, C, N\\) that do not match its type label "
                                    "\\(C-N-C\\).*"));

    BEGIN_HIDE_OUTPUT();
    command("write_data labelmap_write.data nocoeff");
    command("clear");
    command("atom_style full");
    command("read_data labelmap_write.data");
    END_HIDE_OUTPUT();
    platform::unlink("labelmap_write.data");

    atom = lmp->atom;
    ASSERT_NE(atom->lmap, nullptr);
    EXPECT_EQ(atom->natoms, 4);
    EXPECT_EQ(atom->nbonds, 3);
    EXPECT_EQ(atom->nangles, 2);
    EXPECT_THAT(atom->lmap->find_label(1, Atom::ATOM), StrEq("C"));
    EXPECT_THAT(atom->lmap->find_label(2, Atom::ATOM), StrEq("N"));
    EXPECT_THAT(atom->lmap->find_label(1, Atom::BOND), StrEq("C-N"));
    EXPECT_THAT(atom->lmap->find_label(1, Atom::ANGLE), StrEq("C-N-C"));

    // 10 steps of 0.005 at unit velocity along x without any forces

    for (int i = 1; i <= 4; ++i) {
        int j = atom->map(i);
        ASSERT_GE(j, 0);
        EXPECT_EQ(atom->type[j], (i % 2) ? 1 : 2);
        EXPECT_NEAR(atom->x[j][0], i - 1.0 + 0.05, 1.0e-10);
        EXPECT_NEAR(atom->v[j][0], 1.0, 1.0e-10);
    }
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
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();

    // finalize the KOKKOS package explicitly: otherwise Kokkos is torn down by
    // static destructors at program exit, leading to segfaults in some cases
    // same workaround as the force-style and FFT3d test drivers

    lammps_kokkos_finalize();

    MPI_Finalize();
    return rv;
}
