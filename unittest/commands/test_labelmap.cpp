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
    EXPECT_EQ(atom->lmap->find("C1", Atom::ATOM), 1);
    EXPECT_EQ(atom->lmap->find("N2", Atom::ATOM), 2);
    EXPECT_EQ(atom->lmap->find("O#", Atom::ATOM), 3);
    EXPECT_EQ(atom->lmap->find("H", Atom::ATOM), 4);
    EXPECT_EQ(atom->lmap->find("X", Atom::ATOM), -1);
    EXPECT_DOUBLE_EQ(atom->mass[3], 10.0);

    EXPECT_EQ(utils::expand_type(FLERR, "1", Atom::ATOM, lmp), nullptr);
    EXPECT_EQ(utils::expand_type(FLERR, "*3", Atom::ATOM, lmp), nullptr);
    EXPECT_EQ(utils::expand_type(FLERR, "1*2", Atom::ATOM, lmp), nullptr);
    EXPECT_EQ(utils::expand_type(FLERR, "*", Atom::ATOM, lmp), nullptr);
    EXPECT_EQ(utils::expand_type(FLERR, "**", Atom::ATOM, lmp), nullptr);
    EXPECT_EQ(utils::expand_type(FLERR, "1*2*", Atom::ATOM, lmp), nullptr);

    auto expanded = utils::expand_type(FLERR, "C1", Atom::ATOM, lmp);
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
    EXPECT_EQ(atom->lmap->find("C1", Atom::ATOM), 3);
    EXPECT_EQ(atom->lmap->find("N2", Atom::ATOM), 2);

    BEGIN_HIDE_OUTPUT();
    command("labelmap clear");
    command("labelmap atom 1 \"C1'\" 2 'C2\"' 3 \"\"\"C1'-C2\" \"\"\" 4 \"\"\" C2\"-C1'\"\"\"");
    END_HIDE_OUTPUT();
    EXPECT_TRUE(atom->lmap->is_complete(Atom::ATOM));
    EXPECT_EQ(atom->lmap->find("C1'", Atom::ATOM), 1);
    EXPECT_EQ(atom->lmap->find("C2\"", Atom::ATOM), 2);
    EXPECT_EQ(atom->lmap->find("C1'-C2\"", Atom::ATOM), 3);
    EXPECT_EQ(atom->lmap->find("C2\"-C1'", Atom::ATOM), 4);
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
    command("labelmap angle 1 C1-N2-C1 2 \"\"\" N2'-C1\"-N2' \"\"\"");
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
    EXPECT_EQ(atom->lmap->find("C1", Atom::ATOM), 1);
    EXPECT_EQ(atom->lmap->find("N2'", Atom::ATOM), 2);
    EXPECT_EQ(atom->lmap->find("C1-N2", Atom::BOND), 1);
    EXPECT_EQ(atom->lmap->find("[C1][C1]", Atom::BOND), 2);
    EXPECT_EQ(atom->lmap->find("N2=N2", Atom::BOND), 3);
    EXPECT_EQ(atom->lmap->find("C1-N2-C1", Atom::ANGLE), 1);
    EXPECT_EQ(atom->lmap->find("N2'-C1\"-N2'", Atom::ANGLE), 2);
    EXPECT_EQ(atom->lmap->find("C1-N2-C1-N2", Atom::DIHEDRAL), 1);
    EXPECT_EQ(atom->lmap->find("C1-N2-C1-N2", Atom::IMPROPER), 1);
    EXPECT_EQ(atom->lmap->find("X", Atom::ATOM), -1);
    EXPECT_EQ(atom->lmap->find("N2'-C1\"-N2'", Atom::BOND), -1);
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
    EXPECT_EQ(atom->lmap->find("C1", Atom::ATOM), 1);
    EXPECT_EQ(atom->lmap->find("N2'", Atom::ATOM), 2);
    EXPECT_EQ(atom->lmap->find("C1-N2", Atom::BOND), 1);
    EXPECT_EQ(atom->lmap->find("[C1][C1]", Atom::BOND), 2);
    EXPECT_EQ(atom->lmap->find("N2=N2", Atom::BOND), 3);
    EXPECT_EQ(atom->lmap->find("C1-N2-C1", Atom::ANGLE), 1);
    EXPECT_EQ(atom->lmap->find("N2'-C1\"-N2'", Atom::ANGLE), 2);
    EXPECT_EQ(atom->lmap->find("C1-N2-C1-N2", Atom::DIHEDRAL), 1);
    EXPECT_EQ(atom->lmap->find("C1-N2-C1-N2", Atom::IMPROPER), 1);
    EXPECT_EQ(atom->lmap->find("X", Atom::ATOM), -1);
    EXPECT_EQ(atom->lmap->find("N2'-C1\"-N2'", Atom::BOND), -1);
    platform::unlink("labelmap_topology.inc");

    auto expanded = utils::expand_type(FLERR, "N2'", Atom::ATOM, lmp);
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
