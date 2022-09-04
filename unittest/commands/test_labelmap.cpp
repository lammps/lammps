/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
class SetTest : public LAMMPSTest {
protected:
    Atom *atom;
    Domain *domain;
    void SetUp() override
    {
        testbinary = "SetTest";
        args       = {"-log", "none", "-echo", "screen", "-nocite", "-v", "num", "1"};
        LAMMPSTest::SetUp();
        atom   = lmp->atom;
        domain = lmp->domain;
    }

    void TearDown() override { LAMMPSTest::TearDown(); }

    void atomic_system(const std::string &atom_style, const std::string units = "real")
    {
        BEGIN_HIDE_OUTPUT();
        command("atom_style " + atom_style);
        command("atom_modify map array");
        command("units " + units);
        command("lattice sc 1.0 origin 0.125 0.125 0.125");
        command("region box block 0 2 0 2 0 2");
        command("create_box 8 box");
        command("create_atoms 1 box");
        command("mass * 1.0");
        command("region left block 0.0 1.0 INF INF INF INF");
        command("region right block 1.0 2.0 INF INF INF INF");
        command("region top block INF INF 0.0 1.0 INF INF");
        command("region bottom block INF INF 1.0 2.0 INF INF");
        command("region front block INF INF INF INF 0.0 1.0");
        command("region back block INF INF INF 1.0 2.0 INF");
        command("group top region top");
        command("group bottom region bottom");
        END_HIDE_OUTPUT();
    }

    bool molecular_system()
    {
        if (info->has_style("atom", "full")) {
            BEGIN_HIDE_OUTPUT();
            command("variable input_dir index \"" STRINGIFY(TEST_INPUT_FOLDER) "\"");
            command("include \"${input_dir}/in.fourmol\"");
            command("group allwater molecule 3:6");
            command("region half block 0.0 INF INF INF INF INF");
            END_HIDE_OUTPUT();
            return true;
        } else
            return false;
    }
};

TEST_F(SetTest, NoBoxAtoms)
{
    ASSERT_EQ(atom->natoms, 0);
    ASSERT_EQ(domain->box_exist, 0);
    ASSERT_EQ(atom->lmap, nullptr);
    TEST_FAILURE(".*ERROR: Labelmap command before simulation box is.*",
                 command("labelmap atom 3 C1"););

    BEGIN_HIDE_OUTPUT();
    command("region box block 0 2 0 2 0 2");
    command("create_box 4 box");
    command("labelmap atom 2 N1");
    command("labelmap atom 3 O1 4 H1");
    END_HIDE_OUTPUT();
    ASSERT_NE(atom->lmap, nullptr);
    ASSERT_FALSE(atom->lmap->is_complete(Atom::ATOM));

    BEGIN_HIDE_OUTPUT();
    command("labelmap atom 1 C1 2 N2 3 'O#' 1 C1 4 H# 2 N3"); // second '#' starts comment
    END_HIDE_OUTPUT();
    ASSERT_TRUE(atom->lmap->is_complete(Atom::ATOM));
    ASSERT_EQ(atom->lmap->find("C1", Atom::ATOM), 1);
    ASSERT_EQ(atom->lmap->find("N2", Atom::ATOM), 2);
    ASSERT_EQ(atom->lmap->find("O#", Atom::ATOM), 3);
    ASSERT_EQ(atom->lmap->find("H", Atom::ATOM), 4);

    TEST_FAILURE(".*ERROR: Labelmap atom type 0 must be within 1-4.*",
                 command("labelmap atom 0 C1"););
    TEST_FAILURE(".*ERROR: Labelmap atom type 5 must be within 1-4.*",
                 command("labelmap atom 5 C1"););
    TEST_FAILURE(".*ERROR: Label 1C for atom type 1 must not start with a number.*",
                 command("labelmap atom 1 1C"););
    TEST_FAILURE(".*ERROR: Label #C for atom type 1 must not start with a number, a '#', or.*",
                 command("labelmap atom 1 '#C'"););
    TEST_FAILURE(".*ERROR: Label \\*C for atom type 1 must not start with a number, a '#', or.*",
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
    ASSERT_FALSE(atom->lmap->is_complete(Atom::ATOM));
    ASSERT_EQ(atom->lmap->find("C1", Atom::ATOM), 3);
    ASSERT_EQ(atom->lmap->find("N2", Atom::ATOM), 2);

    BEGIN_HIDE_OUTPUT();
    command("labelmap clear");
    command("labelmap atom 1 \"C1'\" 2 'C2\"' 3 \"\"\"C1'-C2\" \"\"\" 4 \"\"\" C2\"-C1'\"\"\"");
    END_HIDE_OUTPUT();
    ASSERT_TRUE(atom->lmap->is_complete(Atom::ATOM));
    ASSERT_EQ(atom->lmap->find("C1'", Atom::ATOM), 1);
    ASSERT_EQ(atom->lmap->find("C2\"", Atom::ATOM), 2);
    ASSERT_EQ(atom->lmap->find("C1'-C2\"", Atom::ATOM), 3);
    ASSERT_EQ(atom->lmap->find("C2\"-C1'", Atom::ATOM), 4);
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
