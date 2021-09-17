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

#include "../testing/core.h"
#include "atom.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "molecule.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <mpi.h>
#include <string>

using namespace LAMMPS_NS;

using testing::MatchesRegex;
using testing::StrEq;

using utils::split_words;

#define test_name test_info_->name()

static void create_molecule_files(const std::string &h2o_filename, const std::string &co2_filename)
{
    // create molecule files
    const char h2o_file[] = "# Water molecule. SPC/E model.\n\n3 atoms\n2 bonds\n1 angles\n\n"
                            "Coords\n\n1 1.12456 0.09298 1.27452\n"
                            "2 1.53683 0.75606 1.89928\n3 0.49482 0.56390 0.65678\n\n"
                            "Types\n\n1 1\n2 2\n3 2\n\n"
                            "Charges\n\n1 -0.8472\n2 0.4236\n3 0.4236\n\n"
                            "Bonds\n\n1 1 1 2\n2 1 1 3\n\n"
                            "Angles\n\n1 1 2 1 3\n\n"
                            "Shake Flags\n\n1 1\n2 1\n3 1\n\n"
                            "Shake Atoms\n\n1 1 2 3\n2 1 2 3\n3 1 2 3\n\n"
                            "Shake Bond Types\n\n1 1 1 1\n2 1 1 1\n3 1 1 1\n\n"
                            "Special Bond Counts\n\n1 2 0 0\n2 1 1 0\n3 1 1 0\n\n"
                            "Special Bonds\n\n1 2 3\n2 1 3\n3 1 2\n\n";
    const char co2_file[] = "# CO2 molecule file. TraPPE model.\n\n"
                            "3 atoms\n2 bonds\n1 angles\n\n"
                            "Coords\n\n1 0.0 0.0 0.0\n2 -1.16 0.0 0.0\n3 1.16 0.0 0.0\n\n"
                            "Types\n\n1 1\n2 2\n3 2\n\n"
                            "Charges\n\n1 0.7\n2 -0.35\n3 -0.35\n\n"
                            "Bonds\n\n1 1 1 2\n2 1 1 3\n\n"
                            "Angles\n\n1 1 2 1 3\n\n"
                            "Special Bond Counts\n\n1 2 0 0\n2 1 1 0\n3 1 1 0\n\n"
                            "Special Bonds\n\n1 2 3\n2 1 3\n3 1 2\n\n";

    FILE *fp = fopen(h2o_filename.c_str(), "w");
    if (fp) {
        fputs(h2o_file, fp);
        fclose(fp);
    }
    fp = fopen(co2_filename.c_str(), "w");
    if (fp) {
        fputs(co2_file, fp);
        fclose(fp);
    }
}

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

class MoleculeFileTest : public LAMMPSTest {
protected:
    static void SetUpTestSuite() { create_molecule_files("moltest.h2o.mol", "moltest.co2.mol"); }

    static void TearDownTestSuite()
    {
        remove("moltest.h2o.mol");
        remove("moltest.co2.mol");
    }

    void SetUp() override
    {
        testbinary = "MoleculeFileTest";
        LAMMPSTest::SetUp();
        ASSERT_NE(lmp, nullptr);
    }

    void TearDown() override { LAMMPSTest::TearDown(); }

    void run_mol_cmd(const std::string &name, const std::string &args, const std::string &content)
    {
        std::string file = fmt::format("moltest_{}.mol", name);
        FILE *fp         = fopen(file.c_str(), "w");
        fputs(content.c_str(), fp);
        fclose(fp);

        command(fmt::format("molecule {} {} {}", name, file, args));
        remove(file.c_str());
    }
};

TEST_F(MoleculeFileTest, nofile)
{
    TEST_FAILURE(".*Cannot open molecule file nofile.mol.*", command("molecule 1 nofile.mol"););
}

TEST_F(MoleculeFileTest, badid)
{
    TEST_FAILURE(".*Molecule template ID must have only "
                 "alphanumeric or underscore characters.*",
                 command("molecule @mol nofile.mol"););
}

TEST_F(MoleculeFileTest, badargs)
{
    TEST_FAILURE(".*Illegal molecule command.*",
                 run_mol_cmd(test_name, "offset 1 2 3 4",
                             "Comment\n1 atoms\n\n Coords\n\n 1 0.0 0.0 0.0\n"););
    TEST_FAILURE(
        ".*Illegal molecule command.*",
        run_mol_cmd(test_name, "toff", "Comment\n1 atoms\n\n Coords\n\n 1 0.0 0.0 0.0\n"););
    TEST_FAILURE(
        ".*Illegal molecule command.*",
        run_mol_cmd(test_name, "boff", "Comment\n1 atoms\n\n Coords\n\n 1 0.0 0.0 0.0\n"););
    TEST_FAILURE(
        ".*Illegal molecule command.*",
        run_mol_cmd(test_name, "aoff", "Comment\n1 atoms\n\n Coords\n\n 1 0.0 0.0 0.0\n"););
    TEST_FAILURE(
        ".*Illegal molecule command.*",
        run_mol_cmd(test_name, "doff", "Comment\n1 atoms\n\n Coords\n\n 1 0.0 0.0 0.0\n"););
    TEST_FAILURE(
        ".*Illegal molecule command.*",
        run_mol_cmd(test_name, "ioff", "Comment\n1 atoms\n\n Coords\n\n 1 0.0 0.0 0.0\n"););
    TEST_FAILURE(
        ".*Illegal molecule command.*",
        run_mol_cmd(test_name, "scale", "Comment\n1 atoms\n\n Coords\n\n 1 0.0 0.0 0.0\n"););
    remove("badargs.mol");
}

TEST_F(MoleculeFileTest, noatom)
{
    TEST_FAILURE(".*ERROR: No atoms or invalid atom count in molecule file.*",
                 run_mol_cmd(test_name, "",
                             "Comment\n0 atoms\n1 bonds\n\n"
                             " Coords\n\nBonds\n\n 1 1 2\n"););
    remove("noatom.mol");
}

TEST_F(MoleculeFileTest, empty)
{
    TEST_FAILURE(".*ERROR: Unexpected end of molecule file.*",
                 run_mol_cmd(test_name, "", "Comment\n\n"););
    remove("empty.mol");
}

TEST_F(MoleculeFileTest, nospecial)
{
    TEST_FAILURE(".*ERROR: Cannot auto-generate special bonds before simulation box is defined.*",
                 run_mol_cmd(test_name, "",
                             "Comment\n3 atoms\n\n2 bonds\n\n"
                             " Coords\n\n 1 1.0 1.0 1.0\n 2 1.0 1.0 0.0\n 3 1.0 0.0 1.0\n"
                             " Bonds\n\n 1 1 1 2\n 2 1 1 3\n"););
    remove("nospecial.mol");
}

TEST_F(MoleculeFileTest, minimal)
{
    BEGIN_CAPTURE_OUTPUT();
    run_mol_cmd(test_name, "", "Comment\n1 atoms\n\n Coords\n\n 1 0.0 0.0 0.0\n");
    auto output = END_CAPTURE_OUTPUT();
    ASSERT_THAT(output, MatchesRegex(".*Read molecule template.*1 molecules.*1 atoms.*0 bonds.*"));
}

TEST_F(MoleculeFileTest, twomols)
{
    BEGIN_CAPTURE_OUTPUT();
    run_mol_cmd(test_name, "",
                "Comment\n2 atoms\n\n"
                " Coords\n\n 1 0.0 0.0 0.0\n 2 0.0 0.0 1.0\n"
                " Molecules\n\n 1 1\n 2 2\n\n Types\n\n 1 1\n 2 2\n\n");
    auto output = END_CAPTURE_OUTPUT();
    ASSERT_THAT(output, MatchesRegex(".*Read molecule template.*2 molecules.*2 atoms "
                                     "with max type 2.*0 bonds.*"));
}

TEST_F(MoleculeFileTest, twofiles)
{
    BEGIN_CAPTURE_OUTPUT();
    command("molecule twomols moltest.h2o.mol moltest.co2.mol offset 2 1 1 0 0");
    auto output = END_CAPTURE_OUTPUT();
    ASSERT_THAT(output, MatchesRegex(".*Read molecule template twomols:.*1 molecules.*3 atoms "
                                     "with max type 2.*2 bonds with max type 1.*"
                                     "1 angles with max type 1.*0 dihedrals.*"
                                     ".*Read molecule template twomols:.*1 molecules.*3 atoms "
                                     "with max type 4.*2 bonds with max type 2.*"
                                     "1 angles with max type 2.*0 dihedrals.*"));
}

TEST_F(MoleculeFileTest, bonds)
{
    if (!LAMMPS::is_installed_pkg("MOLECULE")) GTEST_SKIP();
    BEGIN_CAPTURE_OUTPUT();
    command("atom_style bond");
    command("region box block 0 1 0 1 0 1");
    command("create_box 2 box bond/types 2 extra/bond/per/atom 2 "
            "extra/special/per/atom 4");
    run_mol_cmd(test_name, "",
                "Comment\n"
                "4 atoms\n"
                "2 bonds\n\n"
                " Coords\n\n"
                " 1 1.0 1.0 1.0\n"
                " 2 1.0 1.0 0.0\n"
                " 3 1.0 0.0 1.0\n"
                " 4 1.0 0.0 0.0\n"
                " Types\n\n"
                " 1 1\n"
                " 2 1\n"
                " 3 2\n"
                " 4 2\n\n"
                " Bonds\n\n"
                " 1 1 1 2\n"
                " 2 2 1 3\n\n");
    auto output = END_CAPTURE_OUTPUT();
    ASSERT_THAT(output, MatchesRegex(".*Read molecule template.*1 molecules.*4 atoms.*type.*2.*"
                                     "2 bonds.*type.*2.*0 angles.*"));

    BEGIN_CAPTURE_OUTPUT();
    command("mass * 2.0");
    command("create_atoms 0 single 0.5 0.5 0.5 mol bonds 67235");
    output = END_CAPTURE_OUTPUT();
    ASSERT_THAT(output, MatchesRegex(".*Created 4 atoms.*"));

    BEGIN_HIDE_OUTPUT();
    Molecule *mol = lmp->atom->molecules[0];
    ASSERT_EQ(mol->natoms, 4);
    ASSERT_EQ(lmp->atom->natoms, 4);
    mol->compute_mass();
    mol->compute_com();
    ASSERT_DOUBLE_EQ(mol->masstotal, 8.0);
    EXPECT_DOUBLE_EQ(mol->com[0], 1.0);
    EXPECT_DOUBLE_EQ(mol->com[1], 0.5);
    EXPECT_DOUBLE_EQ(mol->com[2], 0.5);
    EXPECT_DOUBLE_EQ(mol->maxextent, sqrt(2));
    EXPECT_EQ(mol->comatom, 1);
    END_HIDE_OUTPUT();
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (Info::get_mpi_vendor() == "Open MPI" && !LAMMPS_NS::Info::has_exceptions())
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
