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
#include "group.h"
#include "info.h"
#include "input.h"
#include "math_const.h"
#include "region.h"
#include "variable.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstring>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
using ::testing::ContainsRegex;
using ::testing::ExitedWithCode;
using ::testing::StrEq;

class DeleteAtomsTest : public LAMMPSTest {
protected:
    Atom *atom;

    void SetUp() override
    {
        testbinary = "DeleteAtomsTest";
        args       = {"-log", "none", "-echo", "screen", "-nocite", "-v", "num", "1"};
        LAMMPSTest::SetUp();
        atom = lmp->atom;
    }

    void TearDown() override
    {
        LAMMPSTest::TearDown();
        platform::unlink("test_variable.file");
        platform::unlink("test_variable.atomfile");
    }

    void atomic_system()
    {
        BEGIN_HIDE_OUTPUT();
        command("units real");
        command("lattice sc 1.0 origin 0.125 0.125 0.125");
        command("region box block -4 4 -4 4 -4 4");
        command("create_box 8 box");
        command("create_atoms 1 box");
        command("mass * 1.0");
        command("region left block -2.0 -1.0 INF INF INF INF");
        command("region right block 0.5  2.0 INF INF INF INF");
        command("region top block INF INF -2.0 -1.0 INF INF");
        command("region bottom block INF INF 0.0 4.0 INF INF");
        command("set region left type 2");
        command("set region right type 3");
        command("group bottom region bottom");
        command("group top region top");
        END_HIDE_OUTPUT();
    }

    void molecular_system()
    {
        HIDE_OUTPUT([&] {
            command("fix props all property/atom mol rmass q");
        });
        atomic_system();
        BEGIN_HIDE_OUTPUT();
        command("variable molid atom floor(id/4)+1");
        command("variable charge atom 2.0*sin(PI/32*id)");
        command("set atom * mol v_molid");
        command("set atom * charge v_charge");
        command("set type 1 mass 0.5");
        command("set type 2*4 mass 2.0");
        END_HIDE_OUTPUT();
    }
};

TEST_F(DeleteAtomsTest, Simple)
{
    atomic_system();
    ASSERT_EQ(atom->natoms, 512);

    HIDE_OUTPUT([&] {
        command("delete_atoms group top");
    });
    ASSERT_EQ(atom->natoms, 448);

    HIDE_OUTPUT([&] {
        command("delete_atoms region left");
    });
    ASSERT_EQ(atom->natoms, 392);

    HIDE_OUTPUT([&] {
        command("delete_atoms random fraction 0.5 yes all right 43252");
    });
    ASSERT_EQ(atom->natoms, 364);

    HIDE_OUTPUT([&] {
        command("variable checker atom sin(4*PI*x/lx)*sin(4*PI*y/ly)*sin(4*PI*z/lz)>0");
        command("delete_atoms variable checker");
    });
    ASSERT_EQ(atom->natoms, 178);

    HIDE_OUTPUT([&] {
        command("delete_atoms random count 3 no bottom right 443252");
    });
    ASSERT_EQ(atom->natoms, 175);

    HIDE_OUTPUT([&] {
        command("delete_atoms random count 50 no all NULL 434325");
    });
    ASSERT_EQ(atom->natoms, 125);

    HIDE_OUTPUT([&] {
        command("delete_atoms random fraction 0.2 no all NULL 34325");
    });
    ASSERT_EQ(atom->natoms, 104);

    HIDE_OUTPUT([&] {
        command("delete_atoms random count 50 no bottom right 77325");
    });
    ASSERT_EQ(atom->natoms, 102);

    TEST_FAILURE(".*ERROR: Illegal delete_atoms command: missing argument.*",
                 command("delete_atoms"););
    TEST_FAILURE(".*ERROR: Unknown delete_atoms sub-command: xxx.*", command("delete_atoms xxx"););
    TEST_FAILURE(".*ERROR: The delete_atoms 'porosity' keyword has been removed.*",
                 command("delete_atoms porosity 0.5 all right 4325234"););
    TEST_FAILURE(".*ERROR: Illegal delete_atoms random command: missing argument.*",
                 command("delete_atoms random count 50 bottom right 77325"););
    TEST_FAILURE(".*ERROR: Illegal delete_atoms random command: missing argument.*",
                 command("delete_atoms random fraction 0.4 bottom right 77325"););
    TEST_FAILURE(".*ERROR: Delete_atoms random count has invalid value: -5.*",
                 command("delete_atoms random count -5 no bottom right 77325"););
    TEST_FAILURE(".*ERROR: Delete_atoms count of 5 exceeds number of eligible atoms 0.*",
                 command("delete_atoms random count 5 yes bottom right 77325"););
    TEST_FAILURE(".*ERROR: Delete_atoms random fraction has invalid value: -0.4.*",
                 command("delete_atoms random fraction -0.4 no bottom right 77325"););
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
