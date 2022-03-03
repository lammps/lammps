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
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "library.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

class ComputeGlobalTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "ComputeGlobalTest";
        LAMMPSTest::SetUp();
        if (info->has_style("atom", "full")) {
            BEGIN_HIDE_OUTPUT();
            command("variable input_dir index \"" STRINGIFY(TEST_INPUT_FOLDER) "\"");
            command("include \"${input_dir}/in.fourmol\"");
            END_HIDE_OUTPUT();
        }
    }
};

TEST_F(ComputeGlobalTest, Energy)
{
    void *handle = (void *) lmp;
    if (lammps_get_natoms(handle) == 0.0) GTEST_SKIP();

    BEGIN_HIDE_OUTPUT();
    command("group allwater molecule 3:6");
    command("compute ke1 all ke");
    command("compute ke2 allwater ke");
    command("compute pe1 all pe");
    command("compute pe2 all pe bond");
    command("compute pe3 all pe angle dihedral");
    command("pair_style lj/cut 10.0");
    command("pair_coeff * * 0.01 3.0");
    command("bond_style harmonic");
    command("bond_coeff * 100.0 1.5");
    command("thermo_style custom c_ke1 c_ke2 c_pe1 c_pe2 c_pe3");
    command("run 0 post no");
    END_HIDE_OUTPUT();

    double ke1 = *(double *)lammps_extract_compute(handle, "ke1", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR);
    double ke2 = *(double *)lammps_extract_compute(handle, "ke2", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR);
    double pe1 = *(double *)lammps_extract_compute(handle, "pe1", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR);
    double pe2 = *(double *)lammps_extract_compute(handle, "pe2", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR);
    double pe3 = *(double *)lammps_extract_compute(handle, "pe3", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR);

    EXPECT_DOUBLE_EQ(ke1, 2.3405256449146168);
    EXPECT_DOUBLE_EQ(ke2, 1.192924237073665);
    EXPECT_DOUBLE_EQ(pe1, 24280.922367235136);
    EXPECT_DOUBLE_EQ(pe2, 361.37528652881286);
    EXPECT_DOUBLE_EQ(pe3, 0.0);
    

    TEST_FAILURE(".*ERROR: Reuse of compute ID 'pe2'.*", command("compute pe2 all pe"););
    TEST_FAILURE(".*ERROR: Compute pe must use group all.*", command("compute pe allwater pe"););
    TEST_FAILURE(".*ERROR: Illegal compute command.*", command("compute pe potential"););
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (platform::mpi_vendor() == "Open MPI" && !LAMMPS_NS::Info::has_exceptions())
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
