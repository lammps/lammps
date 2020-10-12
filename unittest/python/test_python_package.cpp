/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom.h"
#include "input.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstring>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::Eq;


class PythonPackageTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {"PythonPackageTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        ASSERT_NE(lmp, nullptr);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp->input->one("units real");
        lmp->input->one("dimension 3");
        lmp->input->one("region box block -4 4 -4 4 -4 4");
        lmp->input->one("create_box 1 box");
        lmp->input->one("create_atoms 1 single  0.0  0.0 0.0    units box");
        lmp->input->one("create_atoms 1 single  1.9 -1.9 1.9999 units box");
        lmp->input->one("pair_style zero 2.0");
        lmp->input->one("pair_coeff * *");
        lmp->input->one("mass * 1.0");
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(PythonPackageTest, python_invoke)
{
}
    
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

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
