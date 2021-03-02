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

#include "info.h"
#include "input.h"
#include "lammps.h"
#include "molecule.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <mpi.h>
#include <fstream>
#include <string>

using namespace LAMMPS_NS;

using testing::MatchesRegex;
using testing::StrEq;

using utils::split_words;

#if defined(OMPI_MAJOR_VERSION)
const bool have_openmpi = true;
#else
const bool have_openmpi = false;
#endif

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

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

class MoleculeFileTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {"MoleculeFileTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        ASSERT_NE(lmp, nullptr);
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void run_mol_cmd(const std::string &name, const std::string &args, const std::string &content)
    {
        std::string file = name + ".mol";
        std::ofstream mol(file, std::ofstream::trunc);
        mol << content << std::endl;
        mol.close();

        lmp->input->one(fmt::format("molecule {} {} {}",name,file,args));
        remove(file.c_str());
    }       
};

TEST_F(MoleculeFileTest, minimal)
{
    ::testing::internal::CaptureStdout();
    run_mol_cmd(test_info_->name(),"","Comment\n"
                "1 atoms\n\n"
                " Coords\n\n"
                " 1 0.0 0.0 0.0\n");
    auto output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    ASSERT_THAT(output,MatchesRegex(".*Read molecule template.*1 molecules.*1 atoms.*0 bonds.*"));
}

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
