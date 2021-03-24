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

#include "atom.h"
#include "info.h"
#include "input.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstring>
#include <vector>
#include <functional>

// location of '*.py' files required by tests
#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val
std::string INPUT_FOLDER = STRINGIFY(TEST_INPUT_FOLDER);

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::MatchesRegex;
using ::testing::StrEq;
using ::testing::HasSubstr;

class PythonPackageTest : public ::testing::Test {
protected:
    LAMMPS *lmp;
    Info *info;

    void command(const std::string &line) { lmp->input->one(line.c_str()); }

    void HIDE_OUTPUT(std::function<void()> f) {
        if (!verbose) ::testing::internal::CaptureStdout();
        f();
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    std::string CAPTURE_OUTPUT(std::function<void()> f) {
        ::testing::internal::CaptureStdout();
        f();
        auto output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        return output;
    }

    void SetUp() override
    {
        const char *args[] = {"PythonPackageTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        HIDE_OUTPUT([&] {
            lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        });
        ASSERT_NE(lmp, nullptr);
        info = new Info(lmp);
        if (!info->has_package("PYTHON")) GTEST_SKIP();
        HIDE_OUTPUT([&] {
            command("units real");
            command("dimension 3");
            command("region box block -4 4 -4 4 -4 4");
            command("create_box 1 box");
            command("create_atoms 1 single  0.0  0.0 0.0    units box");
            command("create_atoms 1 single  1.9 -1.9 1.9999 units box");
            command("pair_style zero 2.0");
            command("pair_coeff * *");
            command("mass * 1.0");
            command("variable input_dir index " + INPUT_FOLDER);
        });
    }

    void TearDown() override
    {
        HIDE_OUTPUT([&] {
            delete info;
            delete lmp;
            info = nullptr;
            lmp = nullptr;
        });
    }
};

TEST_F(PythonPackageTest, InvokeFunctionFromFile)
{
    // execute python function from file
    HIDE_OUTPUT([&] {
        command("python printnum file ${input_dir}/func.py");
    });
    
    auto output = CAPTURE_OUTPUT([&]() {
        command("python printnum invoke");
    });
    ASSERT_THAT(output, HasSubstr("2.25\n"));
}

TEST_F(PythonPackageTest, InvokeOtherFunctionFromFile)
{
    // execute another python function from same file
    HIDE_OUTPUT([&] {
        command("python printnum file ${input_dir}/func.py");
        command("python printtxt exists");
    });

    auto output = CAPTURE_OUTPUT([&] {
        command("python printtxt invoke");
    });
    ASSERT_THAT(output, HasSubstr("sometext\n"));
}

TEST_F(PythonPackageTest, InvokeFunctionThatUsesLAMMPSModule)
{
    // execute python function that uses the LAMMPS python module
    HIDE_OUTPUT([&] {
        command("python printnum file ${input_dir}/func.py");
        command("variable idx equal 2.25");
        command("python getidxvar input 1 SELF format p exists");
    });
    auto output = CAPTURE_OUTPUT([&] {
        command("python getidxvar invoke");
    });
    ASSERT_THAT(output, HasSubstr("2.25\n"));
}

TEST_F(PythonPackageTest, python_variable)
{
    HIDE_OUTPUT([&] {
        command("variable sq python square");
        command("variable val index 1.5");
        command("python square input 1 v_val return v_sq format ff file ${input_dir}/func.py");
    });
    std::string output = CAPTURE_OUTPUT([&] {
        command("print \"${sq}\"");
    });
    ASSERT_THAT(output, MatchesRegex("print.*2.25.*"));
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
