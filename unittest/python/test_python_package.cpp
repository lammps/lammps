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

#include "atom.h"
#include "info.h"
#include "input.h"
#include "library.h"
#include "variable.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstring>
#include <functional>
#include <vector>

#include "../testing/core.h"
#include "../testing/systems/melt.h"

#if defined(TEST_HAVE_PYTHON_DEVELOPMENT)
#include <Python.h>
#endif

// location of '*.py' files required by tests
#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val
std::string INPUT_FOLDER = STRINGIFY(TEST_INPUT_FOLDER);

const char *LOREM_IPSUM =
    "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Praesent metus.";
bool verbose = false;

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::ContainsRegex;
using ::testing::Eq;
using ::testing::HasSubstr;
using ::testing::StrEq;

class PythonPackageTest : public LAMMPSTest {
protected:
    void InitSystem() override
    {
        if (!Info::has_package("PYTHON")) GTEST_SKIP();

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
};

class FixPythonInvokeTest : public MeltTest {
protected:
    void InitSystem() override
    {
        if (!Info::has_package("PYTHON")) GTEST_SKIP();

        MeltTest::InitSystem();
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
    ASSERT_THAT(output, HasSubstr("2.25"));
}

#if defined(TEST_HAVE_PYTHON_DEVELOPMENT)
TEST_F(PythonPackageTest, InvokeInitialized)
{
    // execute python function from file
    HIDE_OUTPUT([&] {
        command("python printnum file ${input_dir}/func.py");
    });
    ASSERT_TRUE(Py_IsInitialized());
    HIDE_OUTPUT([&] {
        command("clear");
    });
    ASSERT_TRUE(Py_IsInitialized());
    lammps_python_finalize();
    ASSERT_FALSE(Py_IsInitialized());
}
#endif

TEST_F(PythonPackageTest, InvokeFunctionPassInt)
{
    // execute python function, passing integer as argument
    HIDE_OUTPUT([&] {
        command("variable sq python square");
        command("python square input 1 2 format ii return v_sq file ${input_dir}/func.py");
        command("python square invoke");
    });

    ASSERT_EQ(get_variable_value("sq"), 4.0);
}

TEST_F(PythonPackageTest, InvokeFunctionPassFloat)
{
    // execute python function, passing float as argument
    HIDE_OUTPUT([&] {
        command("variable sq python square");
        command("python square input 1 2.5 format ff return v_sq file ${input_dir}/func.py");
    });

    ASSERT_EQ(get_variable_value("sq"), 6.25);
}

TEST_F(PythonPackageTest, InvokeFunctionPassString)
{
    // execute python function, passing string as argument
    HIDE_OUTPUT([&] {
        command("variable val python bool_to_val");
        command(
            "python bool_to_val input 1 \"true\" format sf return v_val file ${input_dir}/func.py");
    });

    ASSERT_EQ(get_variable_value("val"), 1.0);
}

TEST_F(PythonPackageTest, InvokeFunctionPassStringVariable)
{
    // execute python function, passing string variable as argument
    HIDE_OUTPUT([&] {
        command("variable val python bool_to_val");
        command(
            "python bool_to_val input 1 v_str format sf return v_val file ${input_dir}/func.py");
    });

    HIDE_OUTPUT([&] {
        command("variable str string \"true\"");
    });

    ASSERT_EQ(get_variable_value("val"), 1.0);

    HIDE_OUTPUT([&] {
        command("variable str string \"false\"");
    });

    ASSERT_EQ(get_variable_value("val"), 0.0);
}

TEST_F(PythonPackageTest, InvokeStringFunction)
{
    // execute python function, passing string variable as argument
    HIDE_OUTPUT([&] {
        command("variable str python val_to_bool");
        command(
            "python val_to_bool input 1 v_val format is return v_str file ${input_dir}/func.py");
    });

    HIDE_OUTPUT([&] {
        command("variable val equal 0");
    });

    ASSERT_THAT(get_variable_string("str"), StrEq("False"));

    HIDE_OUTPUT([&] {
        command("variable val equal 1");
    });

    ASSERT_THAT(get_variable_string("str"), StrEq("True"));
}

TEST_F(PythonPackageTest, InvokeLongStringFunction)
{
    // execute python function, passing string variable as argument
    HIDE_OUTPUT([&] {
        command("variable str python longstr");
        command("python longstr format s length 72 return v_str file ${input_dir}/func.py");
    });

    ASSERT_THAT(get_variable_string("str"), StrEq(LOREM_IPSUM));
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
    ASSERT_THAT(output, HasSubstr("sometext"));
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
    ASSERT_THAT(output, HasSubstr("2.25"));
}

TEST_F(PythonPackageTest, python_variable)
{
    // define variable that evaluates a python function
    HIDE_OUTPUT([&] {
        command("variable sq python square");
        command("variable val index 1.5");
        command("python square input 1 v_val return v_sq format ff file ${input_dir}/func.py");
    });
    std::string output = CAPTURE_OUTPUT([&] {
        command("print \"${sq}\"");
    });
    ASSERT_THAT(output, ContainsRegex("print.*\n.*2.25.*"));
}

TEST_F(PythonPackageTest, InlineFunction)
{
    // define variable that evaluates a python function
    HIDE_OUTPUT([&] {
        command("variable fact python factorial");
        command("python factorial input 1 v_n return v_fact format ii here \"\"\"\n"
                "def factorial(n):\n"
                "  if n == 0 or n == 1: return 1\n"
                "  return n*factorial(n-1)\n"
                "\"\"\"");
    });

    HIDE_OUTPUT([&] {
        command("variable n equal 1");
    });

    ASSERT_EQ(get_variable_value("fact"), 1.0);

    HIDE_OUTPUT([&] {
        command("variable n equal 2");
    });

    ASSERT_EQ(get_variable_value("fact"), 2.0);

    HIDE_OUTPUT([&] {
        command("variable n equal 3");
    });

    ASSERT_EQ(get_variable_value("fact"), 6.0);
}

TEST_F(PythonPackageTest, RunSource)
{
    // execute python script from file
    auto output = CAPTURE_OUTPUT([&] {
        command("python source ${input_dir}/run.py");
    });

    ASSERT_THAT(output, HasSubstr(LOREM_IPSUM));
}

TEST_F(PythonPackageTest, RunSourceInline)
{
    // execute inline python script
    auto output = CAPTURE_OUTPUT([&] {
        command("python source here \"\"\"\n"
                "from __future__ import print_function\n"
                "print(2+2)\n"
                "\"\"\"");
    });

    ASSERT_THAT(output, HasSubstr("4"));
}

TEST_F(FixPythonInvokeTest, end_of_step)
{
    HIDE_OUTPUT([&] {
        command("python end_of_step_callback here \"\"\"\n"
                "from __future__ import print_function\n"
                "def end_of_step_callback(ptr):\n"
                "    print(\"PYTHON_END_OF_STEP\")\n"
                "\"\"\"");
        command("fix eos all python/invoke 10 end_of_step end_of_step_callback");
    });

    auto output = CAPTURE_OUTPUT([&] {
        command("run 50");
    });
    fprintf(stderr, "lines: %s\n", output.c_str());
    auto lines = utils::split_lines(output);
    int count  = 0;

    for (auto &line : lines) {
        if (line == "PYTHON_END_OF_STEP") ++count;
    }

    ASSERT_EQ(count, 5);
}

TEST_F(FixPythonInvokeTest, post_force)
{
    HIDE_OUTPUT([&] {
        command("python post_force_callback here \"\"\"\n"
                "from __future__ import print_function\n"
                "def post_force_callback(ptr, vflag):\n"
                "    print(\"PYTHON_POST_FORCE\")\n"
                "\"\"\"");
        command("fix pf all python/invoke 10 post_force post_force_callback");
    });

    auto output = CAPTURE_OUTPUT([&] {
        command("run 50");
    });

    auto lines = utils::split_lines(output);
    int count  = 0;

    for (auto &line : lines) {
        if (line == "PYTHON_POST_FORCE") ++count;
    }

    ASSERT_EQ(count, 5);
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        auto env = split_words(var);
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
