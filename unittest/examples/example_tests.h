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
#ifndef LMP_EXAMPLE_TESTS_H
#define LMP_EXAMPLE_TESTS_H

// Common infrastructure for tests running abbreviated versions of the
// production-like example inputs that the regression tests leave out
// (see EXCLUDED_FOLDERS in tools/regression-tests/run_tests.py).  The
// example inputs define their adjustable run lengths and system sizes
// as index-style variables, so a test can preset shortened values with
// preset() before including the unmodified example input with
// run_input(): index-style variable definitions in the input are
// skipped when the variable already exists.

#include "../testing/core.h"
#include "../testing/utils.h"

#include "output.h"
#include "thermo.h"
#include "update.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mpi.h>
#include <string>
#include <utility>
#include <vector>

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

class ExampleTest : public LAMMPSTest {
protected:
    // report the first style from the list that is not available in the
    // LAMMPS library, or an empty string when all are present.  used by
    // the REQUIRE_STYLES macro below, since GTEST_SKIP() must be issued
    // in the test body itself to skip the whole test
    std::string missing_style(const std::vector<std::pair<std::string, std::string>> &styles)
    {
        for (const auto &[category, name] : styles) {
            if (!info->has_style(category, name))
                return fmt::format("{} style {} is not available", category, name);
        }
        return "";
    }

    // preset an index-style variable to override the default in the input
    void preset(const std::string &name, const std::string &value)
    {
        HIDE_OUTPUT([&] {
            command(fmt::format("variable {} index {}", name, value));
        });
    }

    // run one input script file from the examples folder set at compile
    // time, on a cleared system (index-style variable presets survive the
    // clear command)
    void run_input(const std::string &script)
    {
        BEGIN_HIDE_OUTPUT();
        command("clear");
        command("include \"" STRINGIFY(TEST_EXAMPLES_FOLDER) "/" + script + "\"");
        END_HIDE_OUTPUT();
    }

    // copy a file (may be in a subfolder) from the examples folder to the
    // current working directory, for inputs reading data files or include
    // files with relative paths
    void copy_from_examples(const std::string &name)
    {
        std::ifstream src(std::string(STRINGIFY(TEST_EXAMPLES_FOLDER) "/") + name,
                          std::ios::binary);
        ASSERT_TRUE(src.is_open()) << "cannot read " << name;
        const auto base = platform::path_basename(name);
        std::ofstream dst(base, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(dst.is_open()) << "cannot write " << base;
        dst << src.rdbuf();
    }

    double thermo_value(const std::string &keyword)
    {
        double value = 0.0;
        lmp->output->thermo->evaluate_keyword(keyword, &value);
        return value;
    }
};

#define REQUIRE_STYLES(...)                            \
    {                                                  \
        auto missing = missing_style({__VA_ARGS__});   \
        if (!missing.empty()) GTEST_SKIP() << missing; \
    }

// last output block of a fix ave/time "mode vector" file (2 words in the
// per-block header line: timestep and number-of-rows) or of a fix
// ave/chunk file (3 words: timestep, number-of-chunks, total-count),
// as rows of columns
static std::vector<std::vector<double>> last_vector_block(const std::string &filename,
                                                          std::size_t header_words = 2)
{
    std::vector<std::vector<double>> block;
    std::ifstream data(filename);
    if (!data.is_open()) return block;

    std::string line;
    std::size_t nrows = 0;
    while (std::getline(data, line)) {
        auto words = utils::split_words(line);
        if (words.empty() || (words[0][0] == '#')) continue;
        if (words.size() == header_words) {
            nrows = std::stoul(words[1]);
            block.clear();
        } else {
            std::vector<double> row;
            for (const auto &word : words)
                row.push_back(std::stod(word));
            block.push_back(row);
        }
    }
    EXPECT_EQ(block.size(), nrows);
    return block;
}

} // namespace LAMMPS_NS

// common main() for all example tests
#define EXAMPLE_TEST_MAIN()                                                    \
    int main(int argc, char **argv)                                            \
    {                                                                          \
        MPI_Init(&argc, &argv);                                                \
        ::testing::InitGoogleMock(&argc, argv);                                \
                                                                               \
        /* handle arguments passed via environment variable */                 \
        if (const char *var = getenv("TEST_ARGS")) {                           \
            std::vector<std::string> env = LAMMPS_NS::utils::split_words(var); \
            for (const auto &arg : env) {                                      \
                if (arg == "-v") verbose = true;                               \
            }                                                                  \
        }                                                                      \
        if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;        \
                                                                               \
        int rv = RUN_ALL_TESTS();                                              \
        MPI_Finalize();                                                        \
        return rv;                                                             \
    }

#endif
