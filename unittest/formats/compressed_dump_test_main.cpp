/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "../testing/utils.h"
#include "compressed_dump_test.h"
#include "fmt/format.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>

const char *COMPRESS_SUFFIX    = nullptr;
const char *COMPRESS_EXTENSION = nullptr;
char *COMPRESS_EXECUTABLE      = nullptr;
bool verbose                   = false;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " (gz|zstd)\n\n" << std::endl;
        return 1;
    }

    if (strcmp(argv[1], "gz") == 0) {
        COMPRESS_SUFFIX    = "gz";
        COMPRESS_EXTENSION = "gz";
    } else if (strcmp(argv[1], "zstd") == 0) {
        COMPRESS_SUFFIX    = "zstd";
        COMPRESS_EXTENSION = "zst";
    } else {
        std::cerr << "usage: " << argv[0] << " (gz|zstd)\n\n" << std::endl;
        return 1;
    }

    COMPRESS_EXECUTABLE = getenv("COMPRESS_EXECUTABLE");

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
