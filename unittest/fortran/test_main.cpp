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

#include "test_main.h"
#include "pointers.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <iostream>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (argc < 1) {
        return 1;
    }

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    int iarg = 1;
    while (iarg < argc) {
        if (strcmp(argv[iarg], "-v") == 0) {
            verbose = true;
            ++iarg;
        } else {
            std::cerr << "unknown option: " << argv[iarg] << "\n\n";
            MPI_Finalize();
            return 1;
        }
    }

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
