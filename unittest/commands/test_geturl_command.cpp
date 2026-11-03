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

// The geturl command downloads from the network, which takes long and varies
// with the connection, so this test lives in its own executable carrying the
// "slow" ctest label instead of delaying the SimpleCommands tests.

#include "lammps.h"

#include "info.h"
#include "input.h"
#include "platform.h"
#include "utils.h"

#include "../testing/core.h"
#include "../testing/utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstring>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

class GeturlTest : public LAMMPSTest {};

TEST_F(GeturlTest, Geturl)
{
    if (!Info::has_package("EXTRA-COMMAND")) GTEST_SKIP();
    platform::unlink("index.html");
    platform::unlink("myindex.html");
    if (Info::has_curl_support()) {
        BEGIN_CAPTURE_OUTPUT();
        command("geturl https://github.com/");
        command("geturl https://github.com/ output myindex.html");
        END_CAPTURE_OUTPUT();
        EXPECT_TRUE(platform::file_is_readable("index.html"));
        EXPECT_TRUE(platform::file_is_readable("myindex.html"));
        FILE *fp = fopen("index.html", "wb");
        fputs("just testing\n", fp);
        fclose(fp);
        BEGIN_CAPTURE_OUTPUT();
        command("geturl https://github.com/ overwrite no");
        END_CAPTURE_OUTPUT();
        char checkme[20];
        fp = fopen("index.html", "rb");
        fgets(checkme, 19, fp);
        fclose(fp);
        EXPECT_EQ(strcmp(checkme, "just testing\n"), 0);
        BEGIN_CAPTURE_OUTPUT();
        command("geturl https://github.com/ overwrite yes");
        END_CAPTURE_OUTPUT();
        fp = fopen("index.html", "rb");
        fgets(checkme, 19, fp);
        fclose(fp);
        EXPECT_NE(strcmp(checkme, "just testing\n"), 0);
        platform::unlink("index.html");
        BEGIN_CAPTURE_OUTPUT();
        command("geturl https://github.com");
        END_CAPTURE_OUTPUT();
        EXPECT_TRUE(platform::file_is_readable("index.html"));

        TEST_FAILURE(".*ERROR: Illegal geturl command: missing argument.*", command("geturl "););
        TEST_FAILURE(".*ERROR: URL 'dummy' is not a supported URL.*", command("geturl dummy"););
        TEST_FAILURE(".*ERROR: URL '/tmp' is not a supported URL.*", command("geturl /tmp"););
        TEST_FAILURE(".*ERROR on proc 0: Download of xxx.txt failed.*",
                     command("geturl https://github.com/xxx.txt"););
    } else {
        TEST_FAILURE(".*ERROR: LAMMPS has not been compiled with libcurl support*",
                     command("geturl https:://github.com/"););
    }
    platform::unlink("index.html");
    platform::unlink("myindex.html");
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

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
