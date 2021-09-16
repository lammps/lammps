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
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "library.h"
#include "utils.h"

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::Eq;

class InputConvertTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "InputConvertTest";
        LAMMPSTest::SetUp();
        ASSERT_NE(lmp, nullptr);
    }

    void TearDown() override { LAMMPSTest::TearDown(); }
};

TEST_F(InputConvertTest, logical)
{
    EXPECT_TRUE(utils::logical(FLERR, "yes", false, lmp));
    EXPECT_TRUE(utils::logical(FLERR, "Yes", false, lmp));
    EXPECT_TRUE(utils::logical(FLERR, "YES", false, lmp));
    EXPECT_TRUE(utils::logical(FLERR, "yEs", false, lmp));
    EXPECT_TRUE(utils::logical(FLERR, "true", false, lmp));
    EXPECT_TRUE(utils::logical(FLERR, "True", false, lmp));
    EXPECT_TRUE(utils::logical(FLERR, "TRUE", false, lmp));
    EXPECT_TRUE(utils::logical(FLERR, "on", false, lmp));
    EXPECT_TRUE(utils::logical(FLERR, "On", false, lmp));
    EXPECT_TRUE(utils::logical(FLERR, "ON", false, lmp));
    EXPECT_TRUE(utils::logical(FLERR, "1", false, lmp));
    EXPECT_FALSE(utils::logical(FLERR, "no", false, lmp));
    EXPECT_FALSE(utils::logical(FLERR, "No", false, lmp));
    EXPECT_FALSE(utils::logical(FLERR, "NO", false, lmp));
    EXPECT_FALSE(utils::logical(FLERR, "nO", false, lmp));
    EXPECT_FALSE(utils::logical(FLERR, "off", false, lmp));
    EXPECT_FALSE(utils::logical(FLERR, "Off", false, lmp));
    EXPECT_FALSE(utils::logical(FLERR, "OFF", false, lmp));
    EXPECT_FALSE(utils::logical(FLERR, "OfF", false, lmp));
    EXPECT_FALSE(utils::logical(FLERR, "0", false, lmp));

    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "yay", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "xxx", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "none", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "5", false, lmp););
}

TEST_F(InputConvertTest, numeric)
{
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, "0", false, lmp), 0);

    TEST_FAILURE(".*ERROR: Expected floating point.*", utils::numeric(FLERR, "yay", false, lmp););
}

TEST_F(InputConvertTest, inumeric)
{
    EXPECT_DOUBLE_EQ(utils::inumeric(FLERR, "0", false, lmp), 0);

    TEST_FAILURE(".*ERROR: Expected integer.*", utils::inumeric(FLERR, "yay", false, lmp););
}

TEST_F(InputConvertTest, bnumeric)
{
    EXPECT_DOUBLE_EQ(utils::bnumeric(FLERR, "0", false, lmp), 0);

    TEST_FAILURE(".*ERROR: Expected integer.*", utils::bnumeric(FLERR, "yay", false, lmp););
}
TEST_F(InputConvertTest, tnumeric)
{
    EXPECT_DOUBLE_EQ(utils::tnumeric(FLERR, "0", false, lmp), 0);

    TEST_FAILURE(".*ERROR: Expected integer.*", utils::tnumeric(FLERR, "yay", false, lmp););
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    lammps_mpi_init();
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
    lammps_mpi_finalize();
    return rv;
}
