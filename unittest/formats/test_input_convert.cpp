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
    EXPECT_EQ(utils::logical(FLERR, "yes", false, lmp), 1);
    EXPECT_EQ(utils::logical(FLERR, "true", false, lmp), 1);
    EXPECT_EQ(utils::logical(FLERR, "on", false, lmp), 1);
    EXPECT_EQ(utils::logical(FLERR, "1", false, lmp), 1);
    EXPECT_EQ(utils::logical(FLERR, "no", false, lmp), 0);
    EXPECT_EQ(utils::logical(FLERR, "false", false, lmp), 0);
    EXPECT_EQ(utils::logical(FLERR, "off", false, lmp), 0);
    EXPECT_EQ(utils::logical(FLERR, "0", false, lmp), 0);

    EXPECT_EQ(utils::logical(FLERR, std::string("yes"), false, lmp), 1);
    EXPECT_EQ(utils::logical(FLERR, std::string("true"), false, lmp), 1);
    EXPECT_EQ(utils::logical(FLERR, std::string("on"), false, lmp), 1);
    EXPECT_EQ(utils::logical(FLERR, std::string("1"), false, lmp), 1);
    EXPECT_EQ(utils::logical(FLERR, std::string("no"), false, lmp), 0);
    EXPECT_EQ(utils::logical(FLERR, std::string("false"), false, lmp), 0);
    EXPECT_EQ(utils::logical(FLERR, std::string("off"), false, lmp), 0);
    EXPECT_EQ(utils::logical(FLERR, std::string("0"), false, lmp), 0);

    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "YES", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "Yes", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "On", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "ON", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "TRUE", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "True", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "NO", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "No", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "Off", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "OFF", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "FALSE", false, lmp););
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of.*",
                 utils::logical(FLERR, "False", false, lmp););
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
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, "0.1", false, lmp), 0.1);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, "-.232", false, lmp), -0.232);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, ".2e5", false, lmp), 20000.0);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, "2.5e-10", false, lmp), 2.5e-10);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, "+0.3", false, lmp), 0.3);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, "10000000000", false, lmp), 1e10);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, "2.56E+3", false, lmp), 2560);

    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, std::string("0"), false, lmp), 0);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, std::string("0.1"), false, lmp), 0.1);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, std::string("-.232"), false, lmp), -0.232);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, std::string(".2e5"), false, lmp), 20000.0);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, std::string("2.5e-10"), false, lmp), 2.5e-10);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, std::string("+0.3"), false, lmp), 0.3);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, std::string("10000000000"), false, lmp), 1e10);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, std::string("2.56E+3"), false, lmp), 2560);

    TEST_FAILURE(".*ERROR: Expected floating point.*", utils::numeric(FLERR, "yay", false, lmp););
    TEST_FAILURE(".*ERROR: Expected floating point.*", utils::numeric(FLERR, "", false, lmp););
    TEST_FAILURE(".*ERROR: Expected floating point.*", utils::numeric(FLERR, nullptr, false, lmp););
    TEST_FAILURE(".*ERROR: Expected floating point.*",
                 utils::numeric(FLERR, "2.56D+3", false, lmp););
}

TEST_F(InputConvertTest, inumeric)
{
    EXPECT_EQ(utils::inumeric(FLERR, "0", false, lmp), 0);
    EXPECT_EQ(utils::inumeric(FLERR, "-1", false, lmp), -1);
    EXPECT_EQ(utils::inumeric(FLERR, "10000", false, lmp), 10000);
    EXPECT_EQ(utils::inumeric(FLERR, "-532410", false, lmp), -532410);
    EXPECT_EQ(utils::inumeric(FLERR, "-0", false, lmp), 0);
    EXPECT_EQ(utils::inumeric(FLERR, "0100", false, lmp), 100);

    EXPECT_EQ(utils::inumeric(FLERR, std::string("0"), false, lmp), 0);
    EXPECT_EQ(utils::inumeric(FLERR, std::string("-1"), false, lmp), -1);
    EXPECT_EQ(utils::inumeric(FLERR, std::string("10000"), false, lmp), 10000);
    EXPECT_EQ(utils::inumeric(FLERR, std::string("-532410"), false, lmp), -532410);
    EXPECT_EQ(utils::inumeric(FLERR, std::string("-0"), false, lmp), 0);
    EXPECT_EQ(utils::inumeric(FLERR, std::string("0100"), false, lmp), 100);

    TEST_FAILURE(".*ERROR: Expected integer.*", utils::inumeric(FLERR, "yay", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::inumeric(FLERR, "0.1", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::inumeric(FLERR, "1.1", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::inumeric(FLERR, "1e5", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::inumeric(FLERR, "0x05", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::inumeric(FLERR, "", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::inumeric(FLERR, nullptr, false, lmp););
}

TEST_F(InputConvertTest, bnumeric)
{
    EXPECT_EQ(utils::bnumeric(FLERR, "0", false, lmp), 0);
    EXPECT_EQ(utils::bnumeric(FLERR, "-1", false, lmp), -1);
    EXPECT_EQ(utils::bnumeric(FLERR, "10000", false, lmp), 10000);
    EXPECT_EQ(utils::bnumeric(FLERR, "-532410", false, lmp), -532410);
    EXPECT_EQ(utils::bnumeric(FLERR, "-0", false, lmp), 0);
    EXPECT_EQ(utils::bnumeric(FLERR, "0100", false, lmp), 100);

    EXPECT_EQ(utils::bnumeric(FLERR, std::string("0"), false, lmp), 0);
    EXPECT_EQ(utils::bnumeric(FLERR, std::string("-1"), false, lmp), -1);
    EXPECT_EQ(utils::bnumeric(FLERR, std::string("10000"), false, lmp), 10000);
    EXPECT_EQ(utils::bnumeric(FLERR, std::string("-532410"), false, lmp), -532410);
    EXPECT_EQ(utils::bnumeric(FLERR, std::string("-0"), false, lmp), 0);
    EXPECT_EQ(utils::bnumeric(FLERR, std::string("0100"), false, lmp), 100);

    TEST_FAILURE(".*ERROR: Expected integer.*", utils::bnumeric(FLERR, "yay", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::bnumeric(FLERR, "0.1", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::bnumeric(FLERR, "1.1", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::bnumeric(FLERR, "1e5", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::bnumeric(FLERR, "0x05", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::bnumeric(FLERR, "", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::bnumeric(FLERR, nullptr, false, lmp););
}

TEST_F(InputConvertTest, tnumeric)
{
    EXPECT_EQ(utils::tnumeric(FLERR, "0", false, lmp), 0);
    EXPECT_EQ(utils::tnumeric(FLERR, "-1", false, lmp), -1);
    EXPECT_EQ(utils::tnumeric(FLERR, "10000", false, lmp), 10000);
    EXPECT_EQ(utils::tnumeric(FLERR, "-532410", false, lmp), -532410);
    EXPECT_EQ(utils::tnumeric(FLERR, "-0", false, lmp), 0);
    EXPECT_EQ(utils::tnumeric(FLERR, "0100", false, lmp), 100);

    EXPECT_EQ(utils::tnumeric(FLERR, std::string("0"), false, lmp), 0);
    EXPECT_EQ(utils::tnumeric(FLERR, std::string("-1"), false, lmp), -1);
    EXPECT_EQ(utils::tnumeric(FLERR, std::string("10000"), false, lmp), 10000);
    EXPECT_EQ(utils::tnumeric(FLERR, std::string("-532410"), false, lmp), -532410);
    EXPECT_EQ(utils::tnumeric(FLERR, std::string("-0"), false, lmp), 0);
    EXPECT_EQ(utils::tnumeric(FLERR, std::string("0100"), false, lmp), 100);

    TEST_FAILURE(".*ERROR: Expected integer.*", utils::tnumeric(FLERR, "yay", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::tnumeric(FLERR, "0.1", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::tnumeric(FLERR, "1.1", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::tnumeric(FLERR, "1e5", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::tnumeric(FLERR, "0x05", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::tnumeric(FLERR, "", false, lmp););
    TEST_FAILURE(".*ERROR: Expected integer.*", utils::tnumeric(FLERR, nullptr, false, lmp););
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
