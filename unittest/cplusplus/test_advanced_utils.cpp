// unit tests for utils:: functions requiring a LAMMPS instance

#include "exceptions.h"
#include "info.h"
#include "input.h"
#include "memory.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstring>

using ::testing::StrEq;

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
class Error;

class AdvancedUtils : public LAMMPSTest {
protected:
    Error *error;

    void SetUp() override
    {
        testbinary = "AdvancedUtils";
        LAMMPSTest::SetUp();
        error = lmp->error;
    }

    void atomic_system()
    {
        BEGIN_HIDE_OUTPUT();
        command("units real");
        command("lattice sc 1.0 origin 0.125 0.125 0.125");
        command("region box block -2 2 -2 2 -2 2");
        command("create_box 8 box");
        command("create_atoms 1 box");
        command("pair_style zero 3.5");
        command("pair_coeff * *");
        command("mass * 1.0");
        command("region left block -2.0 -1.0 INF INF INF INF");
        command("region right block 0.5  2.0 INF INF INF INF");
        command("region top block INF INF -2.0 -1.0 INF INF");
        command("set region left type 2");
        command("set region right type 3");
        END_HIDE_OUTPUT();
    }
};

TEST_F(AdvancedUtils, missing_cmd_args)
{
    auto output = CAPTURE_OUTPUT([&] {
        utils::missing_cmd_args(FLERR, "dummy", nullptr);
    });
    EXPECT_EQ(output, "");

    TEST_FAILURE("ERROR: Illegal dummy command: missing argument",
                 utils::missing_cmd_args(FLERR, "dummy", error););
};

TEST_F(AdvancedUtils, logmesg)
{
    auto output = CAPTURE_OUTPUT([&] {
        utils::logmesg(lmp, "test message");
    });
    EXPECT_EQ(output, "test message");

    output = CAPTURE_OUTPUT([&] {
        utils::logmesg(lmp, "test message from test {}", testbinary);
    });
    EXPECT_EQ(output, "test message from test " + testbinary);
};

// death tests only. the other cases are tested in the basic utils unit tester
TEST_F(AdvancedUtils, bounds_int_fail)
{
    int nlo, nhi;
    TEST_FAILURE("ERROR: Invalid range string: 1x ",
                 utils::bounds(FLERR, "1x", 1, 10, nlo, nhi, error););
    TEST_FAILURE("ERROR: Invalid range string: -1 ",
                 utils::bounds(FLERR, "-1", 1, 10, nlo, nhi, error););
    TEST_FAILURE("ERROR: Invalid range string: \\+1 ",
                 utils::bounds(FLERR, "+1", 1, 10, nlo, nhi, error););
    TEST_FAILURE("ERROR: Invalid range string: 1:3 ",
                 utils::bounds(FLERR, "1:3", 1, 10, nlo, nhi, error););
    TEST_FAILURE("ERROR: Invalid range string: \\? ",
                 utils::bounds(FLERR, "?", -10, 5, nlo, nhi, error););
    TEST_FAILURE("ERROR: Invalid range string: 3\\*:2 ",
                 utils::bounds(FLERR, "3*:2", -10, 5, nlo, nhi, error););
    TEST_FAILURE("ERROR: Numeric index 1 is out of bounds \\(5-6\\).*",
                 utils::bounds(FLERR, "1*4", 5, 6, nlo, nhi, error););
    TEST_FAILURE("ERROR: Numeric index 4 is out of bounds \\(1-3\\).*",
                 utils::bounds(FLERR, "1*4", 1, 3, nlo, nhi, error););
}

TEST_F(AdvancedUtils, bounds_bigint_fail)
{
    bigint nlo, nhi;
    TEST_FAILURE("ERROR: Invalid range string: 1x ",
                 utils::bounds(FLERR, "1x", 1, 10, nlo, nhi, error););
    TEST_FAILURE("ERROR: Invalid range string: -1 ",
                 utils::bounds(FLERR, "-1", 1, 10, nlo, nhi, error););
    TEST_FAILURE("ERROR: Invalid range string: \\+1 ",
                 utils::bounds(FLERR, "+1", 1, 10, nlo, nhi, error););
    TEST_FAILURE("ERROR: Invalid range string: 1:3 ",
                 utils::bounds(FLERR, "1:3", 1, 10, nlo, nhi, error););
    TEST_FAILURE("ERROR: Invalid range string: \\? ",
                 utils::bounds(FLERR, "?", -10, 5, nlo, nhi, error););
    TEST_FAILURE("ERROR: Invalid range string: 3\\*:2 ",
                 utils::bounds(FLERR, "3*:2", -10, 5, nlo, nhi, error););
    TEST_FAILURE("ERROR: Numeric index 1 is out of bounds \\(5-6\\).*",
                 utils::bounds(FLERR, "1*4", 5, 6, nlo, nhi, error););
    TEST_FAILURE("ERROR: Numeric index 4 is out of bounds \\(1-3\\).*",
                 utils::bounds(FLERR, "1*4", 1, 3, nlo, nhi, error););
}

TEST_F(AdvancedUtils, expand_args)
{
    atomic_system();
    BEGIN_CAPTURE_OUTPUT();
    try {
        command("compute temp all temp");
        command("variable temp vector c_temp");
        command("variable step equal step");
        command("variable pe equal pe");
        command("variable pe equal pe");
        command("variable epair equal epair");
        command("compute gofr all rdf 20 1 1 1 2");
        command("fix 1 all ave/time 1 1 1 v_step v_pe v_epair");
        command("fix 2 all nve");
        command("run 1 post no");
    } catch (LAMMPSAbortException &ae) {
        fprintf(stderr, "LAMMPS Error: %s\n", ae.what());
        exit(2);
    } catch (LAMMPSException &e) {
        fprintf(stderr, "LAMMPS Error: %s\n", e.what());
        exit(3);
    } catch (fmt::format_error &fe) {
        fprintf(stderr, "fmt::format_error: %s\n", fe.what());
        exit(4);
    } catch (std::exception &e) {
        fprintf(stderr, "General exception: %s\n", e.what());
        exit(5);
    }

    auto output = END_CAPTURE_OUTPUT();
    if (verbose) std::cout << output << std::endl;

    char **args, **earg;
    constexpr int oarg = 9;
    args               = new char *[oarg];
    args[0]            = utils::strdup("v_step");
    args[1]            = utils::strdup("c_temp");
    args[2]            = utils::strdup("f_1[*]");
    args[3]            = utils::strdup("c_temp[2*4]");
    args[4]            = utils::strdup("v_temp[*4]");
    args[5]            = utils::strdup("c_gofr[3*]");
    args[6]            = utils::strdup("c_gofr[1][*]");
    args[7]            = utils::strdup("c_gofr[*2][2]");
    args[8]            = utils::strdup("c_gofr[*][*]");

    // disable use of input->command and input->arg which point to the last run command right now
    lmp->input->command = nullptr;
    lmp->input->arg     = nullptr;

    auto narg = utils::expand_args(FLERR, oarg, args, 0, earg, lmp);
    EXPECT_EQ(narg, 16);
    EXPECT_STREQ(earg[0], "v_step");
    EXPECT_STREQ(earg[1], "c_temp");
    EXPECT_STREQ(earg[2], "f_1[1]");
    EXPECT_STREQ(earg[3], "f_1[2]");
    EXPECT_STREQ(earg[4], "f_1[3]");
    EXPECT_STREQ(earg[5], "c_temp[2]");
    EXPECT_STREQ(earg[6], "c_temp[3]");
    EXPECT_STREQ(earg[7], "c_temp[4]");
    EXPECT_STREQ(earg[8], "v_temp[1]");
    EXPECT_STREQ(earg[9], "v_temp[2]");
    EXPECT_STREQ(earg[10], "v_temp[3]");
    EXPECT_STREQ(earg[11], "v_temp[4]");
    EXPECT_STREQ(earg[12], "c_gofr[3*]");
    EXPECT_STREQ(earg[13], "c_gofr[1][*]");
    EXPECT_STREQ(earg[14], "c_gofr[*2][2]");
    EXPECT_STREQ(earg[15], "c_gofr[*][*]");

    for (int i = 0; i < narg; ++i)
        delete[] earg[i];
    lmp->memory->sfree(earg);

    narg = utils::expand_args(FLERR, oarg, args, 1, earg, lmp);
    EXPECT_EQ(narg, 16);
    EXPECT_NE(args, earg);
    EXPECT_STREQ(earg[0], "v_step");
    EXPECT_STREQ(earg[1], "c_temp");
    EXPECT_STREQ(earg[2], "f_1[*]");
    EXPECT_STREQ(earg[3], "c_temp[2*4]");
    EXPECT_STREQ(earg[4], "v_temp[*4]");
    EXPECT_STREQ(earg[5], "c_gofr[3]");
    EXPECT_STREQ(earg[6], "c_gofr[4]");
    EXPECT_STREQ(earg[7], "c_gofr[5]");
    EXPECT_STREQ(earg[8], "c_gofr[1][*]");
    EXPECT_STREQ(earg[9], "c_gofr[1][2]");
    EXPECT_STREQ(earg[10], "c_gofr[2][2]");
    EXPECT_STREQ(earg[11], "c_gofr[1][*]");
    EXPECT_STREQ(earg[12], "c_gofr[2][*]");
    EXPECT_STREQ(earg[13], "c_gofr[3][*]");
    EXPECT_STREQ(earg[14], "c_gofr[4][*]");
    EXPECT_STREQ(earg[15], "c_gofr[5][*]");

    for (int i = 0; i < narg; ++i)
        delete[] earg[i];
    lmp->memory->sfree(earg);

    args[3][9] = '9';
    TEST_FAILURE("ERROR: Numeric index 9 is out of bounds \\(1-6\\).*",
                 utils::expand_args(FLERR, oarg, args, 0, earg, lmp););

    args[3][9] = '4';
    args[5][7] = '9';
    TEST_FAILURE("ERROR: Numeric index 9 is out of bounds \\(1-5\\).*",
                 utils::expand_args(FLERR, oarg, args, 1, earg, lmp););

    args[5][7] = '3';
    delete[] args[4];
    args[4] = utils::strdup("v_temp[2*]");
    TEST_FAILURE("ERROR: Upper bound required to expand vector style variable temp.*",
                 utils::expand_args(FLERR, oarg, args, 0, earg, lmp););

    delete[] args[4];
    args[4] = utils::strdup("v_temp[*2]");
    narg    = utils::expand_args(FLERR, oarg, args, 0, earg, lmp);
    EXPECT_EQ(narg, 14);
    EXPECT_STREQ(earg[0], "v_step");
    EXPECT_STREQ(earg[1], "c_temp");
    EXPECT_STREQ(earg[2], "f_1[1]");
    EXPECT_STREQ(earg[3], "f_1[2]");
    EXPECT_STREQ(earg[4], "f_1[3]");
    EXPECT_STREQ(earg[5], "c_temp[2]");
    EXPECT_STREQ(earg[6], "c_temp[3]");
    EXPECT_STREQ(earg[7], "c_temp[4]");
    EXPECT_STREQ(earg[8], "v_temp[1]");
    EXPECT_STREQ(earg[9], "v_temp[2]");
    EXPECT_STREQ(earg[10], "c_gofr[3*]");
    EXPECT_STREQ(earg[11], "c_gofr[1][*]");
    EXPECT_STREQ(earg[12], "c_gofr[*2][2]");
    EXPECT_STREQ(earg[13], "c_gofr[*][*]");

    for (int i = 0; i < narg; ++i)
        delete[] earg[i];
    lmp->memory->sfree(earg);
    for (int i = 0; i < oarg; ++i)
        delete[] args[i];
    delete[] args;
}

TEST_F(AdvancedUtils, check_packages_for_style)
{
    auto mesg = utils::check_packages_for_style("pair", "unknown", lmp);
    EXPECT_THAT(mesg, StrEq("Unrecognized pair style 'unknown'"));

    // try styles from multiple packages that are rarely installed in
    // the hope that one of them triggers the error message about the missing package
    if (!Info::has_package("ADIOS")) {
        auto mesg = utils::check_packages_for_style("dump", "atom/adios", lmp);
        EXPECT_THAT(mesg, StrEq("Unrecognized dump style 'atom/adios' is part of the ADIOS "
                                "package which is not enabled in this LAMMPS binary.\n"
                                "For more information see https://docs.lammps.org/err0010"));
    }

    if (!Info::has_package("SCAFACOS")) {
        auto mesg = utils::check_packages_for_style("kspace", "scafacos", lmp);
        EXPECT_THAT(mesg, StrEq("Unrecognized kspace style 'scafacos' is part of the SCAFACOS "
                                "package which is not enabled in this LAMMPS binary.\n"
                                "For more information see https://docs.lammps.org/err0010"));
    }

    // try styles from multiple packages that are commonly installed in
    // the hope that one of them triggers the error message about the dependency
    if (Info::has_package("MANYBODY")) {
        auto mesg = utils::check_packages_for_style("pair", "tersoff", lmp);
        EXPECT_THAT(mesg, StrEq("Unrecognized pair style 'tersoff' is part of the MANYBODY "
                                "package, but seems to be missing because of a dependency"));
    }
    if (Info::has_package("MOLECULE")) {
        auto mesg = utils::check_packages_for_style("bond", "harmonic", lmp);
        EXPECT_THAT(mesg, StrEq("Unrecognized bond style 'harmonic' is part of the MOLECULE "
                                "package, but seems to be missing because of a dependency"));
    }
}

TEST_F(AdvancedUtils, logical)
{
    bool caught = false;
    EXPECT_EQ(utils::logical(FLERR, "yes", false, lmp), true);
    EXPECT_EQ(utils::logical(FLERR, "no", false, lmp), false);
    EXPECT_EQ(utils::logical(FLERR, "on", false, lmp), true);
    EXPECT_EQ(utils::logical(FLERR, "off", false, lmp), false);
    EXPECT_EQ(utils::logical(FLERR, "1", false, lmp), true);
    EXPECT_EQ(utils::logical(FLERR, "0", false, lmp), false);
    EXPECT_EQ(utils::logical(FLERR, "true", false, lmp), true);
    EXPECT_EQ(utils::logical(FLERR, "false", false, lmp), false);

    try {
        BEGIN_HIDE_OUTPUT();
        utils::logical(FLERR, "YES", true, lmp);
    } catch (LAMMPSAbortException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("ERROR on proc 0: Expected boolean parameter "
                                            "instead of 'YES' in input script or data file "));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
    }
    ASSERT_TRUE(caught);

    try {
        BEGIN_HIDE_OUTPUT();
        utils::logical(FLERR, "1.0e-2", false, lmp);
    } catch (LAMMPSException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("ERROR: Expected boolean parameter instead "
                                            "of '1.0e-2' in input script or data file "));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
    }
    ASSERT_TRUE(caught);
};

TEST_F(AdvancedUtils, numeric)
{
    bool caught = false;
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, "0.0", false, lmp), 0.0);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, "-1", false, lmp), -1.0);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, "2.1e-1", false, lmp), 0.21);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, "2.23e-308", false, lmp), 2.23e-308);
    EXPECT_DOUBLE_EQ(utils::numeric(FLERR, "1.79e308", false, lmp), 1.79e308);

    try {
        BEGIN_HIDE_OUTPUT();
        utils::numeric(FLERR, "text", true, lmp);
    } catch (LAMMPSAbortException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("ERROR on proc 0: Expected floating point parameter "
                                            "instead of 'text' in input script or data file "));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
    }
    ASSERT_TRUE(caught);

    try {
        BEGIN_HIDE_OUTPUT();
        utils::numeric(FLERR, "1.0d-2", false, lmp);
    } catch (LAMMPSException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("ERROR: Expected floating point parameter instead "
                                            "of '1.0d-2' in input script or data file "));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
    }
    ASSERT_TRUE(caught);
};

TEST_F(AdvancedUtils, inumeric)
{
    bool caught = false;
    EXPECT_EQ(utils::inumeric(FLERR, "0", false, lmp), 0);
    EXPECT_EQ(utils::inumeric(FLERR, "-1", false, lmp), -1);
    EXPECT_EQ(utils::inumeric(FLERR, "+11", false, lmp), 11);
    EXPECT_EQ(utils::inumeric(FLERR, "2147483647", false, lmp), 2147483647);
    EXPECT_EQ(utils::inumeric(FLERR, "-2147483648", false, lmp), -2147483648);

    try {
        BEGIN_HIDE_OUTPUT();
        utils::inumeric(FLERR, "2147483648", true, lmp);
    } catch (LAMMPSAbortException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("ERROR on proc 0: Integer 2147483648 in input "
                                            "script or data file is out of range"));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
    }
    ASSERT_TRUE(caught);

    try {
        BEGIN_HIDE_OUTPUT();
        utils::inumeric(FLERR, "-2147483649", false, lmp);
    } catch (LAMMPSException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("Integer -2147483649 in input script "
                                            "or data file is out of range"));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
    }
    ASSERT_TRUE(caught);
    try {
        BEGIN_HIDE_OUTPUT();
        utils::inumeric(FLERR, "--10", true, lmp);
    } catch (LAMMPSAbortException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("ERROR on proc 0: Expected integer parameter "
                                            "instead of '--10' in input script or data file "));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
    }
    ASSERT_TRUE(caught);

    try {
        BEGIN_HIDE_OUTPUT();
        utils::inumeric(FLERR, "1.0", false, lmp);
    } catch (LAMMPSException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("ERROR: Expected integer parameter instead "
                                            "of '1.0' in input script or data file "));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
    }
    ASSERT_TRUE(caught);
};

// bigint is *always* 64-bit

TEST_F(AdvancedUtils, bnumeric)
{
    bool caught = false;
    EXPECT_EQ(utils::bnumeric(FLERR, "0", false, lmp), 0);
    EXPECT_EQ(utils::bnumeric(FLERR, "-1", false, lmp), -1);
    EXPECT_EQ(utils::bnumeric(FLERR, "+11", false, lmp), 11);
    EXPECT_EQ(utils::bnumeric(FLERR, "2147483647000", false, lmp), 2147483647000);
    EXPECT_EQ(utils::bnumeric(FLERR, "-2147483648000", false, lmp), -2147483648000);

    try {
        BEGIN_HIDE_OUTPUT();
        utils::bnumeric(FLERR, "9223372036854775808", true, lmp);
    } catch (LAMMPSAbortException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("ERROR on proc 0: Integer 9223372036854775808 in "
                                            "input script or data file is out of range"));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
    }
    ASSERT_TRUE(caught);

    try {
        BEGIN_HIDE_OUTPUT();
        utils::bnumeric(FLERR, "-9223372036854775809", false, lmp);
    } catch (LAMMPSException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("Integer -9223372036854775809 in input script "
                                            "or data file is out of range"));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
    }
    ASSERT_TRUE(caught);
    try {
        BEGIN_HIDE_OUTPUT();
        utils::bnumeric(FLERR, "text", true, lmp);
    } catch (LAMMPSAbortException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("ERROR on proc 0: Expected integer parameter "
                                            "instead of 'text' in input script or data file "));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
    }
    ASSERT_TRUE(caught);

    try {
        BEGIN_HIDE_OUTPUT();
        utils::bnumeric(FLERR, "1.0", false, lmp);
    } catch (LAMMPSException &e) {
        END_HIDE_OUTPUT();
        EXPECT_THAT(e.what(), ContainsRegex("ERROR: Expected integer parameter instead "
                                            "of '1.0' in input script or data file "));
        caught = true;
    } catch (std::exception &e) {
        END_HIDE_OUTPUT();
        GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
    }
    ASSERT_TRUE(caught);
};

TEST_F(AdvancedUtils, tnumeric)
{
    if (sizeof(tagint) == 4) {

        bool caught = false;
        EXPECT_EQ(utils::tnumeric(FLERR, "0", false, lmp), 0);
        EXPECT_EQ(utils::tnumeric(FLERR, "-1", false, lmp), -1);
        EXPECT_EQ(utils::tnumeric(FLERR, "+11", false, lmp), 11);
        EXPECT_EQ(utils::tnumeric(FLERR, "2147483647", false, lmp), 2147483647);
        EXPECT_EQ(utils::tnumeric(FLERR, "-2147483648", false, lmp), -2147483648);

        try {
            BEGIN_HIDE_OUTPUT();
            utils::tnumeric(FLERR, "2147483648", true, lmp);
        } catch (LAMMPSAbortException &e) {
            END_HIDE_OUTPUT();
            EXPECT_THAT(e.what(), ContainsRegex("ERROR on proc 0: Integer 2147483648 in input "
                                                "script or data file is out of range"));
            caught = true;
        } catch (std::exception &e) {
            END_HIDE_OUTPUT();
            GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
        }
        ASSERT_TRUE(caught);

        try {
            BEGIN_HIDE_OUTPUT();
            utils::tnumeric(FLERR, "-2147483649", false, lmp);
        } catch (LAMMPSException &e) {
            END_HIDE_OUTPUT();
            EXPECT_THAT(e.what(), ContainsRegex("Integer -2147483649 in input script "
                                                "or data file is out of range"));
            caught = true;
        } catch (std::exception &e) {
            END_HIDE_OUTPUT();
            GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
        }
        ASSERT_TRUE(caught);
        try {
            BEGIN_HIDE_OUTPUT();
            utils::tnumeric(FLERR, "--10", true, lmp);
        } catch (LAMMPSAbortException &e) {
            END_HIDE_OUTPUT();
            EXPECT_THAT(e.what(), ContainsRegex("ERROR on proc 0: Expected integer parameter "
                                                "instead of '--10' in input script or data file "));
            caught = true;
        } catch (std::exception &e) {
            END_HIDE_OUTPUT();
            GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
        }
        ASSERT_TRUE(caught);

        try {
            BEGIN_HIDE_OUTPUT();
            utils::tnumeric(FLERR, "1.0", false, lmp);
        } catch (LAMMPSException &e) {
            END_HIDE_OUTPUT();
            EXPECT_THAT(e.what(), ContainsRegex("ERROR: Expected integer parameter instead "
                                                "of '1.0' in input script or data file "));
            caught = true;
        } catch (std::exception &e) {
            END_HIDE_OUTPUT();
            GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
        }
        ASSERT_TRUE(caught);
    } else if (sizeof(tagint) == 8) {
        bool caught = false;
        EXPECT_EQ(utils::tnumeric(FLERR, "0", false, lmp), 0);
        EXPECT_EQ(utils::tnumeric(FLERR, "-1", false, lmp), -1);
        EXPECT_EQ(utils::tnumeric(FLERR, "+11", false, lmp), 11);
        EXPECT_EQ(utils::tnumeric(FLERR, "2147483647000", false, lmp), 2147483647000);
        EXPECT_EQ(utils::tnumeric(FLERR, "-2147483648000", false, lmp), -2147483648000);

        try {
            BEGIN_HIDE_OUTPUT();
            utils::tnumeric(FLERR, "9223372036854775808", true, lmp);
        } catch (LAMMPSAbortException &e) {
            END_HIDE_OUTPUT();
            EXPECT_THAT(e.what(), ContainsRegex("ERROR on proc 0: Integer 9223372036854775808 in "
                                                "input script or data file is out of range"));
            caught = true;
        } catch (std::exception &e) {
            END_HIDE_OUTPUT();
            GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
        }
        ASSERT_TRUE(caught);

        try {
            BEGIN_HIDE_OUTPUT();
            utils::tnumeric(FLERR, "-9223372036854775809", false, lmp);
        } catch (LAMMPSException &e) {
            END_HIDE_OUTPUT();
            EXPECT_THAT(e.what(), ContainsRegex("Integer -9223372036854775809 in input script "
                                                "or data file is out of range"));
            caught = true;
        } catch (std::exception &e) {
            END_HIDE_OUTPUT();
            GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
        }
        ASSERT_TRUE(caught);
        try {
            BEGIN_HIDE_OUTPUT();
            utils::tnumeric(FLERR, "text", true, lmp);
        } catch (LAMMPSAbortException &e) {
            END_HIDE_OUTPUT();
            EXPECT_THAT(e.what(), ContainsRegex("ERROR on proc 0: Expected integer parameter "
                                                "instead of 'text' in input script or data file "));
            caught = true;
        } catch (std::exception &e) {
            END_HIDE_OUTPUT();
            GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
        }
        ASSERT_TRUE(caught);

        try {
            BEGIN_HIDE_OUTPUT();
            utils::tnumeric(FLERR, "1.0", false, lmp);
        } catch (LAMMPSException &e) {
            END_HIDE_OUTPUT();
            EXPECT_THAT(e.what(), ContainsRegex("ERROR: Expected integer parameter instead "
                                                "of '1.0' in input script or data file "));
            caught = true;
        } catch (std::exception &e) {
            END_HIDE_OUTPUT();
            GTEST_FAIL() << "Incorrect exception: " << e.what() << "\n";
        }
        ASSERT_TRUE(caught);
    }
};

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
