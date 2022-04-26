// unit tests for utils:: functions requiring a LAMMPS

#include "error.h"
#include "input.h"
#include "lammps.h"
#include "utils.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

class Advanced_utils : public LAMMPSTest {
protected:
    Error *error;

    void SetUp() override
    {
        testbinary = "AdvancedUtils";
        LAMMPSTest::SetUp();
        error = lmp->error;
    }
};

TEST_F(Advanced_utils, missing_cmd_args)
{
    auto output = CAPTURE_OUTPUT([&] {
        utils::missing_cmd_args(FLERR, "dummy", nullptr);
    });
    EXPECT_EQ(output, "");

    TEST_FAILURE("ERROR: Illegal dummy command: missing argument",
                 utils::missing_cmd_args(FLERR, "dummy", error););
};

TEST_F(Advanced_utils, logmesg)
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

TEST_F(Advanced_utils, bounds_case1)
{
    int nlo, nhi;

    nlo = nhi = -1;
    utils::bounds(FLERR, "9", 0, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 9);
    ASSERT_EQ(nhi, 9);
    utils::bounds(FLERR, "1", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 1);
    ASSERT_EQ(nhi, 1);
    utils::bounds(FLERR, "1x", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
    utils::bounds(FLERR, "-1", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
    utils::bounds(FLERR, "+1", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
    utils::bounds(FLERR, "1:3", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
}

TEST_F(Advanced_utils, bounds_case2)
{
    int nlo, nhi;

    nlo = nhi = -1;
    utils::bounds(FLERR, "*", 0, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 0);
    ASSERT_EQ(nhi, 10);
    utils::bounds(FLERR, "*", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -10);
    ASSERT_EQ(nhi, 5);
    utils::bounds(FLERR, "?", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
}

TEST_F(Advanced_utils, bounds_case3)
{
    int nlo, nhi;

    nlo = nhi = -1;
    utils::bounds(FLERR, "2*", 0, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 2);
    ASSERT_EQ(nhi, 10);
    utils::bounds(FLERR, "3*", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 3);
    ASSERT_EQ(nhi, 5);
    utils::bounds(FLERR, "3*:2", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
}

TEST_F(Advanced_utils, boundsbig_case1)
{
    bigint nlo, nhi;

    nlo = nhi = -1;
    utils::bounds(FLERR, "9", 0, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 9);
    ASSERT_EQ(nhi, 9);
    utils::bounds(FLERR, "1", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 1);
    ASSERT_EQ(nhi, 1);
    utils::bounds(FLERR, "1x", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
    utils::bounds(FLERR, "-1", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
    utils::bounds(FLERR, "+1", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
    utils::bounds(FLERR, "1:3", 1, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
}

TEST_F(Advanced_utils, boundsbig_case2)
{
    bigint nlo, nhi;

    nlo = nhi = -1;
    utils::bounds(FLERR, "*", 0, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 0);
    ASSERT_EQ(nhi, 10);
    utils::bounds(FLERR, "*", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -10);
    ASSERT_EQ(nhi, 5);
    utils::bounds(FLERR, "?", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
}

TEST_F(Advanced_utils, boundsbig_case3)
{
    bigint nlo, nhi;

    nlo = nhi = -1;
    utils::bounds(FLERR, "2*", 0, 10, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 2);
    ASSERT_EQ(nhi, 10);
    utils::bounds(FLERR, "3*", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, 3);
    ASSERT_EQ(nhi, 5);
    utils::bounds(FLERR, "3*:2", -10, 5, nlo, nhi, nullptr);
    ASSERT_EQ(nlo, -1);
    ASSERT_EQ(nhi, -1);
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (platform::mpi_vendor() == "Open MPI" && !LAMMPS_NS::Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. Death tests will be skipped\n";

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = utils::split_words(var);
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
