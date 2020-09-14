// unit tests for checking LAMMPS configuration settings  through the library interface

#include "library.h"
#include "lammps.h"
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "test_main.h"

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

using ::testing::HasSubstr;
using ::testing::StartsWith;

class LAMMPS_config : public ::testing::Test {
protected:
    void *lmp;
    std::string INPUT_DIR = STRINGIFY(TEST_INPUT_FOLDER);

    LAMMPS_config(){};
    ~LAMMPS_config() override{};

    void SetUp() override
    {
        const char *args[] = {"LAMMPS_test", "-log",      "none",
                              "-echo",       "screen",    "-nocite",
                              "-var",        "input_dir", STRINGIFY(TEST_INPUT_FOLDER)};

        char **argv = (char **)args;
        int argc    = sizeof(args) / sizeof(char *);

        ::testing::internal::CaptureStdout();
        lmp                = lammps_open_no_mpi(argc, argv, NULL);
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        EXPECT_THAT(output, StartsWith("LAMMPS ("));
    }
    void TearDown() override
    {
        ::testing::internal::CaptureStdout();
        lammps_close(lmp);
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, HasSubstr("Total wall time:"));
        if (verbose) std::cout << output;
        lmp = nullptr;
    }
};

TEST(LAMMPS_config, package_count)
{
    EXPECT_EQ(lammps_config_package_count(), NUM_LAMMPS_PACKAGES);
};

TEST(LAMMPS_config, has_package)
{
    EXPECT_EQ(lammps_config_has_package("MANYBODY"), LAMMPS_HAS_MANYBODY);
};

