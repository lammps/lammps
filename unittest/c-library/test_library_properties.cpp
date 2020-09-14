// unit tests for checking and changing simulation properties through the library interface

#include "lammps.h"
#include "library.h"
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "test_main.h"

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

using ::testing::HasSubstr;
using ::testing::StartsWith;

class LibraryProperties : public ::testing::Test {
protected:
    void *lmp;
    std::string INPUT_DIR = STRINGIFY(TEST_INPUT_FOLDER);

    LibraryProperties(){};
    ~LibraryProperties() override{};

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

TEST_F(LibraryProperties, version)
{
    EXPECT_GE(20200824, lammps_version(lmp));
};

TEST_F(LibraryProperties, get_mpi_comm)
{
    int f_comm = lammps_get_mpi_comm(lmp);
    if (lammps_config_has_mpi_support())
        EXPECT_GE(f_comm, 0);
    else
        EXPECT_EQ(f_comm, -1);
};

TEST_F(LibraryProperties, natoms)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = INPUT_DIR + PATH_SEP + "in.fourmol";
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_file(lmp, input.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_get_natoms(lmp), 29);
};

TEST_F(LibraryProperties, thermo)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = INPUT_DIR + PATH_SEP + "in.fourmol";
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_file(lmp, input.c_str());
    lammps_command(lmp, "run 2 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_get_thermo(lmp, "step"), 2);
    EXPECT_EQ(lammps_get_thermo(lmp, "atoms"), 29);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "vol"), 3375.0);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "density"), 0.12211250945013695);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "cellalpha"), 90.0);
};

TEST_F(LibraryProperties, box)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = INPUT_DIR + PATH_SEP + "in.fourmol";
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_file(lmp, input.c_str());
    lammps_command(lmp, "run 2 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    double boxlo[3], boxhi[3], xy, yz, xz;
    int pflags[3], boxflag;
    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, pflags, &boxflag);
    EXPECT_DOUBLE_EQ(boxlo[0], -6.024572);
    EXPECT_DOUBLE_EQ(boxlo[1], -7.692866);
    EXPECT_DOUBLE_EQ(boxlo[2], -8.086924);
    EXPECT_DOUBLE_EQ(boxhi[0], 8.975428);
    EXPECT_DOUBLE_EQ(boxhi[1], 7.307134);
    EXPECT_DOUBLE_EQ(boxhi[2], 6.913076);
    EXPECT_DOUBLE_EQ(xy, 0.0);
    EXPECT_DOUBLE_EQ(yz, 0.0);
    EXPECT_DOUBLE_EQ(xz, 0.0);
    EXPECT_EQ(pflags[0], 1);
    EXPECT_EQ(pflags[1], 1);
    EXPECT_EQ(pflags[2], 1);
    EXPECT_EQ(boxflag, 0);

    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "vol"), 3375.0);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "density"), 0.12211250945013695);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "cellalpha"), 90.0);
};
