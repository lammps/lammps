// unit tests for checking LAMMPS configuration settings  through the library interface

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
using ::testing::StrEq;

class LibraryConfig : public ::testing::Test {
protected:
    void *lmp;
    std::string INPUT_DIR = STRINGIFY(TEST_INPUT_FOLDER);

    LibraryConfig(){};
    ~LibraryConfig() override{};

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

TEST(LAMMPSConfig, package_count)
{
    EXPECT_EQ(lammps_config_package_count(), NUM_LAMMPS_PACKAGES);
};

TEST(LAMMPSConfig, has_package)
{
    EXPECT_EQ(lammps_config_has_package("MANYBODY"), LAMMPS_HAS_MANYBODY);
};

TEST(LAMMPSConfig, package_name)
{
    char buf[128];
    int numpkgs = lammps_config_package_count();
    if (numpkgs > 0) {
        EXPECT_EQ(lammps_config_package_name(0, buf, 128), 1);
        EXPECT_EQ(lammps_config_package_name(numpkgs + 10, buf, 128), 0);
        EXPECT_THAT(buf, StrEq(""));
    } else {
        EXPECT_EQ(lammps_config_package_name(0, buf, 128), 1);
        EXPECT_THAT(buf, StrEq(""));
    }
};

TEST_F(LibraryConfig, has_style)
{
    EXPECT_EQ(lammps_has_style(lmp, "atom", "atomic"), 1);
    EXPECT_EQ(lammps_has_style(lmp, "compute", "temp"), 1);
    EXPECT_EQ(lammps_has_style(lmp, "fix", "nve"), 1);
    EXPECT_EQ(lammps_has_style(lmp, "pair", "lj/cut"), 1);
    EXPECT_EQ(lammps_has_style(lmp, "bond", "zero"), 1);
    EXPECT_EQ(lammps_has_style(lmp, "dump", "custom"), 1);
    EXPECT_EQ(lammps_has_style(lmp, "region", "sphere"), 1);
    EXPECT_EQ(lammps_has_style(lmp, "xxxxx", "lj/cut"), 0);
    EXPECT_EQ(lammps_has_style(lmp, "pair", "xxxxxx"), 0);
#if LAMMPS_HAS_MANYBODY
    EXPECT_EQ(lammps_has_style(lmp, "pair", "sw"), 1);
#else
    EXPECT_EQ(lammps_has_style(lmp, "pair", "sw"), 0);
#endif
};

TEST_F(LibraryConfig, style_count)
{
    EXPECT_GT(lammps_style_count(lmp, "atom"), 1);
    EXPECT_GT(lammps_style_count(lmp, "bond"), 1);
    EXPECT_GT(lammps_style_count(lmp, "angle"), 1);
    EXPECT_GT(lammps_style_count(lmp, "dihedral"), 1);
    EXPECT_GT(lammps_style_count(lmp, "improper"), 1);
    EXPECT_GT(lammps_style_count(lmp, "pair"), 1);
    EXPECT_GT(lammps_style_count(lmp, "kspace"), 1);
    EXPECT_GT(lammps_style_count(lmp, "compute"), 1);
    EXPECT_GT(lammps_style_count(lmp, "fix"), 1);
    EXPECT_GT(lammps_style_count(lmp, "region"), 1);
    EXPECT_GT(lammps_style_count(lmp, "dump"), 1);
    EXPECT_GT(lammps_style_count(lmp, "integrate"), 1);
    EXPECT_GT(lammps_style_count(lmp, "minimize"), 1);
};

TEST_F(LibraryConfig, style_name)
{
    char buf[128];
    int numstyles = lammps_style_count(lmp, "atom");
    EXPECT_EQ(lammps_style_name(lmp, "atom", 0, buf, 128), 1);
    EXPECT_EQ(lammps_style_name(lmp, "atom", numstyles + 10, buf, 128), 0);
    EXPECT_THAT(buf, StrEq(""));
};

TEST(LAMMPSConfig, exceptions)
{
    EXPECT_EQ(lammps_config_has_exceptions(), LAMMPS_HAS_EXCEPTIONS);
};

TEST(LAMMPSConfig, mpi_support)
{
    EXPECT_EQ(lammps_config_has_mpi_support(), LAMMPS_HAS_MPI);
};

TEST(LAMMPSConfig, png_support)
{
    EXPECT_EQ(lammps_config_has_png_support(), LAMMPS_HAS_PNG);
};

TEST(LAMMPSConfig, jpeg_support)
{
    EXPECT_EQ(lammps_config_has_jpeg_support(), LAMMPS_HAS_JPEG);
};

TEST(LAMMPSConfig, gzip_support)
{
    EXPECT_EQ(lammps_config_has_gzip_support(), LAMMPS_HAS_GZIP);
};

TEST(LAMMPSConfig, ffmpeg_support)
{
    EXPECT_EQ(lammps_config_has_ffmpeg_support(), LAMMPS_HAS_FFMPEG);
};
