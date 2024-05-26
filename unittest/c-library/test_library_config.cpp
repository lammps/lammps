// unit tests for checking LAMMPS configuration settings  through the library interface

#include "lammps.h"
#include "library.h"
#include "timer.h"
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

    LibraryConfig()           = default;
    ~LibraryConfig() override = default;

    void SetUp() override
    {
        const char *args[] = {"LAMMPS_test", "-log",      "none",
                              "-echo",       "screen",    "-nocite",
                              "-var",        "input_dir", STRINGIFY(TEST_INPUT_FOLDER),
                              nullptr};

        char **argv = (char **)args;
        int argc    = (sizeof(args) / sizeof(char *)) - 1;

        ::testing::internal::CaptureStdout();
        lmp = lammps_open_no_mpi(argc, argv, nullptr);
        lammps_command(lmp, "fix charge all property/atom q ghost yes");
        lammps_command(lmp, "region box block 0 1 0 1 0 1");
        lammps_command(lmp, "create_box 1 box");
        lammps_command(lmp, "group none empty");
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
        EXPECT_EQ(lammps_config_package_name(0, buf, 128), 0);
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
    EXPECT_GE(lammps_style_count(lmp, "kspace"), 0);
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

TEST_F(LibraryConfig, has_id)
{
    EXPECT_EQ(lammps_has_id(lmp, "compute", "thermo_temp"), 1);
    EXPECT_EQ(lammps_has_id(lmp, "compute", "thermo_press"), 1);
    EXPECT_EQ(lammps_has_id(lmp, "compute", "thermo_pe"), 1);
    EXPECT_EQ(lammps_has_id(lmp, "dump", "xxx"), 0);
    EXPECT_EQ(lammps_has_id(lmp, "fix", "charge"), 1);
    EXPECT_EQ(lammps_has_id(lmp, "fix", "xxx"), 0);
    EXPECT_EQ(lammps_has_id(lmp, "group", "all"), 1);
    EXPECT_EQ(lammps_has_id(lmp, "group", "none"), 1);
    EXPECT_EQ(lammps_has_id(lmp, "group", "xxx"), 0);
    EXPECT_EQ(lammps_has_id(lmp, "molecule", "xxx"), 0);
    EXPECT_EQ(lammps_has_id(lmp, "region", "box"), 1);
    EXPECT_EQ(lammps_has_id(lmp, "region", "xxx"), 0);
    EXPECT_EQ(lammps_has_id(lmp, "variable", "input_dir"), 1);
    EXPECT_EQ(lammps_has_id(lmp, "variable", "xxx"), 0);
};

TEST_F(LibraryConfig, id_count)
{
    EXPECT_EQ(lammps_id_count(lmp, "compute"), 3);
    EXPECT_EQ(lammps_id_count(lmp, "dump"), 0);
    EXPECT_EQ(lammps_id_count(lmp, "fix"), 1);
    EXPECT_EQ(lammps_id_count(lmp, "group"), 2);
    EXPECT_EQ(lammps_id_count(lmp, "molecule"), 0);
    EXPECT_EQ(lammps_id_count(lmp, "region"), 1);
    EXPECT_EQ(lammps_id_count(lmp, "variable"), 1);
};

TEST_F(LibraryConfig, id_name)
{
    const int bufsize = 128;
    char buf[bufsize];
    EXPECT_EQ(lammps_id_name(lmp, "compute", 2, buf, bufsize), 1);
    EXPECT_THAT(buf, StrEq("thermo_pe"));
    EXPECT_EQ(lammps_id_name(lmp, "compute", 10, buf, bufsize), 0);
    EXPECT_THAT(buf, StrEq(""));
    EXPECT_EQ(lammps_id_name(lmp, "fix", 0, buf, bufsize), 1);
    EXPECT_THAT(buf, StrEq("charge"));
    EXPECT_EQ(lammps_id_name(lmp, "fix", 10, buf, bufsize), 0);
    EXPECT_THAT(buf, StrEq(""));
    EXPECT_EQ(lammps_id_name(lmp, "group", 0, buf, bufsize), 1);
    EXPECT_THAT(buf, StrEq("all"));
    EXPECT_EQ(lammps_id_name(lmp, "group", 1, buf, bufsize), 1);
    EXPECT_THAT(buf, StrEq("none"));
    EXPECT_EQ(lammps_id_name(lmp, "group", 10, buf, bufsize), 0);
    EXPECT_THAT(buf, StrEq(""));
    EXPECT_EQ(lammps_id_name(lmp, "region", 0, buf, bufsize), 1);
    EXPECT_THAT(buf, StrEq("box"));
    EXPECT_EQ(lammps_id_name(lmp, "region", 10, buf, bufsize), 0);
    EXPECT_THAT(buf, StrEq(""));
    EXPECT_EQ(lammps_id_name(lmp, "variable", 0, buf, bufsize), 1);
    EXPECT_THAT(buf, StrEq("input_dir"));
    EXPECT_EQ(lammps_id_name(lmp, "variable", 10, buf, bufsize), 0);
    EXPECT_THAT(buf, StrEq(""));
};

TEST_F(LibraryConfig, is_running)
{
    EXPECT_EQ(lammps_is_running(lmp), 0);
}

TEST_F(LibraryConfig, force_timeout)
{
    LAMMPS_NS::Timer *timer = ((LAMMPS_NS::LAMMPS *)lmp)->timer;
    EXPECT_EQ(timer->is_timeout(), false);
    lammps_force_timeout(lmp);
    EXPECT_EQ(timer->is_timeout(), true);
}

TEST(LAMMPSConfig, exceptions)
{
    EXPECT_EQ(lammps_config_has_exceptions(), 1);
};

TEST(LAMMPSConfig, mpi_support)
{
    if (LAMMPS_HAS_MPI)
        EXPECT_GT(lammps_config_has_mpi_support(), 0);
    else
        EXPECT_EQ(lammps_config_has_mpi_support(), 0);
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
