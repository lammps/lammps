// unit tests for getting LAMMPS properties through the Fortran wrapper

#include "info.h"
#include "lammps.h"
#include "library.h"

#include <cstdint>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

// prototypes for fortran reverse wrapper functions
extern "C" {
void *f_lammps_with_args();
void f_lammps_close();
void f_lammps_memory_usage(double *);
int f_lammps_get_mpi_comm();
int f_lammps_extract_setting(const char *);
int f_lammps_has_error();
int f_lammps_get_last_error_message(char *, int);
int f_lammps_get_image_flags_int(int, int, int);
int64_t f_lammps_get_image_flags_bigint(int, int, int);
void f_lammps_decode_image_flags(int, int *);
void f_lammps_decode_image_flags_bigbig(int64_t, int *);
}

namespace LAMMPS_NS {

using ::testing::ContainsRegex;

class LAMMPS_properties : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        ::testing::internal::CaptureStdout();
        lmp = (LAMMPS *)f_lammps_with_args();

        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0, 8).c_str(), "LAMMPS (");
    }

    void TearDown() override
    {
        ::testing::internal::CaptureStdout();
        f_lammps_close();
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0, 16).c_str(), "Total wall time:");
        lmp = nullptr;
    }
};

TEST_F(LAMMPS_properties, memory_usage)
{
    // copied from c-library, with a two-character modification
    double meminfo[3];
    f_lammps_memory_usage(meminfo);
    EXPECT_GT(meminfo[0], 0.0);
#if defined(__linux__) || defined(_WIN32)
    EXPECT_GE(meminfo[1], 0.0);
#endif
#if (defined(__linux__) || defined(__APPLE__) || defined(_WIN32)) && !defined(__INTEL_LLVM_COMPILER)
    EXPECT_GT(meminfo[2], 0.0);
#endif
};

TEST_F(LAMMPS_properties, get_mpi_comm)
{
    int f_comm = f_lammps_get_mpi_comm();
    if (lammps_config_has_mpi_support())
        EXPECT_GE(f_comm, 0);
    else
        EXPECT_EQ(f_comm, -1);
};

TEST_F(LAMMPS_properties, extract_setting)
{
#if defined(LAMMPS_SMALLSMALL)
    EXPECT_EQ(f_lammps_extract_setting("bigint"), 4);
#else
    EXPECT_EQ(f_lammps_extract_setting("bigint"), 8);
#endif
#if defined(LAMMPS_BIGBIG)
    EXPECT_EQ(f_lammps_extract_setting("tagint"), 8);
    EXPECT_EQ(f_lammps_extract_setting("imageint"), 8);
    EXPECT_EQ(f_lammps_extract_setting("IMGMASK"), 2097151);
    EXPECT_EQ(f_lammps_extract_setting("IMGMAX"), 1048576);
    EXPECT_EQ(f_lammps_extract_setting("IMGBITS"), 21);
    EXPECT_EQ(f_lammps_extract_setting("IMG2BITS"), 42);
#else
    EXPECT_EQ(f_lammps_extract_setting("tagint"), 4);
    EXPECT_EQ(f_lammps_extract_setting("imageint"), 4);
    EXPECT_EQ(f_lammps_extract_setting("IMGMASK"), 1023);
    EXPECT_EQ(f_lammps_extract_setting("IMGMAX"), 512);
    EXPECT_EQ(f_lammps_extract_setting("IMGBITS"), 10);
    EXPECT_EQ(f_lammps_extract_setting("IMG2BITS"), 20);
#endif

    EXPECT_EQ(f_lammps_extract_setting("box_exist"), 0);
    EXPECT_EQ(f_lammps_extract_setting("dimension"), 3);
    EXPECT_EQ(f_lammps_extract_setting("world_size"), 1);
    EXPECT_EQ(f_lammps_extract_setting("world_rank"), 0);
    EXPECT_EQ(f_lammps_extract_setting("universe_size"), 1);
    EXPECT_EQ(f_lammps_extract_setting("universe_rank"), 0);
    EXPECT_GT(f_lammps_extract_setting("nthreads"), 0);
    EXPECT_EQ(f_lammps_extract_setting("newton_pair"), 1);
    EXPECT_EQ(f_lammps_extract_setting("newton_bond"), 1);

    EXPECT_EQ(f_lammps_extract_setting("ntypes"), 0);
    EXPECT_EQ(f_lammps_extract_setting("nbondtypes"), 0);
    EXPECT_EQ(f_lammps_extract_setting("nangletypes"), 0);
    EXPECT_EQ(f_lammps_extract_setting("ndihedraltypes"), 0);
    EXPECT_EQ(f_lammps_extract_setting("nimpropertypes"), 0);

    EXPECT_EQ(f_lammps_extract_setting("molecule_flag"), 0);
    EXPECT_EQ(f_lammps_extract_setting("q_flag"), 0);
    EXPECT_EQ(f_lammps_extract_setting("mu_flag"), 0);
    EXPECT_EQ(f_lammps_extract_setting("rmass_flag"), 0);
    EXPECT_EQ(f_lammps_extract_setting("UNKNOWN"), -1);
};

TEST_F(LAMMPS_properties, has_error)
{
    EXPECT_EQ(f_lammps_has_error(), lammps_has_error(lmp));
    EXPECT_EQ(f_lammps_has_error(), 0);

    // trigger an error, but hide output
    ::testing::internal::CaptureStdout();
    lammps_command(lmp, "this_is_not_a_known_command");
    ::testing::internal::GetCapturedStdout();

    EXPECT_EQ(f_lammps_has_error(), lammps_has_error(lmp));
    EXPECT_EQ(f_lammps_has_error(), 1);

    // retrieve error message
    char errmsg[1024];
    int err = f_lammps_get_last_error_message(errmsg, 1023);
    EXPECT_EQ(err, 1);
    EXPECT_THAT(errmsg, ContainsRegex(".*ERROR: Unknown command: this_is_not_a_known_command.*"));

    // retrieving the error message clear the error status
    EXPECT_EQ(f_lammps_has_error(), 0);
    err = f_lammps_get_last_error_message(errmsg, 1023);
    EXPECT_EQ(err, 0);
    EXPECT_THAT(errmsg, ContainsRegex("                                                  "));
};

TEST_F(LAMMPS_properties, get_image_flags)
{
#ifdef LAMMPS_BIGBIG
    int64_t image  = f_lammps_get_image_flags_bigint(0, 0, 0);
    int64_t Cimage = lammps_encode_image_flags(0, 0, 0);
    EXPECT_EQ(image, Cimage);
    image  = f_lammps_get_image_flags_bigint(1, -1, 1);
    Cimage = lammps_encode_image_flags(1, -1, 1);
    EXPECT_EQ(image, Cimage);
#else
    int image  = f_lammps_get_image_flags_int(0, 0, 0);
    int Cimage = lammps_encode_image_flags(0, 0, 0);
    EXPECT_EQ(image, Cimage);
    image  = f_lammps_get_image_flags_int(1, -1, 1);
    Cimage = lammps_encode_image_flags(1, -1, 1);
    EXPECT_EQ(image, Cimage);
#endif
}

TEST_F(LAMMPS_properties, decode_image_flags)
{
    int flag[3];
#ifdef LAMMPS_BIGBIG
    int64_t image = f_lammps_get_image_flags_bigint(1, 3, -2);
    f_lammps_decode_image_flags_bigbig(image, flag);
#else
    int image = f_lammps_get_image_flags_int(1, 3, -2);
    f_lammps_decode_image_flags(image, flag);
#endif
    EXPECT_EQ(flag[0], 1);
    EXPECT_EQ(flag[1], 3);
    EXPECT_EQ(flag[2], -2);
};
} // namespace LAMMPS_NS
