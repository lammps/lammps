// unit tests for getting LAMMPS properties through the Fortran wrapper

#include "lammps.h"
#include "library.h"

#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

// prototypes for fortran reverse wrapper functions
extern "C" {
void *f_lammps_with_args();
void f_lammps_close();
int f_lammps_version();
void f_lammps_memory_usage(double *);
int f_lammps_get_mpi_comm();
int f_lammps_extract_setting(const char *);
int f_lammps_has_error();
int f_lammps_get_last_error_message(char *, int);
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

TEST_F(LAMMPS_properties, version)
{
    EXPECT_LT(20200917, f_lammps_version());
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
#else
    EXPECT_EQ(f_lammps_extract_setting("tagint"), 4);
    EXPECT_EQ(f_lammps_extract_setting("imageint"), 4);
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
    // need errors to throw exceptions to be able to intercept them.
    if (!lammps_config_has_exceptions()) GTEST_SKIP();

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
    EXPECT_THAT(errmsg, ContainsRegex(".*ERRORx: Unknown command: this_is_not_a_known_command.*"));
};
} // namespace LAMMPS_NS
