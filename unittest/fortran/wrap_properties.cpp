// unit tests for getting LAMMPS properties through the Fortran wrapper

#include "lammps.h"
//#include <cstdio> // for stdin, stdout
#include "library.h"
#include <mpi.h>
#include <string>

#include "gtest/gtest.h"

// prototypes for fortran reverse wrapper functions
extern "C" {
void *f_lammps_with_args();
void f_lammps_close();
int f_lammps_version();
void f_lammps_memory_usage(double*);
int f_lammps_get_mpi_comm();
int f_lammps_extract_setting(const char*);
}

class LAMMPS_properties : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_properties()           = default;
    ~LAMMPS_properties() override = default;

    void SetUp() override
    {
        ::testing::internal::CaptureStdout();
        lmp                = (LAMMPS_NS::LAMMPS *)f_lammps_with_args();
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
   if ( lammps_config_has_mpi_support() )
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
