
// unit tests for gathering and scattering data from a LAMMPS instance through
// the Fortran wrapper

#include "lammps.h"
#include "library.h"
#include <cstdint>
#include <cstdlib>
#include <mpi.h>
#include <string>

#include "gtest/gtest.h"

// prototypes for Fortran reverse wrapper functions
extern "C" {
void *f_lammps_with_args();
void f_lammps_close();
void f_lammps_setup_fix_external();
void f_lammps_set_fix_external_callbacks();
void f_lammps_get_force(int, double*);
void f_lammps_reverse_direction();
}

using namespace LAMMPS_NS;

class LAMMPS_fixexternal : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_fixexternal()           = default;
    ~LAMMPS_fixexternal() override = default;

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

TEST_F(LAMMPS_fixexternal, callback)
{
    f_lammps_setup_fix_external();
    f_lammps_set_fix_external_callbacks();
    lammps_command(lmp, "run 0");
    double f[3];
    f_lammps_get_force(1,f);
    EXPECT_DOUBLE_EQ(f[0], 3.0);
    EXPECT_DOUBLE_EQ(f[1], -3.0);
    EXPECT_DOUBLE_EQ(f[2], 3.75);
    f_lammps_get_force(2,f);
    EXPECT_DOUBLE_EQ(f[0], -3.0);
    EXPECT_DOUBLE_EQ(f[1], 3.0);
    EXPECT_DOUBLE_EQ(f[2], -3.75);

    f_lammps_reverse_direction();
    f_lammps_set_fix_external_callbacks();
    lammps_command(lmp, "run 0");
    f_lammps_get_force(1,f);
    EXPECT_DOUBLE_EQ(f[0], -1.0);
    EXPECT_DOUBLE_EQ(f[1], 1.0);
    EXPECT_DOUBLE_EQ(f[2], -1.25);
    f_lammps_get_force(2,f);
    EXPECT_DOUBLE_EQ(f[0], 1.0);
    EXPECT_DOUBLE_EQ(f[1], -1.0);
    EXPECT_DOUBLE_EQ(f[2], 1.25);
};
