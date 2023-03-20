// unit tests for extracting box dimensions fom a LAMMPS instance through the Fortran wrapper

#include "lammps.h"
#include <mpi.h>
#include <string>

#include "gtest/gtest.h"

// prototypes for fortran reverse wrapper functions
extern "C" {
void *f_lammps_with_args();
void f_lammps_close();
void f_lammps_box_setup();
double f_lammps_extract_box_xlo();
double f_lammps_extract_box_xhi();
double f_lammps_extract_box_ylo();
double f_lammps_extract_box_yhi();
double f_lammps_extract_box_zlo();
double f_lammps_extract_box_zhi();
void f_lammps_delete_everything();
void f_lammps_reset_box_2x();
}

class LAMMPS_commands : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_commands()           = default;
    ~LAMMPS_commands() override = default;

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

TEST_F(LAMMPS_commands, get_thermo)
{
    f_lammps_box_setup();
    EXPECT_DOUBLE_EQ(f_lammps_extract_box_xlo(), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_box_xhi(), 2.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_box_ylo(), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_box_yhi(), 2.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_box_zlo(), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_box_zhi(), 2.0);
    f_lammps_delete_everything();
    f_lammps_reset_box_2x();
    EXPECT_DOUBLE_EQ(f_lammps_extract_box_xlo(), -1.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_box_xhi(), 3.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_box_ylo(), -1.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_box_yhi(), 3.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_box_zlo(), -1.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_box_zhi(), 3.0);
};
