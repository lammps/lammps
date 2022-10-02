// unit tests for gathering and scattering data from a LAMMPS instance through
// the Fortran wrapper

#include "lammps.h"
#include "library.h"
#include <mpi.h>
#include <string>
#include <cstdlib>
#include <cstdint>

#include "gtest/gtest.h"

// prototypes for Fortran reverse wrapper functions
extern "C" {
void *f_lammps_with_args();
void f_lammps_close();
void f_lammps_setup_gather_scatter();
int f_lammps_gather_mask(int);
double f_lammps_gather_position(int);
int f_lammps_gather_concat_mask(int);
double f_lammps_gather_concat_position(int,int);
int f_lammps_gather_subset_mask(int);
double f_lammps_gather_subset_position(int,int);
}

class LAMMPS_gather_scatter : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_gather_scatter()           = default;
    ~LAMMPS_gather_scatter() override = default;

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

TEST_F(LAMMPS_gather_scatter, gather_masks)
{
   f_lammps_setup_gather_scatter();
   EXPECT_EQ(f_lammps_gather_mask(1), 1);
   EXPECT_EQ(f_lammps_gather_mask(2), 1);
   EXPECT_EQ(f_lammps_gather_mask(3), 1);
   lammps_command(lmp, "group special id 1");
   lammps_command(lmp, "group other id 2");
   lammps_command(lmp, "group spiffy id 3");
   EXPECT_EQ(f_lammps_gather_mask(1), 3);
   EXPECT_EQ(f_lammps_gather_mask(2), 5);
   EXPECT_EQ(f_lammps_gather_mask(3), 9);
   lammps_command(lmp, "group other id 1");
   EXPECT_EQ(f_lammps_gather_mask(1), 7);
   EXPECT_EQ(f_lammps_gather_mask(2), 5);
   EXPECT_EQ(f_lammps_gather_mask(3), 9);
};

TEST_F(LAMMPS_gather_scatter, gather_positions)
{
   f_lammps_setup_gather_scatter();
   EXPECT_EQ(f_lammps_gather_position(1), 1.0);
   EXPECT_EQ(f_lammps_gather_position(2), 1.0);
   EXPECT_EQ(f_lammps_gather_position(3), 1.5);
   EXPECT_EQ(f_lammps_gather_position(4), 0.2);
   EXPECT_EQ(f_lammps_gather_position(5), 0.1);
   EXPECT_EQ(f_lammps_gather_position(6), 0.1);
   EXPECT_EQ(f_lammps_gather_position(7), 0.5);
   EXPECT_EQ(f_lammps_gather_position(8), 0.5);
   EXPECT_EQ(f_lammps_gather_position(9), 0.5);   
};

TEST_F(LAMMPS_gather_scatter, gather_masks_concat)
{
   f_lammps_setup_gather_scatter();
   EXPECT_EQ(f_lammps_gather_concat_mask(1), 1);
   EXPECT_EQ(f_lammps_gather_concat_mask(2), 1);
   EXPECT_EQ(f_lammps_gather_concat_mask(3), 1);
   lammps_command(lmp, "group special id 1");
   lammps_command(lmp, "group other id 2");
   lammps_command(lmp, "group spiffy id 3");
   EXPECT_EQ(f_lammps_gather_concat_mask(1), 3);
   EXPECT_EQ(f_lammps_gather_concat_mask(2), 5);
   EXPECT_EQ(f_lammps_gather_concat_mask(3), 9);
   lammps_command(lmp, "group other id 1");
   EXPECT_EQ(f_lammps_gather_concat_mask(1), 7);
   EXPECT_EQ(f_lammps_gather_concat_mask(2), 5);
   EXPECT_EQ(f_lammps_gather_concat_mask(3), 9);
};

TEST_F(LAMMPS_gather_scatter, gather_positions_concat)
{
   f_lammps_setup_gather_scatter();
   EXPECT_EQ(f_lammps_gather_concat_position(1,1), 1.0);
   EXPECT_EQ(f_lammps_gather_concat_position(2,1), 1.0);
   EXPECT_EQ(f_lammps_gather_concat_position(3,1), 1.5);
   EXPECT_EQ(f_lammps_gather_concat_position(1,2), 0.2);
   EXPECT_EQ(f_lammps_gather_concat_position(2,2), 0.1);
   EXPECT_EQ(f_lammps_gather_concat_position(3,2), 0.1);
   EXPECT_EQ(f_lammps_gather_concat_position(1,3), 0.5);
   EXPECT_EQ(f_lammps_gather_concat_position(2,3), 0.5);
   EXPECT_EQ(f_lammps_gather_concat_position(3,3), 0.5);
};

TEST_F(LAMMPS_gather_scatter, gather_masks_subset)
{
   f_lammps_setup_gather_scatter();
   EXPECT_EQ(f_lammps_gather_subset_mask(2), 1);
   EXPECT_EQ(f_lammps_gather_subset_mask(3), 1);
   lammps_command(lmp, "group special id 1");
   lammps_command(lmp, "group other id 2");
   lammps_command(lmp, "group spiffy id 3");
   EXPECT_EQ(f_lammps_gather_subset_mask(2), 5);
   EXPECT_EQ(f_lammps_gather_subset_mask(3), 9);
   lammps_command(lmp, "group other id 3");
   EXPECT_EQ(f_lammps_gather_subset_mask(2), 5);
   EXPECT_EQ(f_lammps_gather_subset_mask(3), 13);
};

TEST_F(LAMMPS_gather_scatter, gather_positions_subset)
{
   f_lammps_setup_gather_scatter();
//   EXPECT_EQ(f_lammps_gather_subset_position(1,1), 1.0);
//   EXPECT_EQ(f_lammps_gather_subset_position(2,1), 1.0);
//   EXPECT_EQ(f_lammps_gather_subset_position(3,1), 1.5);
   EXPECT_EQ(f_lammps_gather_subset_position(1,2), 0.2);
   EXPECT_EQ(f_lammps_gather_subset_position(2,2), 0.1);
   EXPECT_EQ(f_lammps_gather_subset_position(3,2), 0.1);
   EXPECT_EQ(f_lammps_gather_subset_position(1,3), 0.5);
   EXPECT_EQ(f_lammps_gather_subset_position(2,3), 0.5);
   EXPECT_EQ(f_lammps_gather_subset_position(3,3), 0.5);
};
