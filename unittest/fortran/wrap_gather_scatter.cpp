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
int f_lammps_gather_atoms_mask(int);
double f_lammps_gather_atoms_position(int);
int f_lammps_gather_atoms_concat_mask(int);
double f_lammps_gather_atoms_concat_position(int,int);
int f_lammps_gather_atoms_subset_mask(int);
double f_lammps_gather_atoms_subset_position(int,int);
void f_lammps_scatter_atoms_masks();
void f_lammps_scatter_atoms_positions();
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

TEST_F(LAMMPS_gather_scatter, gather_atoms_masks)
{
  f_lammps_setup_gather_scatter();
  EXPECT_EQ(f_lammps_gather_atoms_mask(1), 1);
  EXPECT_EQ(f_lammps_gather_atoms_mask(2), 1);
  EXPECT_EQ(f_lammps_gather_atoms_mask(3), 1);
  lammps_command(lmp, "group special id 1");
  lammps_command(lmp, "group other id 2");
  lammps_command(lmp, "group spiffy id 3");
  EXPECT_EQ(f_lammps_gather_atoms_mask(1), 3);
  EXPECT_EQ(f_lammps_gather_atoms_mask(2), 5);
  EXPECT_EQ(f_lammps_gather_atoms_mask(3), 9);
  lammps_command(lmp, "group other id 1");
  EXPECT_EQ(f_lammps_gather_atoms_mask(1), 7);
  EXPECT_EQ(f_lammps_gather_atoms_mask(2), 5);
  EXPECT_EQ(f_lammps_gather_atoms_mask(3), 9);
};

TEST_F(LAMMPS_gather_scatter, gather_atoms_positions)
{
  f_lammps_setup_gather_scatter();
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_position(1), 1.0);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_position(2), 1.0);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_position(3), 1.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_position(4), 0.2);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_position(5), 0.1);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_position(6), 0.1);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_position(7), 0.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_position(8), 0.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_position(9), 0.5);   
};

TEST_F(LAMMPS_gather_scatter, gather_atoms_concat_masks)
{
  f_lammps_setup_gather_scatter();
  EXPECT_EQ(f_lammps_gather_atoms_concat_mask(1), 1);
  EXPECT_EQ(f_lammps_gather_atoms_concat_mask(2), 1);
  EXPECT_EQ(f_lammps_gather_atoms_concat_mask(3), 1);
  lammps_command(lmp, "group special id 1");
  lammps_command(lmp, "group other id 2");
  lammps_command(lmp, "group spiffy id 3");
  EXPECT_EQ(f_lammps_gather_atoms_concat_mask(1), 3);
  EXPECT_EQ(f_lammps_gather_atoms_concat_mask(2), 5);
  EXPECT_EQ(f_lammps_gather_atoms_concat_mask(3), 9);
  lammps_command(lmp, "group other id 1");
  EXPECT_EQ(f_lammps_gather_atoms_concat_mask(1), 7);
  EXPECT_EQ(f_lammps_gather_atoms_concat_mask(2), 5);
  EXPECT_EQ(f_lammps_gather_atoms_concat_mask(3), 9);
};

TEST_F(LAMMPS_gather_scatter, gather_atoms_concat_positions)
{
  f_lammps_setup_gather_scatter();
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1,1), 1.0);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2,1), 1.0);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3,1), 1.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1,2), 0.2);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2,2), 0.1);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3,2), 0.1);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1,3), 0.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2,3), 0.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3,3), 0.5);
};

TEST_F(LAMMPS_gather_scatter, gather_atoms_subset_masks)
{
  f_lammps_setup_gather_scatter();
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_mask(2), 1);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_mask(3), 1);
  lammps_command(lmp, "group special id 1");
  lammps_command(lmp, "group other id 2");
  lammps_command(lmp, "group spiffy id 3");
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_mask(2), 5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_mask(3), 9);
  lammps_command(lmp, "group other id 3");
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_mask(2), 5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_mask(3), 13);
};

TEST_F(LAMMPS_gather_scatter, gather_atoms_subset_positions)
{
  f_lammps_setup_gather_scatter();
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_position(1,2), 0.2);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_position(2,2), 0.1);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_position(3,2), 0.1);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_position(1,3), 0.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_position(2,3), 0.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_position(3,3), 0.5);
};

TEST_F(LAMMPS_gather_scatter, scatter_atoms_masks)
{
  f_lammps_setup_gather_scatter();
  lammps_command(lmp, "group special id 1");
  lammps_command(lmp, "group other id 2");
  lammps_command(lmp, "group spiffy id 3");
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_mask(1), 3);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_mask(2), 5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_mask(3), 9);
  f_lammps_scatter_atoms_masks();
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_mask(1), 9);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_mask(2), 5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_mask(3), 3);
};

TEST_F(LAMMPS_gather_scatter, scatter_atoms_positions)
{
  f_lammps_setup_gather_scatter();
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1,1), 1.0);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2,1), 1.0);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3,1), 1.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1,2), 0.2);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2,2), 0.1);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3,2), 0.1);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1,3), 0.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2,3), 0.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3,3), 0.5); 
  f_lammps_scatter_atoms_positions();
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1,3), 1.0);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2,3), 1.0);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3,3), 1.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1,2), 0.2);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2,2), 0.1);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3,2), 0.1);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1,1), 0.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2,1), 0.5);
  EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3,1), 0.5); 
};

TEST_F(LAMMPS_gather_scatter, scatter_atoms_subset_mask)
{
  f_lammps_setup_gather_scatter();
  EXPECT_EQ(f_lammps_gather_atoms_mask(1), 1);
  EXPECT_EQ(f_lammps_gather_atoms_mask(3), 1);
  lammps_command(lmp, "group special id 1");
  lammps_command(lmp, "group other id 2");
  lammps_command(lmp, "group spiffy id 3");
  EXPECT_EQ(f_lammps_gather_atoms_mask(1), 3);
  EXPECT_EQ(f_lammps_gather_atoms_mask(3), 9);
  f_lammps_scatter_atoms_masks();
  EXPECT_EQ(f_lammps_gather_atoms_mask(1), 9);
  EXPECT_EQ(f_lammps_gather_atoms_mask(3), 3);
};
