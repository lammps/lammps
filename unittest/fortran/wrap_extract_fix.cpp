// unit tests for extracting compute data from a LAMMPS instance through the
// Fortran wrapper
#include <cstdio>

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
void f_lammps_setup_extract_fix();
double f_lammps_extract_fix_global_scalar();
double f_lammps_extract_fix_global_vector(int);
double f_lammps_extract_fix_global_array(int, int);
double f_lammps_extract_fix_peratom_vector(int);
double f_lammps_extract_fix_peratom_array(int, int);
double f_lammps_extract_fix_local_vector(int);
double f_lammps_extract_fix_local_array(int, int);
}

class LAMMPS_extract_fix : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_extract_fix()           = default;
    ~LAMMPS_extract_fix() override = default;

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

TEST_F(LAMMPS_extract_fix, global_scalar)
{
    f_lammps_setup_extract_fix();
    double *scalar =
        (double *)lammps_extract_fix(lmp, "recenter", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, -1, -1);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_global_scalar(), *scalar);
    lammps_free(scalar);
};

TEST_F(LAMMPS_extract_fix, global_vector)
{
    f_lammps_setup_extract_fix();
    double *x =
        (double *)lammps_extract_fix(lmp, "recenter", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, 0, -1);
    double *y =
        (double *)lammps_extract_fix(lmp, "recenter", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, 1, -1);
    double *z =
        (double *)lammps_extract_fix(lmp, "recenter", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, 2, -1);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_global_vector(1), *x);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_global_vector(2), *y);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_global_vector(3), *z);
    lammps_free(x);
    lammps_free(y);
    lammps_free(z);
};

TEST_F(LAMMPS_extract_fix, global_array)
{
    f_lammps_setup_extract_fix();
    double natoms = lammps_get_natoms(lmp);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_global_array(1, 1), natoms);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_global_array(1, 2), natoms);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_global_array(2, 1), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_global_array(2, 2), 1.0);
};

TEST_F(LAMMPS_extract_fix, peratom_vector)
{
    f_lammps_setup_extract_fix();
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_peratom_vector(1), 1.5);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_peratom_vector(2), 0.1);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_peratom_vector(3), 0.5);
};

TEST_F(LAMMPS_extract_fix, peratom_array)
{
    f_lammps_setup_extract_fix();
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_peratom_array(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_peratom_array(2, 1), 1.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_peratom_array(3, 1), 1.5);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_peratom_array(1, 2), 0.2);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_peratom_array(2, 2), 0.1);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_peratom_array(3, 2), 0.1);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_peratom_array(1, 3), 0.5);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_peratom_array(2, 3), 0.5);
    EXPECT_DOUBLE_EQ(f_lammps_extract_fix_peratom_array(3, 3), 0.5);
};
