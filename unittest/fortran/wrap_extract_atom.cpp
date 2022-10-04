// unit tests for extracting Atom class data from a LAMMPS instance through the
// Fortran wrapper

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
void f_lammps_setup_extract_atom();
double f_lammps_extract_atom_mass();
int f_lammps_extract_atom_tag_int(int);
int64_t f_lammps_extract_atom_tag_int64(int64_t);
int f_lammps_extract_atom_type(int);
int f_lammps_extract_atom_mask(int);
void f_lammps_extract_atom_x(int, double *);
void f_lammps_extract_atom_v(int, double *);
}

class LAMMPS_extract_atom : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_extract_atom()           = default;
    ~LAMMPS_extract_atom() override = default;

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

TEST_F(LAMMPS_extract_atom, mass)
{
    f_lammps_setup_extract_atom();
    EXPECT_DOUBLE_EQ(f_lammps_extract_atom_mass(), 2.0);
};

TEST_F(LAMMPS_extract_atom, tag)
{
    f_lammps_setup_extract_atom();
#if defined(LAMMPS_BIGBIG)
    EXPECT_EQ(f_lammps_extract_atom_tag_int64(1l), 1l);
    EXPECT_EQ(f_lammps_extract_atom_tag_int64(2l), 2l);
#else
    EXPECT_EQ(f_lammps_extract_atom_tag_int(1), 1);
    EXPECT_EQ(f_lammps_extract_atom_tag_int(2), 2);
#endif
};

TEST_F(LAMMPS_extract_atom, type)
{
    f_lammps_setup_extract_atom();
    EXPECT_EQ(f_lammps_extract_atom_type(1), 1);
    EXPECT_EQ(f_lammps_extract_atom_type(2), 1);
};

TEST_F(LAMMPS_extract_atom, mask)
{
    f_lammps_setup_extract_atom();
    EXPECT_EQ(f_lammps_extract_atom_mask(1), 1);
    EXPECT_EQ(f_lammps_extract_atom_mask(2), 1);
    lammps_command(lmp, "group 1 id 1");
    lammps_command(lmp, "group 2 id 2");
    EXPECT_EQ(f_lammps_extract_atom_mask(1), 3);
    EXPECT_EQ(f_lammps_extract_atom_mask(2), 5);
};

TEST_F(LAMMPS_extract_atom, x)
{
    f_lammps_setup_extract_atom();
    double x1[3];
    double x2[3];
    f_lammps_extract_atom_x(1, x1);
    EXPECT_DOUBLE_EQ(x1[0], 1.0);
    EXPECT_DOUBLE_EQ(x1[1], 1.0);
    EXPECT_DOUBLE_EQ(x1[2], 1.5);
    f_lammps_extract_atom_x(2, x2);
    EXPECT_DOUBLE_EQ(x2[0], 0.2);
    EXPECT_DOUBLE_EQ(x2[1], 0.1);
    EXPECT_DOUBLE_EQ(x2[2], 0.1);
}

TEST_F(LAMMPS_extract_atom, v)
{
    f_lammps_setup_extract_atom();
    double v1[3];
    double v2[3];
    f_lammps_extract_atom_v(1, v1);
    EXPECT_DOUBLE_EQ(v1[0], 0.0);
    EXPECT_DOUBLE_EQ(v1[1], 0.0);
    EXPECT_DOUBLE_EQ(v1[2], 0.0);
    f_lammps_extract_atom_v(2, v2);
    EXPECT_DOUBLE_EQ(v2[0], 0.0);
    EXPECT_DOUBLE_EQ(v2[1], 0.0);
    EXPECT_DOUBLE_EQ(v2[2], 0.0);
    lammps_command(lmp, "group one id 1");
    lammps_command(lmp, "velocity one set 1 2 3");
    f_lammps_extract_atom_v(1, v1);
    EXPECT_DOUBLE_EQ(v1[0], 1.0);
    EXPECT_DOUBLE_EQ(v1[1], 2.0);
    EXPECT_DOUBLE_EQ(v1[2], 3.0);
}
