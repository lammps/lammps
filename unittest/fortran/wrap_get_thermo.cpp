// unit tests for getting thermodynamic output from a LAMMPS instance through the Fortran wrapper

#include "lammps.h"
#include "lmptype.h"
#include <mpi.h>
#include <string>

#include "gtest/gtest.h"

// prototypes for fortran reverse wrapper functions
extern "C" {
void *f_lammps_with_args();
void f_lammps_close();
void f_lammps_get_thermo_setup();
double f_lammps_get_thermo_natoms();
double f_lammps_get_thermo_dt();
double f_lammps_get_thermo_vol();
double f_lammps_get_thermo_lx();
double f_lammps_get_thermo_ly();
double f_lammps_get_thermo_lz();
double f_lammps_get_thermo_xlo();
double f_lammps_get_thermo_xhi();
double f_lammps_get_thermo_ylo();
double f_lammps_get_thermo_yhi();
double f_lammps_get_thermo_zlo();
double f_lammps_get_thermo_zhi();

void f_lammps_last_thermo_setup();
int f_lammps_last_thermo_step();
int f_lammps_last_thermo_num();
int f_lammps_last_thermo_type(int);
const char *f_lammps_last_thermo_string(int);
int f_lammps_last_thermo_int(int);
double f_lammps_last_thermo_double(int);
}

using LAMMPS_NS::multitype;

class LAMMPS_thermo : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_thermo()           = default;
    ~LAMMPS_thermo() override = default;

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

TEST_F(LAMMPS_thermo, get_thermo)
{
    EXPECT_DOUBLE_EQ(f_lammps_get_thermo_natoms(), 0.0);
    f_lammps_get_thermo_setup();
    EXPECT_DOUBLE_EQ(f_lammps_get_thermo_natoms(), 2.0);
    EXPECT_DOUBLE_EQ(f_lammps_get_thermo_dt(), 0.005);
    EXPECT_DOUBLE_EQ(f_lammps_get_thermo_vol(), 24.0);
    EXPECT_DOUBLE_EQ(f_lammps_get_thermo_lx(), 2.0);
    EXPECT_DOUBLE_EQ(f_lammps_get_thermo_ly(), 3.0);
    EXPECT_DOUBLE_EQ(f_lammps_get_thermo_lz(), 4.0);
    EXPECT_DOUBLE_EQ(f_lammps_get_thermo_xlo(), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_get_thermo_xhi(), 2.0);
    EXPECT_DOUBLE_EQ(f_lammps_get_thermo_ylo(), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_get_thermo_yhi(), 3.0);
    EXPECT_DOUBLE_EQ(f_lammps_get_thermo_zlo(), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_get_thermo_zhi(), 4.0);
};

TEST_F(LAMMPS_thermo, last_thermo)
{
    EXPECT_EQ(f_lammps_last_thermo_step(), -1);
    EXPECT_EQ(f_lammps_last_thermo_type(1), multitype::LAMMPS_NONE);
    EXPECT_EQ(f_lammps_last_thermo_type(2), multitype::LAMMPS_NONE);
    f_lammps_last_thermo_setup();
    EXPECT_EQ(f_lammps_last_thermo_step(), 15);
    EXPECT_EQ(f_lammps_last_thermo_num(), 6);
    EXPECT_STREQ(f_lammps_last_thermo_string(1), "Step");
    EXPECT_STREQ(f_lammps_last_thermo_string(2), "Temp");
    EXPECT_STREQ(f_lammps_last_thermo_string(3), "E_pair");
    EXPECT_STREQ(f_lammps_last_thermo_string(6), "Press");
#if defined(LAMMPS_SMALLSMALL)
    EXPECT_EQ(f_lammps_last_thermo_type(1), multitype::LAMMPS_INT);
#else
    EXPECT_EQ(f_lammps_last_thermo_type(1), multitype::LAMMPS_INT64);
#endif
    EXPECT_EQ(f_lammps_last_thermo_int(1), 15);
    EXPECT_EQ(f_lammps_last_thermo_type(2), multitype::LAMMPS_DOUBLE);
    EXPECT_EQ(f_lammps_last_thermo_type(3), multitype::LAMMPS_DOUBLE);
    EXPECT_EQ(f_lammps_last_thermo_type(4), multitype::LAMMPS_DOUBLE);
    EXPECT_EQ(f_lammps_last_thermo_type(5), multitype::LAMMPS_DOUBLE);
    EXPECT_EQ(f_lammps_last_thermo_type(6), multitype::LAMMPS_DOUBLE);
    EXPECT_DOUBLE_EQ(f_lammps_last_thermo_double(2), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_last_thermo_double(3), -0.13713425198078993);
    EXPECT_DOUBLE_EQ(f_lammps_last_thermo_double(4), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_last_thermo_double(5), -0.13713425198078993);
    EXPECT_DOUBLE_EQ(f_lammps_last_thermo_double(6), -0.022421073321023492);
};
