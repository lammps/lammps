// unit tests for getting thermodynamic output from a LAMMPS instance through the Fortran wrapper

#include "lammps.h"
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
}

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
