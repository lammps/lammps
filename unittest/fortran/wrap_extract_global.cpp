// unit tests for extracting global data from a LAMMPS instance through the
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
void f_lammps_setup_extract_global();
void f_lammps_setup_full_extract_global();
int f_lammps_extract_global_units();
int f_lammps_extract_global_ntimestep();
int64_t f_lammps_extract_global_ntimestep_big();
double f_lammps_extract_global_dt();
void f_lammps_extract_global_boxlo(double[3]);
void f_lammps_extract_global_boxhi(double[3]);
double f_lammps_extract_global_boxxlo();
double f_lammps_extract_global_boxylo();
double f_lammps_extract_global_boxzlo();
double f_lammps_extract_global_boxxhi();
double f_lammps_extract_global_boxyhi();
double f_lammps_extract_global_boxzhi();
void f_lammps_extract_global_periodicity(int[3]);
int f_lammps_extract_global_triclinic();
double f_lammps_extract_global_xy();
double f_lammps_extract_global_yz();
double f_lammps_extract_global_xz();
int f_lammps_extract_global_natoms();
int64_t f_lammps_extract_global_natoms_big();
int f_lammps_extract_global_nbonds();
int64_t f_lammps_extract_global_nbonds_big();
int f_lammps_extract_global_nangles();
int64_t f_lammps_extract_global_nangles_big();
int f_lammps_extract_global_ndihedrals();
int64_t f_lammps_extract_global_ndihedrals_big();
int f_lammps_extract_global_nimpropers();
int64_t f_lammps_extract_global_nimpropers_big();
int f_lammps_extract_global_ntypes();
int f_lammps_extract_global_nlocal();
int f_lammps_extract_global_nghost();
int f_lammps_extract_global_nmax();
double f_lammps_extract_global_boltz();
double f_lammps_extract_global_hplanck();
double f_lammps_extract_global_angstrom();
double f_lammps_extract_global_femtosecond();
}

class LAMMPS_extract_global : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_extract_global()           = default;
    ~LAMMPS_extract_global() override = default;

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

TEST_F(LAMMPS_extract_global, units)
{
    f_lammps_setup_extract_global();
    EXPECT_EQ(f_lammps_extract_global_units(), 1);
};

TEST_F(LAMMPS_extract_global, ntimestep)
{
    f_lammps_setup_extract_global();
#ifdef LAMMPS_SMALLSMALL
    EXPECT_EQ(f_lammps_extract_global_ntimestep(), 0);
#else
    EXPECT_EQ(f_lammps_extract_global_ntimestep_big(), 0l);
#endif
};

TEST_F(LAMMPS_extract_global, dt)
{
    f_lammps_setup_extract_global();
    EXPECT_DOUBLE_EQ(f_lammps_extract_global_dt(), 0.005);
};

TEST_F(LAMMPS_extract_global, boxprops)
{
    f_lammps_setup_extract_global();
    double boxlo[3], boxhi[3];
    f_lammps_extract_global_boxlo(boxlo);
    EXPECT_DOUBLE_EQ(boxlo[0], 0.0);
    EXPECT_DOUBLE_EQ(boxlo[1], 0.0);
    EXPECT_DOUBLE_EQ(boxlo[2], 0.0);
    f_lammps_extract_global_boxhi(boxhi);
    EXPECT_DOUBLE_EQ(boxhi[0], 2.0);
    EXPECT_DOUBLE_EQ(boxhi[1], 3.0);
    EXPECT_DOUBLE_EQ(boxhi[2], 4.0);

    EXPECT_DOUBLE_EQ(f_lammps_extract_global_boxxlo(), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_global_boxxhi(), 2.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_global_boxylo(), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_global_boxyhi(), 3.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_global_boxzlo(), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_global_boxzhi(), 4.0);

    int periodicity[3];
    f_lammps_extract_global_periodicity(periodicity);
    EXPECT_EQ(periodicity[0], 1);
    EXPECT_EQ(periodicity[1], 1);
    EXPECT_EQ(periodicity[2], 1);

    EXPECT_EQ(f_lammps_extract_global_triclinic(), 0);

    EXPECT_DOUBLE_EQ(f_lammps_extract_global_xy(), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_global_yz(), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_global_xz(), 0.0);
};

TEST_F(LAMMPS_extract_global, atomprops)
{
    f_lammps_setup_extract_global();
#ifdef LAMMPS_SMALLSMALL
    EXPECT_EQ(f_lammps_extract_global_natoms(), 2);
    EXPECT_EQ(f_lammps_extract_global_nbonds(), 0);
    EXPECT_EQ(f_lammps_extract_global_nangles(), 0);
    EXPECT_EQ(f_lammps_extract_global_ndihedrals(), 0);
#else
    EXPECT_EQ(f_lammps_extract_global_natoms_big(), 2l);
    EXPECT_EQ(f_lammps_extract_global_nbonds_big(), 0l);
    EXPECT_EQ(f_lammps_extract_global_nangles_big(), 0l);
    EXPECT_EQ(f_lammps_extract_global_ndihedrals_big(), 0l);
#endif

    EXPECT_EQ(f_lammps_extract_global_ntypes(), 1);
    EXPECT_EQ(f_lammps_extract_global_nlocal(), 2);
    EXPECT_EQ(f_lammps_extract_global_nghost(), 41);
    EXPECT_EQ(f_lammps_extract_global_nmax(), 16384);

    EXPECT_DOUBLE_EQ(f_lammps_extract_global_boltz(), 1.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_global_hplanck(), 1.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_global_angstrom(), 1.0);

    EXPECT_DOUBLE_EQ(f_lammps_extract_global_femtosecond(), 1.0);
};

TEST_F(LAMMPS_extract_global, fullprops)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    // This is not currently the world's most convincing test....
    f_lammps_setup_full_extract_global();
#ifdef LAMMPS_SMALLSMALL
    EXPECT_EQ(f_lammps_extract_global_natoms(), 2);
    EXPECT_EQ(f_lammps_extract_global_nbonds(), 0);
    EXPECT_EQ(f_lammps_extract_global_nangles(), 0);
    EXPECT_EQ(f_lammps_extract_global_ndihedrals(), 0);
#else
    EXPECT_EQ(f_lammps_extract_global_natoms_big(), 2l);
    EXPECT_EQ(f_lammps_extract_global_nbonds_big(), 0l);
    EXPECT_EQ(f_lammps_extract_global_nangles_big(), 0l);
    EXPECT_EQ(f_lammps_extract_global_ndihedrals_big(), 0l);
#endif
}
