
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
void f_lammps_setup_fix_external_callback();
void f_lammps_setup_fix_external_array();
void f_lammps_set_fix_external_callbacks();
void f_lammps_get_force(int, double*);
void f_lammps_reverse_direction();
void f_lammps_find_forces();
void f_lammps_add_energy();
void f_lammps_set_virial();
double f_lammps_find_peratom_energy(int);
void f_lammps_find_peratom_virial(double[6], int);
void f_lammps_fixexternal_set_vector();
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
    f_lammps_setup_fix_external_callback();
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

TEST_F(LAMMPS_fixexternal, array)
{
    f_lammps_setup_fix_external_array();
    double **f;
    f = (double**) lammps_extract_atom(lmp, "f");
    f_lammps_find_forces();
    lammps_command(lmp, "run 0");
    EXPECT_DOUBLE_EQ(f[0][0], 14.0);
    EXPECT_DOUBLE_EQ(f[0][1], -14.0);
    EXPECT_DOUBLE_EQ(f[0][2], 18.0);
    EXPECT_DOUBLE_EQ(f[1][0], 16.0);
    EXPECT_DOUBLE_EQ(f[1][1], -16.0);
    EXPECT_DOUBLE_EQ(f[1][2], 20.0);
};

TEST_F(LAMMPS_fixexternal, energy_global)
{
    f_lammps_setup_fix_external_array();
    double energy;
    f_lammps_add_energy();
    lammps_command(lmp, "run 0");
    energy = lammps_get_thermo(lmp, "etotal");
    EXPECT_DOUBLE_EQ(energy, -20.2);
};

TEST_F(LAMMPS_fixexternal, virial_global)
{
    f_lammps_setup_fix_external_array();
    double virial[6], volume;
    f_lammps_set_virial();
    lammps_command(lmp, "run 0");
    volume = lammps_get_thermo(lmp, "vol");
    virial[0] = lammps_get_thermo(lmp, "pxx");
    virial[1] = lammps_get_thermo(lmp, "pyy");
    virial[2] = lammps_get_thermo(lmp, "pzz");
    virial[3] = lammps_get_thermo(lmp, "pxy");
    virial[4] = lammps_get_thermo(lmp, "pxz");
    virial[5] = lammps_get_thermo(lmp, "pyz");
    EXPECT_DOUBLE_EQ(virial[0], 1.0/volume);
    EXPECT_DOUBLE_EQ(virial[1], 2.0/volume);
    EXPECT_DOUBLE_EQ(virial[2], 2.5/volume);
    EXPECT_DOUBLE_EQ(virial[3], -1.0/volume);
    EXPECT_DOUBLE_EQ(virial[4], -2.25/volume);
    EXPECT_DOUBLE_EQ(virial[5], -3.02/volume);
};

TEST_F(LAMMPS_fixexternal, energy_peratom)
{
    f_lammps_setup_fix_external_callback();
    f_lammps_set_fix_external_callbacks();
    lammps_command(lmp, "compute peratom all pe/atom");
    double energy;
    lammps_command(lmp, "run 0");
    int nlocal = lammps_extract_setting(lmp, "nlocal");
    for (int i = 1; i <= nlocal; i++)
    {
        energy = f_lammps_find_peratom_energy(i);
        if (i == 1)
          EXPECT_DOUBLE_EQ(energy, 1.0);
        else
          EXPECT_DOUBLE_EQ(energy, 10.0);
    }
};

TEST_F(LAMMPS_fixexternal, virial_peratom)
{
    f_lammps_setup_fix_external_callback();
    f_lammps_set_fix_external_callbacks();
    lammps_command(lmp, "compute vperatom all stress/atom NULL");
    double virial[6];
    lammps_command(lmp, "run 0");
    int nlocal = lammps_extract_setting(lmp, "nlocal");
    for (int i = 1; i <= nlocal; i++)
    {
        f_lammps_find_peratom_virial(virial, i);
        if (i == 1)
        {
          EXPECT_DOUBLE_EQ(virial[0], -1.0);
          EXPECT_DOUBLE_EQ(virial[1], -2.0);
          EXPECT_DOUBLE_EQ(virial[2], 1.0);
          EXPECT_DOUBLE_EQ(virial[3], 2.0);
          EXPECT_DOUBLE_EQ(virial[4], -3.0);
          EXPECT_DOUBLE_EQ(virial[5], 3.0);
        }
        else
        {
          EXPECT_DOUBLE_EQ(virial[0], -10.0);
          EXPECT_DOUBLE_EQ(virial[1], -20.0);
          EXPECT_DOUBLE_EQ(virial[2], 10.0);
          EXPECT_DOUBLE_EQ(virial[3], 20.0);
          EXPECT_DOUBLE_EQ(virial[4], -30.0);
          EXPECT_DOUBLE_EQ(virial[5], 30.0);
        }
    }
};

TEST_F(LAMMPS_fixexternal, vector)
{
    f_lammps_setup_fix_external_callback();
    f_lammps_set_fix_external_callbacks();
    f_lammps_fixexternal_set_vector();
    lammps_command(lmp, "run 0");
    double *v;
    for (int i = 0; i < 8; i++)
    {
      v = (double*) lammps_extract_fix(lmp, "ext2", LMP_STYLE_GLOBAL,
        LMP_TYPE_VECTOR, i, 1);
      EXPECT_DOUBLE_EQ(i+1, *v);
      std::free(v);
    }
};
