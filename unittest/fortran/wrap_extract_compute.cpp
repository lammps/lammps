// unit tests for extracting compute data from a LAMMPS instance through the
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
void f_lammps_setup_extract_compute();
double f_lammps_extract_compute_peratom_vector(int);
double f_lammps_extract_compute_peratom_array(int, int);
double f_lammps_extract_compute_global_scalar();
double f_lammps_extract_compute_global_vector(int);
double f_lammps_extract_compute_global_array(int, int);
double f_lammps_extract_compute_local_vector(int);
double f_lammps_extract_compute_local_array(int, int);
}

class LAMMPS_extract_compute : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_extract_compute()           = default;
    ~LAMMPS_extract_compute() override = default;

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

TEST_F(LAMMPS_extract_compute, peratom_vector)
{
    f_lammps_setup_extract_compute();
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_vector(1), -0.599703102447981);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_vector(2), 391.817623795857);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_vector(3), 391.430665759871);
};

TEST_F(LAMMPS_extract_compute, peratom_array)
{
    f_lammps_setup_extract_compute();
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_array(1, 1), 0.8837067009319107);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_array(2, 1), 0.3588584939803668);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_array(3, 1), 1.2799807127711049);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_array(4, 1), 0.20477632346642258);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_array(5, 1), 0.400429511840588);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_array(6, 1), 0.673995757699694);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_array(1, 2), -1070.0291234709418);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_array(2, 2), -1903.651817128683);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_array(3, 2), -1903.5121520875714);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_array(4, 2), -1427.867483013);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_array(5, 2), -1427.8560790941347);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_peratom_array(6, 2), -1903.5971655908565);
};

TEST_F(LAMMPS_extract_compute, global_scalar)
{
    f_lammps_setup_extract_compute();
    double *scalar;
    scalar = (double *)lammps_extract_compute(lmp, "totalpe", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR);
    // EXPECT_DOUBLE_EQ(f_lammps_extract_compute_global_scalar(), 782.64858645328);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_global_scalar(), *scalar);
};

TEST_F(LAMMPS_extract_compute, global_vector)
{
    f_lammps_setup_extract_compute();
    double *vector;
    vector = (double *)lammps_extract_compute(lmp, "COM", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_global_vector(1), vector[0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_global_vector(2), vector[1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_global_vector(3), vector[2]);
};

TEST_F(LAMMPS_extract_compute, global_array)
{
    f_lammps_setup_extract_compute();
    double **array;
    array = (double **)lammps_extract_compute(lmp, "RDF", LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_global_array(1, 1), array[0][0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_global_array(2, 1), array[0][1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_global_array(1, 2), array[1][0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_global_array(2, 2), array[1][1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_global_array(1, 3), array[2][0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_global_array(1, 4), array[3][0]);
};
TEST_F(LAMMPS_extract_compute, local_vector)
{
    f_lammps_setup_extract_compute();
    double *vector;
    vector = (double *)lammps_extract_compute(lmp, "pairdist", LMP_STYLE_LOCAL, LMP_TYPE_VECTOR);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_vector(1), vector[0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_vector(2), vector[1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_vector(3), vector[2]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_vector(4), vector[3]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_vector(5), vector[4]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_vector(6), vector[5]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_vector(7), vector[6]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_vector(8), vector[7]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_vector(9), vector[8]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_vector(10), vector[9]);
};

TEST_F(LAMMPS_extract_compute, local_array)
{
    f_lammps_setup_extract_compute();
    double **array;
    array = (double **)lammps_extract_compute(lmp, "pairlocal", LMP_STYLE_LOCAL, LMP_TYPE_ARRAY);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(1, 1), array[0][0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(2, 1), array[0][1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(3, 1), array[0][2]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(4, 1), array[0][3]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(1, 2), array[1][0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(2, 2), array[1][1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(3, 2), array[1][2]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(4, 2), array[1][3]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(1, 3), array[2][0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(2, 3), array[2][1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(3, 3), array[2][2]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(4, 3), array[2][3]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(1, 4), array[3][0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(2, 4), array[3][1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(3, 4), array[3][2]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(4, 4), array[3][3]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(1, 5), array[4][0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(2, 5), array[4][1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(3, 5), array[4][2]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(4, 5), array[4][3]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(1, 6), array[5][0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(2, 6), array[5][1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(3, 6), array[5][2]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(4, 6), array[5][3]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(1, 7), array[6][0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(2, 7), array[6][1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(3, 7), array[6][2]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(4, 7), array[6][3]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(1, 8), array[7][0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(2, 8), array[7][1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(3, 8), array[7][2]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(4, 8), array[7][3]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(1, 9), array[8][0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(2, 9), array[8][1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(3, 9), array[8][2]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(4, 9), array[8][3]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(1, 10), array[9][0]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(2, 10), array[9][1]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(3, 10), array[9][2]);
    EXPECT_DOUBLE_EQ(f_lammps_extract_compute_local_array(4, 10), array[9][3]);
};
