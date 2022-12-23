// unit tests for gathering and scattering data from a LAMMPS instance through
// the Fortran wrapper

#include "lammps.h"
#include "library.h"
#include "atom.h"
#include <cstdint>
#include <cstdlib>
#include <mpi.h>
#include <string>

#include "gtest/gtest.h"

// prototypes for Fortran reverse wrapper functions
extern "C" {
void *f_lammps_with_args();
void f_lammps_close();
void f_lammps_setup_gather_scatter();
int f_lammps_gather_atoms_mask(int);
double f_lammps_gather_atoms_position(int);
int f_lammps_gather_atoms_concat_mask(int);
double f_lammps_gather_atoms_concat_position(int, int);
int f_lammps_gather_atoms_subset_mask(int);
double f_lammps_gather_atoms_subset_position(int, int);
void f_lammps_scatter_atoms_masks();
void f_lammps_scatter_atoms_positions();
void f_lammps_setup_gather_bonds();
int f_lammps_test_gather_bonds_small();
int f_lammps_test_gather_bonds_big();
double f_lammps_gather_pe_atom(int);
double f_lammps_gather_pe_atom_concat(int);
void f_lammps_gather_pe_atom_subset(int*, double*);
void f_lammps_scatter_compute();
void f_lammps_scatter_subset_compute();
}

using namespace LAMMPS_NS;

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
    if (lammps_extract_setting(nullptr, "tagint") == 8) GTEST_SKIP();
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
    if (lammps_extract_setting(nullptr, "tagint") == 8) GTEST_SKIP();
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
    if (lammps_extract_setting(nullptr, "tagint") == 8) GTEST_SKIP();
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
    if (lammps_extract_setting(nullptr, "tagint") == 8) GTEST_SKIP();
    f_lammps_setup_gather_scatter();
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2, 1), 1.0);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3, 1), 1.5);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1, 2), 0.2);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2, 2), 0.1);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3, 2), 0.1);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1, 3), 0.5);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2, 3), 0.5);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3, 3), 0.5);
};

TEST_F(LAMMPS_gather_scatter, gather_atoms_subset_masks)
{
    if (lammps_extract_setting(nullptr, "tagint") == 8) GTEST_SKIP();
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
    if (lammps_extract_setting(nullptr, "tagint") == 8) GTEST_SKIP();
    f_lammps_setup_gather_scatter();
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_position(1, 2), 0.2);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_position(2, 2), 0.1);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_position(3, 2), 0.1);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_position(1, 3), 0.5);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_position(2, 3), 0.5);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_subset_position(3, 3), 0.5);
};

TEST_F(LAMMPS_gather_scatter, scatter_atoms_masks)
{
    if (lammps_extract_setting(nullptr, "tagint") == 8) GTEST_SKIP();
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
    if (lammps_extract_setting(nullptr, "tagint") == 8) GTEST_SKIP();
    f_lammps_setup_gather_scatter();
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2, 1), 1.0);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3, 1), 1.5);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1, 2), 0.2);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2, 2), 0.1);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3, 2), 0.1);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1, 3), 0.5);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2, 3), 0.5);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3, 3), 0.5);
    f_lammps_scatter_atoms_positions();
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1, 3), 1.0);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2, 3), 1.0);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3, 3), 1.5);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1, 2), 0.2);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2, 2), 0.1);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3, 2), 0.1);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(1, 1), 0.5);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(2, 1), 0.5);
    EXPECT_DOUBLE_EQ(f_lammps_gather_atoms_concat_position(3, 1), 0.5);
};

TEST_F(LAMMPS_gather_scatter, scatter_atoms_subset_mask)
{
    if (lammps_extract_setting(nullptr, "tagint") == 8) GTEST_SKIP();
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

TEST_F(LAMMPS_gather_scatter, gather_bonds)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    f_lammps_setup_gather_bonds();
#ifdef LAMMPS_BIGBIG
    EXPECT_EQ(f_lammps_test_gather_bonds_big(), 1);
#else
    EXPECT_EQ(f_lammps_test_gather_bonds_small(), 1);
#endif
};

TEST_F(LAMMPS_gather_scatter, gather_compute)
{
#ifdef LAMMPS_BIGBIG
    GTEST_SKIP();
#else
    f_lammps_setup_gather_scatter();
    lammps_command(lmp, "run 0");
    int natoms = lmp->atom->natoms;
    int *tag = lmp->atom->tag;
    double *pe = (double*) lammps_extract_compute(lmp, "pe", LMP_STYLE_ATOM,
        LMP_TYPE_VECTOR);
    for (int i = 0; i < natoms; i++)
        EXPECT_DOUBLE_EQ(f_lammps_gather_pe_atom(tag[i]), pe[i]);
#endif
};

TEST_F(LAMMPS_gather_scatter, gather_compute_concat)
{
#ifdef LAMMPS_BIGBIG
    GTEST_SKIP();
#else
    f_lammps_setup_gather_scatter();
    lammps_command(lmp, "run 0");
    int natoms = lmp->atom->natoms;
    int *tag = lmp->atom->tag;
    double *pe = (double*) lammps_extract_compute(lmp, "pe", LMP_STYLE_ATOM,
        LMP_TYPE_VECTOR);
    for (int i = 0; i < natoms; i++)
        EXPECT_DOUBLE_EQ(f_lammps_gather_pe_atom(tag[i]), pe[i]);
#endif
};

TEST_F(LAMMPS_gather_scatter, gather_compute_subset)
{
#ifdef LAMMPS_BIGBIG
    GTEST_SKIP();
#else
    f_lammps_setup_gather_scatter();
    lammps_command(lmp, "run 0");
    int ids[2] = {3, 1};
    int *tag = lmp->atom->tag;
    double pe[2] = {0.0, 0.0};
    int nlocal = lammps_extract_setting(lmp, "nlocal");
    double *pa_pe = (double*) lammps_extract_compute(lmp, "pe", LMP_STYLE_ATOM,
      LMP_TYPE_VECTOR);

    for (int i = 0; i < nlocal; i++) {
      if(tag[i] == ids[0]) pe[0] = pa_pe[i];
      if(tag[i] == ids[1]) pe[1] = pa_pe[i];
    }

    double ftn_pe[2];
    f_lammps_gather_pe_atom_subset(ids, ftn_pe);
    EXPECT_DOUBLE_EQ(ftn_pe[0], pe[0]);
    EXPECT_DOUBLE_EQ(ftn_pe[1], pe[1]);
#endif
};

TEST_F(LAMMPS_gather_scatter, scatter_compute)
{
#ifdef LAMMPS_BIGBIG
    GTEST_SKIP();
#else
    f_lammps_setup_gather_scatter();
    int natoms = lmp->atom->natoms;
    double *pe = new double[natoms];
    lammps_command(lmp, "run 0");
    lammps_gather(lmp, "c_pe", 1, 1, pe);
    double *old_pe = new double[natoms];
    for (int i = 0; i < natoms; i++)
      old_pe[i] = pe[i];
    EXPECT_DOUBLE_EQ(pe[0], old_pe[0]);
    EXPECT_DOUBLE_EQ(pe[1], old_pe[1]);
    EXPECT_DOUBLE_EQ(pe[2], old_pe[2]);
    f_lammps_scatter_compute();
    lammps_gather(lmp, "c_pe", 1, 1, pe);
    EXPECT_DOUBLE_EQ(pe[0], old_pe[2]);
    EXPECT_DOUBLE_EQ(pe[1], old_pe[1]);
    EXPECT_DOUBLE_EQ(pe[2], old_pe[0]);
    delete[] old_pe;
    delete[] pe;
#endif
};

TEST_F(LAMMPS_gather_scatter, scatter_subset_compute)
{
#ifdef LAMMPS_BIGBIG
    GTEST_SKIP();
#else
    f_lammps_setup_gather_scatter();
    int natoms = lmp->atom->natoms;
    double *pe = new double[natoms];
    lammps_command(lmp, "run 0");
    lammps_gather(lmp, "c_pe", 1, 1, pe);
    double *old_pe = new double[natoms];
    for (int i = 0; i < natoms; i++)
      old_pe[i] = pe[i];
    EXPECT_DOUBLE_EQ(pe[0], old_pe[0]);
    EXPECT_DOUBLE_EQ(pe[1], old_pe[1]);
    EXPECT_DOUBLE_EQ(pe[2], old_pe[2]);
    f_lammps_scatter_subset_compute();
    lammps_gather(lmp, "c_pe", 1, 1, pe);
    EXPECT_DOUBLE_EQ(pe[0], old_pe[2]);
    EXPECT_DOUBLE_EQ(pe[1], old_pe[1]);
    EXPECT_DOUBLE_EQ(pe[2], old_pe[0]);
    delete[] old_pe;
    delete[] pe;
#endif
};
