// unit tests for creating atoms in a LAMMPS instance through the Fortran wrapper

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
void f_lammps_setup_create_atoms();
void f_lammps_create_three_atoms();
}

class LAMMPS_create_atoms : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_create_atoms()           = default;
    ~LAMMPS_create_atoms() override = default;

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

TEST_F(LAMMPS_create_atoms, create_two)
{
    f_lammps_setup_create_atoms();
#ifdef LAMMPS_BIGBIG
    int64_t *tag, *image;
#else
    int *tag, *image;
#endif
    int *type;
    double **x, **v;
    EXPECT_EQ(lmp->atom->nlocal, 3);
    tag = lmp->atom->tag;
    image = lmp->atom->image;
    x = lmp->atom->x;
    v = lmp->atom->v;
    type = lmp->atom->type;
    f_lammps_create_three_atoms();
    EXPECT_EQ(lmp->atom->nlocal, 6);
    for (int i = 0; i < lmp->atom->nlocal; i++) {
        if (tag[i] == 4) {
            EXPECT_EQ(image[i],0);
            EXPECT_DOUBLE_EQ(x[i][0],1.0);
            EXPECT_DOUBLE_EQ(x[i][1],1.8);
            EXPECT_DOUBLE_EQ(x[i][2],2.718281828);
            EXPECT_DOUBLE_EQ(v[i][0],0.0);
            EXPECT_DOUBLE_EQ(v[i][1],1.0);
            EXPECT_DOUBLE_EQ(v[i][2],-1.0);
        }
        if (tag[i] == 5) {
            EXPECT_EQ(image[i],1);
            EXPECT_DOUBLE_EQ(x[i][0],1.8);
            EXPECT_DOUBLE_EQ(x[i][1],0.1);
            EXPECT_DOUBLE_EQ(x[i][2],1.8);
            EXPECT_DOUBLE_EQ(v[i][0],1.0);
            EXPECT_DOUBLE_EQ(v[i][1],-1.0);
            EXPECT_DOUBLE_EQ(v[i][2],3.0);
        }
        if (tag[i] == 6) {
            EXPECT_EQ(image[i],0);
            EXPECT_DOUBLE_EQ(x[i][0],0.6);
            EXPECT_DOUBLE_EQ(x[i][1],0.8);
            EXPECT_DOUBLE_EQ(x[i][2],2.2);
            EXPECT_DOUBLE_EQ(v[i][0],0.1);
            EXPECT_DOUBLE_EQ(v[i][1],0.2);
            EXPECT_DOUBLE_EQ(v[i][2],-0.2);
        }
    }
};
