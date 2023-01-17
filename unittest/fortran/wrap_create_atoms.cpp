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
void f_lammps_create_two_more();
void f_lammps_create_two_more_small();
void f_lammps_create_two_more_big();
void f_lammps_create_two_more_small2();
void f_lammps_create_two_more_big2();
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

TEST_F(LAMMPS_create_atoms, create_three)
{
    f_lammps_setup_create_atoms();
#ifdef LAMMPS_BIGBIG
    int64_t *tag, *image;
#else
    int *tag, *image;
#endif
    double **x, **v;
    EXPECT_EQ(lmp->atom->nlocal, 3);
    tag = lmp->atom->tag;
    image = lmp->atom->image;
    x = lmp->atom->x;
    v = lmp->atom->v;
    f_lammps_create_three_atoms();
    EXPECT_EQ(lmp->atom->nlocal, 6);
    for (int i = 0; i < lmp->atom->nlocal; i++) {
        if (tag[i] == 4) {
            EXPECT_EQ(image[i],lammps_encode_image_flags(1,-1,3));
            EXPECT_DOUBLE_EQ(x[i][0],1.0);
            EXPECT_DOUBLE_EQ(x[i][1],1.8);
            EXPECT_DOUBLE_EQ(x[i][2],2.718281828);
            EXPECT_DOUBLE_EQ(v[i][0],0.0);
            EXPECT_DOUBLE_EQ(v[i][1],1.0);
            EXPECT_DOUBLE_EQ(v[i][2],-1.0);
        }
        if (tag[i] == 5) {
            EXPECT_EQ(image[i],lammps_encode_image_flags(-2,-2,1));
            EXPECT_DOUBLE_EQ(x[i][0],1.8);
            EXPECT_DOUBLE_EQ(x[i][1],0.1);
            EXPECT_DOUBLE_EQ(x[i][2],1.8);
            EXPECT_DOUBLE_EQ(v[i][0],1.0);
            EXPECT_DOUBLE_EQ(v[i][1],-1.0);
            EXPECT_DOUBLE_EQ(v[i][2],3.0);
        }
        if (tag[i] == 6) {
            EXPECT_EQ(image[i],lammps_encode_image_flags(-2,0,0));
            EXPECT_DOUBLE_EQ(x[i][0],0.6);
            EXPECT_DOUBLE_EQ(x[i][1],0.8);
            EXPECT_DOUBLE_EQ(x[i][2],2.2);
            EXPECT_DOUBLE_EQ(v[i][0],0.1);
            EXPECT_DOUBLE_EQ(v[i][1],0.2);
            EXPECT_DOUBLE_EQ(v[i][2],-0.2);
        }
    }
};

TEST_F(LAMMPS_create_atoms, create_two_more)
{
    f_lammps_setup_create_atoms();
#ifdef LAMMPS_BIGBIG
    int64_t *tag, *image;
#else
    int *tag, *image;
#endif
    double **x, **v;
    f_lammps_create_three_atoms();
    EXPECT_EQ(lmp->atom->nlocal, 6);
    f_lammps_create_two_more();
    EXPECT_EQ(lmp->atom->nlocal, 8);
    tag = lmp->atom->tag;
    image = lmp->atom->image;
    x = lmp->atom->x;
    v = lmp->atom->v;
    for (int i = 0; i < lmp->atom->nlocal; i++) {
        if (tag[i] == 7) {
            EXPECT_EQ(image[i],lammps_encode_image_flags(0,0,0));
            EXPECT_DOUBLE_EQ(x[i][0],0.1);
            EXPECT_DOUBLE_EQ(x[i][1],1.9);
            EXPECT_DOUBLE_EQ(x[i][2],3.8);
            EXPECT_DOUBLE_EQ(v[i][0],0.0);
            EXPECT_DOUBLE_EQ(v[i][1],0.0);
            EXPECT_DOUBLE_EQ(v[i][2],0.0);
        }
        if (tag[i] == 8) {
            EXPECT_EQ(image[i],lammps_encode_image_flags(0,0,0));
            EXPECT_DOUBLE_EQ(x[i][0],1.2);
            EXPECT_DOUBLE_EQ(x[i][1],2.1);
            EXPECT_DOUBLE_EQ(x[i][2],1.25);
            EXPECT_DOUBLE_EQ(v[i][0],0.0);
            EXPECT_DOUBLE_EQ(v[i][1],0.0);
            EXPECT_DOUBLE_EQ(v[i][2],0.0);
        }
   }
};

TEST_F(LAMMPS_create_atoms, create_two_more_bigsmall)
{
    f_lammps_setup_create_atoms();
#ifdef LAMMPS_BIGBIG
    int64_t *tag, *image;
#else
    int *tag, *image;
#endif
    double **x, **v;
    f_lammps_create_three_atoms();
    EXPECT_EQ(lmp->atom->nlocal, 6);
#ifdef LAMMPS_BIGBIG
    f_lammps_create_two_more_big();
#else
    f_lammps_create_two_more_small();
#endif
    EXPECT_EQ(lmp->atom->nlocal, 8);
    tag = lmp->atom->tag;
    image = lmp->atom->image;
    x = lmp->atom->x;
    v = lmp->atom->v;
    for (int i = 0; i < lmp->atom->nlocal; i++) {
        if (tag[i] == 7) {
            EXPECT_EQ(image[i],lammps_encode_image_flags(-1,0,0));
            EXPECT_DOUBLE_EQ(x[i][0],1.2);
            EXPECT_DOUBLE_EQ(x[i][1],2.1);
            EXPECT_DOUBLE_EQ(x[i][2],1.25);
            EXPECT_DOUBLE_EQ(v[i][0],0.0);
            EXPECT_DOUBLE_EQ(v[i][1],0.0);
            EXPECT_DOUBLE_EQ(v[i][2],0.0);
        }
        if (tag[i] == 8) {
            EXPECT_EQ(image[i],lammps_encode_image_flags(1,0,0));
            EXPECT_DOUBLE_EQ(x[i][0],0.1);
            EXPECT_DOUBLE_EQ(x[i][1],1.9);
            EXPECT_DOUBLE_EQ(x[i][2],3.8);
            EXPECT_DOUBLE_EQ(v[i][0],0.0);
            EXPECT_DOUBLE_EQ(v[i][1],0.0);
            EXPECT_DOUBLE_EQ(v[i][2],0.0);
        }
   }
};

TEST_F(LAMMPS_create_atoms, create_two_more_bigsmall2)
{
    f_lammps_setup_create_atoms();
#ifdef LAMMPS_BIGBIG
    int64_t *tag, *image;
#else
    int *tag, *image;
#endif
    double **x, **v;
    f_lammps_create_three_atoms();
    EXPECT_EQ(lmp->atom->nlocal, 6);
#ifdef LAMMPS_BIGBIG
    f_lammps_create_two_more_big2();
#else
    f_lammps_create_two_more_small2();
#endif
    EXPECT_EQ(lmp->atom->nlocal, 8);
    tag = lmp->atom->tag;
    image = lmp->atom->image;
    x = lmp->atom->x;
    v = lmp->atom->v;
    for (int i = 0; i < lmp->atom->nlocal; i++) {
        if (tag[i] == 7) {
            EXPECT_EQ(image[i],lammps_encode_image_flags(0,0,0));
            EXPECT_DOUBLE_EQ(x[i][0],1.2);
            EXPECT_DOUBLE_EQ(x[i][1],2.1);
            EXPECT_DOUBLE_EQ(x[i][2],1.25);
            EXPECT_DOUBLE_EQ(v[i][0],0.0);
            EXPECT_DOUBLE_EQ(v[i][1],0.0);
            EXPECT_DOUBLE_EQ(v[i][2],0.0);
        }
        if (tag[i] == 8) {
            EXPECT_EQ(image[i],lammps_encode_image_flags(0,0,0));
            EXPECT_DOUBLE_EQ(x[i][0],0.1);
            EXPECT_DOUBLE_EQ(x[i][1],1.9);
            EXPECT_DOUBLE_EQ(x[i][2],3.8);
            EXPECT_DOUBLE_EQ(v[i][0],0.0);
            EXPECT_DOUBLE_EQ(v[i][1],0.0);
            EXPECT_DOUBLE_EQ(v[i][2],0.0);
        }
   }
};
