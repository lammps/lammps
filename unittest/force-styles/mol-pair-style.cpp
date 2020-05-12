// unit tests for pair styles intended for molecular systems

#include "lammps.h"
#include "atom.h"
#include "input.h"
#include <mpi.h>
#include <cstring>

#include "gtest/gtest.h"

class MolPairStyle : public ::testing::Test
{
protected:
    LAMMPS_NS::LAMMPS *lmp;
    MolPairStyle() {
        const char *args[] = {"LAMMPS_test",
                              "-log", "none",
                              "-echo", "screen",
                              "-nocite" };
        char **argv = (char **)args;
        int argc = sizeof(args)/sizeof(char *);

        int flag;
        MPI_Initialized(&flag);
        if (!flag) MPI_Init(&argc,&argv);

        ::testing::internal::CaptureStdout();
        lmp = new LAMMPS_NS::LAMMPS(argc, argv, MPI_COMM_WORLD);
        ::testing::internal::GetCapturedStdout();
    }
    ~MolPairStyle() override {}

    void SetUp() override {
        lmp->input->one("clear");
        lmp->input->file("in.fourmol");
    }
    void TearDown() override {
    }
};

TEST_F(MolPairStyle, initial) {
    double **f=lmp->atom->f;
    double **v=lmp->atom->v;
    const int nlocal = lmp->atom->nlocal;

    // abort if running in parallel
    ASSERT_EQ(lmp->atom->natoms,nlocal);

    lmp->input->one("velocity all set 0.0 0.0 0.0");
    lmp->input->one("run 0 post no");
    for (int i=0; i < nlocal; ++i) {
        EXPECT_DOUBLE_EQ(f[i][0],0.0);
        EXPECT_DOUBLE_EQ(f[i][1],0.0);
        EXPECT_DOUBLE_EQ(f[i][2],0.0);
        EXPECT_DOUBLE_EQ(v[i][0],0.0);
        EXPECT_DOUBLE_EQ(v[i][1],0.0);
        EXPECT_DOUBLE_EQ(v[i][2],0.0);
    }
};

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
