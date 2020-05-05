
#include "lammps.h"
#include <mpi.h>
#include <cstdio>

#include "gtest/gtest.h"

namespace LAMMPS_NS 
{
    // test fixture for regular tests
    class LAMMPS_test_plain : public ::testing::Test {
    protected:
        LAMMPS *lmp;
        LAMMPS_test_plain() : lmp(nullptr) {
            const char *args[] = {"LAMMPS_test"};
            char **argv = (char **)args;
            int argc = sizeof(args)/sizeof(char *);

            int flag;
            MPI_Initialized(&flag);
            if (!flag) MPI_Init(&argc,&argv);
        }

        ~LAMMPS_test_plain() override {
            lmp = nullptr;
        }

        void SetUp() override {
            const char *args[] = {"LAMMPS_test",
                                  "-log", "none",
                                  "-echo", "screen"};
            char **argv = (char **)args;
            int argc = sizeof(args)/sizeof(char *);
            lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        }

        void TearDown() override {
            delete lmp;
        }
    };

    TEST_F(LAMMPS_test_plain, InitMembers) 
    {
        EXPECT_NE(lmp->memory, nullptr);
        EXPECT_NE(lmp->error, nullptr);
        EXPECT_NE(lmp->universe, nullptr);
        EXPECT_NE(lmp->input, nullptr);

        EXPECT_NE(lmp->atom, nullptr);
        EXPECT_NE(lmp->update, nullptr);
        EXPECT_NE(lmp->neighbor, nullptr);
        EXPECT_NE(lmp->comm, nullptr);
        EXPECT_NE(lmp->domain, nullptr);
        EXPECT_NE(lmp->force, nullptr);
        EXPECT_NE(lmp->modify, nullptr);
        EXPECT_NE(lmp->group, nullptr);
        EXPECT_NE(lmp->output, nullptr);
        EXPECT_NE(lmp->timer, nullptr);

        EXPECT_EQ(lmp->world, MPI_COMM_WORLD);
        EXPECT_EQ(lmp->infile, stdin);
        EXPECT_EQ(lmp->screen, stdout);
        EXPECT_EQ(lmp->logfile, nullptr);
        EXPECT_NE(lmp->initclock, 0.0);

        EXPECT_EQ(lmp->suffix_enable, 0);
        EXPECT_EQ(lmp->suffix, nullptr);
        EXPECT_EQ(lmp->suffix2, nullptr);

        EXPECT_STREQ(lmp->exename, "LAMMPS_test");
    }

    // test fixture for OpenMP with 2 threads
    class LAMMPS_test_omp : public ::testing::Test {
    protected:
        LAMMPS *lmp;
        LAMMPS_test_omp() : lmp(nullptr) {
            const char *args[] = {"LAMMPS_test"};
            char **argv = (char **)args;
            int argc = sizeof(args)/sizeof(char *);

            int flag;
            MPI_Initialized(&flag);
            if (!flag) MPI_Init(&argc,&argv);
        }

        ~LAMMPS_test_omp() override {
            lmp = nullptr;
        }

        void SetUp() override {
            const char *args[] = {"LAMMPS_test",
                                  "-log", "none",
                                  "-echo", "screen",
                                  "-pk", "omp","2", "neigh", "yes",
                                  "-sf", "omp"
            };
            char **argv = (char **)args;
            int argc = sizeof(args)/sizeof(char *);

            // only run this test fixture with omp suffix if USER-OMP is installed
            
            if (LAMMPS::is_installed_pkg("USER-OMP"))
                lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        }

        void TearDown() override {
            delete lmp;
        }
    };

    TEST_F(LAMMPS_test_omp, InitMembers) 
    {
        // skip tests if base class is not available
        if (lmp == nullptr) return;
        EXPECT_NE(lmp->memory, nullptr);
        EXPECT_NE(lmp->error, nullptr);
        EXPECT_NE(lmp->universe, nullptr);
        EXPECT_NE(lmp->input, nullptr);

        EXPECT_NE(lmp->atom, nullptr);
        EXPECT_NE(lmp->update, nullptr);
        EXPECT_NE(lmp->neighbor, nullptr);
        EXPECT_NE(lmp->comm, nullptr);
        EXPECT_NE(lmp->domain, nullptr);
        EXPECT_NE(lmp->force, nullptr);
        EXPECT_NE(lmp->modify, nullptr);
        EXPECT_NE(lmp->group, nullptr);
        EXPECT_NE(lmp->output, nullptr);
        EXPECT_NE(lmp->timer, nullptr);

        EXPECT_EQ(lmp->world, MPI_COMM_WORLD);
        EXPECT_EQ(lmp->infile, stdin);
        EXPECT_EQ(lmp->screen, stdout);
        EXPECT_EQ(lmp->logfile, nullptr);
        EXPECT_GT(lmp->initclock, 0.0);

        EXPECT_EQ(lmp->suffix_enable, 1);
        EXPECT_STREQ(lmp->suffix, "omp");
        EXPECT_EQ(lmp->suffix2, nullptr);

        EXPECT_STREQ(lmp->exename, "LAMMPS_test");
    }

    // test fixture for Kokkos tests
    class LAMMPS_test_kokkos : public ::testing::Test {
    protected:
        LAMMPS *lmp;
        LAMMPS_test_kokkos() : lmp(nullptr) {
            const char *args[] = {"LAMMPS_test"};
            char **argv = (char **)args;
            int argc = sizeof(args)/sizeof(char *);

            int flag;
            MPI_Initialized(&flag);
            if (!flag) MPI_Init(&argc,&argv);
        }

        ~LAMMPS_test_kokkos() override {
            lmp = nullptr;
        }

        void SetUp() override {
            const char *args[] = {"LAMMPS_test",
                                  "-log", "none",
                                  "-echo", "screen",
                                  "-k", "on","t", "2",
                                  "-sf", "kk"
            };
            char **argv = (char **)args;
            int argc = sizeof(args)/sizeof(char *);

            // only run this test fixture with omp suffix if USER-OMP is installed

            if (LAMMPS::is_installed_pkg("KOKKOS"))
                lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        }

        void TearDown() override {
            delete lmp;
        }
    };

    TEST_F(LAMMPS_test_kokkos, InitMembers)
    {
        // skip tests if base class is not available
        if (lmp == nullptr) return;
        EXPECT_NE(lmp->memory, nullptr);
        EXPECT_NE(lmp->error, nullptr);
        EXPECT_NE(lmp->universe, nullptr);
        EXPECT_NE(lmp->input, nullptr);

        EXPECT_NE(lmp->atom, nullptr);
        EXPECT_NE(lmp->update, nullptr);
        EXPECT_NE(lmp->neighbor, nullptr);
        EXPECT_NE(lmp->comm, nullptr);
        EXPECT_NE(lmp->domain, nullptr);
        EXPECT_NE(lmp->force, nullptr);
        EXPECT_NE(lmp->modify, nullptr);
        EXPECT_NE(lmp->group, nullptr);
        EXPECT_NE(lmp->output, nullptr);
        EXPECT_NE(lmp->timer, nullptr);

        EXPECT_EQ(lmp->world, MPI_COMM_WORLD);
        EXPECT_EQ(lmp->infile, stdin);
        EXPECT_EQ(lmp->screen, stdout);
        EXPECT_EQ(lmp->logfile, nullptr);
        EXPECT_GT(lmp->initclock, 0.0);

        EXPECT_EQ(lmp->suffix_enable, 1);
        EXPECT_STREQ(lmp->suffix, "kk");
        EXPECT_EQ(lmp->suffix2, nullptr);

        EXPECT_STREQ(lmp->exename, "LAMMPS_test");
    }
}

int main(int argc, char **argv) 
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
