// unit tests creating LAMMPS instances via the library interface

#include "lammps.h"
#define LAMMPS_LIB_MPI 1
#include "library.h"
#include <cstdio> // for stdin, stdout
#include <mpi.h>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "test_main.h"

using ::testing::HasSubstr;
using ::testing::StartsWith;

TEST(lammps_open, null_args)
{
    ::testing::internal::CaptureStdout();
    void *handle       = lammps_open(0, nullptr, MPI_COMM_WORLD, nullptr);
    std::string output = ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    if (verbose) std::cout << output;
    int mpi_init = 0;
    MPI_Initialized(&mpi_init);
    EXPECT_GT(mpi_init, 0);
    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *)handle;
    EXPECT_EQ(lmp->world, MPI_COMM_WORLD);
    EXPECT_EQ(lmp->infile, stdin);
    EXPECT_EQ(lmp->screen, stdout);
    EXPECT_NE(lmp->citeme, nullptr);
    ::testing::internal::CaptureStdout();
    lammps_close(handle);
    output = ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, HasSubstr("Total wall time:"));
    if (verbose) std::cout << output;
}

TEST(lammps_open, with_args)
{
    const char *args[] = {"liblammps", "-log", "none", "-nocite", nullptr};
    char **argv        = (char **)args;
    int argc           = (sizeof(args) / sizeof(char *)) - 1;

    // MPI is already initialized
    MPI_Comm mycomm;
    MPI_Comm_split(MPI_COMM_WORLD, 0, 1, &mycomm);
    ::testing::internal::CaptureStdout();
    void *alt_ptr;
    void *handle       = lammps_open(argc, argv, mycomm, &alt_ptr);
    std::string output = ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    if (verbose) std::cout << output;
    EXPECT_EQ(handle, alt_ptr);
    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *)handle;

    // MPI STUBS uses no real communicators
#if !defined(MPI_STUBS)
    EXPECT_NE(lmp->world, MPI_COMM_WORLD);
#endif

    EXPECT_EQ(lmp->world, mycomm);
    EXPECT_EQ(lmp->infile, stdin);
    EXPECT_EQ(lmp->logfile, nullptr);
    EXPECT_EQ(lmp->citeme, nullptr);
    EXPECT_EQ(lmp->kokkos, nullptr);
    EXPECT_EQ(lmp->atomKK, nullptr);
    EXPECT_EQ(lmp->memoryKK, nullptr);
    ::testing::internal::CaptureStdout();
    lammps_close(handle);
    output = ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, HasSubstr("Total wall time:"));
    if (verbose) std::cout << output;
    MPI_Comm_free(&mycomm);
}

TEST(lammps_open, with_kokkos)
{
    if (!LAMMPS_NS::LAMMPS::is_installed_pkg("KOKKOS")) GTEST_SKIP();
    const char *args[] = {"liblammps", "-k", "on", "t", "2", "-sf", "kk", "-log", "none", nullptr};
    char **argv        = (char **)args;
    int argc           = (sizeof(args) / sizeof(char *)) - 1;

    ::testing::internal::CaptureStdout();
    void *alt_ptr;
    void *handle       = lammps_open(argc, argv, MPI_COMM_WORLD, &alt_ptr);
    std::string output = ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    if (verbose) std::cout << output;
    EXPECT_EQ(handle, alt_ptr);
    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *)handle;

    EXPECT_EQ(lmp->world, MPI_COMM_WORLD);
    EXPECT_EQ(lmp->infile, stdin);
    EXPECT_EQ(lmp->logfile, nullptr);
    EXPECT_NE(lmp->citeme, nullptr);
    EXPECT_EQ(lmp->num_package, 0);
    EXPECT_NE(lmp->kokkos, nullptr);
    EXPECT_NE(lmp->atomKK, nullptr);
    EXPECT_NE(lmp->memoryKK, nullptr);
    ::testing::internal::CaptureStdout();
    lammps_close(handle);
    output = ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, HasSubstr("Total wall time:"));
    if (verbose) std::cout << output;
}

TEST(lammps_open_no_mpi, no_screen)
{
    const char *args[] = {"liblammps", "-log", "none", "-screen", "none", "-nocite", nullptr};
    char **argv        = (char **)args;
    int argc           = (sizeof(args) / sizeof(char *)) - 1;

    ::testing::internal::CaptureStdout();
    void *alt_ptr;
    void *handle       = lammps_open_no_mpi(argc, argv, &alt_ptr);
    std::string output = ::testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.c_str(), "");
    EXPECT_EQ(handle, alt_ptr);
    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *)handle;

    EXPECT_EQ(lmp->world, MPI_COMM_WORLD);
    EXPECT_EQ(lmp->infile, stdin);
    EXPECT_EQ(lmp->screen, nullptr);
    EXPECT_EQ(lmp->logfile, nullptr);
    EXPECT_EQ(lmp->citeme, nullptr);
    EXPECT_EQ(lmp->suffix_enable, 0);

    EXPECT_STREQ(lmp->exename, "liblammps");
    EXPECT_EQ(lmp->num_package, 0);
    ::testing::internal::CaptureStdout();
    lammps_close(handle);
    output = ::testing::internal::GetCapturedStdout();
    EXPECT_STREQ(output.c_str(), "");
}

TEST(lammps_open_no_mpi, with_omp)
{
    if (!LAMMPS_NS::LAMMPS::is_installed_pkg("OPENMP")) GTEST_SKIP();
    const char *args[] = {"liblammps", "-pk", "omp",  "2",    "neigh",  "no",
                          "-sf",       "omp", "-log", "none", "-nocite", nullptr};
    char **argv        = (char **)args;
    int argc           = (sizeof(args) / sizeof(char *)) - 1;

    ::testing::internal::CaptureStdout();
    void *alt_ptr;
    void *handle       = lammps_open_no_mpi(argc, argv, &alt_ptr);
    std::string output = ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    if (verbose) std::cout << output;
    EXPECT_EQ(handle, alt_ptr);
    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *)handle;

    EXPECT_EQ(lmp->world, MPI_COMM_WORLD);
    EXPECT_EQ(lmp->infile, stdin);
    EXPECT_EQ(lmp->logfile, nullptr);
    EXPECT_EQ(lmp->citeme, nullptr);
    EXPECT_EQ(lmp->suffix_enable, 1);
    EXPECT_STREQ(lmp->suffix, "omp");
    EXPECT_EQ(lmp->suffix2, nullptr);
    EXPECT_STREQ(lmp->exename, "liblammps");
    EXPECT_EQ(lmp->num_package, 1);
    ::testing::internal::CaptureStdout();
    lammps_close(handle);
    output = ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, HasSubstr("Total wall time:"));
    if (verbose) std::cout << output;
}

TEST(lammps_open_fortran, no_args)
{
    // MPI is already initialized
    MPI_Comm mycomm;
    MPI_Comm_split(MPI_COMM_WORLD, 0, 1, &mycomm);
    int fcomm = MPI_Comm_c2f(mycomm);
    ::testing::internal::CaptureStdout();
    void *handle       = lammps_open_fortran(0, nullptr, fcomm);
    std::string output = ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, StartsWith("LAMMPS ("));
    if (verbose) std::cout << output;
    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *)handle;

    // MPI STUBS uses no real communicators
#if !defined(MPI_STUBS)
    EXPECT_NE(lmp->world, MPI_COMM_WORLD);
#endif

    EXPECT_EQ(lmp->world, mycomm);
    EXPECT_EQ(lmp->infile, stdin);
    EXPECT_EQ(lmp->screen, stdout);
    EXPECT_NE(lmp->logfile, nullptr);
    EXPECT_NE(lmp->citeme, nullptr);
    ::testing::internal::CaptureStdout();
    lammps_close(handle);
    output = ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, HasSubstr("Total wall time:"));
    if (verbose) std::cout << output;
    MPI_Comm_free(&mycomm);
}

TEST(lammps_open_no_mpi, lammps_error)
{
    const char *args[] = {"liblammps", "-log", "none", "-nocite", nullptr};
    char **argv        = (char **)args;
    int argc           = (sizeof(args) / sizeof(char *)) - 1;

    ::testing::internal::CaptureStdout();
    void *alt_ptr;
    void *handle       = lammps_open_no_mpi(argc, argv, &alt_ptr);
    std::string output = ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(handle, alt_ptr);
    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *)handle;

    EXPECT_EQ(lmp->world, MPI_COMM_WORLD);
    EXPECT_EQ(lmp->infile, stdin);
    EXPECT_NE(lmp->screen, nullptr);
    EXPECT_EQ(lmp->logfile, nullptr);
    EXPECT_EQ(lmp->citeme, nullptr);
    EXPECT_EQ(lmp->suffix_enable, 0);

    EXPECT_STREQ(lmp->exename, "liblammps");
    EXPECT_EQ(lmp->num_package, 0);
    ::testing::internal::CaptureStdout();
    lammps_error(handle, 0, "test_warning");
    output = ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, HasSubstr("WARNING: test_warning"));
}

