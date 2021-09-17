// unit tests for the LAMMPS base class

#include "comm.h"
#include "info.h"
#include "lammps.h"
#include <cstdio>  // for stdin, stdout
#include <cstdlib> // for setenv
#include <mpi.h>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::MatchesRegex;
using ::testing::StartsWith;

namespace LAMMPS_NS {
// test fixture for regular tests
class LAMMPS_plain : public ::testing::Test {
protected:
    LAMMPS *lmp;
    LAMMPS_plain() : lmp(nullptr)
    {
        const char *args[] = {"LAMMPS_test"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);

        int flag;
        MPI_Initialized(&flag);
        if (!flag) MPI_Init(&argc, &argv);
    }

    ~LAMMPS_plain() override { lmp = nullptr; }

    void SetUp() override
    {
        const char *args[] = {"LAMMPS_test", "-log", "none", "-echo", "both", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);

        ::testing::internal::CaptureStdout();
        lmp                = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, StartsWith("LAMMPS ("));
    }

    void TearDown() override
    {
        ::testing::internal::CaptureStdout();
        delete lmp;
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, StartsWith("Total wall time:"));
    }
};

TEST_F(LAMMPS_plain, InitMembers)
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
    EXPECT_GE(lmp->initclock, 0.0);

    EXPECT_EQ(lmp->suffix_enable, 0);
    EXPECT_EQ(lmp->suffix, nullptr);
    EXPECT_EQ(lmp->suffix2, nullptr);

    EXPECT_STREQ(lmp->exename, "LAMMPS_test");
    EXPECT_EQ(lmp->num_package, 0);
    EXPECT_EQ(lmp->clientserver, 0);

    EXPECT_EQ(lmp->kokkos, nullptr);
    EXPECT_EQ(lmp->atomKK, nullptr);
    EXPECT_EQ(lmp->memoryKK, nullptr);
    EXPECT_NE(lmp->python, nullptr);
    EXPECT_EQ(lmp->citeme, nullptr);
    if (LAMMPS::has_git_info) {
        EXPECT_STRNE(LAMMPS::git_commit, "");
        EXPECT_STRNE(LAMMPS::git_branch, "");
        EXPECT_STRNE(LAMMPS::git_descriptor, "");
    } else {
        EXPECT_STREQ(LAMMPS::git_commit, "(unknown)");
        EXPECT_STREQ(LAMMPS::git_branch, "(unknown)");
        EXPECT_STREQ(LAMMPS::git_descriptor, "(unknown)");
    }
}

TEST_F(LAMMPS_plain, TestStyles)
{
    // skip tests if base class is not available
    if (lmp == nullptr) return;
    const char *found;

    const char *atom_styles[] = {"atomic", "body",   "charge", "ellipsoid", "hybrid",
                                 "line",   "sphere", "tri",    NULL};
    for (int i = 0; atom_styles[i] != NULL; ++i) {
        found = lmp->match_style("atom", atom_styles[i]);
        EXPECT_STREQ(found, NULL);
    }

    const char *molecule_atom_styles[] = {"angle", "bond", "full", "molecular", "template", NULL};
    for (int i = 0; molecule_atom_styles[i] != NULL; ++i) {
        found = lmp->match_style("atom", molecule_atom_styles[i]);
        EXPECT_STREQ(found, "MOLECULE");
    }

    const char *kokkos_atom_styles[] = {"angle/kk",     "bond/kk",   "full/kk",
                                        "molecular/kk", "hybrid/kk", NULL};
    for (int i = 0; kokkos_atom_styles[i] != NULL; ++i) {
        found = lmp->match_style("atom", kokkos_atom_styles[i]);
        EXPECT_STREQ(found, "KOKKOS");
    }
    found = lmp->match_style("atom", "dipole");
    EXPECT_STREQ(found, "DIPOLE");
    found = lmp->match_style("atom", "peri");
    EXPECT_STREQ(found, "PERI");
    found = lmp->match_style("atom", "spin");
    EXPECT_STREQ(found, "SPIN");
    found = lmp->match_style("atom", "wavepacket");
    EXPECT_STREQ(found, "AWPMD");
    found = lmp->match_style("atom", "dpd");
    EXPECT_STREQ(found, "DPD-REACT");
    found = lmp->match_style("atom", "edpd");
    EXPECT_STREQ(found, "DPD-MESO");
    found = lmp->match_style("atom", "mdpd");
    EXPECT_STREQ(found, "DPD-MESO");
    found = lmp->match_style("atom", "tdpd");
    EXPECT_STREQ(found, "DPD-MESO");
    found = lmp->match_style("atom", "spin");
    EXPECT_STREQ(found, "SPIN");
    found = lmp->match_style("atom", "smd");
    EXPECT_STREQ(found, "MACHDYN");
    found = lmp->match_style("atom", "sph");
    EXPECT_STREQ(found, "SPH");
    found = lmp->match_style("atom", "i_don't_exist");
    EXPECT_STREQ(found, NULL);
}

// test fixture for OpenMP with 2 threads
class LAMMPS_omp : public ::testing::Test {
protected:
    LAMMPS *lmp;
    LAMMPS_omp() : lmp(nullptr)
    {
        const char *args[] = {"LAMMPS_test"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);

        int flag;
        MPI_Initialized(&flag);
        if (!flag) MPI_Init(&argc, &argv);
    }

    ~LAMMPS_omp() override { lmp = nullptr; }

    void SetUp() override
    {
        const char *args[] = {"LAMMPS_test", "-log", "none", "-screen", "none", "-echo", "screen",
                              "-pk",         "omp",  "2",    "neigh",   "yes",  "-sf",   "omp"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);

        // only run this test fixture with omp suffix if OPENMP package is installed

        if (LAMMPS::is_installed_pkg("OPENMP"))
            lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        else
            GTEST_SKIP();
    }

    void TearDown() override { delete lmp; }
};

TEST_F(LAMMPS_omp, InitMembers)
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
    EXPECT_EQ(lmp->screen, nullptr);
    EXPECT_EQ(lmp->logfile, nullptr);
    EXPECT_GE(lmp->initclock, 0.0);

    EXPECT_EQ(lmp->suffix_enable, 1);
    EXPECT_STREQ(lmp->suffix, "omp");
    EXPECT_EQ(lmp->suffix2, nullptr);

    EXPECT_STREQ(lmp->exename, "LAMMPS_test");
    EXPECT_EQ(lmp->num_package, 1);
    EXPECT_EQ(lmp->clientserver, 0);

    EXPECT_EQ(lmp->kokkos, nullptr);
    EXPECT_EQ(lmp->atomKK, nullptr);
    EXPECT_EQ(lmp->memoryKK, nullptr);
    EXPECT_NE(lmp->python, nullptr);
    EXPECT_NE(lmp->citeme, nullptr);
    if (LAMMPS::has_git_info) {
        EXPECT_STRNE(LAMMPS::git_commit, "");
        EXPECT_STRNE(LAMMPS::git_branch, "");
        EXPECT_STRNE(LAMMPS::git_descriptor, "");
    } else {
        EXPECT_STREQ(LAMMPS::git_commit, "(unknown)");
        EXPECT_STREQ(LAMMPS::git_branch, "(unknown)");
        EXPECT_STREQ(LAMMPS::git_descriptor, "(unknown)");
    }
}

// test fixture for Kokkos tests
class LAMMPS_kokkos : public ::testing::Test {
protected:
    LAMMPS *lmp;
    LAMMPS_kokkos() : lmp(nullptr)
    {
        const char *args[] = {"LAMMPS_test"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);

        int flag;
        MPI_Initialized(&flag);
        if (!flag) MPI_Init(&argc, &argv);
    }

    ~LAMMPS_kokkos() override { lmp = nullptr; }

    void SetUp() override
    {
        const char *args[] = {"LAMMPS_test", "-log", "none", "-echo", "none", "-screen", "none",
                              "-k",          "on",   "t",    "2",     "-sf",  "kk"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);

        // only run this test fixture with kk suffix if KOKKOS package is installed
        // also need to figure out a way to find which parallelizations are enabled

        if (LAMMPS::is_installed_pkg("KOKKOS")) {
            ::testing::internal::CaptureStdout();
            lmp                = new LAMMPS(argc, argv, MPI_COMM_WORLD);
            std::string output = ::testing::internal::GetCapturedStdout();
            EXPECT_THAT(output, StartsWith("Kokkos::OpenMP::"));
        } else
            GTEST_SKIP();
    }

    void TearDown() override { delete lmp; }
};

TEST_F(LAMMPS_kokkos, InitMembers)
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
    EXPECT_EQ(lmp->screen, nullptr);
    EXPECT_EQ(lmp->logfile, nullptr);
    EXPECT_GE(lmp->initclock, 0.0);

    EXPECT_EQ(lmp->suffix_enable, 1);
    EXPECT_STREQ(lmp->suffix, "kk");
    EXPECT_EQ(lmp->suffix2, nullptr);

    EXPECT_STREQ(lmp->exename, "LAMMPS_test");
    EXPECT_EQ(lmp->num_package, 0);
    EXPECT_EQ(lmp->clientserver, 0);

    EXPECT_NE(lmp->kokkos, nullptr);
    EXPECT_NE(lmp->atomKK, nullptr);
    EXPECT_NE(lmp->memoryKK, nullptr);
    EXPECT_NE(lmp->python, nullptr);
    EXPECT_NE(lmp->citeme, nullptr);
    if (LAMMPS::has_git_info) {
        EXPECT_STRNE(LAMMPS::git_commit, "");
        EXPECT_STRNE(LAMMPS::git_branch, "");
        EXPECT_STRNE(LAMMPS::git_descriptor, "");
    } else {
        EXPECT_STREQ(LAMMPS::git_commit, "(unknown)");
        EXPECT_STREQ(LAMMPS::git_branch, "(unknown)");
        EXPECT_STREQ(LAMMPS::git_descriptor, "(unknown)");
    }
}

TEST(LAMMPS_init, OpenMP)
{
    if (!LAMMPS::is_installed_pkg("OPENMP")) GTEST_SKIP();
    if (Info::get_openmp_info() == "OpenMP not enabled") GTEST_SKIP();

    FILE *fp = fopen("in.lammps_empty", "w");
    fputs("\n", fp);
    fclose(fp);

    const char *args[] = {"LAMMPS_init", "-in", "in.lammps_empty", "-log", "none", "-nocite"};
    char **argv        = (char **)args;
    int argc           = sizeof(args) / sizeof(char *);

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp        = new LAMMPS(argc, argv, MPI_COMM_WORLD);
    std::string output = ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, MatchesRegex(".*using 2 OpenMP thread.*per MPI task.*"));

    if (LAMMPS_NS::Info::has_accelerator_feature("OPENMP", "api", "openmp"))
        EXPECT_EQ(lmp->comm->nthreads, 2);
    else
        EXPECT_EQ(lmp->comm->nthreads, 1);
    ::testing::internal::CaptureStdout();
    delete lmp;
    ::testing::internal::GetCapturedStdout();

    remove("in.lammps_empty");
}

// check no OMP_NUM_THREADS warning message printing. this must be the
// last OpenMP related test as threads will be locked to 1 from here on.

TEST(LAMMPS_init, NoOpenMP)
{
    if (!LAMMPS_NS::Info::has_accelerator_feature("OPENMP", "api", "openmp"))
        GTEST_SKIP() << "No threading enabled";

    FILE *fp = fopen("in.lammps_class_noomp", "w");
    fputs("\n", fp);
    fclose(fp);
#if defined(__WIN32)
    _putenv("OMP_NUM_THREADS");
#else
    unsetenv("OMP_NUM_THREADS");
#endif

    const char *args[] = {"LAMMPS_init", "-in", "in.lammps_class_noomp", "-log", "none", "-nocite"};
    char **argv        = (char **)args;
    int argc           = sizeof(args) / sizeof(char *);

    ::testing::internal::CaptureStdout();
    LAMMPS *lmp        = new LAMMPS(argc, argv, MPI_COMM_WORLD);
    std::string output = ::testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,
                MatchesRegex(".*OMP_NUM_THREADS environment is not set.*Defaulting to 1 thread.*"));
    EXPECT_EQ(lmp->comm->nthreads, 1);
    ::testing::internal::CaptureStdout();
    delete lmp;
    ::testing::internal::GetCapturedStdout();
}

} // namespace LAMMPS_NS
