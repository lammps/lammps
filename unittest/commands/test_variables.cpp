/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "lammps.h"

#include "atom.h"
#include "domain.h"
#include "group.h"
#include "info.h"
#include "input.h"
#include "region.h"
#include "variable.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstring>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

#if defined(OMPI_MAJOR_VERSION)
const bool have_openmpi = true;
#else
const bool have_openmpi = false;
#endif

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::ExitedWithCode;
using ::testing::MatchesRegex;
using ::testing::StrEq;

#define TEST_FAILURE(errmsg, ...)                                 \
    if (Info::has_exceptions()) {                                 \
        ::testing::internal::CaptureStdout();                     \
        ASSERT_ANY_THROW({__VA_ARGS__});                          \
        auto mesg = ::testing::internal::GetCapturedStdout();     \
        if (verbose) std::cout << mesg;                           \
        ASSERT_THAT(mesg, MatchesRegex(errmsg));                  \
    } else {                                                      \
        if (!have_openmpi) {                                      \
            ::testing::internal::CaptureStdout();                 \
            ASSERT_DEATH({__VA_ARGS__}, "");                      \
            auto mesg = ::testing::internal::GetCapturedStdout(); \
            if (verbose) std::cout << mesg;                       \
            ASSERT_THAT(mesg, MatchesRegex(errmsg));              \
        }                                                         \
    }

class VariableTest : public ::testing::Test {
protected:
    LAMMPS *lmp;
    Group *group;
    Domain *domain;
    Variable *variable;

    void SetUp() override
    {
        const char *args[] = {"VariableTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        group = lmp->group;
        domain = lmp->domain;
        variable = lmp->input->variable;
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
        std::cout.flush();
    }

    void command(const std::string &cmd) { lmp->input->one(cmd); }

    void atomic_system()
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        command("units real");
        command("lattice sc 1.0 origin 0.125 0.125 0.125");
        command("region box block -2 2 -2 2 -2 2");
        command("create_box 8 box");
        command("create_atoms 1 box");
        command("mass * 1.0");
        command("region left block -2.0 -1.0 INF INF INF INF");
        command("region right block 0.5  2.0 INF INF INF INF");
        command("region top block INF INF -2.0 -1.0 INF INF");
        command("set region left type 2");
        command("set region right type 3");
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void molecular_system()
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        command("fix props all property/atom mol rmass q");
        if (!verbose) ::testing::internal::GetCapturedStdout();
        atomic_system();
        if (!verbose) ::testing::internal::CaptureStdout();
        command("variable molid atom floor(id/4)+1");
        command("variable charge atom 2.0*sin(PI/32*id)");
        command("set atom * mol v_molid");
        command("set atom * charge v_charge");
        command("set type 1 mass 0.5");
        command("set type 2*4 mass 2.0");
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(VariableTest, NoBox)
{
    ASSERT_EQ(variable->nvar, 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    command("variable one   index  1");
    command("variable two   equal  2");
    command("variable three string three");
    command("variable four  loop   4");
    command("variable five  loop   100 pad");
    command("variable six   world  one");
    command("variable seven format two \"%5.2f\"");
    command("variable eight getenv HOME");
    command("variable dummy index  0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(variable->nvar, 9);
    if (!verbose) ::testing::internal::CaptureStdout();
    command("variable dummy delete");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(variable->nvar, 8);
    
    TEST_FAILURE(".*ERROR: Illegal variable command.*", command("variable"););
}

TEST_F(VariableTest, AtomicSystem)
{
    atomic_system();

    if (!verbose) ::testing::internal::CaptureStdout();
    command("variable id atom id");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(variable->nvar, 1);
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (have_openmpi && !LAMMPS_NS::Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. "
                     "Death tests will be skipped\n";

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
