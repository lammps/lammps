// unit tests for testing scatter/gather operations through the library interface

#include "lammps.h"
#include "library.h"
#include "lmptype.h"
#include "platform.h"
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "test_main.h"

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

using ::LAMMPS_NS::bigint;
using ::LAMMPS_NS::tagint;
using ::LAMMPS_NS::platform::path_join;
using ::testing::HasSubstr;
using ::testing::StartsWith;

class GatherProperties : public ::testing::Test {
protected:
    void *lmp;
    std::string INPUT_DIR = STRINGIFY(TEST_INPUT_FOLDER);

    GatherProperties()           = default;
    ~GatherProperties() override = default;

    void SetUp() override
    {
        const char *args[] = {"LAMMPS_test", "-log",      "none",
                              "-echo",       "screen",    "-nocite",
                              "-var",        "input_dir", STRINGIFY(TEST_INPUT_FOLDER)};

        char **argv = (char **)args;
        int argc    = sizeof(args) / sizeof(char *);

        ::testing::internal::CaptureStdout();
        lmp                = lammps_open_no_mpi(argc, argv, nullptr);
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        EXPECT_THAT(output, StartsWith("LAMMPS ("));
    }
    void TearDown() override
    {
        ::testing::internal::CaptureStdout();
        lammps_close(lmp);
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, HasSubstr("Total wall time:"));
        if (verbose) std::cout << output;
        lmp = nullptr;
    }
};

TEST_F(GatherProperties, gather_bonds_newton_on)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = path_join(INPUT_DIR, "in.fourmol");
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "newton on on");
    lammps_file(lmp, input.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    bigint nbonds = *(bigint *)lammps_extract_global(lmp, "nbonds");
    EXPECT_EQ(nbonds, 24);

    tagint *bonds = new tagint[3 * nbonds];
    lammps_gather_bonds(lmp, bonds);

#define CHECK_BOND(idx, type, atom1, atom2)                                 \
    if (((bonds[3 * idx + 1] == atom1) && (bonds[3 * idx + 2] == atom2)) || \
        ((bonds[3 * idx + 1] == atom2) && (bonds[3 * idx + 2] == atom1))) { \
        EXPECT_EQ(bonds[3 * idx], type);                                    \
        ++count;                                                            \
    }

    // check validity of a few bonds by comparing the bond type and counting the matches.
    int count = 0;
    for (bigint i = 0; i < nbonds; ++i) {
        CHECK_BOND(i, 5, 1, 2);
        CHECK_BOND(i, 3, 1, 3);
        CHECK_BOND(i, 2, 3, 4);
        CHECK_BOND(i, 2, 3, 5);
        CHECK_BOND(i, 1, 3, 6);
        CHECK_BOND(i, 3, 6, 8);
        CHECK_BOND(i, 4, 6, 7);
        CHECK_BOND(i, 5, 8, 9);
        CHECK_BOND(i, 5, 27, 28);
        CHECK_BOND(i, 5, 27, 29);
    }
    EXPECT_EQ(count, 10);
    delete[] bonds;
};

TEST_F(GatherProperties, gather_bonds_newton_off)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = path_join(INPUT_DIR, "in.fourmol");
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "newton off off");
    lammps_file(lmp, input.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    bigint nbonds = *(bigint *)lammps_extract_global(lmp, "nbonds");
    EXPECT_EQ(nbonds, 24);

    tagint *bonds = new tagint[3 * nbonds];
    lammps_gather_bonds(lmp, bonds);

#define CHECK_BOND(idx, type, atom1, atom2)                                 \
    if (((bonds[3 * idx + 1] == atom1) && (bonds[3 * idx + 2] == atom2)) || \
        ((bonds[3 * idx + 1] == atom2) && (bonds[3 * idx + 2] == atom1))) { \
        EXPECT_EQ(bonds[3 * idx], type);                                    \
        ++count;                                                            \
    }

    // check validity of a few bonds by comparing the bond type and counting the matches.
    int count = 0;
    for (bigint i = 0; i < nbonds; ++i) {
        CHECK_BOND(i, 5, 1, 2);
        CHECK_BOND(i, 3, 1, 3);
        CHECK_BOND(i, 2, 3, 4);
        CHECK_BOND(i, 2, 3, 5);
        CHECK_BOND(i, 1, 3, 6);
        CHECK_BOND(i, 3, 6, 8);
        CHECK_BOND(i, 4, 6, 7);
        CHECK_BOND(i, 5, 8, 9);
        CHECK_BOND(i, 5, 27, 28);
        CHECK_BOND(i, 5, 27, 29);
    }
    EXPECT_EQ(count, 10);
    delete[] bonds;
};
