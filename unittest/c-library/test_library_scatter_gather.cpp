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
                              "-var",        "input_dir", STRINGIFY(TEST_INPUT_FOLDER),
                              nullptr};

        char **argv = (char **)args;
        int argc    = (sizeof(args) / sizeof(char *)) - 1;

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
#undef CHECK_BOND
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
#undef CHECK_BOND
};

TEST_F(GatherProperties, gather_angles_newton_on)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = path_join(INPUT_DIR, "in.fourmol");
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "newton on on");
    lammps_file(lmp, input.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    bigint nangles = *(bigint *)lammps_extract_global(lmp, "nangles");
    EXPECT_EQ(nangles, 30);

    tagint *angles = new tagint[4 * nangles];
    lammps_gather_angles(lmp, angles);

#define CHECK_ANGLE(idx, type, atom1, atom2, atom3)                          \
    if (((angles[4 * idx + 1] == atom1) && (angles[4 * idx + 2] == atom2) && \
         (angles[4 * idx + 3] == atom3)) ||                                  \
        ((angles[4 * idx + 1] == atom3) && (angles[4 * idx + 2] == atom2) && \
         (angles[4 * idx + 3] == atom1))) {                                  \
        EXPECT_EQ(angles[4 * idx], type);                                    \
        ++count;                                                             \
    }

    // check validity of a few angles by comparing the angle type and counting the matches.
    int count = 0;
    for (bigint i = 0; i < nangles; ++i) {
        CHECK_ANGLE(i, 4, 2, 1, 3);
        CHECK_ANGLE(i, 4, 1, 3, 5);
        CHECK_ANGLE(i, 4, 1, 3, 4);
        CHECK_ANGLE(i, 4, 13, 12, 15);
        CHECK_ANGLE(i, 4, 13, 12, 14);
        CHECK_ANGLE(i, 2, 5, 3, 6);
        CHECK_ANGLE(i, 2, 4, 3, 6);
        CHECK_ANGLE(i, 3, 3, 6, 7);
        CHECK_ANGLE(i, 3, 3, 6, 8);
        CHECK_ANGLE(i, 1, 22, 21, 23);
    }
    EXPECT_EQ(count, 10);
    delete[] angles;
#undef CHECK_ANGLE
};

TEST_F(GatherProperties, gather_angles_newton_off)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = path_join(INPUT_DIR, "in.fourmol");
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "newton off off");
    lammps_file(lmp, input.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    bigint nangles = *(bigint *)lammps_extract_global(lmp, "nangles");
    EXPECT_EQ(nangles, 30);

    tagint *angles = new tagint[4 * nangles];
    lammps_gather_angles(lmp, angles);

#define CHECK_ANGLE(idx, type, atom1, atom2, atom3)                          \
    if (((angles[4 * idx + 1] == atom1) && (angles[4 * idx + 2] == atom2) && \
         (angles[4 * idx + 3] == atom3)) ||                                  \
        ((angles[4 * idx + 1] == atom3) && (angles[4 * idx + 2] == atom2) && \
         (angles[4 * idx + 3] == atom1))) {                                  \
        EXPECT_EQ(angles[4 * idx], type);                                    \
        ++count;                                                             \
    }

    // check validity of a few angles by comparing the angle type and counting the matches.
    int count = 0;
    for (bigint i = 0; i < nangles; ++i) {
        CHECK_ANGLE(i, 4, 2, 1, 3);
        CHECK_ANGLE(i, 4, 1, 3, 5);
        CHECK_ANGLE(i, 4, 1, 3, 4);
        CHECK_ANGLE(i, 4, 13, 12, 15);
        CHECK_ANGLE(i, 4, 13, 12, 14);
        CHECK_ANGLE(i, 2, 5, 3, 6);
        CHECK_ANGLE(i, 2, 4, 3, 6);
        CHECK_ANGLE(i, 3, 3, 6, 7);
        CHECK_ANGLE(i, 3, 3, 6, 8);
        CHECK_ANGLE(i, 1, 22, 21, 23);
    }
    EXPECT_EQ(count, 10);
    delete[] angles;
#undef CHECK_ANGLES
};

TEST_F(GatherProperties, gather_dihedrals_newton_on)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = path_join(INPUT_DIR, "in.fourmol");
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "newton on on");
    lammps_file(lmp, input.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    bigint ndihedrals = *(bigint *)lammps_extract_global(lmp, "ndihedrals");
    EXPECT_EQ(ndihedrals, 31);

    tagint *dihedrals = new tagint[5 * ndihedrals];
    lammps_gather_dihedrals(lmp, dihedrals);

#define CHECK_DIHEDRAL(idx, type, atom1, atom2, atom3, atom4)                     \
    if ((dihedrals[5 * idx + 1] == atom1) && (dihedrals[5 * idx + 2] == atom2) && \
        (dihedrals[5 * idx + 3] == atom3) && (dihedrals[5 * idx + 4] == atom4)) { \
        EXPECT_EQ(dihedrals[5 * idx], type);                                      \
        ++count;                                                                  \
    }

    // check validity of a few dihedrals by comparing the dihedral type and counting the matches.
    int count = 0;
    for (bigint i = 0; i < ndihedrals; ++i) {
        CHECK_DIHEDRAL(i, 2, 2, 1, 3, 6);
        CHECK_DIHEDRAL(i, 2, 2, 1, 3, 4);
        CHECK_DIHEDRAL(i, 3, 2, 1, 3, 5);
        CHECK_DIHEDRAL(i, 1, 1, 3, 6, 8);
        CHECK_DIHEDRAL(i, 1, 1, 3, 6, 7);
        CHECK_DIHEDRAL(i, 5, 4, 3, 6, 8);
        CHECK_DIHEDRAL(i, 5, 4, 3, 6, 7);
        CHECK_DIHEDRAL(i, 5, 16, 10, 12, 13);
        CHECK_DIHEDRAL(i, 5, 16, 10, 12, 14);
        CHECK_DIHEDRAL(i, 5, 16, 10, 12, 15);
    }
    EXPECT_EQ(count, 10);
    delete[] dihedrals;
#undef CHECK_DIHEDRAL
};

TEST_F(GatherProperties, gather_dihedrals_newton_off)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = path_join(INPUT_DIR, "in.fourmol");
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "newton off off");
    lammps_file(lmp, input.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    bigint ndihedrals = *(bigint *)lammps_extract_global(lmp, "ndihedrals");
    EXPECT_EQ(ndihedrals, 31);

    tagint *dihedrals = new tagint[5 * ndihedrals];
    lammps_gather_dihedrals(lmp, dihedrals);

#define CHECK_DIHEDRAL(idx, type, atom1, atom2, atom3, atom4)                     \
    if ((dihedrals[5 * idx + 1] == atom1) && (dihedrals[5 * idx + 2] == atom2) && \
        (dihedrals[5 * idx + 3] == atom3) && (dihedrals[5 * idx + 4] == atom4)) { \
        EXPECT_EQ(dihedrals[5 * idx], type);                                      \
        ++count;                                                                  \
    }
    // check validity of a few dihedrals by comparing the dihedral type and counting the matches.
    int count = 0;
    for (bigint i = 0; i < ndihedrals; ++i) {
        CHECK_DIHEDRAL(i, 2, 2, 1, 3, 6);
        CHECK_DIHEDRAL(i, 2, 2, 1, 3, 4);
        CHECK_DIHEDRAL(i, 3, 2, 1, 3, 5);
        CHECK_DIHEDRAL(i, 1, 1, 3, 6, 8);
        CHECK_DIHEDRAL(i, 1, 1, 3, 6, 7);
        CHECK_DIHEDRAL(i, 5, 4, 3, 6, 8);
        CHECK_DIHEDRAL(i, 5, 4, 3, 6, 7);
        CHECK_DIHEDRAL(i, 5, 16, 10, 12, 13);
        CHECK_DIHEDRAL(i, 5, 16, 10, 12, 14);
        CHECK_DIHEDRAL(i, 5, 16, 10, 12, 15);
    }
    EXPECT_EQ(count, 10);
    delete[] dihedrals;
#undef CHECK_DIHEDRALS
};

TEST_F(GatherProperties, gather_impropers_newton_on)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = path_join(INPUT_DIR, "in.fourmol");
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "newton on on");
    lammps_file(lmp, input.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    bigint nimpropers = *(bigint *)lammps_extract_global(lmp, "nimpropers");
    EXPECT_EQ(nimpropers, 2);

    tagint *impropers = new tagint[5 * nimpropers];
    lammps_gather_impropers(lmp, impropers);

#define CHECK_IMPROPER(idx, type, atom1, atom2, atom3, atom4)                     \
    if ((impropers[5 * idx + 1] == atom1) && (impropers[5 * idx + 2] == atom2) && \
        (impropers[5 * idx + 3] == atom3) && (impropers[5 * idx + 4] == atom4)) { \
        EXPECT_EQ(impropers[5 * idx], type);                                      \
        ++count;                                                                  \
    }

    // check validity of a few impropers by comparing the improper type and counting the matches.
    int count = 0;
    for (bigint i = 0; i < nimpropers; ++i) {
        CHECK_IMPROPER(i, 1, 6, 3, 8, 7);
        CHECK_IMPROPER(i, 2, 8, 6, 10, 9);
    }
    EXPECT_EQ(count, 2);
    delete[] impropers;
#undef CHECK_IMPROPER
};

TEST_F(GatherProperties, gather_impropers_newton_off)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = path_join(INPUT_DIR, "in.fourmol");
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "newton off off");
    lammps_file(lmp, input.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();

    bigint nimpropers = *(bigint *)lammps_extract_global(lmp, "nimpropers");
    EXPECT_EQ(nimpropers, 2);

    tagint *impropers = new tagint[5 * nimpropers];
    lammps_gather_impropers(lmp, impropers);

#define CHECK_IMPROPER(idx, type, atom1, atom2, atom3, atom4)                     \
    if ((impropers[5 * idx + 1] == atom1) && (impropers[5 * idx + 2] == atom2) && \
        (impropers[5 * idx + 3] == atom3) && (impropers[5 * idx + 4] == atom4)) { \
        EXPECT_EQ(impropers[5 * idx], type);                                      \
        ++count;                                                                  \
    }
    // check validity of a few impropers by comparing the improper type and counting the matches.
    int count = 0;
    for (bigint i = 0; i < nimpropers; ++i) {
        CHECK_IMPROPER(i, 1, 6, 3, 8, 7);
        CHECK_IMPROPER(i, 2, 8, 6, 10, 9);
    }
    EXPECT_EQ(count, 2);
    delete[] impropers;
#undef CHECK_IMPROPERS
};
