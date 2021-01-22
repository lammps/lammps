// unit tests for checking and changing simulation properties through the library interface

#include "lammps.h"
#include "library.h"
#include "lmptype.h"
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "test_main.h"

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

using ::LAMMPS_NS::tagint;
using ::testing::HasSubstr;
using ::testing::StartsWith;
using ::testing::StrEq;

class LibraryProperties : public ::testing::Test {
protected:
    void *lmp;
    std::string INPUT_DIR = STRINGIFY(TEST_INPUT_FOLDER);

    LibraryProperties(){};
    ~LibraryProperties() override{};

    void SetUp() override
    {
        const char *args[] = {"LAMMPS_test", "-log",      "none",
                              "-echo",       "screen",    "-nocite",
                              "-var",        "input_dir", STRINGIFY(TEST_INPUT_FOLDER)};

        char **argv = (char **)args;
        int argc    = sizeof(args) / sizeof(char *);

        ::testing::internal::CaptureStdout();
        lmp                = lammps_open_no_mpi(argc, argv, NULL);
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

TEST_F(LibraryProperties, version)
{
    EXPECT_LT(20200917, lammps_version(lmp));
};

TEST_F(LibraryProperties, memory_usage)
{
    double meminfo[3];
    lammps_memory_usage(lmp, meminfo);
    EXPECT_GT(meminfo[0], 0.0);
#if defined(__linux__) || defined(_WIN32)
    EXPECT_GE(meminfo[1], 0.0);
#endif
    EXPECT_GT(meminfo[2], 0.0);
};

TEST_F(LibraryProperties, get_mpi_comm)
{
    int f_comm = lammps_get_mpi_comm(lmp);
    if (lammps_config_has_mpi_support())
        EXPECT_GE(f_comm, 0);
    else
        EXPECT_EQ(f_comm, -1);
};

TEST_F(LibraryProperties, natoms)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = INPUT_DIR + PATH_SEP + "in.fourmol";
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_file(lmp, input.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_get_natoms(lmp), 29);
};

TEST_F(LibraryProperties, thermo)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = INPUT_DIR + PATH_SEP + "in.fourmol";
    ::testing::internal::CaptureStdout();
    lammps_file(lmp, input.c_str());
    lammps_command(lmp, "run 2 post no");
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    EXPECT_EQ(lammps_get_thermo(lmp, "step"), 2);
    EXPECT_EQ(lammps_get_thermo(lmp, "atoms"), 29);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "vol"), 3375.0);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "density"), 0.12211250945013695);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "cellalpha"), 90.0);
};

TEST_F(LibraryProperties, box)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = INPUT_DIR + PATH_SEP + "in.fourmol";
    ::testing::internal::CaptureStdout();
    lammps_file(lmp, input.c_str());
    lammps_command(lmp, "run 2 post no");
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    double boxlo[3], boxhi[3], xy, yz, xz;
    int pflags[3], boxflag;
    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, pflags, &boxflag);
    EXPECT_DOUBLE_EQ(boxlo[0], -6.024572);
    EXPECT_DOUBLE_EQ(boxlo[1], -7.692866);
    EXPECT_DOUBLE_EQ(boxlo[2], -8.086924);
    EXPECT_DOUBLE_EQ(boxhi[0], 8.975428);
    EXPECT_DOUBLE_EQ(boxhi[1], 7.307134);
    EXPECT_DOUBLE_EQ(boxhi[2], 6.913076);
    EXPECT_DOUBLE_EQ(xy, 0.0);
    EXPECT_DOUBLE_EQ(yz, 0.0);
    EXPECT_DOUBLE_EQ(xz, 0.0);
    EXPECT_EQ(pflags[0], 1);
    EXPECT_EQ(pflags[1], 1);
    EXPECT_EQ(pflags[2], 1);
    EXPECT_EQ(boxflag, 0);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "vol"), 3375.0);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "density"), 0.12211250945013695);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "cellalpha"), 90.0);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "change_box all boundary p p f triclinic xy final 0.5");
    lammps_command(lmp, "fix box all box/relax x 0.0 y 0.0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, pflags, &boxflag);
    EXPECT_DOUBLE_EQ(boxlo[0], -6.024572);
    EXPECT_DOUBLE_EQ(boxlo[1], -7.692866);
    EXPECT_DOUBLE_EQ(boxlo[2], -8.086924);
    EXPECT_DOUBLE_EQ(boxhi[0], 8.975428);
    EXPECT_DOUBLE_EQ(boxhi[1], 7.307134);
    EXPECT_DOUBLE_EQ(boxhi[2], 6.913076);
    EXPECT_DOUBLE_EQ(xy, 0.5);
    EXPECT_DOUBLE_EQ(yz, 0.0);
    EXPECT_DOUBLE_EQ(xz, 0.0);
    EXPECT_EQ(pflags[0], 1);
    EXPECT_EQ(pflags[1], 1);
    EXPECT_EQ(pflags[2], 0);
    EXPECT_EQ(boxflag, 1);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "vol"), 3375.0);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "density"), 0.12211250945013695);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "cellalpha"), 90.0);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "cellbeta"), 90.0);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "cellgamma"), 88.090847567003621);

    boxlo[0] = -6.1;
    boxhi[1] = 7.3;
    xy       = 0.1;
    lammps_reset_box(lmp, boxlo, boxhi, xy, yz, xz);
    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, pflags, &boxflag);
    EXPECT_DOUBLE_EQ(boxlo[0], -6.1);
    EXPECT_DOUBLE_EQ(boxlo[1], -7.692866);
    EXPECT_DOUBLE_EQ(boxlo[2], -8.086924);
    EXPECT_DOUBLE_EQ(boxhi[0], 8.975428);
    EXPECT_DOUBLE_EQ(boxhi[1], 7.3);
    EXPECT_DOUBLE_EQ(boxhi[2], 6.913076);
    EXPECT_DOUBLE_EQ(xy, 0.1);
    EXPECT_DOUBLE_EQ(yz, 0.0);
    EXPECT_DOUBLE_EQ(xz, 0.0);
    EXPECT_EQ(pflags[0], 1);
    EXPECT_EQ(pflags[1], 1);
    EXPECT_EQ(pflags[2], 0);
    EXPECT_EQ(boxflag, 1);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "vol"), 3390.3580784497199);
    EXPECT_DOUBLE_EQ(lammps_get_thermo(lmp, "cellgamma"), 89.61785205109274);
};

TEST_F(LibraryProperties, setting)
{
#if defined(LAMMPS_SMALLSMALL)
    EXPECT_EQ(lammps_extract_setting(lmp, "bigint"), 4);
#else
    EXPECT_EQ(lammps_extract_setting(lmp, "bigint"), 8);
#endif
#if defined(LAMMPS_BIGBIG)
    EXPECT_EQ(lammps_extract_setting(lmp, "tagint"), 8);
    EXPECT_EQ(lammps_extract_setting(lmp, "imageint"), 8);
#else
    EXPECT_EQ(lammps_extract_setting(lmp, "tagint"), 4);
    EXPECT_EQ(lammps_extract_setting(lmp, "imageint"), 4);
#endif

    EXPECT_EQ(lammps_extract_setting(lmp, "box_exist"), 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "dimension 2");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_extract_setting(lmp, "dimension"), 2);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "dimension 3");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    EXPECT_EQ(lammps_extract_setting(lmp, "world_size"), 1);
    EXPECT_EQ(lammps_extract_setting(lmp, "world_rank"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "universe_size"), 1);
    EXPECT_EQ(lammps_extract_setting(lmp, "universe_rank"), 0);
    EXPECT_GT(lammps_extract_setting(lmp, "nthreads"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "newton_pair"), 1);
    EXPECT_EQ(lammps_extract_setting(lmp, "newton_bond"), 1);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "newton off");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_extract_setting(lmp, "newton_pair"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "newton_bond"), 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "newton on off");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_extract_setting(lmp, "newton_pair"), 1);
    EXPECT_EQ(lammps_extract_setting(lmp, "newton_bond"), 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "newton off on");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_extract_setting(lmp, "newton_pair"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "newton_bond"), 1);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "newton on");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_extract_setting(lmp, "newton_pair"), 1);
    EXPECT_EQ(lammps_extract_setting(lmp, "newton_bond"), 1);

    EXPECT_EQ(lammps_extract_setting(lmp, "ntypes"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "nbondtypes"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "nangletypes"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "ndihedraltypes"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "nimpropertypes"), 0);

    EXPECT_EQ(lammps_extract_setting(lmp, "molecule_flag"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "q_flag"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "mu_flag"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "rmass_flag"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "UNKNOWN"), -1);

    if (lammps_has_style(lmp, "atom", "full")) {
        std::string input = INPUT_DIR + PATH_SEP + "in.fourmol";
        if (!verbose) ::testing::internal::CaptureStdout();
        lammps_file(lmp, input.c_str());
        lammps_command(lmp, "run 2 post no");
        if (!verbose) ::testing::internal::GetCapturedStdout();
        EXPECT_EQ(lammps_extract_setting(lmp, "triclinic"), 0);
        EXPECT_EQ(lammps_extract_setting(lmp, "box_exist"), 1);
        EXPECT_EQ(lammps_extract_setting(lmp, "dimension"), 3);
        EXPECT_EQ(lammps_extract_setting(lmp, "nlocal"), 29);
        EXPECT_EQ(lammps_extract_setting(lmp, "nghost"), 518);
        EXPECT_EQ(lammps_extract_setting(lmp, "nall"), 547);
        EXPECT_EQ(lammps_extract_setting(lmp, "nmax"), 16384);
        EXPECT_EQ(lammps_extract_setting(lmp, "ntypes"), 5);
        EXPECT_EQ(lammps_extract_setting(lmp, "nbondtypes"), 5);
        EXPECT_EQ(lammps_extract_setting(lmp, "nangletypes"), 4);
        EXPECT_EQ(lammps_extract_setting(lmp, "ndihedraltypes"), 5);
        EXPECT_EQ(lammps_extract_setting(lmp, "nimpropertypes"), 2);

        EXPECT_EQ(lammps_extract_setting(lmp, "molecule_flag"), 1);
        EXPECT_EQ(lammps_extract_setting(lmp, "q_flag"), 1);
        EXPECT_EQ(lammps_extract_setting(lmp, "mu_flag"), 0);
        EXPECT_EQ(lammps_extract_setting(lmp, "rmass_flag"), 0);
        EXPECT_EQ(lammps_extract_setting(lmp, "radius_flag"), 0);
        EXPECT_EQ(lammps_extract_setting(lmp, "sphere_flag"), 0);
        EXPECT_EQ(lammps_extract_setting(lmp, "ellipsoid_flag"), 0);
        EXPECT_EQ(lammps_extract_setting(lmp, "omega_flag"), 0);
        EXPECT_EQ(lammps_extract_setting(lmp, "torque_flag"), 0);
        EXPECT_EQ(lammps_extract_setting(lmp, "angmom_flag"), 0);
        if (!verbose) ::testing::internal::CaptureStdout();
        lammps_command(lmp, "change_box all triclinic");
        lammps_command(lmp, "fix rmass all property/atom rmass ghost yes");
        if (!verbose) ::testing::internal::GetCapturedStdout();
        EXPECT_EQ(lammps_extract_setting(lmp, "triclinic"), 1);
        EXPECT_EQ(lammps_extract_setting(lmp, "rmass_flag"), 1);
    }
};

TEST_F(LibraryProperties, global)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();

    std::string input = INPUT_DIR + PATH_SEP + "in.fourmol";
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_file(lmp, input.c_str());
    lammps_command(lmp, "run 2 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    EXPECT_EQ(lammps_extract_global_datatype(lmp, "UNKNOWN"), -1);
    EXPECT_EQ(lammps_extract_global(lmp, "UNKNOWN"), nullptr);

    EXPECT_EQ(lammps_extract_global_datatype(lmp, "units"), LAMMPS_STRING);
    char *c_ptr = (char *)lammps_extract_global(lmp, "units");
    EXPECT_THAT(c_ptr, StrEq("real"));

#if defined(LAMMPS_SMALLSMALL)
    EXPECT_EQ(lammps_extract_global_datatype(lmp, "ntimestep"), LAMMPS_INT);
    int *i_ptr = (int *)lammps_extract_global(lmp, "ntimestep");
    EXPECT_EQ((*i_ptr), 2);
#else
    EXPECT_EQ(lammps_extract_global_datatype(lmp, "ntimestep"), LAMMPS_INT64);
    int64_t *b_ptr = (int64_t *)lammps_extract_global(lmp, "ntimestep");
    EXPECT_EQ((*b_ptr), 2);
#endif

    EXPECT_EQ(lammps_extract_global_datatype(lmp, "dt"), LAMMPS_DOUBLE);
    double *d_ptr = (double *)lammps_extract_global(lmp, "dt");
    EXPECT_DOUBLE_EQ((*d_ptr), 0.1);
};

class AtomProperties : public ::testing::Test {
protected:
    void *lmp;

    AtomProperties(){};
    ~AtomProperties() override{};

    void SetUp() override
    {
        const char *args[] = {"LAMMPS_test", "-log", "none", "-echo", "screen", "-nocite"};

        char **argv = (char **)args;
        int argc    = sizeof(args) / sizeof(char *);

        ::testing::internal::CaptureStdout();
        lmp                = lammps_open_no_mpi(argc, argv, NULL);
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        EXPECT_THAT(output, StartsWith("LAMMPS ("));
        ::testing::internal::CaptureStdout();
        lammps_command(lmp, "region box block 0 2 0 2 0 2");
        lammps_command(lmp, "create_box 1 box");
        lammps_command(lmp, "mass 1 3.0");
        lammps_command(lmp, "create_atoms 1 single 1.0 1.0 1.5");
        lammps_command(lmp, "create_atoms 1 single 0.2 0.1 0.1");
        output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
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

TEST_F(AtomProperties, invalid)
{
    ASSERT_EQ(lammps_extract_atom(lmp, "UNKNOWN"), nullptr);
}

TEST_F(AtomProperties, mass)
{
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "mass"), LAMMPS_DOUBLE);
    double *mass = (double *)lammps_extract_atom(lmp, "mass");
    ASSERT_NE(mass, nullptr);
    ASSERT_DOUBLE_EQ(mass[1], 3.0);
}

TEST_F(AtomProperties, id)
{
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "id"), LAMMPS_TAGINT);
    tagint *id = (tagint *)lammps_extract_atom(lmp, "id");
    ASSERT_NE(id, nullptr);
    ASSERT_EQ(id[0], 1);
    ASSERT_EQ(id[1], 2);
}

TEST_F(AtomProperties, type)
{
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "type"), LAMMPS_INT);
    int *type = (int *)lammps_extract_atom(lmp, "type");
    ASSERT_NE(type, nullptr);
    ASSERT_EQ(type[0], 1);
    ASSERT_EQ(type[1], 1);
}

TEST_F(AtomProperties, position)
{
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "x"), LAMMPS_DOUBLE_2D);
    double **x = (double **)lammps_extract_atom(lmp, "x");
    ASSERT_NE(x, nullptr);
    EXPECT_DOUBLE_EQ(x[0][0], 1.0);
    EXPECT_DOUBLE_EQ(x[0][1], 1.0);
    EXPECT_DOUBLE_EQ(x[0][2], 1.5);
    EXPECT_DOUBLE_EQ(x[1][0], 0.2);
    EXPECT_DOUBLE_EQ(x[1][1], 0.1);
    EXPECT_DOUBLE_EQ(x[1][2], 0.1);
}
