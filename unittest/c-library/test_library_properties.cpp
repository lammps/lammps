// unit tests for checking and changing simulation properties through the library interface

#include "library.h"

#include "atom.h"
#include "lammps.h"
#include "lmptype.h"
#include "platform.h"
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "test_main.h"

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

using ::LAMMPS_NS::Atom;
using ::LAMMPS_NS::bigint;
using ::LAMMPS_NS::tagint;
using ::LAMMPS_NS::platform::path_join;
using ::testing::HasSubstr;
using ::testing::StartsWith;
using ::testing::StrEq;

class LibraryProperties : public ::testing::Test {
protected:
    void *lmp;
    std::string INPUT_DIR = STRINGIFY(TEST_INPUT_FOLDER);

    LibraryProperties()           = default;
    ~LibraryProperties() override = default;

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
#if (defined(__linux__) || defined(__APPLE__) || defined(_WIN32)) && !defined(__INTEL_LLVM_COMPILER)
    EXPECT_GT(meminfo[2], 0.0);
#endif
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
    std::string input = path_join(INPUT_DIR, "in.fourmol");
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_file(lmp, input.c_str());
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_get_natoms(lmp), 29);
};

TEST_F(LibraryProperties, thermo)
{
    bigint bval = *(bigint *)lammps_last_thermo(lmp, "step", 0);
    EXPECT_EQ(bval, -1);

    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = path_join(INPUT_DIR, "in.fourmol");
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

    bval = *(bigint *)lammps_last_thermo(lmp, "step", 0);
    EXPECT_EQ(bval, 2);
    int ival = *(int *)lammps_last_thermo(lmp, "num", 0);
    EXPECT_EQ(ival, 6);

    const char *key = (const char *)lammps_last_thermo(lmp, "keyword", 0);
    EXPECT_THAT(key, StrEq("Step"));
    ival = *(int *)lammps_last_thermo(lmp, "type", 0);
#if defined(LAMMPS_SMALLSMALL)
    EXPECT_EQ(ival, LAMMPS_INT);
    ival = *(int *)lammps_last_thermo(lmp, "data", 0);
    EXPECT_EQ(ival, 2);
#else
    EXPECT_EQ(ival, LAMMPS_INT64);
    bval = *(bigint *)lammps_last_thermo(lmp, "data", 0);
    EXPECT_EQ(bval, 2);
#endif

    key = (const char *)lammps_last_thermo(lmp, "keyword", 1);
    EXPECT_THAT(key, StrEq("Temp"));
    ival = *(int *)lammps_last_thermo(lmp, "type", 1);
    EXPECT_EQ(ival, LAMMPS_DOUBLE);
    double dval = *(double *)lammps_last_thermo(lmp, "data", 1);
    EXPECT_DOUBLE_EQ(dval, 28.042780385852982);

    key = (const char *)lammps_last_thermo(lmp, "keyword", 2);
    EXPECT_THAT(key, StrEq("E_pair"));
    ival = *(int *)lammps_last_thermo(lmp, "type", 2);
    EXPECT_EQ(ival, LAMMPS_DOUBLE);
    dval = *(double *)lammps_last_thermo(lmp, "data", 2);
    EXPECT_DOUBLE_EQ(dval, 0.0);

    key = (const char *)lammps_last_thermo(lmp, "keyword", 3);
    EXPECT_THAT(key, StrEq("E_mol"));
    ival = *(int *)lammps_last_thermo(lmp, "type", 3);
    EXPECT_EQ(ival, LAMMPS_DOUBLE);
    dval = *(double *)lammps_last_thermo(lmp, "data", 3);
    EXPECT_DOUBLE_EQ(dval, 0.0);

    key = (const char *)lammps_last_thermo(lmp, "keyword", 4);
    EXPECT_THAT(key, StrEq("TotEng"));
    ival = *(int *)lammps_last_thermo(lmp, "type", 4);
    EXPECT_EQ(ival, LAMMPS_DOUBLE);
    dval = *(double *)lammps_last_thermo(lmp, "data", 4);
    EXPECT_DOUBLE_EQ(dval, 2.3405256449146163);

    key = (const char *)lammps_last_thermo(lmp, "keyword", 5);
    EXPECT_THAT(key, StrEq("Press"));
    ival = *(int *)lammps_last_thermo(lmp, "type", 5);
    EXPECT_EQ(ival, LAMMPS_DOUBLE);
    dval = *(double *)lammps_last_thermo(lmp, "data", 5);
    EXPECT_DOUBLE_EQ(dval, 31.700964689115658);
};

TEST_F(LibraryProperties, box)
{
    if (!lammps_has_style(lmp, "atom", "full")) GTEST_SKIP();
    std::string input = path_join(INPUT_DIR, "in.fourmol");
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
    if (!verbose) ::testing::internal::CaptureStdout();
    // lammps_reset_box() may only be called without atoms
    lammps_command(lmp, "delete_atoms group all bond yes");
    lammps_reset_box(lmp, boxlo, boxhi, xy, yz, xz);
    if (!verbose) ::testing::internal::GetCapturedStdout();
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

    EXPECT_EQ(lammps_extract_setting(lmp, "kokkos_active"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "kokkos_nthreads"), 0);
    EXPECT_EQ(lammps_extract_setting(lmp, "kokkos_ngpus"), 0);

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
        std::string input = path_join(INPUT_DIR, "in.fourmol");
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
        EXPECT_EQ(lammps_extract_setting(lmp, "sphere_flag"), -1);
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

    std::string input = path_join(INPUT_DIR, "in.fourmol");
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "special_bonds lj 0.0 0.5 0.8 coul 0.1 0.5 1.0");
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
    auto *b_ptr = (int64_t *)lammps_extract_global(lmp, "ntimestep");
    EXPECT_EQ((*b_ptr), 2);
#endif

    EXPECT_EQ(lammps_extract_global_datatype(lmp, "dt"), LAMMPS_DOUBLE);
    auto *d_ptr = (double *)lammps_extract_global(lmp, "dt");
    EXPECT_DOUBLE_EQ((*d_ptr), 0.1);

    EXPECT_EQ(lammps_extract_global_datatype(lmp, "special_lj"), LAMMPS_DOUBLE);
    EXPECT_EQ(lammps_extract_global_datatype(lmp, "special_coul"), LAMMPS_DOUBLE);
    auto *special_lj   = (double *)lammps_extract_global(lmp, "special_lj");
    auto *special_coul = (double *)lammps_extract_global(lmp, "special_coul");
    EXPECT_DOUBLE_EQ(special_lj[0], 1.0);
    EXPECT_DOUBLE_EQ(special_lj[1], 0.0);
    EXPECT_DOUBLE_EQ(special_lj[2], 0.5);
    EXPECT_DOUBLE_EQ(special_lj[3], 0.8);
    EXPECT_DOUBLE_EQ(special_coul[0], 1.0);
    EXPECT_DOUBLE_EQ(special_coul[1], 0.1);
    EXPECT_DOUBLE_EQ(special_coul[2], 0.5);
    EXPECT_DOUBLE_EQ(special_coul[3], 1.0);
    lammps_command(lmp, "special_bonds lj/coul 1.0 1.0 1.0");
    EXPECT_DOUBLE_EQ(special_lj[1], 1.0);
    EXPECT_DOUBLE_EQ(special_lj[2], 1.0);
    EXPECT_DOUBLE_EQ(special_lj[3], 1.0);
    EXPECT_DOUBLE_EQ(special_coul[1], 1.0);
    EXPECT_DOUBLE_EQ(special_coul[2], 1.0);
    EXPECT_DOUBLE_EQ(special_coul[3], 1.0);

    EXPECT_EQ(lammps_extract_global_datatype(lmp, "map_style"), LAMMPS_INT);
#if defined(LAMMPS_BIGBIG)
    EXPECT_EQ(lammps_extract_global_datatype(lmp, "map_tag_max"), LAMMPS_BIGINT);
#else
    EXPECT_EQ(lammps_extract_global_datatype(lmp, "map_tag_max"), LAMMPS_INT);
#endif
    EXPECT_EQ(lammps_extract_global_datatype(lmp, "sametag"), LAMMPS_INT);
    EXPECT_EQ(lammps_extract_global_datatype(lmp, "sortfreq"), LAMMPS_INT);
    EXPECT_EQ(lammps_extract_global_datatype(lmp, "nextsort"), LAMMPS_BIGINT);
    int *sametag  = (int *)lammps_extract_global(lmp, "sametag");
    int map_style = *(int *)lammps_extract_global(lmp, "map_style");
    EXPECT_EQ(map_style, Atom::MAP_ARRAY);
    EXPECT_NE(sametag, nullptr);

    auto *tags              = (tagint *)lammps_extract_atom(lmp, "id");
    const tagint sometags[] = {1, 5, 10, 15, 20};
    for (const auto &sometag : sometags) {
        int idx = lammps_map_atom(lmp, (const void *)&sometag);
        EXPECT_EQ(sometag, tags[idx]);
        int nextidx = sametag[idx];
        if (nextidx >= 0) {
            EXPECT_EQ(sometag, tags[nextidx]);
        }
    }

    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "clear");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    map_style = *(int *)lammps_extract_global(lmp, "map_style");
    EXPECT_EQ(map_style, Atom::MAP_NONE);
    sametag = (int *)lammps_extract_global(lmp, "sametag");
    EXPECT_EQ(sametag, nullptr);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "atom_modify map yes");
    lammps_command(lmp, "region box block 0 1 0 1 0 1");
    lammps_command(lmp, "create_box 1 box");
    lammps_command(lmp, "mass 1 1.0");
    lammps_command(lmp, "run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    map_style = *(int *)lammps_extract_global(lmp, "map_style");
    EXPECT_EQ(map_style, Atom::MAP_YES);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "clear");
    lammps_command(lmp, "atom_modify map hash");
    lammps_command(lmp, "region box block 0 1 0 1 0 1");
    lammps_command(lmp, "create_box 1 box");
    lammps_command(lmp, "mass 1 1.0");
    lammps_command(lmp, "run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    map_style = *(int *)lammps_extract_global(lmp, "map_style");
    EXPECT_EQ(map_style, Atom::MAP_HASH);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "clear");
    lammps_command(lmp, "atom_modify map array");
    lammps_command(lmp, "region box block 0 1 0 1 0 1");
    lammps_command(lmp, "create_box 1 box");
    lammps_command(lmp, "mass 1 1.0");
    lammps_command(lmp, "run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    map_style = *(int *)lammps_extract_global(lmp, "map_style");
    EXPECT_EQ(map_style, Atom::MAP_ARRAY);

    EXPECT_EQ(lammps_extract_global_datatype(lmp, "xlattice"), LAMMPS_DOUBLE);
    EXPECT_EQ(lammps_extract_global_datatype(lmp, "ylattice"), LAMMPS_DOUBLE);
    EXPECT_EQ(lammps_extract_global_datatype(lmp, "zlattice"), LAMMPS_DOUBLE);
    auto *xlattice = (double *)lammps_extract_global(lmp, "xlattice");
    auto *ylattice = (double *)lammps_extract_global(lmp, "ylattice");
    auto *zlattice = (double *)lammps_extract_global(lmp, "zlattice");
    EXPECT_NE(xlattice, nullptr);
    EXPECT_NE(ylattice, nullptr);
    EXPECT_NE(zlattice, nullptr);
    EXPECT_DOUBLE_EQ(*xlattice, 1.0);
    EXPECT_DOUBLE_EQ(*ylattice, 1.0);
    EXPECT_DOUBLE_EQ(*zlattice, 1.0);
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "clear");
    lammps_command(lmp, "units real");
    lammps_command(lmp, "lattice fcc 2.0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    xlattice = (double *)lammps_extract_global(lmp, "xlattice");
    ylattice = (double *)lammps_extract_global(lmp, "ylattice");
    zlattice = (double *)lammps_extract_global(lmp, "zlattice");
    EXPECT_NE(xlattice, nullptr);
    EXPECT_NE(ylattice, nullptr);
    EXPECT_NE(zlattice, nullptr);
    EXPECT_DOUBLE_EQ(*xlattice, 2.0);
    EXPECT_DOUBLE_EQ(*ylattice, 2.0);
    EXPECT_DOUBLE_EQ(*zlattice, 2.0);
};

TEST_F(LibraryProperties, pair1)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "region box block 0 1 0 1 0 1");
    lammps_command(lmp, "create_box 3 box");
    lammps_command(lmp, "mass * 1.0");
    lammps_command(lmp, "pair_style lj/cut 3.0");
    lammps_command(lmp, "pair_coeff 1 1 1.0 1.0");
    lammps_command(lmp, "pair_coeff 2 2 1.5 2.0");
    lammps_command(lmp, "pair_coeff 3 3 1.0 3.0");
    lammps_command(lmp, "run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_extract_pair_dimension(lmp, "epsilon"), 2);
    EXPECT_EQ(lammps_extract_pair_dimension(lmp, "sigma"), 2);
    EXPECT_EQ(lammps_extract_pair_dimension(lmp, "cut_coul"), -1);
    auto **sigma = (double **)lammps_extract_pair(lmp, "sigma");
    EXPECT_DOUBLE_EQ(sigma[1][1], 1.0);
    EXPECT_DOUBLE_EQ(sigma[2][2], 2.0);
    EXPECT_DOUBLE_EQ(sigma[3][3], 3.0);
    EXPECT_DOUBLE_EQ(sigma[1][2], sqrt(2.0));
};

TEST_F(LibraryProperties, pair2)
{
    if (!lammps_has_style(lmp, "pair", "coul/streitz")) GTEST_SKIP();
    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_command(lmp, "units metal");
    lammps_command(lmp, "atom_style charge");
    lammps_command(lmp, "region box block 0 1 0 1 0 1");
    lammps_command(lmp, "create_box 2 box");
    lammps_command(lmp, "mass * 1.0");
    lammps_command(lmp, "pair_style coul/streitz 12.0 wolf 0.31");
    lammps_command(lmp, "pair_coeff * * AlO.streitz Al O");
    lammps_command(lmp, "run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(lammps_extract_pair_dimension(lmp, "chi"), 1);
    EXPECT_EQ(lammps_extract_pair_dimension(lmp, "scale"), 2);
    EXPECT_EQ(lammps_extract_pair_dimension(lmp, "cut_coul"), 0);
    EXPECT_DOUBLE_EQ(*((double *)lammps_extract_pair(lmp, "cut_coul")), 12.0);
    auto *chi = (double *)lammps_extract_pair(lmp, "chi");
    EXPECT_DOUBLE_EQ(chi[1], 0.0);
    EXPECT_FLOAT_EQ(chi[2], 5.484763);
    auto **scale = (double **)lammps_extract_pair(lmp, "scale");
    EXPECT_DOUBLE_EQ(scale[1][1], 1.0);
    EXPECT_DOUBLE_EQ(scale[1][2], 1.0);
    EXPECT_DOUBLE_EQ(scale[2][2], 1.0);
};

TEST_F(LibraryProperties, neighlist)
{
    if (!lammps_has_style(lmp, "pair", "sw")) GTEST_SKIP();
    const char sysinit[] = "boundary f f f\n"
                           "units real\n"
                           "region box block -5 5 -5 5 -5 5\n"
                           "create_box 2 box\n"
                           "mass 1 1.0\n"
                           "mass 2 1.0\n"
                           "pair_style hybrid/overlay lj/cut 4.0 lj/cut 4.0 morse 4.0 sw\n"
                           "pair_coeff * * sw Si.sw Si NULL\n"
                           "pair_coeff 1 2 morse 0.2 2.0 2.0\n"
                           "pair_coeff 2 2 lj/cut 1 0.1 2.0\n"
                           "pair_coeff * * lj/cut 2 0.01 2.0\n"
                           "compute dist all pair/local dist\n"
                           "fix dist all ave/histo 1 1 1 0.0 3.0 4 c_dist mode vector\n"
                           "thermo_style custom f_dist[*]";

    const double pos[] = {0.0, 0.0, 0.0, -1.1, 0.0, 0.0, 1.0,  0.0, 0.0, 0.0, -1.1,
                          0.0, 0.0, 1.0, 0.0,  0.0, 0.0, -1.1, 0.0, 0.0, 1.0};
    const tagint ids[] = {1, 2, 3, 4, 5, 6, 7};
    const int types[]  = {1, 1, 1, 1, 2, 2, 2};

    const int numatoms = sizeof(ids) / sizeof(tagint);

    if (!verbose) ::testing::internal::CaptureStdout();
    lammps_commands_string(lmp, sysinit);
    lammps_create_atoms(lmp, numatoms, ids, types, pos, nullptr, nullptr, 0);
    lammps_command(lmp, "run 0 post no");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    void *ptr  = lammps_extract_fix(lmp, "dist", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, 0, 0);
    int nhisto = *(double *)ptr;
    lammps_free(ptr);

    ptr       = lammps_extract_fix(lmp, "dist", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, 1, 0);
    int nskip = *(double *)ptr;
    lammps_free(ptr);

    ptr           = lammps_extract_fix(lmp, "dist", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, 2, 0);
    double minval = *(double *)ptr;
    lammps_free(ptr);

    ptr           = lammps_extract_fix(lmp, "dist", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, 3, 0);
    double maxval = *(double *)ptr;
    lammps_free(ptr);

    // 21 pair distances counted, none skipped, smallest 1.0, largest 2.1
    EXPECT_EQ(nhisto, 21);
    EXPECT_EQ(nskip, 0);
    EXPECT_DOUBLE_EQ(minval, 1.0);
    EXPECT_DOUBLE_EQ(maxval, 2.1);

    const int nlocal = lammps_extract_setting(lmp, "nlocal");
    EXPECT_EQ(nlocal, numatoms);
    EXPECT_NE(lammps_find_pair_neighlist(lmp, "sw", 1, 0, 0), -1);
    EXPECT_NE(lammps_find_pair_neighlist(lmp, "morse", 1, 0, 0), -1);
    EXPECT_NE(lammps_find_pair_neighlist(lmp, "lj/cut", 1, 1, 0), -1);
    EXPECT_NE(lammps_find_pair_neighlist(lmp, "lj/cut", 1, 2, 0), -1);
    EXPECT_EQ(lammps_find_pair_neighlist(lmp, "lj/cut", 1, 0, 0), -1);
    EXPECT_EQ(lammps_find_pair_neighlist(lmp, "hybrid/overlay", 1, 0, 0), -1);
    EXPECT_NE(lammps_find_compute_neighlist(lmp, "dist", 0), -1);
    EXPECT_EQ(lammps_find_fix_neighlist(lmp, "dist", 0), -1);
    EXPECT_EQ(lammps_find_compute_neighlist(lmp, "xxx", 0), -1);

    // full neighbor list for 4 type 1 atoms
    // all have 3 type 1 atom neighbors
    int idx = lammps_find_pair_neighlist(lmp, "sw", 1, 0, 0);
    int num = lammps_neighlist_num_elements(lmp, idx);
    EXPECT_EQ(num, 4);
    int iatom, inum, *neighbors;
    for (int i = 0; i < num; ++i) {
        lammps_neighlist_element_neighbors(lmp, idx, i, &iatom, &inum, &neighbors);
        EXPECT_EQ(iatom, i);
        EXPECT_EQ(inum, 3);
        EXPECT_NE(neighbors, nullptr);
    }

    // half neighbor list for all pairs between type 1 and type 2
    // 4 type 1 atoms with 3 type 2 neighbors and 3 type 2 atoms without neighbors
    idx = lammps_find_pair_neighlist(lmp, "morse", 0, 0, 0);
    num = lammps_neighlist_num_elements(lmp, idx);
    EXPECT_EQ(num, nlocal);
    for (int i = 0; i < num; ++i) {
        lammps_neighlist_element_neighbors(lmp, idx, i, &iatom, &inum, &neighbors);
        if (i < 4)
            EXPECT_EQ(inum, 3);
        else
            EXPECT_EQ(inum, 0);
        EXPECT_NE(neighbors, nullptr);
    }

    // half neighbor list between type 2 atoms only
    // 3 pairs with 2, 1, 0 neighbors
    idx = lammps_find_pair_neighlist(lmp, "lj/cut", 1, 1, 0);
    num = lammps_neighlist_num_elements(lmp, idx);
    EXPECT_EQ(num, 3);
    for (int i = 0; i < num; ++i) {
        lammps_neighlist_element_neighbors(lmp, idx, i, &iatom, &inum, &neighbors);
        EXPECT_EQ(inum, 2 - i);
        EXPECT_NE(neighbors, nullptr);
    }

    // half neighbor list between all pairs. same as simple lj/cut case
    idx = lammps_find_pair_neighlist(lmp, "lj/cut", 1, 2, 0);
    num = lammps_neighlist_num_elements(lmp, idx);
    EXPECT_EQ(num, nlocal);
    for (int i = 0; i < num; ++i) {
        lammps_neighlist_element_neighbors(lmp, idx, i, &iatom, &inum, &neighbors);
        EXPECT_EQ(inum, nlocal - 1 - i);
        EXPECT_NE(neighbors, nullptr);
    }

    //  the compute has a half neighbor list
    idx = lammps_find_compute_neighlist(lmp, "dist", 0);
    num = lammps_neighlist_num_elements(lmp, idx);
    EXPECT_EQ(num, nlocal);
    for (int i = 0; i < num; ++i) {
        lammps_neighlist_element_neighbors(lmp, idx, i, &iatom, &inum, &neighbors);
        EXPECT_EQ(inum, nlocal - 1 - i);
        EXPECT_NE(neighbors, nullptr);
    }
};

TEST_F(LibraryProperties, has_error)
{
    EXPECT_EQ(lammps_has_error(lmp), 0);

    // trigger an error, but hide output
    ::testing::internal::CaptureStdout();
    lammps_command(lmp, "this_is_not_a_known_command");
    ::testing::internal::GetCapturedStdout();

    EXPECT_EQ(lammps_has_error(lmp), 1);

    // retrieve error message
    char errmsg[1024];
    int err = lammps_get_last_error_message(lmp, errmsg, 1024);
    EXPECT_EQ(err, 1);
    EXPECT_THAT(errmsg, HasSubstr("ERROR: Unknown command: this_is_not_a_known_command"));

    // retrieving the error message clear the error status
    EXPECT_EQ(lammps_has_error(lmp), 0);
    err = lammps_get_last_error_message(lmp, errmsg, 1024);
    EXPECT_EQ(err, 0);
    EXPECT_THAT(errmsg, StrEq(""));
};

class AtomProperties : public ::testing::Test {
protected:
    void *lmp;
    int ntypes, nlocal, nall;

    AtomProperties()           = default;
    ~AtomProperties() override = default;

    void SetUp() override
    {
        const char *args[] = {"LAMMPS_test", "-log", "none", "-echo", "screen", "-nocite", nullptr};

        char **argv = (char **)args;
        int argc    = (sizeof(args) / sizeof(char *)) - 1;

        ::testing::internal::CaptureStdout();
        lmp                = lammps_open_no_mpi(argc, argv, nullptr);
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        EXPECT_THAT(output, StartsWith("LAMMPS ("));
        ::testing::internal::CaptureStdout();
        lammps_command(lmp, "fix props all property/atom i_one i2_two 2 d_three d2_four 2");
        lammps_command(lmp, "fix rmass all property/atom mol q rmass ghost yes");
        lammps_command(lmp, "region box block 0 2 0 2 0 2");
        lammps_command(lmp, "create_box 1 box");
        lammps_command(lmp, "mass 1 3.0");
        lammps_command(lmp, "create_atoms 1 single 1.0 1.0 1.5");
        lammps_command(lmp, "create_atoms 1 single 0.2 0.1 0.1");
        lammps_command(lmp, "set group all mass 2.0");
        lammps_command(lmp, "set atom 1 charge -1");
        lammps_command(lmp, "set atom 2 charge  1");
        lammps_command(lmp, "set atom 1 mol 2");
        lammps_command(lmp, "set atom 2 mol 1");
        lammps_command(lmp, "set atom 1 i_one -3");
        lammps_command(lmp, "set atom 2 i_one  3");
        lammps_command(lmp, "set atom 1 d_three -1.3");
        lammps_command(lmp, "set atom 2 d_three  3.5");
        lammps_command(lmp, "set atom 1 i_two[1] -3");
        lammps_command(lmp, "set atom 2 i_two[2]  3");
        lammps_command(lmp, "set atom * d_four[1] -1.3");
        lammps_command(lmp, "set atom * d_four[2]  3.5");
        ntypes = lammps_extract_setting(lmp, "ntypes");
        nlocal = lammps_extract_setting(lmp, "nlocal");
        nall   = lammps_extract_setting(lmp, "nall");

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
    EXPECT_EQ(lammps_extract_atom_size(lmp, "mass", 0), ntypes + 1);
    auto *mass = (double *)lammps_extract_atom(lmp, "mass");
    ASSERT_NE(mass, nullptr);
    ASSERT_DOUBLE_EQ(mass[1], 3.0);
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "rmass"), LAMMPS_DOUBLE);
    EXPECT_EQ(lammps_extract_atom_size(lmp, "rmass", 0), nall);
    mass = (double *)lammps_extract_atom(lmp, "rmass");
    ASSERT_NE(mass, nullptr);
    ASSERT_DOUBLE_EQ(mass[0], 2.0);
    ASSERT_DOUBLE_EQ(mass[1], 2.0);
}

TEST_F(AtomProperties, charge)
{
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "q"), LAMMPS_DOUBLE);
    EXPECT_EQ(lammps_extract_atom_size(lmp, "rmass", 0), nall);
    auto *charge = (double *)lammps_extract_atom(lmp, "q");
    ASSERT_NE(charge, nullptr);
    ASSERT_DOUBLE_EQ(charge[0], -1.0);
    ASSERT_DOUBLE_EQ(charge[1], 1.0);
}

TEST_F(AtomProperties, molecule)
{
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "molecule"), LAMMPS_TAGINT);
    EXPECT_EQ(lammps_extract_atom_size(lmp, "molecule", 0), nall);
    auto *molecule = (tagint *)lammps_extract_atom(lmp, "molecule");
    ASSERT_NE(molecule, nullptr);
    ASSERT_EQ(molecule[0], 2);
    ASSERT_EQ(molecule[1], 1);
}

TEST_F(AtomProperties, id)
{
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "id"), LAMMPS_TAGINT);
    EXPECT_EQ(lammps_extract_atom_size(lmp, "id", 0), nall);
    auto *id = (tagint *)lammps_extract_atom(lmp, "id");
    ASSERT_NE(id, nullptr);
    ASSERT_EQ(id[0], 1);
    ASSERT_EQ(id[1], 2);
}

TEST_F(AtomProperties, type)
{
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "type"), LAMMPS_INT);
    EXPECT_EQ(lammps_extract_atom_size(lmp, "type", 0), nall);
    int *type = (int *)lammps_extract_atom(lmp, "type");
    ASSERT_NE(type, nullptr);
    ASSERT_EQ(type[0], 1);
    ASSERT_EQ(type[1], 1);
}

TEST_F(AtomProperties, position)
{
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "x"), LAMMPS_DOUBLE_2D);
    EXPECT_EQ(lammps_extract_atom_size(lmp, "x", LMP_SIZE_ROWS), nall);
    EXPECT_EQ(lammps_extract_atom_size(lmp, "x", LMP_SIZE_COLS), 3);
    auto **x = (double **)lammps_extract_atom(lmp, "x");
    ASSERT_NE(x, nullptr);
    EXPECT_DOUBLE_EQ(x[0][0], 1.0);
    EXPECT_DOUBLE_EQ(x[0][1], 1.0);
    EXPECT_DOUBLE_EQ(x[0][2], 1.5);
    EXPECT_DOUBLE_EQ(x[1][0], 0.2);
    EXPECT_DOUBLE_EQ(x[1][1], 0.1);
    EXPECT_DOUBLE_EQ(x[1][2], 0.1);
}

TEST_F(AtomProperties, custom)
{
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "i_one"), LAMMPS_INT);
    EXPECT_EQ(lammps_extract_atom_size(lmp, "i_one", 0), nlocal);
    auto *one = (int *)lammps_extract_atom(lmp, "i_one");
    ASSERT_NE(one, nullptr);
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "i2_two"), LAMMPS_INT_2D);
    EXPECT_EQ(lammps_extract_atom_size(lmp, "i2_two", LMP_SIZE_ROWS), nlocal);
    EXPECT_EQ(lammps_extract_atom_size(lmp, "i2_two", LMP_SIZE_COLS), 2);
    auto **two = (int **)lammps_extract_atom(lmp, "i2_two");
    ASSERT_NE(two, nullptr);
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "d_three"), LAMMPS_DOUBLE);
    EXPECT_EQ(lammps_extract_atom_size(lmp, "d_three", 0), nlocal);
    auto *three = (double *)lammps_extract_atom(lmp, "d_three");
    ASSERT_NE(three, nullptr);
    EXPECT_EQ(lammps_extract_atom_datatype(lmp, "d2_four"), LAMMPS_DOUBLE_2D);
    EXPECT_EQ(lammps_extract_atom_size(lmp, "d2_four", LMP_SIZE_ROWS), nlocal);
    EXPECT_EQ(lammps_extract_atom_size(lmp, "d2_four", LMP_SIZE_COLS), 2);
    auto **four = (double **)lammps_extract_atom(lmp, "d2_four");
    ASSERT_NE(four, nullptr);

    EXPECT_EQ(one[0], -3);
    EXPECT_EQ(one[1], 3);
    EXPECT_EQ(two[0][0], -3);
    EXPECT_EQ(two[0][1], 0);
    EXPECT_EQ(two[1][0], 0);
    EXPECT_EQ(two[1][1], 3);
    EXPECT_DOUBLE_EQ(three[0], -1.3);
    EXPECT_DOUBLE_EQ(three[1], 3.5);
    EXPECT_DOUBLE_EQ(four[0][0], -1.3);
    EXPECT_DOUBLE_EQ(four[0][1], 3.5);
    EXPECT_DOUBLE_EQ(four[1][0], -1.3);
    EXPECT_DOUBLE_EQ(four[1][1], 3.5);
}

TEST(SystemSettings, kokkos)
{
    if (!lammps_config_has_package("KOKKOS")) GTEST_SKIP();
    std::vector<char *> args = {(char *)"lammps", (char *)"-log",   (char *)"none",
                                (char *)"-echo",  (char *)"screen", (char *)"-nocite",
                                (char *)"-sf",    (char *)"kk"};

    char *one  = (char *)"1";
    char *four = (char *)"4";
    char *tee  = (char *)"t";
    char *gee  = (char *)"g";
    char *kay  = (char *)"-k";
    char *yes  = (char *)"on";

    args.push_back(kay);
    args.push_back(yes);

    bool has_gpu     = false;
    bool has_threads = false;

    // when GPU support is enabled in KOKKOS, it *must* be used
    if (lammps_config_accelerator("KOKKOS", "api", "hip") ||
        lammps_config_accelerator("KOKKOS", "api", "cuda") ||
        lammps_config_accelerator("KOKKOS", "api", "sycl")) {
        has_gpu = true;
        args.push_back(gee);
        args.push_back(one);
    }

    // use threads or serial
    args.push_back(tee);
    if (lammps_config_accelerator("KOKKOS", "api", "openmp")) {
        has_threads = true;
        args.push_back(four);
    } else if (lammps_config_accelerator("KOKKOS", "api", "pthreads")) {
        has_threads = true;
        args.push_back(four);
    } else {
        args.push_back(one);
    }
    int argc    = args.size();
    char **argv = args.data();

    ::testing::internal::CaptureStdout();
    void *lmp          = lammps_open_no_mpi(argc, argv, nullptr);
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    EXPECT_THAT(output, StartsWith("LAMMPS ("));

    EXPECT_EQ(lammps_extract_setting(lmp, "kokkos_active"), 1);
    if (has_threads)
        EXPECT_EQ(lammps_extract_setting(lmp, "kokkos_nthreads"), 4);
    else
        EXPECT_EQ(lammps_extract_setting(lmp, "kokkos_nthreads"), 1);
    if (has_gpu)
        EXPECT_EQ(lammps_extract_setting(lmp, "kokkos_ngpus"), 1);
    else
        EXPECT_EQ(lammps_extract_setting(lmp, "kokkos_ngpus"), 0);
    lammps_close(lmp);
}
