/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "../testing/core.h"
#include "info.h"
#include "fmt/format.h"
#include "input.h"
#include "lammps.h"
#include "library.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstring>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

// a bcc crystal of 6x6x6 unit cells that exactly matches the reference lattice,
// so a defect free crystal must be reported as free of defects

class ComputeFrenkelTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "ComputeFrenkelTest";
        LAMMPSTest::SetUp();
        if (!Info::has_package("EXTRA-COMPUTE")) GTEST_SKIP();
        BEGIN_HIDE_OUTPUT();
        command("units metal");
        command("atom_style atomic");
        command("atom_modify map array");
        command("boundary p p p");
        command("lattice bcc 2.8553");
        command("region box block 0 6 0 6 0 6");
        command("create_box 1 box");
        command("create_atoms 1 box");
        command("mass 1 55.845");
        command("pair_style zero 5.0");
        command("pair_coeff * *");
        END_HIDE_OUTPUT();
    }

    // remove the single atom sitting at the given position in lattice units
    void remove_site(const std::string &pos, const std::string &id)
    {
        command(fmt::format("region {} sphere {} 0.2 units lattice", id, pos));
        command(fmt::format("delete_atoms region {} compress no", id));
    }

    void setup_compute(const std::string &args)
    {
        BEGIN_HIDE_OUTPUT();
        command("compute fr all frenkel " + args);
        command("run 0 post no");
        END_HIDE_OUTPUT();
    }

    double *get_vector()
    {
        return (double *)lammps_extract_compute(lmp, "fr", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
    }

    double **get_array()
    {
        return (double **)lammps_extract_compute(lmp, "fr", LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY);
    }
};

TEST_F(ComputeFrenkelTest, defect_free_crystal)
{
    setup_compute("");
    auto *vec = get_vector();
    ASSERT_NE(vec, nullptr);
    EXPECT_DOUBLE_EQ(vec[0], 0.0);
    EXPECT_DOUBLE_EQ(vec[1], 0.0);
    EXPECT_DOUBLE_EQ(vec[2], 0.0);
}

TEST_F(ComputeFrenkelTest, frenkel_pair)
{
    BEGIN_HIDE_OUTPUT();
    remove_site("3.0 3.0 3.0", "hole");
    // a second atom close to an occupied site makes that site an interstitial
    command("create_atoms 1 single 1.5 1.5 1.6 units lattice");
    END_HIDE_OUTPUT();
    setup_compute("");
    auto *vec = get_vector();
    EXPECT_DOUBLE_EQ(vec[0], 1.0);    // one vacancy
    EXPECT_DOUBLE_EQ(vec[1], 1.0);    // one interstitial
    EXPECT_DOUBLE_EQ(vec[2], 0.0);    // no site with more than two atoms
}

TEST_F(ComputeFrenkelTest, irregular_site)
{
    BEGIN_HIDE_OUTPUT();
    command("create_atoms 1 single 1.5 1.5 1.6 units lattice");
    command("create_atoms 1 single 1.5 1.6 1.5 units lattice");
    END_HIDE_OUTPUT();
    setup_compute("");
    auto *vec = get_vector();
    EXPECT_DOUBLE_EQ(vec[0], 0.0);
    EXPECT_DOUBLE_EQ(vec[1], 1.0);    // three atoms on one site are one interstitial
    EXPECT_DOUBLE_EQ(vec[2], 1.0);    // and one irregular site
}

// three pairs of vacancies, at first (0.866a), second (1.0a) and third
// (1.414a) neighbor distance.  The default drvac connects the first two pairs
// but not the third, so there must be two clusters of size 2 and two of size 1.

class ComputeFrenkelClusterTest : public ComputeFrenkelTest {
protected:
    void SetUp() override
    {
        ComputeFrenkelTest::SetUp();
        if (IsSkipped()) return;
        BEGIN_HIDE_OUTPUT();
        remove_site("1.0 1.0 1.0", "a1");
        remove_site("1.5 1.5 1.5", "a2");
        remove_site("3.0 1.0 1.0", "b1");
        remove_site("4.0 1.0 1.0", "b2");
        remove_site("1.0 4.0 1.0", "c1");
        remove_site("2.0 5.0 1.0", "c2");
        END_HIDE_OUTPUT();
    }
};

TEST_F(ComputeFrenkelClusterTest, default_connection_distance)
{
    setup_compute("");
    auto *vec = get_vector();
    auto **arr = get_array();
    ASSERT_NE(arr, nullptr);
    EXPECT_DOUBLE_EQ(vec[0], 6.0);       // six vacancies
    EXPECT_DOUBLE_EQ(arr[0][0], 2.0);    // the third neighbor pair stays apart
    EXPECT_DOUBLE_EQ(arr[0][1], 2.0);    // first and second neighbor pairs merge
}

TEST_F(ComputeFrenkelClusterTest, larger_connection_distance)
{
    // 1.7 nearest neighbor distances also reaches the third neighbor shell of a
    // bcc lattice, which is at 1.633 nearest neighbor distances
    setup_compute("drvac 1.7");
    auto *vec = get_vector();
    auto **arr = get_array();
    EXPECT_DOUBLE_EQ(vec[0], 6.0);
    EXPECT_DOUBLE_EQ(arr[0][0], 0.0);
    EXPECT_DOUBLE_EQ(arr[0][1], 3.0);    // now all three pairs are clusters
}

TEST_F(ComputeFrenkelTest, compute_group)
{
    // atoms outside the group of the compute do not count as occupants, so
    // excluding a single atom from the group leaves its site empty
    BEGIN_HIDE_OUTPUT();
    command("group one id 1");
    command("group rest subtract all one");
    command("compute fr rest frenkel");
    command("run 0 post no");
    END_HIDE_OUTPUT();
    auto *vec = get_vector();
    EXPECT_DOUBLE_EQ(vec[0], 1.0);
    EXPECT_DOUBLE_EQ(vec[1], 0.0);
}

TEST_F(ComputeFrenkelTest, region_restriction)
{
    BEGIN_HIDE_OUTPUT();
    command("region inner block 1 5 1 5 1 5 units lattice");
    END_HIDE_OUTPUT();
    setup_compute("drvac 1.6 drint 1.9 region inner rescale no site_file none");
    auto *vec = get_vector();
    EXPECT_DOUBLE_EQ(vec[0], 0.0);
    EXPECT_DOUBLE_EQ(vec[1], 0.0);
}

TEST_F(ComputeFrenkelTest, unsupported_settings)
{
    TEST_FAILURE(".*ERROR: Unknown compute frenkel keyword: bogus.*",
                 command("compute fr all frenkel bogus 1.0"););
    TEST_FAILURE(".*ERROR: Illegal compute frenkel drvac command: missing argument.*",
                 command("compute fr all frenkel drvac"););
    TEST_FAILURE(".*ERROR: Compute frenkel drvac value must be > 0.0.*",
                 command("compute fr all frenkel drvac -1.0"););
    TEST_FAILURE(".*ERROR: Compute frenkel drint value must be > 0.0.*",
                 command("compute fr all frenkel drint 0.0"););
    TEST_FAILURE(".*ERROR: Region nosuchregion for compute frenkel.*",
                 command("compute fr all frenkel region nosuchregion"););
    TEST_FAILURE(".*ERROR: Compute frenkel site file .* is not readable.*",
                 command("compute fr all frenkel site_file no_such_file.txt"););

    // the settings are keywords of the compute command, not of compute_modify
    BEGIN_HIDE_OUTPUT();
    command("compute fr all frenkel");
    END_HIDE_OUTPUT();
    TEST_FAILURE(".*ERROR: Compute fr frenkel does not support compute_modify drvac command.*",
                 command("compute_modify fr drvac 1.2"););
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();

    // finalize the KOKKOS package explicitly: otherwise Kokkos is torn down by
    // static destructors at program exit, leading to segfaults in some cases
    // same workaround as the force-style and FFT3d test drivers

    lammps_kokkos_finalize();

    MPI_Finalize();
    return rv;
}
