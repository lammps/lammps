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
#include "../testing/utils.h"

#include "fix.h"
#include "info.h"
#include "lammps.h"
#include "library.h"
#include "modify.h"
#include "utils.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mpi.h>
#include <string>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

// when a charge equilibration fix is applied to a subset of the atoms, the
// atoms outside the fix group keep their charges and contribute to the
// electrostatic potential felt by the atoms inside the group.  the tests below
// use a small cluster of six atoms in a non-periodic box: atoms 1-3 are in the
// fix group, atoms 4-6 carry fixed charges.  the cluster is small enough that
// the charges can be obtained by solving the linear equations of each method
// directly; the reference values below were computed that way (with numpy),
// independent of the LAMMPS implementation.  each fix style is run in its
// plain version and, if available, with the OPENMP suffix.  the
// QEqSubgroupKokkos ctest entry runs the same tests with
// LAMMPS_ACCELERATOR_ARGS="-k on t 1 -sf kk", so that the KOKKOS versions are
// checked against the same reference data.

class QEqSubgroupTest : public LAMMPSTest {
protected:
    static constexpr int nfixed            = 3;
    static constexpr double qfixed[nfixed] = {0.7, -0.4, 0.25};

    // fix qeq/reaxff and fix qeq/shielded with the parameters of param_file
    static constexpr double ref_shielded[3] = {-0.00539768877713178, -0.0818737261107374,
                                               0.0872714148878691};
    // fix qeq/point: matrix element 1/r within the cutoff
    static constexpr double ref_point[3] = {-0.0201606016070045, 0.0364449747429133,
                                            -0.0162843731359089};
    // fix qeq/reaxff with equal electronegativities for both types (equal_file).
    // fix qtpie/reaxff then has a vanishing effective electronegativity, which
    // is equivalent, since a constant shift of all electronegativities does not
    // change the charges.
    static constexpr double ref_equalchi[3] = {0.0487818013561279, -0.154416759025177,
                                               0.105634957669049};
    // fix acks2/reaxff with bond softness 5.0 and bond cutoff 4.0 for both types.
    // atoms 2 and 4 are within the bond cutoff: the X matrix only couples atoms
    // in the fix group.
    static constexpr double ref_acks2[3] = {-0.00183602052443018, -0.00943151375913459,
                                            0.0112675342835647};

    const std::string param_file = "qeq_subgroup_param.txt";
    const std::string equal_file = "qeq_subgroup_equal.txt";
    const std::string acks2_file = "qeq_subgroup_acks2.txt";
    const std::string gauss_file = "qeq_subgroup_gauss.txt";
    const std::string qeq_file   = "qeq_subgroup_qeq.txt";

    void SetUp() override
    {
        testbinary = "QEqSubgroupTest";
        LAMMPSTest::SetUp();

        std::ofstream(param_file) << "1 5.3200 14.8732 1.0206\n2 4.5000 13.2000 0.8000\n";
        std::ofstream(equal_file) << "1 5.0000 14.8732 1.0206\n2 5.0000 13.2000 0.8000\n";
        std::ofstream(acks2_file) << "5.0\n1 5.3200 14.8732 1.0206 4.0\n"
                                     "2 4.5000 13.2000 0.8000 4.0\n";
        std::ofstream(gauss_file) << "1 0.2240\n2 0.5434\n";
        // the QEQ package fixes read two more columns (zeta and core charge)
        std::ofstream(qeq_file) << "1 5.3200 14.8732 1.0206 0.0 0.0\n"
                                   "2 4.5000 13.2000 0.8000 0.0 0.0\n";

        // the OPENMP versions of the fixes need the package command
        if (!using_accelerator() && Info::has_package("OPENMP"))
            setup_system("package omp 2");
        else
            setup_system("");
    }

    // (re-)create the cluster, the package command must come before the box
    void setup_system(const std::string &package)
    {
        BEGIN_HIDE_OUTPUT();
        command("clear");
        if (!package.empty()) command(package);
        command("units real");
        command("atom_style charge");
        command("atom_modify map array");
        command("boundary f f f");
        command("region box block -20 20 -20 20 -20 20");
        command("create_box 2 box");
        command("create_atoms 1 single  0.0  0.0  0.0 units box");
        command("create_atoms 2 single  1.6  0.3 -0.2 units box");
        command("create_atoms 1 single -1.1  1.4  0.5 units box");
        command("create_atoms 2 single  3.9 -2.1  1.0 units box");
        command("create_atoms 1 single -4.0 -1.5  2.2 units box");
        command("create_atoms 2 single  2.0  4.5 -3.0 units box");
        command("mass * 1.0");
        for (int i = 0; i < nfixed; ++i)
            command(fmt::format("set atom {} charge {}", i + 4, qfixed[i]));
        command("group fixed id 4 5 6");
        command("group active subtract all fixed");
        command("pair_style zero 12.0");
        command("pair_coeff * *");
        END_HIDE_OUTPUT();
    }

    void TearDown() override
    {
        LAMMPSTest::TearDown();
        for (const auto &file : {param_file, equal_file, acks2_file, gauss_file, qeq_file})
            delete_file(file);
    }

    // suffixes to test a fix style with: the OPENMP version in addition to the
    // plain one, if it exists and the test does not run with an accelerator package
    std::vector<std::string> suffixes(const std::string &style)
    {
        std::vector<std::string> list = {""};
        if (!using_accelerator() && info->has_style("fix", style + "/omp"))
            list.emplace_back("omp");
        return list;
    }

    // run charge equilibration with the given fix command for the atoms in the
    // group and compare the charges with the reference.  the charges of the atoms
    // outside the group must not change.

    void check_charges(const std::string &fixargs, const double *ref, double tol,
                       const std::string &suffix)
    {
        BEGIN_HIDE_OUTPUT();
        command("set group active charge 0.0");
        if (!suffix.empty()) command("suffix " + suffix);
        command("fix qeq active " + fixargs);
        if (!suffix.empty()) command("suffix off");
        command("run 0 post no");
        END_HIDE_OUTPUT();

        // make certain the intended variant of the fix was tested
        const std::string style = lmp->modify->get_fix_by_id("qeq")->style;
        const std::string base  = utils::trim(fixargs.substr(0, fixargs.find(' ')));
        if (!suffix.empty()) {
            EXPECT_EQ(style, base + "/" + suffix);
        } else if (using_accelerator() && info->has_style("fix", base + "/kk")) {
            EXPECT_EQ(style, base + "/kk");
        } else {
            EXPECT_EQ(style, base);
        }

        const auto *tag  = (tagint *)lammps_extract_atom(lmp, "id");
        const auto *q    = (double *)lammps_extract_atom(lmp, "q");
        const int nlocal = lammps_extract_setting(lmp, "nlocal");
        ASSERT_EQ(nlocal, 6);

        const std::string info = fmt::format("fix {} ({})", fixargs, style);
        double qsum            = 0.0;
        for (int i = 0; i < nlocal; ++i) {
            const int idx = tag[i] - 1;
            if (idx < 3) {
                EXPECT_NEAR(q[i], ref[idx], prec_tol(ref[idx], tol)) << info << " atom " << tag[i];
                qsum += q[i];
            } else {
                EXPECT_DOUBLE_EQ(q[i], qfixed[idx - 3]) << info << " atom " << tag[i];
            }
        }
        // the fix group is constrained to be charge neutral
        EXPECT_NEAR(qsum, 0.0, prec_tol(1.0, tol)) << info;

        BEGIN_HIDE_OUTPUT();
        command("unfix qeq");
        END_HIDE_OUTPUT();
    }
};

TEST_F(QEqSubgroupTest, qeq_reaxff)
{
    if (!info->has_style("fix", "qeq/reaxff")) GTEST_SKIP() << "fix qeq/reaxff is not available";

    for (const auto &sfx : suffixes("qeq/reaxff"))
        check_charges("qeq/reaxff 1 0.0 10.0 1.0e-12 " + param_file + " maxiter 1000", ref_shielded,
                      1.0e-10, sfx);

    // also check the dual CG solver of fix qeq/reaxff/omp
    if (!using_accelerator() && Info::has_package("OPENMP"))
        check_charges("qeq/reaxff 1 0.0 10.0 1.0e-12 " + param_file + " maxiter 1000 dual",
                      ref_shielded, 1.0e-10, "omp");

    // check both neighbor list settings and the matrix-free option of fix qeq/reaxff/kk
    if (using_accelerator() && Info::has_package("KOKKOS")) {
        for (const auto &neigh : {"half", "full"}) {
            setup_system(fmt::format("package kokkos neigh/qeq {}", neigh));
            check_charges("qeq/reaxff 1 0.0 10.0 1.0e-12 " + param_file + " maxiter 1000",
                          ref_shielded, 1.0e-10, "");
            check_charges("qeq/reaxff 1 0.0 10.0 1.0e-12 " + param_file + " maxiter 1000 matfree",
                          ref_shielded, 1.0e-10, "");
        }
    }
}

TEST_F(QEqSubgroupTest, qeq_rel_reaxff)
{
    if (!info->has_style("fix", "qeq/rel/reaxff"))
        GTEST_SKIP() << "fix qeq/rel/reaxff is not available";

    // without an external field fix qeq/rel/reaxff is equivalent to fix qeq/reaxff
    for (const auto &sfx : suffixes("qeq/rel/reaxff"))
        check_charges("qeq/rel/reaxff 1 0.0 10.0 1.0e-12 " + param_file + " " + gauss_file +
                          " maxiter 1000",
                      ref_shielded, 1.0e-10, sfx);
}

TEST_F(QEqSubgroupTest, qtpie_reaxff)
{
    if (!info->has_style("fix", "qtpie/reaxff"))
        GTEST_SKIP() << "fix qtpie/reaxff is not available";

    for (const auto &sfx : suffixes("qeq/reaxff"))
        check_charges("qeq/reaxff 1 0.0 10.0 1.0e-12 " + equal_file + " maxiter 1000", ref_equalchi,
                      1.0e-10, sfx);
    for (const auto &sfx : suffixes("qtpie/reaxff"))
        check_charges("qtpie/reaxff 1 0.0 10.0 1.0e-12 " + equal_file + " " + gauss_file +
                          " maxiter 1000",
                      ref_equalchi, 1.0e-10, sfx);
}

TEST_F(QEqSubgroupTest, acks2_reaxff)
{
    if (!info->has_style("fix", "acks2/reaxff"))
        GTEST_SKIP() << "fix acks2/reaxff is not available";
    // fix acks2/reaxff/kk requires a double precision KOKKOS package
    if (kokkos_reduced_precision()) GTEST_SKIP() << "fix acks2/reaxff/kk requires double precision";

    // the BiCGStab solver stops when the squared residual is below the tolerance
    for (const auto &sfx : suffixes("acks2/reaxff"))
        check_charges("acks2/reaxff 1 0.0 10.0 1.0e-12 " + acks2_file + " maxiter 1000", ref_acks2,
                      1.0e-7, sfx);
}

TEST_F(QEqSubgroupTest, qeq_shielded)
{
    if (!info->has_style("fix", "qeq/shielded"))
        GTEST_SKIP() << "fix qeq/shielded is not available";

    // fix qeq/shielded is equivalent to fix qeq/reaxff
    for (const auto &sfx : suffixes("qeq/shielded")) {
        check_charges("qeq/shielded 1 10.0 1.0e-12 1000 " + qeq_file, ref_shielded, 1.0e-10, sfx);
        check_charges("qeq/shielded/xlmd 1 10.0 1.0e-12 1000 " + qeq_file, ref_shielded, 1.0e-10,
                      sfx);
    }
}

TEST_F(QEqSubgroupTest, qeq_point)
{
    if (!info->has_style("fix", "qeq/point")) GTEST_SKIP() << "fix qeq/point is not available";

    for (const auto &sfx : suffixes("qeq/point")) {
        check_charges("qeq/point 1 10.0 1.0e-12 1000 " + qeq_file, ref_point, 1.0e-10, sfx);
        check_charges("qeq/point/xlmd 1 10.0 1.0e-12 1000 " + qeq_file, ref_point, 1.0e-10, sfx);
    }
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
            if (arg == "-v") verbose = true;
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();

    // finalize the KOKKOS package explicitly: otherwise Kokkos is torn down by
    // static destructors at program exit, leading to segfaults in some cases

    lammps_kokkos_finalize();

    MPI_Finalize();
    return rv;
}
