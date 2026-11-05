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
#include <algorithm>
#include <cmath>
#include "../testing/utils.h"

#include "comm.h"
#include "info.h"
#include "lammps.h"
#include "library.h"
#include "utils.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <mpi.h>
#include <string>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

// the dynamical_matrix and third_order commands of the PHONON package compute
// the second and third derivatives of the potential energy by finite
// differences of the forces and write them to a file.  the tests below use a
// two atom group of a perfect fcc lattice of Lennard-Jones particles, which
// keeps the output small, and are run a second time with the KOKKOS versions
// of the two commands through the LAMMPS_ACCELERATOR_ARGS setting of the
// PhononCommandsKokkos ctest entry.

class PhononCommandsTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "PhononCommandsTest";
        LAMMPSTest::SetUp();
        if (info->has_style("command", "dynamical_matrix")) {
            BEGIN_HIDE_OUTPUT();
            command("units lj");
            command("atom_style atomic");
            command("atom_modify map array");
            command("boundary p p p");
            command("lattice fcc 0.8442");
            command("region box block 0 2 0 2 0 2");
            command("create_box 1 box");
            command("create_atoms 1 box");
            command("mass 1 1.0");
            command("pair_style lj/cut 1.2");
            command("pair_coeff 1 1 1.0 1.0");
            command("neighbor 0.3 bin");
            command("neigh_modify delay 0 every 1 check no");
            command("group two id 1 2");
            command("run 0 post no");
            END_HIDE_OUTPUT();
        }
    }

    // split a line of numbers into doubles
    static std::vector<double> numbers(const std::string &line)
    {
        std::vector<double> values;
        for (const auto &word : utils::split_words(line))
            values.push_back(std::stod(word));
        return values;
    }
};

TEST_F(PhononCommandsTest, DynamicalMatrix)
{
    // the dynamical matrix is a finite difference of the forces, so the error
    // of single precision forces is amplified by the inverse of the step size
    // and reaches a few percent, far beyond a tolerance that stays meaningful
    if (kokkos_reduced_precision())
        GTEST_SKIP() << "finite differences are too inaccurate with single precision forces";

    if (!info->has_style("command", "dynamical_matrix")) GTEST_SKIP();

    const std::string outfile = "test_dynmat.dat";
    delete_file(outfile);

    BEGIN_HIDE_OUTPUT();
    command("dynamical_matrix two regular 1.0e-6 file " + outfile);
    END_HIDE_OUTPUT();

    ASSERT_FILE_EXISTS(outfile);
    auto lines = read_lines(outfile);

    // one line of three columns for each degree of freedom of the group (3*2)
    // and each atom of the group (2)
    ASSERT_EQ((int)lines.size(), 12);

    // assemble the 6x6 force constant matrix of the group
    double matrix[6][6];
    for (std::size_t n = 0; n < lines.size(); ++n) {
        auto values = numbers(lines[n]);
        ASSERT_EQ((int)values.size(), 3);
        const int row = (int)n / 2;
        const int col = ((int)n % 2) * 3;
        for (int k = 0; k < 3; ++k)
            matrix[row][col + k] = values[k];
    }

    // the force constant matrix must be symmetric
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            EXPECT_NEAR(matrix[i][j], matrix[j][i], prec_tol(matrix[j][i], 1.0e-6));

    // reference values of the perfect fcc lattice
    EXPECT_NEAR(matrix[0][0], 68.86210274, prec_tol(68.86210274, 1.0e-6));
    EXPECT_NEAR(matrix[1][1], 68.86210274, prec_tol(68.86210274, 1.0e-6));
    EXPECT_NEAR(matrix[2][2], 68.86210274, prec_tol(68.86210274, 1.0e-6));
    EXPECT_NEAR(matrix[0][3], -7.73672374, prec_tol(-7.73672374, 1.0e-6));
    EXPECT_NEAR(matrix[0][4], -5.99464554, prec_tol(-5.99464554, 1.0e-6));
    EXPECT_NEAR(matrix[2][5], -1.74207820, prec_tol(-1.74207820, 1.0e-6));
    EXPECT_NEAR(matrix[0][5], 0.0, prec_tol(0.0, 1.0e-6));

    delete_file(outfile);
}

TEST_F(PhononCommandsTest, ThirdOrder)
{
    if (!info->has_style("command", "third_order")) GTEST_SKIP();

    const std::string outfile = "test_thirdorder.dat";
    delete_file(outfile);

    BEGIN_HIDE_OUTPUT();
    command("third_order two regular 1.0e-6 file " + outfile);
    END_HIDE_OUTPUT();

    ASSERT_FILE_EXISTS(outfile);
    auto lines = read_lines(outfile);

    // atom i, direction alpha, atom j, direction beta, atom k plus three values
    // for the three directions of atom k, for all group atoms and directions
    // entries below the write threshold of the command are omitted, and which
    // ones fall below it differs once the forces are single precision
    if (kokkos_reduced_precision() && ((int)lines.size() != 72))
        GTEST_SKIP() << "third order output of a reduced precision KOKKOS build differs";
    ASSERT_EQ((int)lines.size(), 72);

    std::map<std::string, std::vector<double>> entries;
    for (const auto &line : lines) {
        auto words = utils::split_words(line);
        ASSERT_EQ((int)words.size(), 8);
        const std::string key = words[0] + " " + words[1] + " " + words[2] + " " + words[3] + " " +
            words[4];
        entries[key] = {std::stod(words[5]), std::stod(words[6]), std::stod(words[7])};
    }
    ASSERT_EQ((int)entries.size(), 72);

    // the third derivatives of a Lennard-Jones fcc lattice are large along the
    // lattice directions that connect the two atoms of the group.  the finite
    // difference noise of the vanishing components is of order 1e-3, so the
    // reference values are compared with a matching tolerance.
    EXPECT_NEAR(entries["1 1 1 1 2"][0], 122.32675983, 1.0e-2);
    EXPECT_NEAR(entries["1 1 1 1 2"][1], 136.60317322, 1.0e-2);
    EXPECT_NEAR(entries["1 1 1 2 2"][0], 136.60367282, 1.0e-2);
    EXPECT_NEAR(entries["1 1 1 3 2"][2], -7.13828996, 1.0e-2);
    EXPECT_NEAR(entries["1 1 2 1 2"][0], -122.32675983, 1.0e-2);
    EXPECT_NEAR(entries["1 1 2 2 1"][0], 136.60383935, 1.0e-2);

    delete_file(outfile);
}

// two atoms that interact only through kspace and are farther apart than the ghost
// cutoff.  the PhononCommandsMPI ctest entry runs these tests on two MPI ranks, where
// the rank that owns one atom holds no copy of the other.

class PhononLongRangeTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "PhononLongRangeTest";
        LAMMPSTest::SetUp();
        // there is no KOKKOS version of kspace style ewald
        if (!info->has_style("kspace", "ewald") || using_accelerator())
            GTEST_SKIP() << "needs kspace style ewald without accelerator";

        // NaCl point charges, 33.8 x 11.3 x 11.3 Angstrom, split along x.  atom 13 is at
        // (9.87, 1.41, 1.41), atom 36 at (23.97, 4.23, 4.23), and the ghost cutoff is 6.0
        BEGIN_HIDE_OUTPUT();
        command("units metal");
        command("atom_style charge");
        command("atom_modify map array");
        command("processors * 1 1");
        command("lattice custom 5.64 a1 1 0 0 a2 0 1 0 a3 0 0 1 origin 0.25 0.25 0.25 "
                "basis 0 0 0 basis 0.5 0.5 0 basis 0.5 0 0.5 basis 0 0.5 0.5 "
                "basis 0.5 0 0 basis 0 0.5 0 basis 0 0 0.5 basis 0.5 0.5 0.5");
        command("region box block 0 6 0 2 0 2");
        command("create_box 2 box");
        command("create_atoms 1 box basis 5 2 basis 6 2 basis 7 2 basis 8 2");
        command("mass 1 22.98977");
        command("mass 2 35.453");
        command("set type 1 charge 1.0");
        command("set type 2 charge -1.0");
        command("pair_style coul/long 5.0");
        command("pair_coeff * *");
        command("kspace_style ewald 1.0e-6");
        command("neighbor 1.0 bin");
        command("reset_atoms id sort yes");
        command("group pair id 13 36");
        END_HIDE_OUTPUT();
    }

    // run cmd and compare elements {index, serial reference} of its binary output of
    // "size" values.  only rank 0 touches the file, whose name contains the number of
    // ranks, since the serial and the MPI tests may run at once in the same directory
    void check(const std::string &cmd, int size, const std::vector<std::pair<int, double>> &ref,
               double reltol)
    {
        const auto file = "test_phonon_long_range_" + std::to_string(lmp->comm->nprocs) + ".bin";
        if (lmp->comm->me == 0) delete_file(file);
        MPI_Barrier(MPI_COMM_WORLD);
        BEGIN_HIDE_OUTPUT();
        command(cmd + " file " + file + " binary yes");
        END_HIDE_OUTPUT();
        if (lmp->comm->me == 0) {
            std::vector<double> data(size + 1, 0.0);
            std::ifstream in(file, std::ios::binary);
            in.read(reinterpret_cast<char *>(data.data()), sizeof(double) * (size + 1));
            EXPECT_EQ(in.gcount(), (std::streamsize)(sizeof(double) * size));
            in.close();
            delete_file(file);
            for (const auto &r : ref) {
                EXPECT_NE(data[r.first], 0.0);
                EXPECT_NEAR(data[r.first], r.second, reltol * std::fabs(r.second));
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
};

// rows: displaced atom (13, 36) and direction, columns: atom and direction of the force.
// D(13x,36x), D(13x,36y), D(36x,13x), D(36y,13x)
TEST_F(PhononLongRangeTest, DynamicalMatrix)
{
    check("dynamical_matrix pair regular 0.05", 36,
          {{3, 1.471983646960e-03}, {4, 1.032002541766e-05}, {18, 1.471983646960e-03},
           {24, 1.031732983749e-05}},
          1.0e-8);
}

// the atoms are not neighbors, so rows are (atom i, alpha, atom i, beta), columns (atom k,
// gamma).  Phi(13x,13x,36y), Phi(13y,13y,36y), Phi(36x,36x,13y), Phi(36y,36y,13y).  third
// order finite differences amplify the round-off of the forces more, hence the tolerance
TEST_F(PhononLongRangeTest, ThirdOrder)
{
    check("third_order pair regular 0.05", 108,
          {{4, -1.781198082636e-04}, {28, 1.798481285130e-04}, {55, 1.781198081960e-04},
           {79, -1.798481285398e-04}},
          1.0e-7);
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
