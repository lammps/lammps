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

// Unit tests for nbody MPI communication with large values
// Tests that MPI_Allreduce with MPI_LMP_BIGINT works correctly for values > INT_MAX

#include "lmptype.h"
#include "gtest/gtest.h"

#include <climits>
#include <mpi.h>

using namespace LAMMPS_NS;

namespace {

// Test MPI_Allreduce with MPI_LMP_BIGINT (validates fix_rigid_small.cpp:447)
TEST(NbodyMPI, AllreduceSumBasic)
{
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Each rank contributes 1 (simulates counting local bodies)
    bigint local_count = 1;
    bigint global_count = 0;

    MPI_Allreduce(&local_count, &global_count, 1,
                  MPI_LMP_BIGINT, MPI_SUM, MPI_COMM_WORLD);

    EXPECT_EQ(global_count, nprocs);
}

// Test with large values that would overflow int
TEST(NbodyMPI, AllreduceLargeValues)
{
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Each rank contributes a value that when summed exceeds INT_MAX
    // Use INT_MAX/nprocs + 1 to ensure total > INT_MAX
    bigint local_nbody = static_cast<bigint>(INT_MAX) / nprocs + 1;
    bigint global_nbody = 0;

    MPI_Allreduce(&local_nbody, &global_nbody, 1,
                  MPI_LMP_BIGINT, MPI_SUM, MPI_COMM_WORLD);

    // Verify sum exceeds INT_MAX
    EXPECT_GT(global_nbody, INT_MAX);

    // Verify sum is correct
    bigint expected = (static_cast<bigint>(INT_MAX) / nprocs + 1) * nprocs;
    EXPECT_EQ(global_nbody, expected);
}

// Test boundary case at INT_MAX
TEST(NbodyMPI, AllreduceBoundary)
{
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Rank 0 contributes INT_MAX, others contribute 1
    bigint local_nbody = (rank == 0) ? INT_MAX : 1;
    bigint global_nbody = 0;

    MPI_Allreduce(&local_nbody, &global_nbody, 1,
                  MPI_LMP_BIGINT, MPI_SUM, MPI_COMM_WORLD);

    bigint expected = static_cast<bigint>(INT_MAX) + (nprocs - 1);
    EXPECT_EQ(global_nbody, expected);
    EXPECT_GT(global_nbody, INT_MAX);
}

// Test very large values (billions)
TEST(NbodyMPI, AllreduceVeryLarge)
{
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Each rank contributes 1 billion
    bigint local_nbody = 1000000000LL;
    bigint global_nbody = 0;

    MPI_Allreduce(&local_nbody, &global_nbody, 1,
                  MPI_LMP_BIGINT, MPI_SUM, MPI_COMM_WORLD);

    bigint expected = 1000000000LL * nprocs;
    EXPECT_EQ(global_nbody, expected);

    // With 4 ranks, total should be 4 billion
    if (nprocs >= 4) {
        EXPECT_GE(global_nbody, 4000000000LL);
    }
}

// Test that the actual pattern from fix_rigid_small.cpp works
// This mimics lines 441-447 of the fixed code
TEST(NbodyMPI, ActualCodePattern)
{
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Simulate the actual code pattern
    bigint one = 0;  // Changed from int to bigint (our fix!)

    // Simulate counting bodies (each rank has some)
    int nlocal_body = rank + 1;  // Rank 0 has 1, rank 1 has 2, etc.
    for (int i = 0; i < nlocal_body; i++) {
        one++;
    }

    bigint nbody = 0;
    MPI_Allreduce(&one, &nbody, 1, MPI_LMP_BIGINT, MPI_SUM, MPI_COMM_WORLD);

    // Total should be sum of 1+2+3+...+nprocs = nprocs*(nprocs+1)/2
    bigint expected = static_cast<bigint>(nprocs) * (nprocs + 1) / 2;
    EXPECT_EQ(nbody, expected);
}

// Test zero values (edge case)
TEST(NbodyMPI, AllreduceZero)
{
    bigint local_nbody = 0;
    bigint global_nbody = 0;

    MPI_Allreduce(&local_nbody, &global_nbody, 1,
                  MPI_LMP_BIGINT, MPI_SUM, MPI_COMM_WORLD);

    EXPECT_EQ(global_nbody, 0);
}

// Test asymmetric distribution (realistic scenario)
TEST(NbodyMPI, AllreduceAsymmetric)
{
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Rank 0 has most bodies, others have few
    bigint local_nbody;
    if (rank == 0) {
        local_nbody = static_cast<bigint>(INT_MAX) + 1000;
    } else {
        local_nbody = 100;
    }

    bigint global_nbody = 0;
    MPI_Allreduce(&local_nbody, &global_nbody, 1,
                  MPI_LMP_BIGINT, MPI_SUM, MPI_COMM_WORLD);

    bigint expected = static_cast<bigint>(INT_MAX) + 1000 + 100 * (nprocs - 1);
    EXPECT_EQ(global_nbody, expected);
    EXPECT_GT(global_nbody, INT_MAX);
}

}  // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
