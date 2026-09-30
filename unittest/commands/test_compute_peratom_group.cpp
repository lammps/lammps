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

// Tests that per-atom computes report zeros for atoms outside the compute
// group.  Several computes only assigned values to the atoms in the group,
// so the entries of all other atoms kept whatever the per-atom array held:
// uninitialized memory after an allocation, or the values of a previous
// invocation.  The latter does not depend on the memory allocator, so each
// test evaluates the compute twice and swaps the lower and upper half of
// the system in and out of the compute group in between.  The atoms that
// left the group must then report zeros instead of their previous values.
// The same executable is registered for the plain, OPENMP and KOKKOS
// versions of the computes and as a test with 4 MPI ranks.

#include "../testing/core.h"
#include "../testing/mpitesting.h"

#include "atom.h"
#include "compute.h"
#include "fmt/format.h"
#include "group.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "library.h"
#include "modify.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include <string>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

class ComputePeratomGroupTest : public LAMMPSTest {
protected:
    // the compute group "sub" starts out as the lower half of the box
    void define_groups(double zsplit)
    {
        HIDE_OUTPUT([&] {
            command(fmt::format("region lower block INF INF INF INF INF {} units box", zsplit));
            command(fmt::format("region upper block INF INF INF INF {} INF units box", zsplit));
            command("group sub region lower");
        });
    }

    void swap_groups()
    {
        HIDE_OUTPUT([&] {
            command("group sub clear");
            command("group sub region upper");
        });
    }

    // a slightly distorted fcc lattice of Lennard-Jones atoms
    void lj_lattice()
    {
        HIDE_OUTPUT([&] {
            command("units lj");
            command("atom_style atomic");
            command("atom_modify map array sort 0 0.0");
            command("lattice fcc 0.8442");
            command("region box block 0 4 0 4 0 4");
            command("create_box 1 box");
            command("create_atoms 1 box");
            command("mass 1 1.0");
            command("displace_atoms all random 0.05 0.05 0.05 12345 units box");
            command("pair_style lj/cut 2.5");
            command("pair_coeff 1 1 1.0 1.0");
            command("neighbor 1.0 bin");
        });
        define_groups(3.3);
    }

    // evaluate the compute and check its per-atom data: atoms outside the
    // group "sub" must report exactly zero in all columns.  atoms in the
    // group must report some non-zero data, since otherwise values left
    // over from the previous evaluation could not be detected.  with
    // multiple MPI ranks, the group atoms may all be owned by other ranks.
    void check_peratom(const std::string &id)
    {
        HIDE_OUTPUT([&] {
            command("run 0 post no");
        });
        auto *compute = lmp->modify->get_compute_by_id(id);
        ASSERT_NE(compute, nullptr);
        compute->compute_peratom();

        const int groupbit = lmp->group->get_bitmask_by_id(FLERR, "sub", "test");
        const int ncols    = compute->size_peratom_cols;
        const int *mask    = lmp->atom->mask;
        const auto *tag    = lmp->atom->tag;
        int nonzero        = 0;
        for (int i = 0; i < lmp->atom->nlocal; ++i) {
            for (int m = 0; m < ((ncols == 0) ? 1 : ncols); ++m) {
                const double value = (ncols == 0) ? compute->vector_atom[i]
                                                  : compute->array_atom[i][m];
                if (mask[i] & groupbit) {
                    if (value != 0.0) ++nonzero;
                } else {
                    EXPECT_EQ(value, 0.0) << "compute " << compute->style << " atom " << tag[i]
                                          << " column " << m + 1;
                }
            }
        }
        int allnonzero = 0;
        MPI_Allreduce(&nonzero, &allnonzero, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        EXPECT_GT(allnonzero, 0) << "compute " << compute->style;
    }

    // check before and after swapping the atoms in and out of the group
    void check_swap(const std::string &id)
    {
        check_peratom(id);
        swap_groups();
        check_peratom(id);
    }
};

TEST_F(ComputePeratomGroupTest, centro_atom_axes)
{
    lj_lattice();
    HIDE_OUTPUT([&] {
        command("compute test sub centro/atom fcc axes yes");
    });
    check_swap("test");
}

TEST_F(ComputePeratomGroupTest, entropy_atom)
{
    if (!Info::has_package("EXTRA-COMPUTE")) GTEST_SKIP();
    lj_lattice();
    HIDE_OUTPUT([&] {
        command("compute test sub entropy/atom 0.25 2.0");
    });
    check_swap("test");
}

TEST_F(ComputePeratomGroupTest, entropy_atom_avg)
{
    if (!Info::has_package("EXTRA-COMPUTE")) GTEST_SKIP();
    lj_lattice();
    HIDE_OUTPUT([&] {
        command("compute test sub entropy/atom 0.25 2.0 avg yes 1.3");
    });
    check_swap("test");
}

// the average over the neighbors of an atom in the group includes all its
// neighbors, so the values must be the same as for the group "all"

TEST_F(ComputePeratomGroupTest, entropy_atom_avg_neighbors)
{
    if (!Info::has_package("EXTRA-COMPUTE")) GTEST_SKIP();
    lj_lattice();
    HIDE_OUTPUT([&] {
        command("compute eall all entropy/atom 0.25 2.0 avg yes 1.3");
        command("compute esub sub entropy/atom 0.25 2.0 avg yes 1.3");
        command("run 0 post no");
    });
    auto *eall = lmp->modify->get_compute_by_id("eall");
    auto *esub = lmp->modify->get_compute_by_id("esub");
    ASSERT_NE(eall, nullptr);
    ASSERT_NE(esub, nullptr);
    eall->compute_peratom();
    esub->compute_peratom();

    const int groupbit = lmp->group->get_bitmask_by_id(FLERR, "sub", "test");
    const int *mask    = lmp->atom->mask;
    int ningroup       = 0;
    int allingroup     = 0;
    for (int i = 0; i < lmp->atom->nlocal; ++i) {
        if (mask[i] & groupbit) {
            ++ningroup;
            EXPECT_DOUBLE_EQ(esub->vector_atom[i], eall->vector_atom[i])
                << "atom " << lmp->atom->tag[i];
        } else {
            EXPECT_EQ(esub->vector_atom[i], 0.0) << "atom " << lmp->atom->tag[i];
        }
    }
    MPI_Allreduce(&ningroup, &allingroup, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_GT(allingroup, 0);
}

TEST_F(ComputePeratomGroupTest, dpd_atom)
{
    if (!Info::has_package("DPD-REACT")) GTEST_SKIP();
    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style dpd");
        command("atom_modify map array sort 0 0.0");
        command("lattice sc 1.0");
        command("region box block 0 4 0 4 0 4");
        command("create_box 1 box");
        command("create_atoms 1 box");
        command("mass 1 1.0");
        command("set group all dpd/theta 1.5");
        command("pair_style zero 1.0");
        command("pair_coeff * *");
    });
    define_groups(1.9);
    HIDE_OUTPUT([&] {
        command("compute test sub dpd/atom");
    });
    check_swap("test");
}

// a peridynamic bar that is stretched along z until the bonds yield, so
// that the dilatation and the plastic stretch are non-zero

TEST_F(ComputePeratomGroupTest, peri_atom)
{
    if (!Info::has_package("PERI")) GTEST_SKIP();
    if (lmp->kokkos) GTEST_SKIP() << "atom style peri has no KOKKOS version";
    HIDE_OUTPUT([&] {
        command("units lj");
        command("boundary s s s");
        command("atom_style peri");
        command("atom_modify map array sort 0 0.0");
        command("lattice sc 1.0");
        command("region box block 0 3 0 3 0 10 units box");
        command("create_box 1 box");
        command("create_atoms 1 box");
        command("pair_style peri/eps");
        command("pair_coeff * * 1.0 1.0 3.01 0.5 0.25 1.0e-4");
        command("set group all density 1.0");
        command("set group all volume 1.0");
        command("fix 1 all nve");
        command("timestep 1.0e-4");
        command("run 0 post no");
        command("displace_atoms all ramp z -0.1 0.1 z 0.0 10.0 units box");
        command("run 1 post no");
    });
    define_groups(4.5);
    HIDE_OUTPUT([&] {
        command("compute dil sub dilatation/atom");
        command("compute pla sub plasticity/atom");
    });
    check_peratom("dil");
    check_peratom("pla");
    swap_groups();
    check_peratom("dil");
    check_peratom("pla");
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

    // print the test results only once with multiple MPI ranks

    auto &listeners       = UnitTest::GetInstance()->listeners();
    auto default_listener = listeners.Release(listeners.default_result_printer());
    listeners.Append(new MPIPrinter(default_listener));

    int rv = RUN_ALL_TESTS();

    // finalize the KOKKOS package explicitly: otherwise Kokkos is torn down by
    // static destructors at program exit, leading to segfaults in some cases
    // same workaround as the force-style and FFT3d test drivers

    lammps_kokkos_finalize();

    MPI_Finalize();
    return rv;
}
