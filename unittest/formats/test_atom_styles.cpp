/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom.h"
#include "input.h"
#include "lammps.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstring>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
using ::testing::Eq;

class AtomStyleTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {"SimpleCommandsTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(AtomStyleTest, atomic)
{
    ASSERT_THAT(std::string(lmp->atom->atom_style), Eq("atomic"));
    ASSERT_NE(lmp->atom->avec, nullptr);
    ASSERT_EQ(lmp->atom->natoms, 0);
    ASSERT_EQ(lmp->atom->nlocal, 0);
    ASSERT_EQ(lmp->atom->nghost, 0);
    ASSERT_EQ(lmp->atom->nmax, 1);
    ASSERT_EQ(lmp->atom->tag_enable, 1);
    ASSERT_EQ(lmp->atom->molecular, 0);
    ASSERT_EQ(lmp->atom->nellipsoids, 0);
    ASSERT_EQ(lmp->atom->nlines, 0);
    ASSERT_EQ(lmp->atom->ntris, 0);
    ASSERT_EQ(lmp->atom->nbodies, 0);
    ASSERT_EQ(lmp->atom->nbonds, 0);
    ASSERT_EQ(lmp->atom->nangles, 0);
    ASSERT_EQ(lmp->atom->ndihedrals, 0);
    ASSERT_EQ(lmp->atom->nimpropers, 0);
    ASSERT_EQ(lmp->atom->ntypes, 0);
    ASSERT_EQ(lmp->atom->nbondtypes, 0);
    ASSERT_EQ(lmp->atom->nangletypes, 0);
    ASSERT_EQ(lmp->atom->ndihedraltypes, 0);
    ASSERT_EQ(lmp->atom->nimpropertypes, 0);
    ASSERT_EQ(lmp->atom->bond_per_atom, 0);
    ASSERT_EQ(lmp->atom->angle_per_atom, 0);
    ASSERT_EQ(lmp->atom->dihedral_per_atom, 0);
    ASSERT_EQ(lmp->atom->improper_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_bond_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_angle_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_dihedral_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_improper_per_atom, 0);

    ASSERT_EQ(lmp->atom->tag, nullptr);
    ASSERT_EQ(lmp->atom->type, nullptr);
    ASSERT_EQ(lmp->atom->mask, nullptr);
    ASSERT_EQ(lmp->atom->image, nullptr);
    ASSERT_EQ(lmp->atom->x, nullptr);
    ASSERT_EQ(lmp->atom->v, nullptr);
    ASSERT_EQ(lmp->atom->f, nullptr);
    ASSERT_EQ(lmp->atom->rmass, nullptr);
    ASSERT_EQ(lmp->atom->q, nullptr);
    ASSERT_EQ(lmp->atom->mu, nullptr);
    //    ASSERT_NE(lmp->atom->,nullptr);
    //    ASSERT_EQ(lmp->atom->firstgroup, -1);
    //    ASSERT_EQ(lmp->atom->firstgroupname, nullptr);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("atom_style charge");
    lmp->input->one("atom_style atomic");
    if (!verbose) ::testing::internal::GetCapturedStdout();
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);
    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;
    return RUN_ALL_TESTS();
}
