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

#include "lammps.h"

#include "atom.h"
#include "comm.h"
#include "group.h"
#include "info.h"
#include "input.h"
#include "math_const.h"
#include "region.h"
#include "variable.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstring>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
using ::testing::ContainsRegex;
using ::testing::ExitedWithCode;
using ::testing::StrEq;

class DeleteAtomsTest : public LAMMPSTest {
protected:
    Atom *atom;

    void SetUp() override
    {
        testbinary = "DeleteAtomsTest";
        args       = {"-log", "none", "-echo", "screen", "-nocite", "-v", "num", "1"};
        LAMMPSTest::SetUp();
        atom = lmp->atom;
    }

    void TearDown() override
    {
        LAMMPSTest::TearDown();
        platform::unlink("test_variable.file");
        platform::unlink("test_variable.atomfile");
    }

    void atomic_system()
    {
        BEGIN_HIDE_OUTPUT();
        command("units real");
        command("lattice sc 1.0 origin 0.125 0.125 0.125");
        command("region box block -4 4 -4 4 -4 4");
        command("create_box 8 box");
        command("create_atoms 1 box");
        command("mass * 1.0");
        command("region left block -2.0 -1.0 INF INF INF INF");
        command("region right block 0.5  2.0 INF INF INF INF");
        command("region top block INF INF -2.0 -1.0 INF INF");
        command("region bottom block INF INF 0.0 4.0 INF INF");
        command("set region left type 2");
        command("set region right type 3");
        command("group bottom region bottom");
        command("group top region top");
        END_HIDE_OUTPUT();
    }

    void molecular_system()
    {
        HIDE_OUTPUT([&] {
            command("fix props all property/atom mol rmass q");
        });
        atomic_system();
        BEGIN_HIDE_OUTPUT();
        command("variable molid atom floor(id/4)+1");
        command("variable charge atom 2.0*sin(PI/32*id)");
        command("set atom * mol v_molid");
        command("set atom * charge v_charge");
        command("set type 1 mass 0.5");
        command("set type 2*4 mass 2.0");
        END_HIDE_OUTPUT();
    }
};

TEST_F(DeleteAtomsTest, Simple)
{
    ASSERT_EQ(atom->map_tag_max, -1);
    HIDE_OUTPUT([&] {
        command("atom_modify map yes");
    });
    atomic_system();
    ASSERT_EQ(atom->natoms, 512);
    ASSERT_EQ(atom->map_tag_max, 512);

    HIDE_OUTPUT([&] {
        command("delete_atoms group top compress no");
    });
    ASSERT_EQ(atom->natoms, 448);
    ASSERT_EQ(atom->map_tag_max, 512);

    HIDE_OUTPUT([&] {
        command("delete_atoms region left");
    });
    ASSERT_EQ(atom->natoms, 392);
    ASSERT_EQ(atom->map_tag_max, 392);

    HIDE_OUTPUT([&] {
        command("delete_atoms random fraction 0.5 yes all right 43252");
    });
    ASSERT_EQ(atom->natoms, 364);

    HIDE_OUTPUT([&] {
        command("variable checker atom sin(4*PI*x/lx)*sin(4*PI*y/ly)*sin(4*PI*z/lz)>0");
        command("delete_atoms variable checker");
    });
    ASSERT_EQ(atom->natoms, 178);

    HIDE_OUTPUT([&] {
        command("delete_atoms random count 3 no bottom right 443252");
    });
    ASSERT_EQ(atom->natoms, 175);

    HIDE_OUTPUT([&] {
        command("delete_atoms random count 50 no all NULL 434325");
    });
    ASSERT_EQ(atom->natoms, 125);

    HIDE_OUTPUT([&] {
        command("delete_atoms random fraction 0.2 no all NULL 34325");
    });
    ASSERT_EQ(atom->natoms, 104);

    HIDE_OUTPUT([&] {
        command("delete_atoms random count 50 no bottom right 77325");
    });
    ASSERT_EQ(atom->natoms, 102);

    TEST_FAILURE(".*ERROR: Illegal delete_atoms command: missing argument.*",
                 command("delete_atoms"););
    TEST_FAILURE(".*ERROR: Unknown delete_atoms sub-command: xxx.*", command("delete_atoms xxx"););
    TEST_FAILURE(".*ERROR: The delete_atoms 'porosity' keyword has been removed.*",
                 command("delete_atoms porosity 0.5 all right 4325234"););
    TEST_FAILURE(".*ERROR: Illegal delete_atoms random command: missing argument.*",
                 command("delete_atoms random count 50 bottom right 77325"););
    TEST_FAILURE(".*ERROR: Illegal delete_atoms random command: missing argument.*",
                 command("delete_atoms random fraction 0.4 bottom right 77325"););
    TEST_FAILURE(".*ERROR: Delete_atoms random count has invalid value: -5.*",
                 command("delete_atoms random count -5 no bottom right 77325"););
    TEST_FAILURE(".*ERROR: Delete_atoms count of 5 exceeds number of eligible atoms 0.*",
                 command("delete_atoms random count 5 yes bottom right 77325"););
    TEST_FAILURE(".*ERROR: Delete_atoms random fraction has invalid value: -0.4.*",
                 command("delete_atoms random fraction -0.4 no bottom right 77325"););
}

TEST_F(DeleteAtomsTest, CondenseAtomic)
{
    HIDE_OUTPUT([&] {
        command("atom_modify id yes map yes ");
        command("region box block -4 4 -4 4 -4 4");
        command("create_box 1 box");
        command("mass 1 1.0");
        command("fix oldid all property/atom i_oldid ghost yes");
        command("create_atoms 1 random 8 9648523 box");
        command("variable oldid atom id");
        command("set atom * i_oldid v_oldid");
        command("group odd id 1:8:2");
    });
    ASSERT_EQ(atom->natoms, 8);
    ASSERT_EQ(atom->map_tag_max, 8);
    int flag = -1, cols = -1;
    auto index_custom = atom->find_custom("oldid", flag, cols);
    ASSERT_TRUE(index_custom >= 0);
    ASSERT_EQ(flag, 0);
    ASSERT_EQ(cols, 0);

    HIDE_OUTPUT([&] {
        command("delete_atoms group odd compress no");
    });
    ASSERT_EQ(atom->natoms, 4);
    ASSERT_EQ(atom->map_tag_max, 8);
    for (int i = 0; i < 4; ++i)
        ASSERT_EQ(atom->ivector[index_custom][i], atom->tag[i]);

    HIDE_OUTPUT([&] {
        command("delete_atoms group all");
        command("create_atoms 1 random 8 9648523 box");
        command("variable oldid atom id");
        command("set atom * i_oldid v_oldid");
        command("group odd delete");
        command("group odd id 1:8:2");
    });

    HIDE_OUTPUT([&] {
        command("delete_atoms group odd condense yes");
    });
    ASSERT_EQ(atom->natoms, 4);
    ASSERT_EQ(atom->map_tag_max, 4);
    for (int i = 0; i < 4; ++i)
        ASSERT_EQ(atom->ivector[index_custom][i]/2, atom->tag[i]);
}

TEST_F(DeleteAtomsTest, CondenseMolecular)
{
    // Create a 3-atom linear molecule template
    if (lmp->comm->me == 0) {
        FILE *fp = fopen("tri.mol", "w");
        fprintf(fp, "# Linear Triatomic\n");
        fprintf(fp, "3 atoms\n2 bonds\n1 angles\n");
        fprintf(fp, "\nCoords\n\n1 0 0 0\n2 1 0 0\n3 2 0 0\n");
        fprintf(fp, "\nTypes\n\n1 1\n2 1\n3 1\n");
        fprintf(fp, "\nBonds\n\n1 1 1 2\n2 1 2 3\n");
        fprintf(fp, "\nAngles\n\n1 1 1 2 3\n");
        fclose(fp);
    }
    MPI_Barrier(lmp->world);

    HIDE_OUTPUT([&] {
        command("units lj");
        command("atom_style full");
        command("atom_modify map yes");
        command("region box block 0 20 0 20 0 20");
        command("create_box 1 box bond/types 1 angle/types 1 extra/bond/per/atom 2 extra/angle/per/atom 1 extra/special/per/atom 2");
        command("mass 1 1.0");
        command("bond_style harmonic");
        command("bond_coeff 1 100.0 1.0");
        command("angle_style harmonic");
        command("angle_coeff 1 100.0 180.0");

        command("molecule tri tri.mol");

        // Create 4 molecules (12 atoms total)
        // IDs: 1-3, 4-6, 7-9, 10-12
        command("create_atoms 0 single 2 2 2 mol tri 12345");
        command("create_atoms 0 single 5 2 2 mol tri 12346");
        command("create_atoms 0 single 8 2 2 mol tri 12347");
        command("create_atoms 0 single 11 2 2 mol tri 12348");
    });

    ASSERT_EQ(atom->natoms, 12);
    ASSERT_EQ(atom->nbonds, 8);
    ASSERT_EQ(atom->nangles, 4);

    // Delete middle molecules (IDs 4-9)
    HIDE_OUTPUT([&] {
        command("group del id 4:9");
        command("delete_atoms group del compress no");
    });

    // Remaining IDs: 1,2,3, 10,11,12
    ASSERT_EQ(atom->natoms, 6);

    // Condense IDs -> Should become 1,2,3, 4,5,6
    HIDE_OUTPUT([&] {
        command("delete_atoms group all condense yes");
    });

    // Verify Topology Mapping
    // Old ID 11 (center of last mol) should become New ID 5.
    // It should be bonded to Old 10 (New 4) and Old 12 (New 6).
    int nlocal = atom->nlocal;
    for(int i=0; i<nlocal; i++) {
        tagint id = atom->tag[i];
        ASSERT_LE(id, 6);
        ASSERT_GE(id, 1);

        if (id == 5) {
            int found4=0, found6=0;
            for(int m=0; m<atom->num_bond[i]; m++) {
                if (atom->bond_atom[i][m] == 4) found4=1;
                if (atom->bond_atom[i][m] == 6) found6=1;
            }
            EXPECT_EQ(found4, 1);
            EXPECT_EQ(found6, 1);
        }
    }

    if (lmp->comm->me == 0) platform::unlink("tri.mol");
}

TEST_F(DeleteAtomsTest, CondenseFixes)
{

#if defined(LMP_MOLECULE)

    // ------------------------------------------------------------------
    // 3. FIX CMAP
    // ------------------------------------------------------------------
    // Note: molecule command does not support 'Crossterms' section natively
    // for read-in, so we must use read_data to create a system with CMAP.

    if (lmp->comm->me == 0) {
        // Write a small Data file with CMAP section
        FILE *fp = fopen("data.cmap", "w");
        fprintf(fp, "LAMMPS CMAP test\n\n");
        fprintf(fp, "8 atoms\n6 bonds\n4 angles\n2 dihedrals\n1 crossterms\n\n");
        fprintf(fp, "1 atom types\n1 bond types\n1 angle types\n1 dihedral types\n1 crossterm types\n\n");
        fprintf(fp, "0 20 xlo xhi\n0 20 ylo yhi\n0 20 zlo zhi\n\n");

        fprintf(fp, "Masses\n\n1 1.0\n\n");

        // Two linear chains of 4 atoms: 1-2-3-4 and 5-6-7-8
        fprintf(fp, "Atoms\n\n");
        fprintf(fp, "1 1 1 0.0 0 0 0\n");
        fprintf(fp, "2 1 1 0.0 1 0 0\n");
        fprintf(fp, "3 1 1 0.0 2 0 0\n");
        fprintf(fp, "4 1 1 0.0 3 0 0\n");
        fprintf(fp, "5 2 1 0.0 5 0 0\n");
        fprintf(fp, "6 2 1 0.0 6 0 0\n");
        fprintf(fp, "7 2 1 0.0 7 0 0\n");
        fprintf(fp, "8 2 1 0.0 8 0 0\n");

        fprintf(fp, "Bonds\n\n");
        fprintf(fp, "1 1 1 2\n2 1 2 3\n3 1 3 4\n");
        fprintf(fp, "4 1 5 6\n5 1 6 7\n6 1 7 8\n");

        fprintf(fp, "Angles\n\n");
        fprintf(fp, "1 1 1 2 3\n2 1 2 3 4\n");
        fprintf(fp, "3 1 5 6 7\n4 1 6 7 8\n");

        fprintf(fp, "Dihedrals\n\n");
        fprintf(fp, "1 1 1 2 3 4\n");
        fprintf(fp, "2 1 5 6 7 8\n");

        // One CMAP term involving the second molecule (5-6-7-8)
        // If we delete atoms 1-4, atoms 5-8 become 1-4.
        // Format: ID Type A B C D E
        fprintf(fp, "Crossterms\n\n");
        fprintf(fp, "1 1 5 6 7 8 5\n"); // 5-body term
        fclose(fp);

        // Write minimal dummy CMAP coefficient file
        fp = fopen("cmap.coef", "w");
        fprintf(fp, "cmap\n");
        fprintf(fp, "1 1 1\n"); // Ntypes
        fprintf(fp, "1 24 24\n"); // Type GridDims
        for(int k=0; k<24*24; k++) fprintf(fp, "0.0 ");
        fprintf(fp, "\n");
        fclose(fp);
    }
    MPI_Barrier(lmp->world);

    HIDE_OUTPUT([&] {
        command("clear");
        command("units real");
        command("atom_style full");
        command("atom_modify map yes");

        // Define styles supporting CMAP
        command("bond_style harmonic");
        command("angle_style harmonic");
        command("dihedral_style charmm");
        command("improper_style harmonic");
        command("crossterm_style cmap");

        command("read_data data.cmap");

        command("bond_coeff 1 100.0 1.0");
        command("angle_coeff 1 50.0 180.0");
        command("dihedral_coeff 1 1.0 1 180.0 0.0");

        // Fix CMAP using the dummy file
        command("fix 3 all cmap charmm cmap.coef");
    });

    // Delete the first molecule (atoms 1-4)
    HIDE_OUTPUT([&] {
        command("group firstmol id 1 2 3 4");
        command("delete_atoms group firstmol compress no");
    });

    ASSERT_EQ(atom->natoms, 4); // 5,6,7,8 remain

    // Condense -> 5,6,7,8 become 1,2,3,4
    HIDE_OUTPUT([&] {
        command("delete_atoms group all condense yes");
    });

    // Run 0 to ensure FixCMAP::post_force or setup logic runs with new IDs
    // without crashing (which checks if atoms exist)
    HIDE_OUTPUT([&] {
        command("run 0");
    });

    // Cleanup
    command("unfix 3");

    if (lmp->comm->me == 0) {
        platform::unlink("tri_shake.mol");
        platform::unlink("data.cmap");
        platform::unlink("cmap.coef");
    }

#endif // defined(LMP_MOLECULE)

#if defined(LMP_RIGID)
    // ------------------------------------------------------------------
    // 1. FIX SHAKE
    // ------------------------------------------------------------------
    if (lmp->comm->me == 0) {
        FILE *fp = fopen("tri_shake.mol", "w");
        fprintf(fp, "# Linear Triatomic\n");
        fprintf(fp, "3 atoms\n2 bonds\n1 angles\n");
        fprintf(fp, "\nCoords\n\n1 0 0 0\n2 1 0 0\n3 2 0 0\n");
        fprintf(fp, "\nTypes\n\n1 1\n2 1\n3 1\n");
        fprintf(fp, "\nBonds\n\n1 1 1 2\n2 1 2 3\n");
        fprintf(fp, "\nAngles\n\n1 1 1 2 3\n");
        fprintf(fp, "\nShake Flags\n\n1 1\n2 1\n3 1\n");
        fprintf(fp, "\nShake Atoms\n\n1 1 2 3\n2 1 2 3\n3 1 2 3\n");
        fprintf(fp, "\nShake Bond Types\n\n1 1 1 1\n2 1 1 1\n3 1 1 1\n");
        fclose(fp);
    }
    MPI_Barrier(lmp->world);

    HIDE_OUTPUT([&] {
        command("clear");
        command("units lj");
        command("atom_style full");
        command("atom_modify map yes");
        command("region box block 0 20 0 20 0 20");
        command("create_box 1 box bond/types 1 angle/types 1 extra/bond/per/atom 2 extra/angle/per/atom 1 extra/special/per/atom 2");
        command("mass 1 1.0");
        command("bond_style harmonic");
        command("bond_coeff 1 100.0 1.0");
        command("angle_style harmonic");
        command("angle_coeff 1 100.0 180.0");

        command("molecule tri_s tri_shake.mol");

        command("create_atoms 0 single 2 2 2 mol tri_s 12345");
        command("create_atoms 0 single 5 2 2 mol tri_s 12346");
        command("create_atoms 0 single 8 2 2 mol tri_s 12347");

        command("fix 1 all shake 0.0001 20 0 a 1");
    });

    HIDE_OUTPUT([&] {
        command("group del id 4:6");
        command("delete_atoms group del compress no");
        command("delete_atoms group all condense yes");
        command("run 0"); // Trigger fix setup/checks
    });

    command("unfix 1");

    // ------------------------------------------------------------------
    // 2. FIX RIGID/SMALL
    // ------------------------------------------------------------------
    HIDE_OUTPUT([&] {
        command("delete_atoms group all");
        command("create_atoms 0 single 2 2 2 mol tri_s 22345");
        command("create_atoms 0 single 5 2 2 mol tri_s 22346");
        command("create_atoms 0 single 8 2 2 mol tri_s 22347");

        command("fix 2 all rigid/small molecule");
        command("group del2 id 4:6");
    });

    HIDE_OUTPUT([&] {
        command("delete_atoms group del2 compress no");
        command("delete_atoms group all condense yes");
        command("run 0");
    });

    command("unfix 2");

#endif // defined(LMP_RIGID)

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
    MPI_Finalize();
    return rv;
}

