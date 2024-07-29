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
#include "atom.h"
#include "fmt/format.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "library.h"
#include "output.h"
#include "update.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <random>
#include <vector>
#include <mpi.h>

using ::testing::Eq;
using testing::HasSubstr;

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

#define GETIDX(i) lmp->atom->map(i)

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

class PsfTest : public LAMMPSTest {
protected:

  void create_fourmol_psf()
  {
    FILE *fp = fopen("fourmol.psf", "w");

    fmt::print(fp,"PSF EXT XPLOR\n\n         5 !NTITLE\n");

    fmt::print(fp,"* LAMMPS psf file via write_psf, version {}, timestep = {}, units = {}\n",
             lmp->version, lmp->update->ntimestep, lmp->update->unit_style);

    fmt::print(fp,"* II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I)\n");
    fmt::print(fp,"* expanded format EXT:\n");
    fmt::print(fp,"* (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8) XPLOR\n");
    fmt::print(fp,"* [https://userguide.mdanalysis.org/stable/formats/reference/psf.html]\n");

    std::vector<const char *> segments = {"PROA","DNAA","SOLV","IONS"};
    std::vector<const char *> residues = {"GLY","ASP","GLU"};
    std::vector<const char *> names = {"FOO","BAR","QUUX"};
    std::vector<const char *> types = {"C","HE","N","OA","OB"};

    std::default_random_engine generator;
    std::uniform_int_distribution<int> rnd_segment(0,3);
    std::uniform_int_distribution<int> rnd_residue(0,2);
    std::uniform_int_distribution<int> rnd_name(0,2);

    fmt::print(fp,"\n {:8} !NATOM\n", 29);

    for( int i=0; i<29 ; i++ ) {
      auto atom_tag = lmp->atom->tag[i];
      auto atom_type = lmp->atom->type[GETIDX(atom_tag)];
      //utils::logmesg(lmp, "create_fourmol_psf()... atom_tag={} atom_type={}\n", atom_tag, atom_type);
      fmt::print(fp, "{:10} ", atom_tag);
      fmt::print(fp, "{0:<8} ", segments[rnd_segment(generator)] );
      fmt::print(fp, "{0:<8} ", lmp->atom->molecule[i] );
      fmt::print(fp, "{0:<8} ", residues[rnd_residue(generator)] );
      fmt::print(fp, "{0:<8} ", names[rnd_name(generator)] );
      fmt::print(fp, "{0:<4} ", types[atom_type-1] );
      fmt::print(fp, "{:12.6F}      ", lmp->atom->q[i] );
      fmt::print(fp, "{:8g}           0\n", lmp->atom->mass[atom_type] );
    }

    fclose(fp);
  }

  void SetUp() override
  {
    testbinary = "PsfTest";
    LAMMPSTest::SetUp();
    if (info->has_style("atom", "full")) {
      BEGIN_HIDE_OUTPUT();
      command("variable input_dir index \"" STRINGIFY(TEST_INPUT_FOLDER) "\"");
      command("include \"${input_dir}/in.fourmol\"");
      create_fourmol_psf();
      END_HIDE_OUTPUT();
    }
  }

  void TearDown() override
  {
    platform::unlink("fourmol.psf");
    platform::unlink("fourmol-all.psf");
    platform::unlink("fourmol-oxygen.psf");
  }

};

TEST_F(PsfTest, FailWritePsfBeforeReadPsf)
{
    TEST_FAILURE(".*ERROR: write_psf command needs psf_segment_residue_name per-atom custom array created by read_psf command.", command("write_psf all fourmol.psf"););
}

TEST_F(PsfTest, ReadWritePsfGroupAll)
{
    BEGIN_HIDE_OUTPUT();
    command("read_psf fourmol.psf");
    command("write_psf all fourmol-all.psf");
    END_HIDE_OUTPUT();

    ASSERT_FILE_EXISTS("fourmol-all.psf");
    ASSERT_FILE_EQUAL("fourmol.psf", "fourmol-all.psf");
}

TEST_F(PsfTest, ReadWritePsfGroupSubset)
{
    BEGIN_HIDE_OUTPUT();
    command("read_psf fourmol.psf");
    command("group oxygen type OA OB");
    command("write_psf oxygen fourmol-oxygen.psf");
    END_HIDE_OUTPUT();

    ASSERT_FILE_EXISTS("fourmol-oxygen.psf");
    auto lines = read_lines("fourmol-oxygen.psf");
    ASSERT_EQ(lines.size(), 57);
    ASSERT_THAT(lines[9], Eq("        6 !NATOM"));

    for( int i=10 ; i<=15 ; i++ )
        ASSERT_THAT(lines[i], AnyOf(HasSubstr("OA"),HasSubstr("OB")));
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
