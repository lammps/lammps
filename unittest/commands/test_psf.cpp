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
#include "error.h"
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

/* ----------------------------------------------------------------------
   comparison function invoked by merge_sort()
------------------------------------------------------------------------- */

int compare_tags(const int i, const int j, void *ptr)
{
  tagint *tags = (tagint *)ptr;
  if (tags[i] < tags[j]) return -1;
  else return 1;
}

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

    tagint *tags = new tagint[29];
    std::vector<std::string> psf_lines;

    for( int i=0; i<29 ; i++ ) {
      auto atom_tag = lmp->atom->tag[i];
      auto atom_type = lmp->atom->type[GETIDX(atom_tag)];
      //utils::logmesg(lmp, "create_fourmol_psf()... atom_tag={} atom_type={}\n", atom_tag, atom_type);
      tags[i] = atom_tag;
      psf_lines.push_back(fmt::format("{:10} {:<8} {:<8} {:<8} {:<8} {:<4} {:12.6F}      {:8g}           0\n", atom_tag,
        segments[rnd_segment(generator)],
        lmp->atom->molecule[i],
        residues[rnd_residue(generator)],
        names[rnd_name(generator)],
        types[atom_type-1],
        lmp->atom->q[i],
        lmp->atom->mass[atom_type] ));

    }

    int *order = new int[29];
    for (int i = 0; i < 29; i++) order[i] = i;
    utils::merge_sort(order, 29, (void *)tags, compare_tags);

    for( int i=0; i<29 ; i++ )
      fmt::print(fp, psf_lines[order[i]] );

    fmt::print(fp, R"(
        24 !NBOND: bonds
        10        11        10        12        10        16        12        13
        12        14        12        15         3         4         3         5
         3         6         6         8         6         7        18        19
        18        20         8         9         8        10        16        17
         1         2         1         3        21        22        21        23
        24        25        24        26        27        28        27        29

        30 !NTHETA: angles
         8        10        11         8        10        16        11        10        12
        12        10        16         8        10        12        11        10        16
        10        12        15        10        12        14        10        12        13
        13        12        15        13        12        14        14        12        15
         1         3         5         1         3         4         1         3         6
         4         3         5         5         3         6         4         3         6
         3         6         7         3         6         8         7         6         8
        19        18        20         6         8         9         9         8        10
         6         8        10        10        16        17         2         1         3
        22        21        23        25        24        26        28        27        29

        31 !NPHI: dihedrals
         8        10        12        13         8        10        12        14
         8        10        12        15         8        10        16        17
        11        10        12        13        11        10        12        14
        11        10        12        15        11        10        16        17
        12        10        16        17        16        10        12        13
        16        10        12        14        16        10        12        15
         1         3         6         8         1         3         6         7
         4         3         6         8         4         3         6         7
         5         3         6         8         5         3         6         7
         3         6         8         9         3         6         8        10
         7         6         8         9         7         6         8        10
         6         8        10        12         6         8        10        16
         6         8        10        11         9         8        10        12
         9         8        10        16         9         8        10        11
         2         1         3         6         2         1         3         4
         2         1         3         5

         2 !NIMPHI: impropers
         6         3         8         7         8         6        10         9
)");

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
