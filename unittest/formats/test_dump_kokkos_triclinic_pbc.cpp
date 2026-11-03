/* ----------------------------------------------------------------------
   Regression test for GitHub issues #4923 and #4940:
   "dump_modify pbc yes" hands plain host copies of x, v and image to
   Domain::pbc() and, for a triclinic box, to Domain::x2lamda() and
   Domain::lamda2x().  The KOKKOS versions of those three functions work on
   the Kokkos views instead, so they remapped the real atom data and left the
   copies that are written to the dump file untouched.  The copies were then
   wrapped against the lamda bounds [0,1) rather than the box bounds, which
   moves atoms by up to a box length and increments the wrong image flags.

   Wrapping an atom into the box must not move it, so the unwrapped position
   is the same with and without "pbc yes".  Two dumps of the same timestep,
   one with the remapping and one without, therefore have to write the same
   unwrapped coordinates.  Comparing the two inside a single run keeps the
   test independent of the KOKKOS backend and of the floating point precision
   the package was compiled with, both of which the CI varies.  The atoms have
   to have left the box before the dump for the difference to show up, so a
   short run precedes it.
------------------------------------------------------------------------- */

#include "../testing/core.h"
#include "../testing/utils.h"
#include "info.h"
#include "library.h"
#include "utils.h"

#include "fmt/format.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <mpi.h>

#include <algorithm>
#include <string>
#include <vector>

using LAMMPS_NS::Info;
using LAMMPS_NS::LAMMPS;
using LAMMPS_NS::utils::split_words;

namespace {

// one line of the dump file: atom id and its unwrapped coordinates

struct DumpAtom {
    int id;
    double x, y, z;
};

// read the atom block of a dump file, sorted by atom id

std::vector<DumpAtom> read_dump(const std::string &path)
{
    std::vector<DumpAtom> atoms;
    auto lines = read_lines(path);
    std::size_t i = 0;
    for (; i < lines.size(); ++i)
        if (lines[i].rfind("ITEM: ATOMS", 0) == 0) break;
    for (++i; i < lines.size(); ++i) {
        auto words = split_words(lines[i]);
        if (words.size() < 4) continue;
        atoms.push_back({std::stoi(words[0]), std::stod(words[1]), std::stod(words[2]),
                         std::stod(words[3])});
    }
    std::sort(atoms.begin(), atoms.end(),
              [](const DumpAtom &a, const DumpAtom &b) { return a.id < b.id; });
    return atoms;
}

}    // namespace

TEST(DumpKokkosTriclinicPbc, pbc_yes_keeps_unwrapped_coords)
{
    if (!Info::has_package("KOKKOS")) GTEST_SKIP() << "KOKKOS package not available";
    if (Info::has_accelerator_feature("KOKKOS", "api", "cuda") ||
        Info::has_accelerator_feature("KOKKOS", "api", "hip") ||
        Info::has_accelerator_feature("KOKKOS", "api", "sycl"))
        GTEST_SKIP() << "KOKKOS GPU build needs a GPU, use a CPU-only preset";

    // the serial backend refuses to start with more than one thread

    const bool threaded = (Info::has_accelerator_feature("KOKKOS", "api", "openmp") ||
                           Info::has_accelerator_feature("KOKKOS", "api", "pthreads"));
    const char *nthreads = threaded ? "2" : "1";

    const std::string plainfile = "dump_tric_pbc_plain.melt";
    const std::string remapfile = "dump_tric_pbc_remap.melt";
    delete_file(plainfile);
    delete_file(remapfile);

    LAMMPS::argv args = {"LAMMPS_test", "-log", "none", "-echo",    "none",
                         "-screen",     "none", "-k",   "on",       "t",
                         nthreads,      "-sf",  "kk"};

    ::testing::internal::CaptureStdout();
    auto *lmp = new LAMMPS(args, MPI_COMM_WORLD);
    lmp->input->one("units lj");
    lmp->input->one("atom_style atomic");
    lmp->input->one("boundary p p p");
    lmp->input->one("lattice fcc 0.8442");
    lmp->input->one("region box prism 0 4 0 4 0 4 1.2 0.8 0.5");
    lmp->input->one("create_box 1 box");
    lmp->input->one("create_atoms 1 box");
    lmp->input->one("mass 1 1.0");
    lmp->input->one("pair_style lj/cut 2.5");
    lmp->input->one("pair_coeff 1 1 1.0 1.0 2.5");
    lmp->input->one("neighbor 0.3 bin");
    lmp->input->one("velocity all create 2.0 12345");
    lmp->input->one("fix 1 all nve");

    // let the atoms drift out of the box, then dump the same timestep twice

    lmp->input->one("run 20 post no");
    lmp->input->one(fmt::format("dump plain all custom 1 {} id xu yu zu", plainfile));
    lmp->input->one(fmt::format("dump remap all custom 1 {} id xu yu zu", remapfile));
    lmp->input->one("dump_modify remap pbc yes");
    lmp->input->one("run 0 post no");
    lmp->input->one("undump plain");
    lmp->input->one("undump remap");
    delete lmp;
    (void) ::testing::internal::GetCapturedStdout();

    auto plain = read_dump(plainfile);
    auto remap = read_dump(remapfile);

    // the remapping is done on host copies in double precision, so the
    // unwrapped coordinates have to agree to within rounding of the box
    // lengths that were added and subtracted again

    const double tol = 1.0e-9;
    EXPECT_FALSE(plain.empty());
    ASSERT_EQ(plain.size(), remap.size());
    for (std::size_t i = 0; i < plain.size(); ++i) {
        ASSERT_EQ(plain[i].id, remap[i].id);
        EXPECT_NEAR(plain[i].x, remap[i].x, tol) << "atom " << plain[i].id;
        EXPECT_NEAR(plain[i].y, remap[i].y, tol) << "atom " << plain[i].id;
        EXPECT_NEAR(plain[i].z, remap[i].z, tol) << "atom " << plain[i].id;
    }

    delete_file(plainfile);
    delete_file(remapfile);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);
    const int rv = RUN_ALL_TESTS();

    // finalize the KOKKOS package explicitly: otherwise Kokkos is torn down by
    // static destructors at program exit, leading to segfaults in some cases
    // same workaround as the force-style and FFT3d test drivers

    lammps_kokkos_finalize();

    MPI_Finalize();
    return rv;
}
