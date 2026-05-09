/* ----------------------------------------------------------------------
   Regression for GitHub issues #4923 and #4940:
   dump_modify pbc yes with KOKKOS must remap triclinic coords on the same
   host buffers used for packing (xpbc), not only on device k_x.
------------------------------------------------------------------------- */

#include "../testing/core.h"
#include "../testing/utils.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "domain.h"
#include "fmt/format.h"
#include "info.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <mpi.h>

#include <cstdlib>
#include <exception>
#include <cmath>
#include <string>
#include <vector>

using LAMMPS_NS::Info;

using namespace LAMMPS_NS;

namespace {

void parse_dump_xyz(const std::string &path, double &x, double &y, double &z)
{
  auto lines = read_lines(path);
  ASSERT_FALSE(lines.empty());
  std::size_t i = 0;
  for (; i < lines.size(); ++i) {
    if (lines[i].rfind("ITEM: ATOMS", 0) == 0) break;
  }
  ASSERT_LT(i + 1, lines.size());
  auto w = utils::split_words(lines[i + 1]);
  ASSERT_GE(w.size(), 4u);
  x = std::stod(w[1]);
  y = std::stod(w[2]);
  z = std::stod(w[3]);
}

}    // namespace

class DumpKokkosTriclinicPbcTest : public ::testing::Test {
protected:
  LAMMPS *lmp = nullptr;

  void SetUp() override
  {
    if (!Info::has_package("KOKKOS")) GTEST_SKIP();

    if (Info::has_accelerator_feature("KOKKOS", "api", "cuda") ||
        Info::has_accelerator_feature("KOKKOS", "api", "hip") ||
        Info::has_accelerator_feature("KOKKOS", "api", "sycl"))
      GTEST_SKIP() << "GPU Kokkos build requires GPU; use CPU-only CI preset";

    LAMMPS::argv args = {"LAMMPS_test", "-log", "none", "-echo", "none", "-screen", "none",
                         "-k", "on", "t", "1", "-sf", "kk"};
    if (Info::has_accelerator_feature("KOKKOS", "api", "openmp")) args[10] = "2";

    ::testing::internal::CaptureStdout();
    lmp = new LAMMPS(args, MPI_COMM_WORLD);
    (void) ::testing::internal::GetCapturedStdout();

    // General triclinic primitive cell (same lattice as examples/triclinic/in.fcc.primitive)
    ::testing::internal::CaptureStdout();
    lmp->input->one("units lj");
    lmp->input->one("atom_style atomic");
    lmp->input->one(
        "lattice custom 1.1 a2 0.0 0.5 0.5 a3 0.5 0.0 0.5 a1 0.5 0.5 0.0 basis 0.0 0.0 0.0 "
        "triclinic/general");
    lmp->input->one("create_box 1 NULL 0 1 0 1 0 1");
    lmp->input->one("create_atoms 1 box");
    lmp->input->one("mass * 1.0");
    lmp->input->one("pair_style lj/cut 1.2");
    lmp->input->one("pair_coeff * * 1.0 1.0");
    lmp->input->one("neighbor 0.0 bin");
    lmp->input->one("run 0 post no");
    (void) ::testing::internal::GetCapturedStdout();
  }

  void TearDown() override
  {
    if (lmp) {
      ::testing::internal::CaptureStdout();
      delete lmp;
      (void) ::testing::internal::GetCapturedStdout();
      lmp = nullptr;
    }
  }
};

TEST_F(DumpKokkosTriclinicPbcTest, pbc_yes_matches_reference_general_triclinic)
{
  auto *atomkk = dynamic_cast<LAMMPS_NS::AtomKokkos *>(lmp->atom);
  ASSERT_NE(atomkk, nullptr);

  atomkk->sync(LAMMPS_NS::Host, X_MASK);
  double xref[3] = {lmp->atom->x[0][0], lmp->atom->x[0][1], lmp->atom->x[0][2]};
  double xexpect[3] = {xref[0], xref[1], xref[2]};
  lmp->domain->restricted_to_general_coords(xexpect);

  const std::string dumpfile = "dump_kk_tric_pbc_test.melt";
  delete_file(dumpfile);

  ::testing::internal::CaptureStdout();
  lmp->input->one(fmt::format("dump d0 all custom 1 {} id x y z", dumpfile));
  lmp->input->one("dump_modify d0 pbc yes triclinic/general yes units yes");
  lmp->input->one("run 0 post no");
  lmp->input->one("undump d0");
  (void) ::testing::internal::GetCapturedStdout();

  double xd, yd, zd;
  ASSERT_NO_FATAL_FAILURE(parse_dump_xyz(dumpfile, xd, yd, zd));

  const double tol = 1.0e-10;
  EXPECT_NEAR(xd, xexpect[0], tol);
  EXPECT_NEAR(yd, xexpect[1], tol);
  EXPECT_NEAR(zd, xexpect[2], tol);

  delete_file(dumpfile);
}

// Kokkos OpenMP teardown on macOS throws std::system_error from a static
// destructor after pthreads has partially cleaned up.  Intercept it here so
// the test binary exits with the real Google-Test result code.
// (Same workaround as unittest/force-styles/test_main.cpp.)
static int g_test_result = 1;
static void kokkos_omp_teardown_terminate()
{
  _Exit(g_test_result);
}

int main(int argc, char **argv)
{
  std::set_terminate(kokkos_omp_teardown_terminate);
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleMock(&argc, argv);
  const int rv = RUN_ALL_TESTS();
  g_test_result = rv;
  MPI_Finalize();
  return rv;
}
