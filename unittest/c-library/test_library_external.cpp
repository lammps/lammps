// unit tests creating LAMMPS instances via the library interface

#include "library.h"

#include <cinttypes>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "test_main.h"

using ::testing::HasSubstr;
using ::testing::StartsWith;

extern "C" {
#ifdef LAMMPS_SMALLSMALL
typedef int32_t step_t;
typedef int32_t tag_t;
#elif LAMMPS_SMALLBIG
typedef int64_t step_t;
typedef int32_t tag_t;
#else
typedef int64_t step_t;
typedef int64_t tag_t;
#endif
static void callback_one(void *lmp, step_t timestep, int nlocal, tag_t *ids, double **x, double **f)
{
    for (int i = 0; i < nlocal; ++i)
        f[i][0] = f[i][1] = f[i][2] = (double)timestep;
}
}

TEST(lammps_external_pf, null_args)
{
    const char *args[] = {"liblammps", "-log", "none", "-nocite"};
    char **argv        = (char **)args;
    int argc           = sizeof(args) / sizeof(char *);

    ::testing::internal::CaptureStdout();
    void *handle       = lammps_open_no_mpi(argc, argv, NULL);
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    ::testing::internal::CaptureStdout();
    lammps_commands_string(handle, "lattice sc 1.0\n"
                                   "region box block -1 1 -1 1 -1 1\n"
                                   "create_box 1 box\n"
                                   "create_atoms 1 box\n"
                                   "mass 1 1.0\n"
                                   "pair_style zero 0.1\n"
                                   "pair_coeff 1 1\n"
                                   "velocity all set 0.1 0.0 -0.1\n"
                                   "thermo 5\n"
                                   "fix 1 all nve\n"
                                   "fix ext all external pf/callback 5 1\n");

    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    ::testing::internal::CaptureStdout();
    lammps_set_fix_external_callback(handle, (char *)"ext", &callback_one, handle);
    lammps_command(handle, "run 10 post no");
    double temp = lammps_get_thermo(handle,"temp");
    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    EXPECT_DOUBLE_EQ(temp,1.0/30.0);

    ::testing::internal::CaptureStdout();
    lammps_close(handle);
    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
}
