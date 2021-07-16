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
static void callback_one(void *handle, step_t timestep, int nlocal, tag_t *, double **, double **f)
{
    for (int i = 0; i < nlocal; ++i)
        f[i][0] = f[i][1] = f[i][2] = (double)timestep;
    lammps_fix_external_set_energy_global(handle, "ext", 1.0);
    double v[6] = {1.0,1.0,1.0,0.0,0.0,0.0 };
    lammps_fix_external_set_virial_global(handle, "ext", v);
}
}

TEST(lammps_external, callback)
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
                                   "fix ext all external pf/callback 5 1\n"
                                   "fix_modify ext energy yes virial yes\n");

    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    ::testing::internal::CaptureStdout();
    lammps_set_fix_external_callback(handle, "ext", &callback_one, handle);
    lammps_command(handle, "run 10 post no");
    double temp = lammps_get_thermo(handle, "temp");
    double pe = lammps_get_thermo(handle, "pe");
    double press = lammps_get_thermo(handle, "press");
    output      = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    EXPECT_DOUBLE_EQ(temp, 1.0 / 30.0);
    EXPECT_DOUBLE_EQ(pe, 1.0 / 8.0);
    EXPECT_DOUBLE_EQ(press, 0.15416666666666667);

    ::testing::internal::CaptureStdout();
    lammps_close(handle);
    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
}

TEST(lammps_external, array)
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
                                   "fix ext all external pf/array 1\n");

    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    ::testing::internal::CaptureStdout();
    double **force = lammps_fix_external_get_force(handle, "ext");
    int nlocal     = lammps_extract_setting(handle, "nlocal");
    for (int i = 0; i < nlocal; ++i)
        force[i][0] = force[i][1] = force[i][2] = 0.0;
    lammps_command(handle, "run 5 post no");
    double temp = lammps_get_thermo(handle, "temp");
    output      = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    EXPECT_DOUBLE_EQ(temp, 4.0 / 525.0);

    ::testing::internal::CaptureStdout();
    nlocal = lammps_extract_setting(handle, "nlocal");
    force  = lammps_fix_external_get_force(handle, "ext");
    for (int i = 0; i < nlocal; ++i)
        force[i][0] = force[i][1] = force[i][2] = 6.0;
    lammps_command(handle, "run 5 post no");
    temp   = lammps_get_thermo(handle, "temp");
    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    EXPECT_DOUBLE_EQ(temp, 1.0 / 30.0);

    ::testing::internal::CaptureStdout();
    lammps_close(handle);
    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
}
