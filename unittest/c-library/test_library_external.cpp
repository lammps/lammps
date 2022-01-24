// unit tests for interfacing with fix external via the library interface

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
static void callback(void *handle, step_t timestep, int nlocal, tag_t *, double **, double **f)
{
    for (int i = 0; i < nlocal; ++i)
        f[i][0] = f[i][1] = f[i][2] = (double)timestep;

    double v[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
    lammps_fix_external_set_virial_global(handle, "ext", v);
    if (timestep < 10) {
        lammps_fix_external_set_energy_global(handle, "ext", 0.5);
        lammps_fix_external_set_vector(handle, "ext", 1, timestep);
        lammps_fix_external_set_vector(handle, "ext", 3, 1.0);
        lammps_fix_external_set_vector(handle, "ext", 4, -0.25);
    } else {
        lammps_fix_external_set_energy_global(handle, "ext", 1.0);
        lammps_fix_external_set_vector(handle, "ext", 2, timestep);
        lammps_fix_external_set_vector(handle, "ext", 5, -1.0);
        lammps_fix_external_set_vector(handle, "ext", 6, 0.25);
    }
    double *eatom  = new double[nlocal];
    double **vatom = new double *[nlocal];
    vatom[0]       = new double[nlocal * 6];
    eatom[0]       = 0.0;
    vatom[0][0] = vatom[0][1] = vatom[0][2] = vatom[0][3] = vatom[0][4] = vatom[0][5] = 0.0;

    for (int i = 1; i < nlocal; ++i) {
        eatom[i]    = 0.1 * i;
        vatom[i]    = vatom[0] + 6 * i;
        vatom[i][0] = vatom[i][1] = vatom[i][2] = 0.1;
        vatom[i][3] = vatom[i][4] = vatom[i][5] = -0.2;
    }
    lammps_fix_external_set_energy_peratom(handle, "ext", eatom);
    lammps_fix_external_set_virial_peratom(handle, "ext", vatom);
    delete[] eatom;
    delete[] vatom[0];
    delete[] vatom;
}
}

TEST(lammps_external, callback)
{
    const char *args[] = {"liblammps", "-log", "none", "-nocite"};
    char **argv        = (char **)args;
    int argc           = sizeof(args) / sizeof(char *);

    ::testing::internal::CaptureStdout();
    void *handle       = lammps_open_no_mpi(argc, argv, nullptr);
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
                                   "fix 1 all nve\n"
                                   "fix ext all external pf/callback 5 1\n"
                                   "compute eatm all pe/atom fix\n"
                                   "compute vatm all stress/atom NULL fix\n"
                                   "compute sum all reduce sum c_eatm c_vatm[*]\n"
                                   "thermo_style custom step temp pe ke etotal press c_sum[*]\n"
                                   "thermo 5\n"
                                   "fix_modify ext energy yes virial yes\n");
    lammps_fix_external_set_vector_length(handle, "ext", 6);
    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    ::testing::internal::CaptureStdout();
    lammps_set_fix_external_callback(handle, "ext", &callback, handle);
    lammps_command(handle, "run 10 post no");
    double temp  = lammps_get_thermo(handle, "temp");
    double pe    = lammps_get_thermo(handle, "pe");
    double press = lammps_get_thermo(handle, "press");
    double val   = 0.0;
    double *valp;
    for (int i = 0; i < 6; ++i) {
        valp = (double *)lammps_extract_fix(handle, "ext", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, i, 0);
        val += *valp;
        lammps_free(valp);
    }
    double *reduce =
        (double *)lammps_extract_compute(handle, "sum", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR);
    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    EXPECT_DOUBLE_EQ(temp, 1.0 / 30.0);
    EXPECT_DOUBLE_EQ(pe, 1.0 / 8.0);
    EXPECT_DOUBLE_EQ(press, 0.15416666666666667);
    EXPECT_DOUBLE_EQ(val, 15);
    EXPECT_DOUBLE_EQ(reduce[0], 2.8);
    EXPECT_DOUBLE_EQ(reduce[1], -0.7);
    EXPECT_DOUBLE_EQ(reduce[2], -0.7);
    EXPECT_DOUBLE_EQ(reduce[3], -0.7);
    EXPECT_DOUBLE_EQ(reduce[4], 1.4);
    EXPECT_DOUBLE_EQ(reduce[5], 1.4);
    EXPECT_DOUBLE_EQ(reduce[6], 1.4);

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
    void *handle       = lammps_open_no_mpi(argc, argv, nullptr);
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
                                   "fix 1 all nve\n"
                                   "fix ext all external pf/array 1\n"
                                   "thermo 5\n");

    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    ::testing::internal::CaptureStdout();
    double **force = lammps_fix_external_get_force(handle, "ext");
    int nlocal     = lammps_extract_setting(handle, "nlocal");
    for (int i = 0; i < nlocal; ++i)
        force[i][0] = force[i][1] = force[i][2] = 0.0;
    lammps_fix_external_set_energy_global(handle, "ext", 0.5);
    double v[6] = {0.5, 0.5, 0.5, 0.0, 0.0, 0.0};
    lammps_fix_external_set_virial_global(handle, "ext", v);
    lammps_command(handle, "run 5 post no");
    double temp  = lammps_get_thermo(handle, "temp");
    double pe    = lammps_get_thermo(handle, "pe");
    double press = lammps_get_thermo(handle, "press");
    output       = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    EXPECT_DOUBLE_EQ(temp, 4.0 / 525.0);
    EXPECT_DOUBLE_EQ(pe, 1.0 / 16.0);
    EXPECT_DOUBLE_EQ(press, 0.069166666666666668);

    ::testing::internal::CaptureStdout();
    nlocal = lammps_extract_setting(handle, "nlocal");
    force  = lammps_fix_external_get_force(handle, "ext");
    for (int i = 0; i < nlocal; ++i)
        force[i][0] = force[i][1] = force[i][2] = 6.0;
    lammps_fix_external_set_energy_global(handle, "ext", 1.0);
    v[0] = v[1] = v[2] = 1.0;
    v[3] = v[4] = v[5] = 0.0;
    lammps_fix_external_set_virial_global(handle, "ext", v);
    lammps_command(handle, "run 5 post no");
    temp   = lammps_get_thermo(handle, "temp");
    pe     = lammps_get_thermo(handle, "pe");
    press  = lammps_get_thermo(handle, "press");
    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
    EXPECT_DOUBLE_EQ(temp, 1.0 / 30.0);
    EXPECT_DOUBLE_EQ(pe, 1.0 / 8.0);
    EXPECT_DOUBLE_EQ(press, 0.15416666666666667);

    ::testing::internal::CaptureStdout();
    lammps_close(handle);
    output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;
}
