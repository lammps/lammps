// unit tests for extracting compute data from a LAMMPS instance through the
// Fortran wrapper

#include "lammps.h"
#include "library.h"
#include "platform.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <string>
#include <thread>

#include "gtest/gtest.h"

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

// prototypes for Fortran reverse wrapper functions
extern "C" {
void *f_lammps_with_c_args(int, char **);
void f_lammps_close();
void f_lammps_setup_extract_variable();
int f_lammps_extract_variable_index_1();
int f_lammps_extract_variable_index_2();
int f_lammps_extract_variable_loop();
char *f_lammps_extract_variable_loop_pad();
char *f_lammps_extract_variable_world();
char *f_lammps_extract_variable_universe();
int f_lammps_extract_variable_uloop();
char *f_lammps_extract_variable_string();
char *f_lammps_extract_variable_format();
char *f_lammps_extract_variable_format_pad();
char *f_lammps_extract_variable_getenv();
char *f_lammps_extract_variable_file();
double f_lammps_extract_variable_atomfile(int);
double f_lammps_extract_variable_python();
double f_lammps_extract_variable_timer();
double f_lammps_extract_variable_internal();
double f_lammps_extract_variable_equal();
double f_lammps_extract_variable_atom(int);
double f_lammps_extract_variable_vector(int);
void f_lammps_set_variable_string();
char *c_path_join(const char *, const char *);
}

char *c_path_join(const char *a, const char *b)
{
    std::string A = a;
    std::string B = b;
    std::string C = LAMMPS_NS::platform::path_join(A, B);
    size_t length = C.length() + 1;
    char *retval  = (char *)malloc(length * sizeof(char));
    C.copy(retval, length);
    retval[length - 1] = '\0';
    return retval;
}

constexpr char input_dir[] = STRINGIFY(TEST_INPUT_FOLDER);
class LAMMPS_extract_variable : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_extract_variable()           = default;
    ~LAMMPS_extract_variable() override = default;

    void SetUp() override
    {
        // clang-format off
        const char *args[] =
            { "LAMMPS_Fortran_test", "-l", "none", "-echo", "screen", "-nocite",
              "-var", "input_dir", input_dir, "-var", "zpos", "1.5", "-var", "x", "2" };
        // clang-format on
        char **argv = (char **)args;
        int argc    = sizeof(args) / sizeof(const char *);
        ::testing::internal::CaptureStdout();
        lmp = (LAMMPS_NS::LAMMPS *)f_lammps_with_c_args(argc, argv);

        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0, 8).c_str(), "LAMMPS (");
    }

    void TearDown() override
    {
        ::testing::internal::CaptureStdout();
        f_lammps_close();
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0, 16).c_str(), "Total wall time:");
        lmp = nullptr;
    }
};

TEST_F(LAMMPS_extract_variable, index)
{
    f_lammps_setup_extract_variable();
    EXPECT_EQ(f_lammps_extract_variable_index_1(), 1);
    EXPECT_EQ(f_lammps_extract_variable_index_2(), 0);
    lammps_command(lmp, "next idx");
    EXPECT_EQ(f_lammps_extract_variable_index_1(), 0);
    EXPECT_EQ(f_lammps_extract_variable_index_2(), 1);
};

TEST_F(LAMMPS_extract_variable, loop)
{
    f_lammps_setup_extract_variable();
    int i;
    for (i = 1; i <= 10; i++) {
        EXPECT_EQ(f_lammps_extract_variable_loop(), i);
        lammps_command(lmp, "next lp");
    }
};

TEST_F(LAMMPS_extract_variable, loop_pad)
{
    f_lammps_setup_extract_variable();
    int i;
    char str[10];
    char *fstr;
    for (i = 1; i <= 10; i++) {
        std::sprintf(str, "%02d", i);
        fstr = f_lammps_extract_variable_loop_pad();
        EXPECT_STREQ(fstr, str);
        std::free(fstr);
        lammps_command(lmp, "next lp_pad");
    }
};

TEST_F(LAMMPS_extract_variable, world)
{
    f_lammps_setup_extract_variable();
    char *fstr = f_lammps_extract_variable_world();
    EXPECT_STREQ(fstr, "group1");
    std::free(fstr);
};

TEST_F(LAMMPS_extract_variable, universe)
{
    f_lammps_setup_extract_variable();
    char *fstr = f_lammps_extract_variable_universe();
    EXPECT_STREQ(fstr, "universe1");
    std::free(fstr);
};

TEST_F(LAMMPS_extract_variable, uloop)
{
    f_lammps_setup_extract_variable();
    EXPECT_EQ(f_lammps_extract_variable_uloop(), 1);
};

TEST_F(LAMMPS_extract_variable, string)
{
    f_lammps_setup_extract_variable();
    char *fstr = f_lammps_extract_variable_string();
    EXPECT_STREQ(fstr, "this is a string");
    std::free(fstr);
    f_lammps_set_variable_string();
    fstr = f_lammps_extract_variable_string();
    EXPECT_STREQ(fstr, "this is the new string");
    std::free(fstr);
};

TEST_F(LAMMPS_extract_variable, format)
{
    f_lammps_setup_extract_variable();
    int i;
    char str[16];
    char *fstr;
    for (i = 1; i <= 10; i++) {
        std::sprintf(str, "%.6G", std::exp(i));
        fstr = f_lammps_extract_variable_format();
        EXPECT_STREQ(fstr, str);
        std::free(fstr);
        lammps_command(lmp, "next lp");
    }
};

TEST_F(LAMMPS_extract_variable, format_pad)
{
    f_lammps_setup_extract_variable();
    int i;
    char str[16];
    char *fstr;
    for (i = 1; i <= 10; i++) {
        std::sprintf(str, "%08.6G", std::exp(i));
        fstr = f_lammps_extract_variable_format_pad();
        EXPECT_STREQ(fstr, str);
        std::free(fstr);
        lammps_command(lmp, "next lp");
    }
};

TEST_F(LAMMPS_extract_variable, getenv)
{
    LAMMPS_NS::platform::putenv("FORTRAN_USER=myuser");
    f_lammps_setup_extract_variable();
    char *env  = std::getenv("FORTRAN_USER");
    char *fenv = f_lammps_extract_variable_getenv();
    EXPECT_STREQ(fenv, env);
    std::free(fenv);
};

TEST_F(LAMMPS_extract_variable, file)
{
    f_lammps_setup_extract_variable();
    int i;
    const char *str[9] = {"hello",      "god_dag", "hola",  "bonjour",  "guten_Tag",
                          "konnichiwa", "shalom",  "salve", "goedendag"};
    char *fstr;
    for (i = 0; i < 9; i++) {
        fstr = f_lammps_extract_variable_file();
        EXPECT_STREQ(fstr, str[i]);
        std::free(fstr);
        lammps_command(lmp, "next greeting");
    }
};

TEST_F(LAMMPS_extract_variable, atomfile)
{
    f_lammps_setup_extract_variable();
    EXPECT_DOUBLE_EQ(f_lammps_extract_variable_atomfile(1), 5.2);
    EXPECT_DOUBLE_EQ(f_lammps_extract_variable_atomfile(2), 1.6);
    EXPECT_DOUBLE_EQ(f_lammps_extract_variable_atomfile(3), -1.4);
    lammps_command(lmp, "next atfile");
    EXPECT_DOUBLE_EQ(f_lammps_extract_variable_atomfile(1), -1.1);
    EXPECT_DOUBLE_EQ(f_lammps_extract_variable_atomfile(2), 0.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_variable_atomfile(3), 2.5);
};

TEST_F(LAMMPS_extract_variable, python)
{
    if (lammps_config_has_package("PYTHON")) {
        f_lammps_setup_extract_variable();
        for (int i = 1; i <= 10; i++) {
            EXPECT_DOUBLE_EQ(f_lammps_extract_variable_python(), static_cast<double>(i * i));
            lammps_command(lmp, "next lp");
        }
    }
};

TEST_F(LAMMPS_extract_variable, timer)
{
    f_lammps_setup_extract_variable();
    double initial_t, final_t;
    initial_t = f_lammps_extract_variable_timer();
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    lammps_command(lmp, "variable time timer"); // update the time
    final_t = f_lammps_extract_variable_timer();
    EXPECT_GT(final_t, initial_t + 0.1);
};

TEST_F(LAMMPS_extract_variable, internal)
{
    f_lammps_setup_extract_variable();
    EXPECT_DOUBLE_EQ(f_lammps_extract_variable_internal(), 4.0);
};

TEST_F(LAMMPS_extract_variable, equal)
{
    f_lammps_setup_extract_variable();
    int i;
    for (i = 1; i <= 10; i++) {
        EXPECT_DOUBLE_EQ(f_lammps_extract_variable_equal(), std::exp(static_cast<double>(i)));
        lammps_command(lmp, "next lp");
    }
};

TEST_F(LAMMPS_extract_variable, atom)
{
    f_lammps_setup_extract_variable();
    EXPECT_DOUBLE_EQ(f_lammps_extract_variable_atom(1), 1.5);
    EXPECT_DOUBLE_EQ(f_lammps_extract_variable_atom(2), 0.1);
    EXPECT_DOUBLE_EQ(f_lammps_extract_variable_atom(3), 0.5);
};

TEST_F(LAMMPS_extract_variable, vector)
{
    f_lammps_setup_extract_variable();
    EXPECT_DOUBLE_EQ(f_lammps_extract_variable_vector(1), (1 + 0.2 + 0.5) / 3.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_variable_vector(2), (1 + 0.1 + 0.5) / 3.0);
    EXPECT_DOUBLE_EQ(f_lammps_extract_variable_vector(3), (1.5 + 0.1 + 0.5) / 3.0);
};
