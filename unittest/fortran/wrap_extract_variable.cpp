// unit tests for extracting compute data from a LAMMPS instance through the
// Fortran wrapper
#include <cstdio>

#include "lammps.h"
#include "library.h"
#include <mpi.h>
#include <string>
#include <cstdlib>
#include <cstdint>

#include "gtest/gtest.h"

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

// prototypes for Fortran reverse wrapper functions
extern "C" {
void *f_lammps_with_c_args(int,char**);
void f_lammps_close();
void f_lammps_setup_extract_variable();
int f_lammps_extract_variable_index_1();
int f_lammps_extract_variable_index_2();
int f_lammps_extract_variable_loop();
double f_lammps_extract_variable_loop_pad();
void f_lammps_setup_extract_variable_world();
void f_lammps_setup_extract_variable_universe();
int f_lammps_setup_extract_variable_uloop();
void f_lammps_setup_extract_variable_string();
void f_lammps_setup_extract_variable_format();
void f_lammps_setup_extract_variable_getenv();
void f_lammps_setup_extract_variable_file();
void f_lammps_setup_extract_variable_atomfile();
double f_lammps_setup_extract_variable_python();
double f_lammps_setup_extract_variable_timer();
double f_lammps_setup_extract_variable_internal();
double f_lammps_extract_variable_equal_natoms();
double f_lammps_extract_variable_equal_dt();
double f_lammps_extract_variable_vector(int);
double f_lammps_extract_variable_atom(int);
}

class LAMMPS_extract_variable : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_extract_variable()           = default;
    ~LAMMPS_extract_variable() override = default;

    void SetUp() override
    {
        const char *args[] = {"LAMMPS_Fortran_test", "-l", "none",
                              "-echo", "screen", "-nocite", "-var",
                              "input_dir", STRINGIFY(TEST_INPUT_FOLDER),
                              "-var", "zpos", "1.5", "-var", "x", "2"};
        char** argv = (char**) args;
        int argc = sizeof(args) / sizeof(const char*);
        ::testing::internal::CaptureStdout();
std::fprintf(stderr,"THIS IS A TEST\n");
        lmp = (LAMMPS_NS::LAMMPS*)f_lammps_with_c_args(argc, argv);
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
   EXPECT_EQ(f_lammps_extract_variable_loop(), 1);
};
