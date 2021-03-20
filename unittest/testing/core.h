/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#ifndef TESTING_CORE__H
#define TESTING_CORE__H

#include "info.h"
#include "input.h"
#include "lammps.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace LAMMPS_NS;

using ::testing::MatchesRegex;

#define TEST_FAILURE(errmsg, ...)                                 \
    if (Info::has_exceptions()) {                                 \
        ::testing::internal::CaptureStdout();                     \
        ASSERT_ANY_THROW({__VA_ARGS__});                          \
        auto mesg = ::testing::internal::GetCapturedStdout();     \
        ASSERT_THAT(mesg, MatchesRegex(errmsg));                  \
    } else {                                                      \
        if (Info::get_mpi_vendor() != "Open MPI") {               \
            ::testing::internal::CaptureStdout();                 \
            ASSERT_DEATH({__VA_ARGS__}, "");                      \
            auto mesg = ::testing::internal::GetCapturedStdout(); \
            ASSERT_THAT(mesg, MatchesRegex(errmsg));              \
        }                                                         \
    }

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

class LAMMPSTest : public ::testing::Test {
public:
    void command(const std::string &line) { lmp->input->one(line.c_str()); }

protected:
    const char *testbinary = "LAMMPSTest";
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {testbinary, "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        InitSystem();
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    virtual void InitSystem() {}

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        lmp = nullptr;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

#endif
