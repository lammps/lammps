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
#include "variable.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "exceptions.h"

#include <functional>

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
extern bool verbose;

class LAMMPSTest : public ::testing::Test {
public:
    void command(const std::string &line) { lmp->input->one(line.c_str()); }

    void BEGIN_HIDE_OUTPUT() {
        if (!verbose) ::testing::internal::CaptureStdout();
    }

    void END_HIDE_OUTPUT() {
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void BEGIN_CAPTURE_OUTPUT() {
        if (!verbose) ::testing::internal::CaptureStdout();
    }

    std::string END_CAPTURE_OUTPUT() {
        auto output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        return output;
    }

    void HIDE_OUTPUT(std::function<void()> f) {
        if (!verbose) ::testing::internal::CaptureStdout();
        try {
            f();
        } catch(LAMMPSException & e) {
            if (!verbose) std::cout << ::testing::internal::GetCapturedStdout();
            throw e;
        }
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    std::string CAPTURE_OUTPUT(std::function<void()> f) {
        ::testing::internal::CaptureStdout();
        try {
            f();
        } catch(LAMMPSException & e) {
            if (verbose) std::cout << ::testing::internal::GetCapturedStdout();
            throw e;
        }
        auto output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        return output;
    }

    double get_variable_value(const std::string & name) {
        char * str = utils::strdup(fmt::format("v_{}", name));
        double value = lmp->input->variable->compute_equal(str);
        delete [] str;
        return value;
    }

    std::string get_variable_string(const std::string & name) {
        return lmp->input->variable->retrieve(name.c_str());
    }

protected:
    const char *testbinary = "LAMMPSTest";
    LAMMPS *lmp;
    Info *info;

    void SetUp() override
    {
        const char *args[] = {testbinary, "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        HIDE_OUTPUT([&] {
            lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
            info = new Info(lmp);
        });
        InitSystem();
    }

    virtual void InitSystem() {}

    void TearDown() override
    {
        HIDE_OUTPUT([&] {
            delete info;
            delete lmp;
            info = nullptr;
            lmp = nullptr;
        });
        std::cout.flush();
    }
};

#endif
