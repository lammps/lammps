/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#ifndef TESTING_CORE__H
#define TESTING_CORE__H

#include "exceptions.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "platform.h"
#include "utils.h"
#include "variable.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdlib>
#include <functional>
#include <string>
#include <vector>

using LAMMPS_NS::Info;
using LAMMPS_NS::LAMMPS;
using LAMMPS_NS::LAMMPSException;

using ::testing::ContainsRegex;

#if defined(LAMMPS_SKIP_DEATH_TESTS)
#define TEST_FAILURE(errmsg, ...)                             \
    {                                                         \
        ;                                                     \
    }
#else
#define TEST_FAILURE(errmsg, ...)                             \
    {                                                         \
        ::testing::internal::CaptureStdout();                 \
        ASSERT_ANY_THROW({__VA_ARGS__});                      \
        auto mesg = ::testing::internal::GetCapturedStdout(); \
        ASSERT_THAT(mesg, ContainsRegex(errmsg));             \
    }
#endif

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
extern bool verbose;

class LAMMPSTest : public ::testing::Test {
public:
    void command(const std::string &line) { lmp->input->one(line); }

    // GoogleTest supports only one active stdout capture at a time.  When an
    // error is thrown between BEGIN_HIDE_OUTPUT() and END_HIDE_OUTPUT(), the
    // capture stays active and the next one aborts the entire test program
    // with "Only one stdout capturer can exist at a time", which hides the
    // original error message.  Track the state so a capture left behind can
    // be dropped when the test ends.

    void RESET_OUTPUT()
    {
        if (capturing) {
            ::testing::internal::GetCapturedStdout();
            capturing = false;
        }
    }

    void BEGIN_HIDE_OUTPUT()
    {
        if (!verbose) {
            ::testing::internal::CaptureStdout();
            capturing = true;
        }
    }

    void END_HIDE_OUTPUT()
    {
        if (!verbose) RESET_OUTPUT();
    }

    void BEGIN_CAPTURE_OUTPUT()
    {
        ::testing::internal::CaptureStdout();
        capturing = true;
    }

    std::string END_CAPTURE_OUTPUT()
    {
        auto output = ::testing::internal::GetCapturedStdout();
        capturing   = false;
        if (verbose) std::cout << output;
        return output;
    }

    void HIDE_OUTPUT(std::function<void()> f)
    {
        BEGIN_HIDE_OUTPUT();
        try {
            f();
        } catch (LAMMPSException &e) {
            if (!verbose) {
                capturing = false;
                std::cout << ::testing::internal::GetCapturedStdout();
            }
            throw e;
        }
        END_HIDE_OUTPUT();
    }

    std::string CAPTURE_OUTPUT(std::function<void()> f)
    {
        BEGIN_CAPTURE_OUTPUT();
        try {
            f();
        } catch (LAMMPSException &e) {
            capturing = false;
            auto mesg = ::testing::internal::GetCapturedStdout();
            if (verbose) std::cout << mesg;
            throw e;
        }
        return END_CAPTURE_OUTPUT();
    }

    double get_variable_value(const std::string &name)
    {
        char *str    = LAMMPS_NS::utils::strdup(fmt::format("v_{}", name));
        double value = lmp->input->variable->compute_equal(str);
        delete[] str;
        return value;
    }

    std::string get_variable_string(const std::string &name)
    {
        return lmp->input->variable->retrieve(name.c_str());
    }

protected:
    std::string testbinary = "LAMMPSTest";
    bool capturing         = false;
    LAMMPS::argv args      = {"-log", "none", "-echo", "screen", "-nocite"};
    LAMMPS *lmp;
    Info *info;

    void SetUp() override
    {
        LAMMPS::argv full_args = {testbinary};
        full_args.insert(full_args.end(), args.begin(), args.end());

        // append the accelerator command line flags given in the environment.
        // this lets CTest run the very same fixtures a second time with, e.g.,
        // LAMMPS_ACCELERATOR_ARGS="-k on t 1 -sf kk" so that the KOKKOS styles
        // and commands are exercised without duplicating any test bodies.
        const char *accel_args = std::getenv("LAMMPS_ACCELERATOR_ARGS");
        if (accel_args && (accel_args[0] != '\0')) {
            auto accel = LAMMPS_NS::utils::split_words(accel_args);
            full_args.insert(full_args.end(), accel.begin(), accel.end());
        }

        HIDE_OUTPUT([&] {
            lmp  = new LAMMPS(full_args, MPI_COMM_WORLD);
            info = new Info(lmp);
        });
        InitSystem();
    }

    virtual void InitSystem() {}

    void TearDown() override
    {
        // an error thrown inside a captured block may have left a capture
        // active.  drop it, so the deletion below does not abort the program
        RESET_OUTPUT();
        HIDE_OUTPUT([&] {
            delete info;
            delete lmp;
            info = nullptr;
            lmp  = nullptr;
        });
        std::cout.flush();
    }
};

#endif
