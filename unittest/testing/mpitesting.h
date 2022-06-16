/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <deque>
#include <mpi.h>

using ::testing::TestCase;
using ::testing::TestEventListener;
using ::testing::TestInfo;
using ::testing::TestPartResult;
using ::testing::TestSuite;
using ::testing::UnitTest;

class MPIPrinter : public TestEventListener {
    MPI_Comm comm;
    TestEventListener *default_listener;
    int me;
    int nprocs;
    char *buffer;
    size_t buffer_size;
    std::deque<TestPartResult> results;
    bool finalize_test;

public:
    MPIPrinter(TestEventListener *default_listener) : default_listener(default_listener)
    {
        comm = MPI_COMM_WORLD;
        MPI_Comm_rank(comm, &me);
        MPI_Comm_size(comm, &nprocs);
        buffer_size   = 1024;
        buffer        = new char[buffer_size];
        finalize_test = false;
    }

    ~MPIPrinter() override
    {
        delete default_listener;
        default_listener = nullptr;

        delete[] buffer;
        buffer      = nullptr;
        buffer_size = 0;
    }

    virtual void OnTestProgramStart(const UnitTest &unit_test) override
    {
        if (me == 0) default_listener->OnTestProgramStart(unit_test);
    }

    virtual void OnTestIterationStart(const UnitTest &unit_test, int iteration) override
    {
        if (me == 0) default_listener->OnTestIterationStart(unit_test, iteration);
    }

    virtual void OnEnvironmentsSetUpStart(const UnitTest &unit_test) override
    {
        if (me == 0) default_listener->OnEnvironmentsSetUpStart(unit_test);
    }

    virtual void OnEnvironmentsSetUpEnd(const UnitTest &unit_test) override
    {
        if (me == 0) default_listener->OnEnvironmentsSetUpEnd(unit_test);
    }

    virtual void OnTestSuiteStart(const TestSuite &test_suite) override
    {
        if (me == 0) default_listener->OnTestSuiteStart(test_suite);
    }

    //  Legacy API is deprecated but still available
#ifndef GTEST_REMOVE_LEGACY_TEST_CASEAPI_
    virtual void OnTestCaseStart(const TestCase &test_case) override
    {
        if (me == 0) default_listener->OnTestSuiteStart(test_case);
    }
#endif //  GTEST_REMOVE_LEGACY_TEST_CASEAPI_

    virtual void OnTestStart(const TestInfo &test_info) override
    {
        // Called before a test starts.
        if (me == 0) default_listener->OnTestStart(test_info);
        results.clear();
        finalize_test = false;
    }

    virtual void OnTestPartResult(const TestPartResult &test_part_result) override
    {
        // Called after a failed assertion or a SUCCESS().
        // test_part_result()

        if (me == 0 && finalize_test) {
            default_listener->OnTestPartResult(test_part_result);
        } else {
            std::stringstream proc_message;
            std::istringstream msg(test_part_result.message());
            std::string line;

            while (std::getline(msg, line)) {
                proc_message << "[Rank " << me << "] " << line << std::endl;
            }

            results.push_back(TestPartResult(test_part_result.type(), test_part_result.file_name(),
                                             test_part_result.line_number(),
                                             proc_message.str().c_str()));
        }
    }

    virtual void OnTestEnd(const TestInfo &test_info) override
    {
        // Called after a test ends.
        MPI_Barrier(comm);

        // other procs send their test part results
        if (me != 0) {
            int nresults = results.size();
            MPI_Send(&nresults, 1, MPI_INT, 0, 0, comm);

            for (auto &test_part_result : results) {

                int type = test_part_result.type();
                MPI_Send(&type, 1, MPI_INT, 0, 0, comm);

                const char *str = test_part_result.file_name();
                int length      = 0;
                if (str) length = strlen(str) + 1;
                MPI_Send(&length, 1, MPI_INT, 0, 0, comm);
                if (str) MPI_Send(str, length, MPI_CHAR, 0, 0, comm);

                int lineno = test_part_result.line_number();
                MPI_Send(&lineno, 1, MPI_INT, 0, 0, comm);

                str    = test_part_result.message();
                length = 0;
                if (str) length = strlen(str) + 1;
                MPI_Send(&length, 1, MPI_INT, 0, 0, comm);
                if (str) MPI_Send(str, length, MPI_CHAR, 0, 0, comm);
            }
        }

        if (me == 0) {
            // collect results from other procs
            for (int p = 1; p < nprocs; p++) {
                int nresults = 0;
                MPI_Recv(&nresults, 1, MPI_INT, p, 0, comm, MPI_STATUS_IGNORE);

                for (int r = 0; r < nresults; r++) {

                    int type;
                    MPI_Recv(&type, 1, MPI_INT, p, 0, comm, MPI_STATUS_IGNORE);

                    int length = 0;
                    MPI_Recv(&length, 1, MPI_INT, p, 0, comm, MPI_STATUS_IGNORE);
                    std::string file_name;

                    if (length > 0) {
                        if (length > buffer_size) {
                            delete[] buffer;
                            buffer      = new char[length];
                            buffer_size = length;
                        }
                        MPI_Recv(buffer, length, MPI_CHAR, p, 0, comm, MPI_STATUS_IGNORE);
                        file_name = buffer;
                    }

                    int lineno;
                    MPI_Recv(&lineno, 1, MPI_INT, p, 0, comm, MPI_STATUS_IGNORE);

                    MPI_Recv(&length, 1, MPI_INT, p, 0, comm, MPI_STATUS_IGNORE);
                    std::string message;

                    if (length > 0) {
                        if (length > buffer_size) {
                            delete[] buffer;
                            buffer      = new char[length];
                            buffer_size = length;
                        }
                        MPI_Recv(buffer, length, MPI_CHAR, p, 0, comm, MPI_STATUS_IGNORE);
                        message = std::string(buffer);
                    }

                    results.push_back(TestPartResult((TestPartResult::Type)type, file_name.c_str(),
                                                     lineno, message.c_str()));
                }
            }

            // ensure failures are reported
            finalize_test = true;

            // add all failures
            while (!results.empty()) {
                auto result = results.front();
                if (result.failed()) {
                    ADD_FAILURE_AT(result.file_name(), result.line_number()) << result.message();
                } else {
                    default_listener->OnTestPartResult(result);
                }
                results.pop_front();
            }

            default_listener->OnTestEnd(test_info);
        }
    }

    virtual void OnTestSuiteEnd(const TestSuite &test_suite) override
    {
        if (me == 0) default_listener->OnTestSuiteEnd(test_suite);
    }

#ifndef GTEST_REMOVE_LEGACY_TEST_CASEAPI_
    virtual void OnTestCaseEnd(const TestCase &test_case) override
    {
        if (me == 0) default_listener->OnTestCaseEnd(test_case);
    }
#endif //  GTEST_REMOVE_LEGACY_TEST_CASEAPI_

    virtual void OnEnvironmentsTearDownStart(const UnitTest &unit_test) override
    {
        if (me == 0) default_listener->OnEnvironmentsTearDownStart(unit_test);
    }

    virtual void OnEnvironmentsTearDownEnd(const UnitTest &unit_test) override
    {
        if (me == 0) default_listener->OnEnvironmentsTearDownEnd(unit_test);
    }

    virtual void OnTestIterationEnd(const UnitTest &unit_test, int iteration) override
    {
        if (me == 0) default_listener->OnTestIterationEnd(unit_test, iteration);
    }

    virtual void OnTestProgramEnd(const UnitTest &unit_test) override
    {
        if (me == 0) default_listener->OnTestProgramEnd(unit_test);
    }
};
