// unit tests for checking LAMMPS configuration settings  through the library interface

#include "lammps.h"
#include "library.h"
#include "timer.h"
#include <string>
#include <deque>

#include <mpi.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::ExitedWithCode;
using ::testing::HasSubstr;
using ::testing::StartsWith;
using ::testing::StrEq;

using ::testing::TestEventListener;
using ::testing::TestCase;
using ::testing::TestSuite;
using ::testing::UnitTest;
using ::testing::TestPartResult;
using ::testing::TestInfo;

class MPIPrinter : public TestEventListener {
    MPI_Comm comm;
    TestEventListener * default_listener;
    int me;
    int nprocs;
    char * buffer;
    size_t buffer_size;
    std::deque<TestPartResult> results;
    bool finalize_test;
public:
    MPIPrinter(TestEventListener * default_listener) : default_listener(default_listener) {
        comm = MPI_COMM_WORLD;
        MPI_Comm_rank(comm, &me);
        MPI_Comm_size(comm, &nprocs);
        buffer_size = 1024;
        buffer = new char[buffer_size];
        finalize_test = false;
    }

    ~MPIPrinter() override {
        delete default_listener;
        default_listener = nullptr;

        delete [] buffer;
        buffer = nullptr;
        buffer_size = 0;
    }

    virtual void OnTestProgramStart(const UnitTest& unit_test) override {
        if(me == 0) default_listener->OnTestProgramStart(unit_test);
    }

    virtual void OnTestIterationStart(const UnitTest& unit_test, int iteration) override {
        if(me == 0) default_listener->OnTestIterationStart(unit_test, iteration);
    }

    virtual void OnEnvironmentsSetUpStart(const UnitTest& unit_test) override {
        if(me == 0) default_listener->OnEnvironmentsSetUpStart(unit_test);
    }

    virtual void OnEnvironmentsSetUpEnd(const UnitTest& unit_test) override {
        if(me == 0) default_listener->OnEnvironmentsSetUpEnd(unit_test);
    }

    virtual void OnTestSuiteStart(const TestSuite& test_suite) override {
        if(me == 0) default_listener->OnTestSuiteStart(test_suite);
    }

    //  Legacy API is deprecated but still available
#ifndef GTEST_REMOVE_LEGACY_TEST_CASEAPI_
    virtual void OnTestCaseStart(const TestCase& test_case) override {
        if(me == 0) default_listener->OnTestSuiteStart(test_case);
    }
#endif  //  GTEST_REMOVE_LEGACY_TEST_CASEAPI_


    virtual void OnTestStart(const TestInfo& test_info) override {
        // Called before a test starts.
        if(me == 0) default_listener->OnTestStart(test_info);
        results.clear();
        finalize_test = false;
    }


    virtual void OnTestPartResult(const TestPartResult& test_part_result) override {
        // Called after a failed assertion or a SUCCESS().
        // test_part_result()

        if (me == 0 && finalize_test) {
            default_listener->OnTestPartResult(test_part_result);
        } else {
            std::stringstream proc_message;
            std::istringstream msg(test_part_result.message());
            std::string line;

            while(std::getline(msg, line)) {
                proc_message << "[Rank " << me << "] " << line << std::endl;
            }

            results.push_back(TestPartResult(test_part_result.type(), test_part_result.file_name(), test_part_result.line_number(), proc_message.str().c_str()));
        }
    }

    virtual void OnTestEnd(const TestInfo& test_info) override {
        // Called after a test ends.
        MPI_Barrier(comm);

        // other procs send their test part results
        if(me != 0) {
            int nresults = results.size();
            MPI_Send(&nresults, 1, MPI_INT, 0, 0, comm);

            for(auto& test_part_result : results) {

                int type = test_part_result.type();
                MPI_Send(&type, 1, MPI_INT, 0, 0, comm);

                const char * str = test_part_result.file_name();
                int length = 0;
                if(str) length = strlen(str)+1;
                MPI_Send(&length, 1, MPI_INT, 0, 0, comm);
                if(str) MPI_Send(str, length, MPI_CHAR, 0, 0, comm);

                int lineno = test_part_result.line_number();
                MPI_Send(&lineno, 1, MPI_INT, 0, 0, comm);

                str = test_part_result.message();
                length = 0;
                if(str) length = strlen(str)+1;
                MPI_Send(&length, 1, MPI_INT, 0, 0, comm);
                if(str) MPI_Send(str, length, MPI_CHAR, 0, 0, comm);
            }
        }

        if(me == 0) {
            // collect results from other procs
            for(int p = 1; p < nprocs; p++) {
                int nresults = 0;
                MPI_Recv(&nresults, 1, MPI_INT, p, 0, comm, MPI_STATUS_IGNORE);

                for(int r = 0; r < nresults; r++) {

                    int type;
                    MPI_Recv(&type, 1, MPI_INT, p, 0, comm, MPI_STATUS_IGNORE);

                    int length = 0;
                    MPI_Recv(&length, 1, MPI_INT, p, 0, comm, MPI_STATUS_IGNORE);
                    std::string file_name;

                    if (length > 0) {
                        if (length > buffer_size) {
                            delete [] buffer;
                            buffer = new char[length];
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
                            delete [] buffer;
                            buffer = new char[length];
                            buffer_size = length;
                        }
                        MPI_Recv(buffer, length, MPI_CHAR, p, 0, comm, MPI_STATUS_IGNORE);
                        message = std::string(buffer);
                    }

                    results.push_back(TestPartResult((TestPartResult::Type)type, file_name.c_str(), lineno, message.c_str()));
                }
            }

            // ensure failures are reported
            finalize_test = true;

            // add all failures
            while(!results.empty()) {
                auto result = results.front();
                if(result.failed()) {
                    ADD_FAILURE_AT(result.file_name(), result.line_number()) << result.message();
                } else {
                    default_listener->OnTestPartResult(result);
                }
                results.pop_front();
            }

            default_listener->OnTestEnd(test_info);
        }
    }

    virtual void OnTestSuiteEnd(const TestSuite& test_suite) override {
        if(me == 0) default_listener->OnTestSuiteEnd(test_suite);
    }

#ifndef GTEST_REMOVE_LEGACY_TEST_CASEAPI_
    virtual void OnTestCaseEnd(const TestCase& test_case) override {
        if(me == 0) default_listener->OnTestCaseEnd(test_case);
    }
#endif  //  GTEST_REMOVE_LEGACY_TEST_CASEAPI_

    virtual void OnEnvironmentsTearDownStart(const UnitTest& unit_test) override {
        if(me == 0) default_listener->OnEnvironmentsTearDownStart(unit_test);
    }

    virtual void OnEnvironmentsTearDownEnd(const UnitTest& unit_test) override {
        if(me == 0) default_listener->OnEnvironmentsTearDownEnd(unit_test);
    }

    virtual void OnTestIterationEnd(const UnitTest& unit_test, int iteration) override {
        if(me == 0) default_listener->OnTestIterationEnd(unit_test, iteration);
    }

    virtual void OnTestProgramEnd(const UnitTest& unit_test) override {
        if(me == 0) default_listener->OnTestProgramEnd(unit_test);
    }
};


TEST(MPI, global_box)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    EXPECT_EQ(nprocs, 4);
    EXPECT_GT(me, -1);
    EXPECT_LT(me, 5);

    double boxlo[3];
    double boxhi[3];
    double xy  = 0.0;
    double yz  = 0.0;
    double xz  = 0.0;
    int pflags[3];
    int boxflag;

    ::testing::internal::CaptureStdout();
    const char *args[] = {"LAMMPS_test", "-log",      "none",
                          "-echo",       "screen",    "-nocite"};
    char **argv = (char **)args;
    int argc    = sizeof(args) / sizeof(char *);
    void * lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
    lammps_command(lmp, "units           lj");
    lammps_command(lmp, "atom_style      atomic");
    lammps_command(lmp, "region          box block 0 2 0 2 0 2");
    lammps_command(lmp, "create_box      1 box");

    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, pflags, &boxflag);
    ::testing::internal::GetCapturedStdout();

    EXPECT_EQ(boxlo[0], 0.0);
    EXPECT_EQ(boxlo[1], 0.0);
    EXPECT_EQ(boxlo[2], 0.0);

    EXPECT_EQ(boxhi[0], 2.0);
    EXPECT_EQ(boxhi[1], 2.0);
    EXPECT_EQ(boxhi[2], 2.0);

    ::testing::internal::CaptureStdout();
    lammps_close(lmp);
    ::testing::internal::GetCapturedStdout();
};

TEST(MPI, sub_box)
{
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    EXPECT_EQ(nprocs, 4);
    EXPECT_GT(me, -1);
    EXPECT_LT(me, 5);

    double boxlo[3];
    double boxhi[3];
    double xy  = 0.0;
    double yz  = 0.0;
    double xz  = 0.0;
    int pflags[3];
    int boxflag;

    ::testing::internal::CaptureStdout();
    const char *args[] = {"LAMMPS_test", "-log",      "none",
                          "-echo",       "screen",    "-nocite"};
    char **argv = (char **)args;
    int argc    = sizeof(args) / sizeof(char *);
    void * lmp = lammps_open(argc, argv, MPI_COMM_WORLD, nullptr);
    lammps_command(lmp, "units           lj");
    lammps_command(lmp, "atom_style      atomic");
    lammps_command(lmp, "region          box block 0 2 0 2 0 2");
    lammps_command(lmp, "create_box      1 box");

    lammps_extract_box(lmp, boxlo, boxhi, &xy, &yz, &xz, pflags, &boxflag);
    ::testing::internal::GetCapturedStdout();

    EXPECT_EQ(boxlo[0], 0.0);
    EXPECT_EQ(boxlo[1], 0.0);
    EXPECT_EQ(boxlo[2], 0.0);

    EXPECT_EQ(boxhi[0], 2.0);
    EXPECT_EQ(boxhi[1], 2.0);
    EXPECT_EQ(boxhi[2], 2.0);

    double * sublo = (double*)lammps_extract_global(lmp, "sublo");
    double * subhi = (double*)lammps_extract_global(lmp, "subhi");

    ASSERT_NE(sublo, nullptr);
    ASSERT_NE(subhi, nullptr);

    EXPECT_GE(sublo[0], boxlo[0]);
    EXPECT_GE(sublo[1], boxlo[1]);
    EXPECT_GE(sublo[2], boxlo[2]);
    EXPECT_LE(subhi[0], boxhi[0]);
    EXPECT_LE(subhi[1], boxhi[1]);
    EXPECT_LE(subhi[2], boxhi[2]);

    ::testing::internal::CaptureStdout();
    lammps_close(lmp);
    ::testing::internal::GetCapturedStdout();
};


bool verbose = false;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (argc < 1) {
        return 1;
    }

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    int iarg = 1;
    while (iarg < argc) {
        if (strcmp(argv[iarg], "-v") == 0) {
            verbose = true;
            ++iarg;
        } else {
            std::cerr << "unknown option: " << argv[iarg] << "\n\n";
            MPI_Finalize();
            return 1;
        }
    }

    auto & listeners = UnitTest::GetInstance()->listeners();

    // Remove default listener
    auto default_listener = listeners.Release(listeners.default_result_printer());

    // Adds a listener to the end.  googletest takes the ownership.
    listeners.Append(new MPIPrinter(default_listener));

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
