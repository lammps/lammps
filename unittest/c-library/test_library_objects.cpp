// unit tests for checking and changing simulation properties through the library interface

#include "library.h"

#include "atom.h"
#include "input.h"
#include "lammps.h"
#include "lmptype.h"
#include "platform.h"
#include "variable.h"

#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "test_main.h"

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

using ::LAMMPS_NS::Atom;
using ::LAMMPS_NS::bigint;
using ::LAMMPS_NS::Input;
using ::LAMMPS_NS::tagint;
using ::LAMMPS_NS::Variable;
using ::testing::HasSubstr;
using ::testing::StartsWith;
using ::testing::StrEq;

class LibraryObjects : public ::testing::Test {
protected:
    void *lmp;
    Variable *variable;
    std::string INPUT_DIR = STRINGIFY(TEST_INPUT_FOLDER);

    LibraryObjects()           = default;
    ~LibraryObjects() override = default;

    void SetUp() override
    {
        const char *args[] = {"LAMMPS_test", "-log",      "none",
                              "-echo",       "screen",    "-nocite",
                              "-var",        "input_dir", STRINGIFY(TEST_INPUT_FOLDER),
                              nullptr};

        char **argv = (char **)args;
        int argc    = (sizeof(args) / sizeof(char *)) - 1;

        ::testing::internal::CaptureStdout();
        lmp                = lammps_open_no_mpi(argc, argv, nullptr);
        variable           = ((LAMMPS_NS::LAMMPS *)lmp)->input->variable;
        std::string output = ::testing::internal::GetCapturedStdout();
        if (verbose) std::cout << output;
        EXPECT_THAT(output, StartsWith("LAMMPS ("));
    }

    void TearDown() override
    {
        ::testing::internal::CaptureStdout();
        lammps_close(lmp);
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, HasSubstr("Total wall time:"));
        if (verbose) std::cout << output;
        lmp = nullptr;
    }
};

TEST_F(LibraryObjects, variables)
{
    FILE *fp = fopen("test_variable.file", "w");
    fputs("# test file for file style variable\n\n\none\n  two  \n\n"
          "three  # with comment\nfour   ! with non-comment\n"
          "# comments only\n    five\n#END\n",
          fp);
    fclose(fp);

    ::testing::internal::CaptureStdout();
    lammps_command(lmp, "region box block 0 2 0 2 0 2");
    lammps_command(lmp, "create_box 1 box");
    lammps_command(lmp, "mass 1 3.0");
    lammps_command(lmp, "create_atoms 1 single 1.0 1.0 1.5");
    lammps_command(lmp, "create_atoms 1 single 0.2 0.1 0.1");
    lammps_command(lmp, "shell putenv TEST_VARIABLE=simpletest2");
    lammps_command(lmp, "shell putenv TEST_VARIABLE2=simpletest OTHER_VARIABLE=2");
    lammps_command(lmp, "variable one    index     1 2 3 4");
    lammps_command(lmp, "variable two    equal     1");
    lammps_command(lmp, "variable two    equal     2");
    lammps_command(lmp, "variable three  string    four");
    lammps_command(lmp, "variable three  string    three");
    lammps_command(lmp, "variable four1  loop      4");
    lammps_command(lmp, "variable four2  loop      2 4");
    lammps_command(lmp, "variable five1  loop      100 pad");
    lammps_command(lmp, "variable five2  loop      10 200 pad");
    lammps_command(lmp, "variable six    world     one");
    lammps_command(lmp, "variable seven  format    two \"%5.2f\"");
    lammps_command(lmp, "variable eight  getenv    TEST_VARIABLE2");
    lammps_command(lmp, "variable eight  getenv    XXX");
    lammps_command(lmp, "variable nine   file      test_variable.file");
    lammps_command(lmp, "variable ten    internal  1.0");
    lammps_command(lmp, "variable ten    internal  10.0");
    lammps_command(lmp, "variable ten1   universe  1 2 3 4");
    lammps_command(lmp, "variable ten2   uloop     4");
    lammps_command(lmp, "variable ten3   uloop     4 pad");
    lammps_command(lmp, "variable ten4   vector    [0,1,2,3,5,7,11]");
    lammps_command(lmp, "variable ten5   vector    [0.5,1.25]");
    lammps_command(lmp, "variable dummy  index     0");
    lammps_command(lmp, "variable file   equal     is_file(MYFILE)");
    lammps_command(lmp, "variable iswin  equal     is_os(^Windows)");
    lammps_command(lmp, "variable islin  equal     is_os(^Linux)");
    std::string output = ::testing::internal::GetCapturedStdout();
    if (verbose) std::cout << output;

    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "unknown"), -1);
    void *ptr = lammps_extract_variable(lmp, "unknown", NULL);
    EXPECT_EQ(ptr, nullptr);
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "one"), LMP_VAR_STRING);
    ptr = lammps_extract_variable(lmp, "one", NULL);
    EXPECT_NE(ptr, nullptr);
    EXPECT_THAT((char *)ptr, StrEq("1"));
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "two"), LMP_VAR_EQUAL);
    ptr = lammps_extract_variable(lmp, "two", NULL);
    EXPECT_NE(ptr, nullptr);
    EXPECT_THAT(*(double *)ptr, 2.0);
    lammps_free(ptr);
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "three"), LMP_VAR_STRING);
    ptr = lammps_extract_variable(lmp, "three", NULL);
    EXPECT_THAT((char *)ptr, StrEq("three"));
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "four1"), LMP_VAR_STRING);
    ptr = lammps_extract_variable(lmp, "four1", NULL);
    EXPECT_THAT((char *)ptr, StrEq("1"));
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "four2"), LMP_VAR_STRING);
    ptr = lammps_extract_variable(lmp, "four2", NULL);
    EXPECT_THAT((char *)ptr, StrEq("2"));
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "five1"), LMP_VAR_STRING);
    ptr = lammps_extract_variable(lmp, "five1", NULL);
    EXPECT_THAT((char *)ptr, StrEq("001"));
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "five2"), LMP_VAR_STRING);
    ptr = lammps_extract_variable(lmp, "five2", NULL);
    EXPECT_THAT((char *)ptr, StrEq("010"));
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "six"), LMP_VAR_STRING);
    ptr = lammps_extract_variable(lmp, "six", NULL);
    EXPECT_THAT((char *)ptr, StrEq("one"));
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "seven"), LMP_VAR_STRING);
    ptr = lammps_extract_variable(lmp, "seven", NULL);
    EXPECT_THAT((char *)ptr, StrEq(" 2.00"));
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "eight"), LMP_VAR_STRING);
    ptr = lammps_extract_variable(lmp, "eight", NULL);
    EXPECT_THAT((char *)ptr, StrEq(""));
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "nine"), LMP_VAR_STRING);
    ptr = lammps_extract_variable(lmp, "nine", NULL);
    EXPECT_THAT((char *)ptr, StrEq("one"));

    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "ten"), LMP_VAR_EQUAL);
    ptr = lammps_extract_variable(lmp, "ten", NULL);
    EXPECT_THAT(*(double *)ptr, 1.0);
    lammps_free(ptr);
    variable->internal_set(variable->find("ten"), 2.5);
    ptr = lammps_extract_variable(lmp, "ten", NULL);
    EXPECT_THAT(*(double *)ptr, 2.5);
    lammps_free(ptr);

    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "ten1"), LMP_VAR_STRING);
    ptr = lammps_extract_variable(lmp, "ten1", NULL);
    EXPECT_THAT((char *)ptr, StrEq("1"));
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "ten2"), LMP_VAR_STRING);
    ptr = lammps_extract_variable(lmp, "ten2", NULL);
    EXPECT_THAT((char *)ptr, StrEq("1"));
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "ten3"), LMP_VAR_STRING);
    ptr = lammps_extract_variable(lmp, "ten3", NULL);
    EXPECT_THAT((char *)ptr, StrEq("1"));

    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "ten4"), LMP_VAR_VECTOR);
    ptr = lammps_extract_variable(lmp, "ten4", (const char *)1);
    double *dptr = (double *)lammps_extract_variable(lmp, "ten4", NULL);
    EXPECT_EQ((*(int *)ptr), 7);
    lammps_free(ptr);
    EXPECT_EQ(dptr[0], 0);
    EXPECT_EQ(dptr[4], 5);
    EXPECT_EQ(dptr[6], 11);
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "ten5"), LMP_VAR_VECTOR);
    ptr = lammps_extract_variable(lmp, "ten5", (const char *)1);
    dptr = (double *)lammps_extract_variable(lmp, "ten5", NULL);
    EXPECT_EQ((*(int *)ptr), 2);
    lammps_free(ptr);
    EXPECT_EQ(dptr[0], 0.5);
    EXPECT_EQ(dptr[1], 1.25);

    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "iswin"), LMP_VAR_EQUAL);
    EXPECT_EQ(lammps_extract_variable_datatype(lmp, "islin"), LMP_VAR_EQUAL);
#if defined(_WIN32)
    ptr = lammps_extract_variable(lmp, "iswin", NULL);
    EXPECT_THAT(*(double *)ptr, 1.0);
    lammps_free(ptr);
    ptr = lammps_extract_variable(lmp, "islin", NULL);
    EXPECT_THAT(*(double *)ptr, 0.0);
    lammps_free(ptr);
#elif defined(__linux__)
    ptr = lammps_extract_variable(lmp, "iswin", NULL);
    EXPECT_THAT(*(double *)ptr, 0.0);
    lammps_free(ptr);
    ptr = lammps_extract_variable(lmp, "islin", NULL);
    EXPECT_THAT(*(double *)ptr, 1.0);
    lammps_free(ptr);
#else
    ptr = lammps_extract_variable(lmp, "iswin", NULL);
    EXPECT_THAT(*(double *)ptr, 0.0);
    lammps_free(ptr);
    ptr = lammps_extract_variable(lmp, "islin", NULL);
    EXPECT_THAT(*(double *)ptr, 0.0);
    lammps_free(ptr);
#endif

    LAMMPS_NS::platform::unlink("test_variable.file");
}
