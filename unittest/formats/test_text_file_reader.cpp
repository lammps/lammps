/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "info.h"
#include "input.h"
#include "text_file_reader.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../testing/core.h"

#include <cstring>
#include <iostream>
#include <mpi.h>
#include <vector>

using namespace LAMMPS_NS;
using LAMMPS_NS::utils::split_words;

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

class TextFileReaderTest : public ::testing::Test {

protected:
    void TearDown() override
    {
        unlink("text_reader_one.file");
        unlink("text_reader_two.file");
    }

    void test_files()
    {
        FILE *fp = fopen("text_reader_one.file", "w");
        fputs("# test file 1 for text file reader\n\n\none\n  two  \n\n"
              "three  # with comment\nfour   ! with non-comment\n"
              "# comments only\n        five\n#END\n",
              fp);
        fclose(fp);

        fp = fopen("text_reader_two.file", "w");
        fputs("# test file for atomfile style variable\n\n"
              "4  # four lines\n4 0.5   #with comment\n"
              "2 -0.5         \n3 1.5\n1 -1.5\n\n"
              "2\n10 1.0 # test\n13 1.0\n\n######\n"
              "4\n1 4.0 # test\n2 3.0\n3 2.0\n4 1.0\n#END\n",
              fp);
        fclose(fp);
    }
};

TEST_F(TextFileReaderTest, nofile)
{
    ASSERT_THROW({ TextFileReader reader("text_reader_noexist.file", "test"); },
                 FileReaderException);
}

TEST_F(TextFileReaderTest, permissions)
{
    FILE *fp = fopen("text_reader_noperms.file", "w");
    fputs("word\n", fp);
    fclose(fp);
    chmod("text_reader_noperms.file", 0);
    ASSERT_THROW({ TextFileReader reader("text_reader_noperms.file", "test"); },
                 FileReaderException);
    unlink("text_reader_noperms.file");
}

TEST_F(TextFileReaderTest, nofp)
{
    ASSERT_THROW({ TextFileReader reader(nullptr, "test"); }, FileReaderException);
}

TEST_F(TextFileReaderTest, usefp)
{
    test_files();
    FILE *fp = fopen("text_reader_two.file", "r");
    ASSERT_NE(fp, nullptr);

    auto reader = new TextFileReader(fp, "test");
    auto line   = reader->next_line();
    ASSERT_STREQ(line, "4  ");
    line = reader->next_line(1);
    ASSERT_STREQ(line, "4 0.5   ");
    ASSERT_NO_THROW({ reader->skip_line(); });
    auto values = reader->next_values(1);
    ASSERT_EQ(values.count(), 2);
    ASSERT_EQ(values.next_int(), 3);
    ASSERT_STREQ(values.next_string().c_str(), "1.5");
    ASSERT_NE(reader->next_line(), nullptr);
    double data[20];
    ASSERT_THROW({ reader->next_dvector(data, 20); }, FileReaderException);
    ASSERT_THROW({ reader->skip_line(); }, EOFException);
    ASSERT_EQ(reader->next_line(), nullptr);
    delete reader;

    // check that we reached EOF and the destructor didn't close the file.
    ASSERT_EQ(feof(fp), 1);
    ASSERT_EQ(fclose(fp), 0);
}

TEST_F(TextFileReaderTest, comments)
{
    test_files();
    TextFileReader reader("text_reader_two.file", "test");
    reader.ignore_comments = true;
    auto line              = reader.next_line();
    ASSERT_STREQ(line, "4  ");
    line = reader.next_line(1);
    ASSERT_STREQ(line, "4 0.5   ");
    ASSERT_NO_THROW({ reader.skip_line(); });
    auto values = reader.next_values(1);
    ASSERT_EQ(values.count(), 2);
    ASSERT_EQ(values.next_int(), 3);
    ASSERT_STREQ(values.next_string().c_str(), "1.5");
    ASSERT_NE(reader.next_line(), nullptr);
    double data[20];
    ASSERT_THROW({ reader.next_dvector(data, 20); }, FileReaderException);
    ASSERT_THROW({ reader.skip_line(); }, EOFException);
    ASSERT_EQ(reader.next_line(), nullptr);
}

TEST_F(TextFileReaderTest, nocomments)
{
    test_files();
    TextFileReader reader("text_reader_one.file", "test");
    reader.ignore_comments = false;
    auto line              = reader.next_line();
    ASSERT_STREQ(line, "# test file 1 for text file reader\n");
    line = reader.next_line(1);
    ASSERT_STREQ(line, "one\n");
    ASSERT_NO_THROW({ reader.skip_line(); });
    auto values = reader.next_values(4);
    ASSERT_EQ(values.count(), 4);
    ASSERT_STREQ(values.next_string().c_str(), "three");
    ASSERT_STREQ(values.next_string().c_str(), "#");
    ASSERT_STREQ(values.next_string().c_str(), "with");
    try {
        reader.next_values(100);
        FAIL() << "No exception thrown\n";
    } catch (EOFException &e) {
        ASSERT_STREQ(e.what(), "Incorrect format in test file! 9/100 parameters");
    }
    ASSERT_THROW({ reader.skip_line(); }, EOFException);
    ASSERT_EQ(reader.next_line(), nullptr);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (Info::get_mpi_vendor() == "Open MPI" && !LAMMPS_NS::Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. "
                     "Death tests will be skipped\n";

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }
    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
