
#include "platform.h"
#include "utils.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>

using namespace LAMMPS_NS;
using testing::EndsWith;
using testing::Eq;
using testing::IsEmpty;
using testing::StartsWith;
using testing::StrEq;

TEST(Platform, clock)
{
    const double wt_start = platform::walltime();
    const double ct_start = platform::cputime();

    // spend some time computing pi
    constexpr double known_pi = 3.141592653589793238462643;
    constexpr int n           = 100000000;
    constexpr double h        = 1.0 / (double)n;
    double my_pi              = 0.0, x;
    for (int i = 0; i < n; ++i) {
        x = h * ((double)i + 0.5);
        my_pi += 4.0 / (1.0 + x * x);
    }
    my_pi *= h;
    const double wt_used = platform::walltime() - wt_start;
    const double ct_used = platform::cputime() - ct_start;

    ASSERT_NEAR(my_pi, known_pi, 1e-12);
    ASSERT_GT(wt_used, 1e-4);
    // windows sometimes incorrectly reports used CPU time as 0.0
    if (ct_used != 0.0) ASSERT_GT(ct_used, 1e-4);
}

TEST(Platform, putenv_unsetenv)
{
    const char *var = getenv("UNITTEST_VAR1");
    ASSERT_EQ(var, nullptr);
    int rv = platform::putenv("UNITTEST_VAR1");
    var    = getenv("UNITTEST_VAR1");
    ASSERT_EQ(rv, 0);
    ASSERT_NE(var, nullptr);
    // we cannot set environment variables without a value on windows with _putenv()
#if defined(_WIN32)
    ASSERT_THAT(var, StrEq("1"));
#else
    ASSERT_THAT(var, StrEq(""));
#endif

    rv  = platform::putenv("UNITTEST_VAR1=one");
    var = getenv("UNITTEST_VAR1");
    ASSERT_EQ(rv, 0);
    ASSERT_NE(var, nullptr);
    ASSERT_THAT(var, StrEq("one"));

    rv  = platform::putenv("UNITTEST_VAR1=one=two");
    var = getenv("UNITTEST_VAR1");
    ASSERT_EQ(rv, 0);
    ASSERT_NE(var, nullptr);
    ASSERT_THAT(var, StrEq("one=two"));

    ASSERT_EQ(platform::putenv(""), -1);

    ASSERT_EQ(platform::unsetenv(""), -1);
    ASSERT_EQ(platform::unsetenv("UNITTEST_VAR3=two"), -1);
    var = getenv("UNITTEST_VAR1");
    ASSERT_NE(var, nullptr);
    ASSERT_EQ(platform::unsetenv("UNITTEST_VAR1"), 0);
    var = getenv("UNITTEST_VAR1");
    ASSERT_EQ(var, nullptr);
}

TEST(Platform, list_pathenv)
{
    auto dirs = platform::list_pathenv("PATH");
    ASSERT_GT(dirs.size(), 1);
}

TEST(Platform, find_cmd_path)
{
#if defined(_WIN32)
    ASSERT_THAT(platform::find_exe_path("notepad"), EndsWith("\\notepad.exe"));
    ASSERT_THAT(platform::find_exe_path("cmd"), EndsWith("\\cmd.exe"));
    ASSERT_THAT(platform::find_exe_path("some_bogus_command"), IsEmpty());
#else
    ASSERT_THAT(platform::find_exe_path("ls"), EndsWith("bin/ls"));
    ASSERT_THAT(platform::find_exe_path("sh"), EndsWith("bin/sh"));
    ASSERT_THAT(platform::find_exe_path("some_bogus_command"), IsEmpty());
#endif
}

#if defined(TEST_SHARED_LOAD)
#define stringify(s) mkstring(s)
#define mkstring(s) #s
TEST(Platform, sharedload)
{
    const std::vector<std::string> objs = {stringify(TEST_SHARED_OBJ), stringify(TEST_SHARED_LIB)};
    const int *intvar;
    const double *doublevar;
    void *handle;
    int (*intfunc)(int);
    double (*doublefunc)(double, int);

    for (const auto &obj : objs) {
        handle = platform::dlopen(obj);
        EXPECT_NE(handle, nullptr);
        intvar = (int *)platform::dlsym(handle, "some_int_val");
        EXPECT_NE(intvar, nullptr);
        EXPECT_EQ(*intvar, 12345);
        doublevar = (double *)platform::dlsym(handle, "some_double_val");
        EXPECT_NE(doublevar, nullptr);
        EXPECT_DOUBLE_EQ(*doublevar, 6.78e-9);
        intfunc = (int (*)(int))platform::dlsym(handle, "some_int_function");
        EXPECT_NE(intfunc, nullptr);
        EXPECT_EQ((*intfunc)(12), 144);
        doublefunc = (double (*)(double, int))platform::dlsym(handle, "some_double_function");
        EXPECT_NE(doublefunc, nullptr);
        EXPECT_DOUBLE_EQ((*doublefunc)(0.5, 6), 3.0);
        EXPECT_EQ(platform::dlsym(handle, "some_nonexisting_symbol"), nullptr);
        EXPECT_EQ(platform::dlclose(handle), 0);
    }
}
#undef stringify
#undef mkstring
#endif

TEST(Platform, guesspath)
{
    char buf[256];
    FILE *fp = fopen("test_guesspath.txt", "w");
#if defined(__linux__) || defined(__APPLE__) || defined(_WIN32)
    const char *path = platform::guesspath(fp, buf, sizeof(buf));
    ASSERT_THAT(path, EndsWith("test_guesspath.txt"));
#else
    const char *path = platform::guesspath(fp, buf, sizeof(buf));
    ASSERT_THAT(path, EndsWith("(unknown)"));
#endif
    fclose(fp);
    platform::unlink("test_guesspath.txt");
}

TEST(Platform, unlink)
{
    const char test[] = "12345678901234567890";
    platform::unlink("unlink.dat");
    ASSERT_EQ(platform::unlink("dummy.dat"), -1);
    FILE *fp = fopen("unlink.dat", "w");
    fwrite(test, sizeof(test), 1, fp);
    fclose(fp);
    ASSERT_EQ(platform::unlink("unlink.dat"), 0);
    ASSERT_EQ(platform::unlink("unlink.dat"), -1);
    fp = fopen("unlink.dat", "r");
    ASSERT_EQ(fp, nullptr);

    platform::mkdir("unlink.dir");
    ASSERT_EQ(platform::unlink("unlink.dir"), -1);
    platform::rmdir("unlink.dir");
}

TEST(Platform, fseek_ftell)
{
    const char test[] = "12345678901234567890";
    platform::unlink("seek_tell.dat");
    FILE *fp = fopen("seek_tell.dat", "w");
    fwrite(test, sizeof(test), 1, fp);
    fflush(fp);
    ASSERT_EQ(platform::ftell(fp), sizeof(test));
    fclose(fp);
    fp = fopen("seek_tell.dat", "r+");
    ASSERT_EQ(fgetc(fp), '1');
    ASSERT_EQ(fgetc(fp), '2');
    ASSERT_EQ(platform::ftell(fp), 2);
    ASSERT_EQ(platform::fseek(fp, 15), 0);
    ASSERT_EQ(fgetc(fp), '6');
    fflush(fp);
    ASSERT_EQ(platform::fseek(fp, platform::END_OF_FILE), 0);
    ASSERT_EQ(platform::ftell(fp), 21);
    ASSERT_EQ(platform::fseek(fp, 20), 0);
    ASSERT_EQ(fgetc(fp), 0);
    ASSERT_EQ(platform::ftell(fp), 21);
    fclose(fp);
    platform::unlink("seek_tell.dat");
}

TEST(Platform, ftruncate)
{
    platform::unlink("truncate.dat");
    FILE *fp = fopen("truncate.dat", "w");
    fputs("header one\n", fp);
    fputs("header two\n", fp);
    fflush(fp);
    bigint filepos = platform::ftell(fp);
    fputs("line one\n", fp);
    fputs("line two\n", fp);
    fputs("line three\n", fp);
    fflush(fp);
    ASSERT_EQ(platform::ftruncate(fp, filepos), 0);
    fputs("line four\n", fp);
    ASSERT_GT(platform::ftell(fp), filepos);
    fputs("line five\n", fp);
    fflush(fp);
    fclose(fp);

    // check file
    fp = fopen("truncate.dat", "r");
    char buf[128];
    char *ptr = fgets(buf, 127, fp);
    ASSERT_THAT(ptr, StartsWith("header one"));
    ptr = fgets(buf, 127, fp);
    ASSERT_THAT(ptr, StartsWith("header two"));
    ptr = fgets(buf, 127, fp);
    ASSERT_THAT(ptr, StartsWith("line four"));
    ptr = fgets(buf, 127, fp);
    ASSERT_THAT(ptr, StartsWith("line five"));
    ptr = fgets(buf, 127, fp);
    ASSERT_EQ(ptr, nullptr);
    fclose(fp);
    platform::unlink("truncate.dat");
}

TEST(Platform, path_basename)
{
#if defined(_WIN32)
    ASSERT_THAT(platform::path_basename("c:\\parent\\folder\\filename"), Eq("filename"));
    ASSERT_THAT(platform::path_basename("folder\\"), Eq(""));
    ASSERT_THAT(platform::path_basename("c:/parent/folder/filename"), Eq("filename"));
#else
    ASSERT_THAT(platform::path_basename("/parent/folder/filename"), Eq("filename"));
    ASSERT_THAT(platform::path_basename("/parent/folder/"), Eq(""));
#endif
}

TEST(Platform, path_dirname)
{
#if defined(_WIN32)
    ASSERT_THAT(platform::path_dirname("c:/parent/folder/filename"), Eq("c:/parent/folder"));
    ASSERT_THAT(platform::path_dirname("c:\\parent\\folder\\filename"), Eq("c:\\parent\\folder"));
    ASSERT_THAT(platform::path_dirname("c:filename"), Eq("."));
#else
    ASSERT_THAT(platform::path_dirname("/parent/folder/filename"), Eq("/parent/folder"));
#endif
    ASSERT_THAT(platform::path_dirname("filename"), Eq("."));
}

TEST(Platform, path_join)
{
#if defined(_WIN32)
    ASSERT_THAT(platform::path_join("c:\\folder", "filename"), Eq("c:\\folder\\filename"));
    ASSERT_THAT(platform::path_join("c:\\folder\\", "filename"), Eq("c:\\folder\\filename"));
    ASSERT_THAT(platform::path_join("c:\\folder", "\\filename"), Eq("c:\\folder\\filename"));
    ASSERT_THAT(platform::path_join("c:\\folder\\", "\\filename"), Eq("c:\\folder\\filename"));
    ASSERT_THAT(platform::path_join("c:\\folder", "/filename"), Eq("c:\\folder\\filename"));
    ASSERT_THAT(platform::path_join("c:\\folder\\\\", "\\filename"), Eq("c:\\folder\\filename"));
    ASSERT_THAT(platform::path_join("c:\\folder\\", "\\\\filename"), Eq("c:\\folder\\filename"));
    ASSERT_THAT(platform::path_join("c:\\folder/\\", "/\\filename"), Eq("c:\\folder\\filename"));
    ASSERT_THAT(platform::path_join("c:\\folder\\/", "\\/filename"), Eq("c:\\folder\\filename"));
    ASSERT_THAT(platform::path_join("c:\\folder", ""), Eq("c:\\folder"));
    ASSERT_THAT(platform::path_join("", "\\/filename"), Eq("\\/filename"));
#else
    ASSERT_THAT(platform::path_join("/parent/folder", "filename"), Eq("/parent/folder/filename"));
    ASSERT_THAT(platform::path_join("/parent/folder/", "filename"), Eq("/parent/folder/filename"));
    ASSERT_THAT(platform::path_join("/parent/folder", "/filename"), Eq("/parent/folder/filename"));
    ASSERT_THAT(platform::path_join("/parent/folder/", "/filename"), Eq("/parent/folder/filename"));
    ASSERT_THAT(platform::path_join("/parent/folder//", "filename"), Eq("/parent/folder/filename"));
    ASSERT_THAT(platform::path_join("/parent/folder", "//filename"), Eq("/parent/folder/filename"));
    ASSERT_THAT(platform::path_join("/parent/folder///", "/filename"),
                Eq("/parent/folder/filename"));
    ASSERT_THAT(platform::path_join("/parent/folder/", "///filename"),
                Eq("/parent/folder/filename"));
    ASSERT_THAT(platform::path_join("/parent/folder/", ""), Eq("/parent/folder/"));
    ASSERT_THAT(platform::path_join("", "\\/filename"), Eq("\\/filename"));
#endif
}

TEST(Platform, is_console)
{
    platform::unlink("file_is_no_console.txt");
    FILE *fp = fopen("file_is_no_console.txt", "w");
    fputs("some text\n", fp);
    EXPECT_FALSE(platform::is_console(fp));
    fclose(fp);
    platform::unlink("file_is_no_console.txt");
}

TEST(Platform, path_and_directory)
{
    platform::unlink("path_is_directory");
    platform::rmdir("path_is_directory");
    platform::unlink("path_is_file");
    platform::mkdir("path_is_directory");
    FILE *fp = fopen("path_is_file", "w");
    fputs("some text\n", fp);
    fclose(fp);

    ASSERT_TRUE(platform::path_is_directory("path_is_directory"));
    ASSERT_FALSE(platform::path_is_directory("path_is_file"));
    ASSERT_FALSE(platform::path_is_directory("path_does_not_exist"));
    platform::unlink("path_is_file");

#if defined(_WIN32)
    fp = fopen("path_is_directory\\path_is_file", "w");
#else
    fp = fopen("path_is_directory/path_is_file", "w");
#endif
    fputs("some text\n", fp);
    fclose(fp);
#if defined(_WIN32)
    platform::mkdir("path_is_directory\\path_is_directory");
    fp = fopen("path_is_directory\\path_is_other_file", "w");
#else
    platform::mkdir("path_is_directory/path_is_directory");
    fp = fopen("path_is_directory/path_is_other_file", "w");
#endif
    fputs("some text\n", fp);
    fclose(fp);
    auto dirs = platform::list_directory("path_is_directory");
    ASSERT_EQ(dirs.size(), 3);
    platform::rmdir("path_is_directory");
    ASSERT_FALSE(platform::path_is_directory("path_is_directory"));

#if defined(_WIN32)
    ASSERT_EQ(platform::mkdir("path_is_directory\\path_is_directory"), 0);
    ASSERT_TRUE(platform::path_is_directory("path_is_directory\\path_is_directory"));
#else
    ASSERT_EQ(platform::mkdir("path_is_directory/path_is_directory"), 0);
    ASSERT_TRUE(platform::path_is_directory("path_is_directory/path_is_directory"));
#endif
    platform::rmdir("path_is_directory");
}

TEST(Platform, get_change_directory)
{
    platform::unlink("working_directory");
    platform::rmdir("working_directory");

    auto cwd = platform::current_directory();
    ASSERT_GT(cwd.size(), 0);

    platform::mkdir("working_directory");
    ASSERT_EQ(platform::chdir("working_directory"), 0);
    ASSERT_THAT(platform::current_directory(), EndsWith("working_directory"));

    ASSERT_EQ(platform::chdir(".."), 0);
    ASSERT_THAT(platform::current_directory(), StrEq(cwd));
    platform::rmdir("working_directory");
}

TEST(Platform, file_is_readable)
{
    platform::unlink("file_is_readable.txt");
    FILE *fp = fopen("file_is_readable.txt", "w");
    fputs("some text\n", fp);
    fclose(fp);

    ASSERT_TRUE(platform::file_is_readable("file_is_readable.txt"));
    ASSERT_FALSE(platform::file_is_readable("file_does_not_exist.txt"));
    platform::unlink("file_is_readable.txt");

    // windows does not have permission flags
#if !defined(_WIN32)
    platform::unlink("file_is_not_readable.txt");
    fp = fopen("file_is_not_readable.txt", "w");
    fputs("some text\n", fp);
    fclose(fp);
    chmod("file_is_not_readable.txt", 0);
    ASSERT_FALSE(platform::file_is_readable("file_is_not_readable.txt"));
    platform::unlink("file_is_not_readable.txt");
#endif
}

TEST(Platform, has_compress_extension)
{
    ASSERT_FALSE(platform::has_compress_extension("dummy"));
    ASSERT_FALSE(platform::has_compress_extension("dum.my"));
    ASSERT_TRUE(platform::has_compress_extension("dummy.gz"));
    ASSERT_TRUE(platform::has_compress_extension("dummy.bz2"));
    ASSERT_TRUE(platform::has_compress_extension("dummy.zst"));
    ASSERT_TRUE(platform::has_compress_extension("dummy.xz"));
    ASSERT_TRUE(platform::has_compress_extension("dummy.lzma"));
    ASSERT_TRUE(platform::has_compress_extension("dummy.lz4"));
}

TEST(Platform, compress_read_write)
{
    const std::vector<std::string> test_files = {"zip_test.zip", "zip_test.gz",  "zip_test.bz2",
                                                 "zip_test.zst", "zip_test.xz",  "zip_test.lzma",
                                                 "zip_test.lz4", "zip_test.unk", "zip test.gz"};
    for (const auto &file : test_files) {
        platform::unlink(file);
        FILE *fp = platform::compressed_write(file);
        if (!fp) {
            platform::unlink(file);
            continue;
        }

        clearerr(fp);
        fputs("line one\n", fp);
        fputs("line two\n", fp);
        ASSERT_EQ(ferror(fp), 0);
        fflush(fp);
        platform::pclose(fp);

        fp = platform::compressed_read(file);
        ASSERT_NE(fp, nullptr);
        char buf[128];
        char *ptr = fgets(buf, 128, fp);
        EXPECT_THAT(ptr, StartsWith("line one"));
        ptr = fgets(buf, 128, fp);
        EXPECT_THAT(ptr, StartsWith("line two"));
        ASSERT_EQ(ferror(fp), 0);
        platform::pclose(fp);
        platform::unlink(file);
    }
}
