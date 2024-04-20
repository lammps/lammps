/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/** \file platform.cpp
 * This file provides abstractions for a variety of platform specific
 * functionality in a namespace "platform".  This is a companion to
 * the "utils" namespace with convenience and utility functions. */

#include "platform.h"

#include "fmt/format.h"
#include "text_file_reader.h"
#include "utils.h"

#include <deque>
#include <exception>
#include <mpi.h>

////////////////////////////////////////////////////////////////////////
// include system headers and tweak system settings
#if defined(_WIN32)

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#if defined(_WIN32_WINNT)
#undef _WIN32_WINNT
#endif

// target Windows version is windows 7 and later
#define _WIN32_WINNT _WIN32_WINNT_WIN7
#define PSAPI_VERSION 2

#include <direct.h>
#include <io.h>    // for _get_osfhandle()
#include <sys/stat.h>
#include <windows.h>

#else    // not Windows  ///////////////////////////////////////////////

#include <dirent.h>
#include <dlfcn.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>
#endif

#if defined(__APPLE__)
#include <fcntl.h>
#include <sys/syslimits.h>
#endif

// for disk_free()
#if defined(__linux__) || defined(__APPLE__) || defined(__FreeBSD__) || defined(__DragonFly__) || \
    defined(__OpenBSD__) || defined(__NetBSD__)
#include <sys/statvfs.h>
#endif

////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cstring>
#include <thread>

/* ------------------------------------------------------------------ */

/// Struct for listing on-the-fly compression/decompression commands
struct compress_info {
  /// identifier for the different compression algorithms
  enum styles { NONE, GZIP, BZIP2, ZSTD, XZ, LZMA, LZ4 };
  const std::string extension;          ///< filename extension for the current algorithm
  const std::string command;            ///< command to perform compression or decompression
  const std::string compressflags;      ///< flags to append to compress from stdin to stdout
  const std::string uncompressflags;    ///< flags to decompress file to stdout
  const int style;                      ///< compression style flag
};

// clang-format off
static const std::vector<compress_info> compress_styles = {
    {"",     "",      "",       "",        compress_info::NONE},
    {"gz",   "gzip",  " > ",    " -cdf ",  compress_info::GZIP},
    {"bz2",  "bzip2", " > ",    " -cdf ",  compress_info::BZIP2},
    {"zst",  "zstd",  " -q > ", " -cdf ",  compress_info::ZSTD},
    {"xz",   "xz",    " > ",    " -cdf ",  compress_info::XZ},
    {"lzma", "xz", " --format=lzma > ", " --format=lzma -cdf ", compress_info::LZMA},
    {"lz4",  "lz4",   " > ",    " -cdf ",  compress_info::LZ4},
};
// clang-format on

/* ------------------------------------------------------------------ */

static const compress_info &find_compress_type(const std::string &file)
{
  std::size_t dot = file.find_last_of('.');
  if (dot != std::string::npos) {
    const std::string ext = file.substr(dot + 1);
    for (const auto &i : compress_styles) {
      if (i.extension == ext) return i;
    }
  }
  return compress_styles[0];
}

/* ------------------------------------------------------------------ */

// set reference time stamp during executable/library init.
// should provide better resolution than using epoch, if the system clock supports it.
static auto initial_time = std::chrono::steady_clock::now();

using namespace LAMMPS_NS;

// get CPU time

// clang-format off
// clang compilers are optimizing this function too aggressively returning always 0
#if defined(__clang__)
[[clang::optnone]]
#elif defined(_MSC_VER)
#pragma optimize("",off)
#endif
double platform::cputime()
// clang-format on
{
  double rv = 0.0;

#ifdef _WIN32

  // from MSD docs.
  FILETIME ct, et, kt, ut;
  union {
    FILETIME ft;
    uint64_t ui;
  } cpu;
  if (GetProcessTimes(GetCurrentProcess(), &ct, &et, &kt, &ut)) {
    cpu.ft = ut;
    rv = cpu.ui * 0.0000001;
  }

#else /* ! _WIN32 */

  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
    rv = (double) ru.ru_utime.tv_sec;
    rv += (double) ru.ru_utime.tv_usec * 0.000001;
  }

#endif

  return rv;
}
#if defined(__clang__)
#elif defined(_MSC_VER)
#pragma optimize("", on)
#endif

/* ----------------------------------------------------------------------
   get wall time
------------------------------------------------------------------------ */
double platform::walltime()
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - initial_time).count();
}

/* ----------------------------------------------------------------------
   sleep with microsecond resolution
------------------------------------------------------------------------ */
void platform::usleep(int usec)
{
  return std::this_thread::sleep_for(std::chrono::microseconds(usec));
}

/* ----------------------------------------------------------------------
   get Operating system and version info
------------------------------------------------------------------------- */

std::string platform::os_info()
{
  std::string buf;

#if defined(_WIN32)

  // Get Windows Edition name from registry
  char value[1024];
  DWORD value_length = 1024;
  const char *subkey = "SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion";
  const char *entry = "CurrentBuild";
  RegGetValue(HKEY_LOCAL_MACHINE, subkey, entry, RRF_RT_REG_SZ, nullptr, &value,
              (LPDWORD) &value_length);
  // enforce zero termination
  value[1023] = '\0';
  auto build = std::string(value);

  if (build == "6002") {
    buf = "Windows Vista";
  } else if (build == "6003") {
    buf = "Windows Server 2008";
  } else if (build == "7601") {
    buf = "Windows 7";
  } else if (build == "9200") {
    buf = "Windows 8";
  } else if (build == "9600") {
    buf = "Windows 8.1";
  } else if (build == "10240") {
    buf = "Windows 10 1507";
  } else if (build == "10586") {
    buf = "Windows 10 1511";
  } else if (build == "14393") {
    buf = "Windows 10 1607";
  } else if (build == "15063") {
    buf = "Windows 10 1703";
  } else if (build == "16299") {
    buf = "Windows 10 1709";
  } else if (build == "17134") {
    buf = "Windows 10 1803";
  } else if (build == "17763") {
    buf = "Windows 10 1809";
  } else if (build == "18362") {
    buf = "Windows 10 1903";
  } else if (build == "18363") {
    buf = "Windows 10 1909";
  } else if (build == "19041") {
    buf = "Windows 10 2004";
  } else if (build == "19042") {
    buf = "Windows 10 20H2";
  } else if (build == "19043") {
    buf = "Windows 10 21H1";
  } else if (build == "19044") {
    buf = "Windows 10 21H2";
  } else if (build == "19045") {
    buf = "Windows 10 22H2";
  } else if (build == "20348") {
    buf = "Windows Server 2022";
  } else if (build == "22000") {
    buf = "Windows 11 21H2";
  } else if (build == "22621") {
    buf = "Windows 11 22H2";
  } else if (build == "22631") {
    buf = "Windows 11 23H2";
  } else {
    buf = "Windows Build " + build;
  }
  DWORD fullversion, majorv, minorv, buildv = 0;
  fullversion = GetVersion();
  majorv = (DWORD) (LOBYTE(LOWORD(fullversion)));
  minorv = (DWORD) (HIBYTE(LOWORD(fullversion)));
  if (fullversion < 0x80000000) buildv = (DWORD) (HIWORD(fullversion));

  buf += ", Windows ABI " + std::to_string(majorv) + "." + std::to_string(minorv) + " (" +
      std::to_string(buildv) + ") on ";

  SYSTEM_INFO si;
  GetSystemInfo(&si);

  switch (si.wProcessorArchitecture) {
    case PROCESSOR_ARCHITECTURE_AMD64:
      buf += "x86_64";
      break;
    case PROCESSOR_ARCHITECTURE_ARM:
      buf += "arm";
      break;
    case PROCESSOR_ARCHITECTURE_IA64:
      buf += "ia64";
      break;
    case PROCESSOR_ARCHITECTURE_INTEL:
      buf += "i386";
      break;
    default:
      buf += "(unknown)";
  }
#else
  struct utsname ut;
  uname(&ut);

  // try to get OS distribution name, if available
  buf = ut.sysname;

  if (platform::file_is_readable("/etc/os-release")) {
    try {
      TextFileReader reader("/etc/os-release", "");
      while (true) {
        auto words = reader.next_values(0, "=");
        if ((words.count() > 1) && (words.next_string() == "PRETTY_NAME")) {
          buf += " " + utils::trim(words.next_string());
          break;
        }
      }
    } catch (std::exception &e) {
      ;    // EOF but keyword not found
    }
  }

  buf += std::string(" ") + ut.release + " " + ut.machine;
#endif
  return buf;
}

/* ----------------------------------------------------------------------
   identify C++ standard version
------------------------------------------------------------------------- */

std::string platform::cxx_standard()
{
#if __cplusplus > 202002L
  return "newer than C++20";
#elif __cplusplus == 202002L
  return "C++20";
#elif __cplusplus == 201703L
  return "C++17";
#elif __cplusplus == 201402L
  return "C++14";
#elif __cplusplus == 201103L
  return "C++11";
#elif __cplusplus == 199711L
  return "C++98";
#else
  return "unknown";
#endif
}

/* ----------------------------------------------------------------------
   identify compiler and its version
------------------------------------------------------------------------- */

std::string platform::compiler_info()
{
  std::string buf = "(Unknown)";
#if defined(__INTEL_LLVM_COMPILER)
  double version = static_cast<double>(__INTEL_LLVM_COMPILER) * 0.01;
  buf = fmt::format("Intel LLVM C++ {:.1f} / {}", version, __VERSION__);
#elif defined(__ibmxl__)
  buf = fmt::format("IBM XL C/C++ (Clang) {}.{}.{}", __ibmxl_version__, __ibmxl_release__,
                    __ibmxl_modification__);
#elif defined(__clang__)
  buf = fmt::format("Clang C++ {}", __VERSION__);
#elif defined(__PGI)
  buf = fmt::format("PGI C++ {}.{}", __PGIC__, __PGIC_MINOR__);
#elif defined(__INTEL_COMPILER)
#if !defined(__VERSION__)
#define __VERSION__ __INTEL_COMPILER_BUILD_DATE
#endif
  double version = static_cast<double>(__INTEL_COMPILER) * 0.01;
  buf = fmt::format("Intel Classic C++ {:.2f}.{} / {}", version, __INTEL_COMPILER_UPDATE,
                    __VERSION__);
#elif defined(__MINGW64__)
  buf = fmt::format("MinGW-w64 64bit {}.{} / GNU C++ {}", __MINGW64_VERSION_MAJOR,
                    __MINGW64_VERSION_MINOR, __VERSION__);
#elif defined(__MINGW32__)
  buf = fmt::format("MinGW-w64 32bit {}.{} / GNU C++ {}", __MINGW32_MAJOR_VERSION,
                    __MINGW32_MINOR_VERSION, __VERSION__);
#elif defined(__GNUC__)
  buf = fmt::format("GNU C++ {}", __VERSION__);
#elif defined(_MSC_VER) && (_MSC_VER >= 1920) && (_MSC_VER < 1930)
  constexpr int major = _MSC_VER / 100;
  constexpr int minor = _MSC_VER - major * 100;
  constexpr int patch = minor - 20;
  buf = fmt::format("Microsoft Visual Studio 2019 Version 16.{}, C/C++ {}.{}", patch, major - 5,
                    minor);
#elif defined(_MSC_VER) && (_MSC_VER >= 1930) && (_MSC_VER < 2000)
  constexpr int major = _MSC_VER / 100;
  constexpr int minor = _MSC_VER - major * 100;
  constexpr int patch = minor - 30;
  buf = fmt::format("Microsoft Visual Studio 2022 Version 17.{}, C/C++ {}.{}", patch, major - 5,
                    minor);
#else
  buf = "(Unknown)";
#endif
  return buf;
}

/* ----------------------------------------------------------------------
   detect OpenMP standard
------------------------------------------------------------------------- */

std::string platform::openmp_standard()
{

#if !defined(_OPENMP)
  return "OpenMP not enabled";
#else

  // Supported OpenMP version corresponds to the release date of the
  // specifications as posted at https://www.openmp.org/specifications/

#if _OPENMP > 202411
  return "OpenMP newer than version 6.0";
#elif _OPENMP == 202411
  return "OpenMP 6.0";
#elif _OPENMP == 202311
  return "OpenMP 6.0 preview 2";
#elif _OPENMP == 202211
  return "OpenMP 6.0 preview 1";
#elif _OPENMP == 202111
  return "OpenMP 5.2";
#elif _OPENMP == 202011
  return "OpenMP 5.1";
#elif _OPENMP == 201811
  return "OpenMP 5.0";
#elif _OPENMP == 201611
  return "OpenMP 5.0 preview 1";
#elif _OPENMP == 201511
  return "OpenMP 4.5";
#elif _OPENMP == 201307
  return "OpenMP 4.0";
#elif _OPENMP == 201107
  return "OpenMP 3.1";
#elif _OPENMP == 200805
  return "OpenMP 3.0";
#elif _OPENMP == 200505
  return "OpenMP 2.5";
#elif _OPENMP == 200203
  return "OpenMP 2.0";
#else
  return "unknown OpenMP version";
#endif

#endif
}

/* ----------------------------------------------------------------------
   identify MPI vendor from defines in the mpi.h file.
------------------------------------------------------------------------- */

std::string platform::mpi_vendor()
{
#if defined(MPI_STUBS)
  return "MPI STUBS";
#elif defined(OPEN_MPI)
  return "Open MPI";
#elif defined(MPICH_NAME)
  return "MPICH";
#elif defined(I_MPI_VERSION)
  return "Intel MPI";
#elif defined(PLATFORM_MPI)
  return "Platform MPI";
#elif defined(HP_MPI)
  return "HP MPI";
#elif defined(MSMPI_VER)
  // Get Microsoft MPI version from registry
  char value[1024];
  DWORD value_length = 1024;
  const char *subkey = "SOFTWARE\\Microsoft\\MPI";
  const char *entry = "Version";
  auto rv = RegGetValueA(HKEY_LOCAL_MACHINE, subkey, entry, RRF_RT_REG_SZ, nullptr, &value,
                         (LPDWORD) &value_length);
  std::string buf = "Microsoft MPI";
  if (rv == ERROR_SUCCESS) buf += std::string(" v") + value;
  return buf;
#else
  return "Unknown MPI implementation";
#endif
}

/* ----------------------------------------------------------------------
   detect MPI version info
------------------------------------------------------------------------- */

std::string platform::mpi_info(int &major, int &minor)
{
#if (defined(MPI_VERSION) && (MPI_VERSION > 2)) || defined(MPI_STUBS)
  int len = 0;
  static char version[MPI_MAX_LIBRARY_VERSION_STRING];
  MPI_Get_library_version(version, &len);
  if (len > 80) {
    char *ptr = strchr(version + 80, '\n');
    if (ptr) *ptr = '\0';
  }
#else
  constexpr int MAX_VERSION_STRING = 32;
  static char version[MAX_VERSION_STRING];
  strncpy(version, mpi_vendor().c_str(), MAX_VERSION_STRING);
#endif

#if defined(MPI_VERSION)
  MPI_Get_version(&major, &minor);
#else
  major = 1;
  minor = 0;
#endif
  return {version};
}

/* ----------------------------------------------------------------------
   collect available compression tool info
------------------------------------------------------------------------- */

std::string platform::compress_info()
{
  std::string buf = "Available compression formats:\n\n";
  bool none_found = true;
  for (const auto &cmpi : compress_styles) {
    if (cmpi.style == ::compress_info::NONE) continue;
    if (find_exe_path(cmpi.command).size()) {
      none_found = false;
      buf += fmt::format("Extension: .{:6} Command: {}\n", cmpi.extension, cmpi.command);
    }
  }
  if (none_found) buf += "None\n";
  return buf;
}
/* ----------------------------------------------------------------------
   set environment variable
------------------------------------------------------------------------- */

int platform::putenv(const std::string &vardef)
{
  if (vardef.size() == 0) return -1;

  auto found = vardef.find_first_of('=');
#ifdef _WIN32
  // must assign a value to variable with _putenv_s()
  if (found == std::string::npos)
    return _putenv_s(vardef.c_str(), "1");
  else
    return _putenv_s(vardef.substr(0, found).c_str(), vardef.substr(found + 1).c_str());
#else
  if (found == std::string::npos)
    return setenv(vardef.c_str(), "", 1);
  else
    return setenv(vardef.substr(0, found).c_str(), vardef.substr(found + 1).c_str(), 1);
#endif
  return -1;
}

/* ----------------------------------------------------------------------
   unset environment variable
------------------------------------------------------------------------- */

int platform::unsetenv(const std::string &variable)
{
  if (variable.size() == 0) return -1;
#ifdef _WIN32
  // emulate POSIX semantics by returning -1 on trying to unset non-existing variable
  const char *ptr = getenv(variable.c_str());
  if (!ptr) return -1;
  // empty _putenv_s() definition deletes variable
  return _putenv_s(variable.c_str(), "");
#else
  return ::unsetenv(variable.c_str());
#endif
}

/* ----------------------------------------------------------------------
   split a "path" environment variable into a list
------------------------------------------------------------------------- */

std::vector<std::string> platform::list_pathenv(const std::string &var)
{
  std::vector<std::string> dirs;
  const char *ptr = getenv(var.c_str());
  if (ptr == nullptr) return dirs;

  std::string pathvar = ptr;
  std::size_t first = 0, next;
  while (true) {
    next = pathvar.find_first_of(pathvarsep, first);
    if (next == std::string::npos) {
      dirs.push_back(pathvar.substr(first));
      break;
    } else {
      dirs.push_back(pathvar.substr(first, next - first));
      first = next + 1;
    }
  }
  return dirs;
}

/* ----------------------------------------------------------------------
   find the full path name of an executable
------------------------------------------------------------------------- */

std::string platform::find_exe_path(const std::string &cmd)
{
  if (cmd.size() == 0) return "";
  auto pathdirs = list_pathenv("PATH");
#ifdef _WIN32
  // windows always looks in "." and does it first
  pathdirs.insert(pathdirs.begin(), ".");
#else
  struct stat info;
#endif
  for (const auto &dir : pathdirs) {
    std::string exe = path_join(dir, cmd);
#ifdef _WIN32
    const char *extensions[] = {".exe", ".com", ".bat", nullptr};
    for (auto ext = extensions; *ext != nullptr; ++ext) {
      auto exe_path = exe + *ext;
      if (file_is_readable(exe_path)) return exe_path;
    }
#else
    memset(&info, 0, sizeof(info));
    if (stat(exe.c_str(), &info) != 0) continue;
    if ((info.st_mode & (S_IXOTH | S_IXGRP | S_IXUSR)) != 0) return exe;
#endif
  }
  return "";
}

/* ----------------------------------------------------------------------
   wrapper functions for loading shared objects and libraries
------------------------------------------------------------------------- */

#ifdef _WIN32

// open a shared object file
void *platform::dlopen(const std::string &fname)
{
  return (void *) LoadLibrary(fname.c_str());
}

// return dynamic linker error string

std::string platform::dlerror()
{
  return "";
}

// close a shared object
int platform::dlclose(void *handle)
{
  /* FreeLibrary returns nonzero on success unlike dlclose() */
  return (FreeLibrary((HINSTANCE) handle) == 0);
}

// resolve a symbol in shared object
void *platform::dlsym(void *handle, const std::string &symbol)
{
  return (void *) GetProcAddress((HINSTANCE) handle, symbol.c_str());
}

#else

// open a shared object file
void *platform::dlopen(const std::string &fname)
{
  return ::dlopen(fname.c_str(), RTLD_NOW | RTLD_GLOBAL);
}

// return dynamic linker error string

std::string platform::dlerror()
{
  const char *errmesg = ::dlerror();
  if (errmesg)
    return {errmesg};
  else
    return {""};
}

// close a shared object
int platform::dlclose(void *handle)
{
  return ::dlclose(handle);
}

// resolve a symbol in shared object
void *platform::dlsym(void *handle, const std::string &symbol)
{
  return ::dlsym(handle, symbol.c_str());
}
#endif

/* ---------------------------------------------------------------------- */

/** On Linux the folder /proc/self/fd holds symbolic links to the actual
 * pathnames associated with each open file descriptor of the current process.
 * On MacOS the same kind of information can be obtained using ``fcntl(fd,F_GETPATH,buf)``.
 * On Windows we use ``GetFinalPathNameByHandleA()`` which is available with
 * Windows Vista and later. If the buffer is too small (< 16 bytes) a null pointer is returned.
 *
 * This function is used to provide a filename with error messages in functions
 * where the filename is not passed as an argument, but the FILE * pointer.  */

const char *platform::guesspath(FILE *fp, char *buf, int len)
{
  // no point in guessing a path with a short buffer or NULL pointer as buffer
  if ((buf == nullptr) || (len < 16)) return nullptr;

  // zero buffer and reserve last character in buffer for terminating '\0'
  memset(buf, 0, len);
  len--;

#if defined(__linux__)

  int fd = fileno(fp);
  // get pathname from /proc or copy (unknown)
  if (readlink((std::string("/proc/self/fd/") + std::to_string(fd)).c_str(), buf, len) <= 0)
    strncpy(buf, "(unknown)", len);

#elif defined(__APPLE__)

  int fd = fileno(fp);
  char filepath[PATH_MAX];
  if (fcntl(fd, F_GETPATH, filepath) != -1)
    strncpy(buf, filepath, len);
  else
    strncpy(buf, "(unknown)", len);

#elif defined(_WIN32)

  char filepath[MAX_PATH];
  HANDLE h = (HANDLE) _get_osfhandle(_fileno(fp));
  if (GetFinalPathNameByHandleA(h, filepath, MAX_PATH, FILE_NAME_NORMALIZED) > 0)
    strncpy(buf, filepath, len);
  else
    strncpy(buf, "(unknown)", len);

#else    // unsupported OS

  strncpy(buf, "(unknown)", len);

#endif

  return buf;
}

/* ----------------------------------------------------------------------
   detect terminal, e.g. for using a pager automatically
------------------------------------------------------------------------- */

bool platform::is_console(FILE *fp)
{
  if (!fp) return false;
#if defined(_WIN32)
  return (_isatty(_fileno(fp)) == 1);
#else
  return (isatty(fileno(fp)) == 1);
#endif
}

/* ----------------------------------------------------------------------
   Get string with path to the current directory
   PATH_MAX may not be a compile time constant, so we must allocate and delete a buffer.
------------------------------------------------------------------------- */

std::string platform::current_directory()
{
  std::string cwd;

#if defined(_WIN32)
  char *buf = new char[MAX_PATH];
  if (_getcwd(buf, MAX_PATH)) { cwd = buf; }
  delete[] buf;
#else
  auto buf = new char[PATH_MAX];
  if (::getcwd(buf, PATH_MAX)) { cwd = buf; }
  delete[] buf;
#endif
  return cwd;
}

/* ----------------------------------------------------------------------
   check if a path is a directory
------------------------------------------------------------------------- */

bool platform::path_is_directory(const std::string &path)
{
#if defined(_WIN32)
  struct _stat info;
  memset(&info, 0, sizeof(info));
  if (_stat(path.c_str(), &info) != 0) return false;
#else
  struct stat info;
  memset(&info, 0, sizeof(info));
  if (stat(path.c_str(), &info) != 0) return false;
#endif
  return ((info.st_mode & S_IFDIR) != 0);
}

/* ----------------------------------------------------------------------
   get directory listing in string vector
------------------------------------------------------------------------- */

std::vector<std::string> platform::list_directory(const std::string &dir)
{
  std::vector<std::string> files;
  if (!path_is_directory(dir)) return files;

#if defined(_WIN32)
  HANDLE handle;
  WIN32_FIND_DATA fd;
  std::string searchname = dir + filepathsep[0] + "*";
  handle = FindFirstFile(searchname.c_str(), &fd);
  if (handle == ((HANDLE) -1)) return files;
  while (FindNextFile(handle, &fd)) {
    std::string entry(fd.cFileName);
    if ((entry == "..") || (entry == ".")) continue;
    files.push_back(entry);
  }
  FindClose(handle);
#else
  std::string dirname = dir + filepathsep[0];
  DIR *handle = opendir(dirname.c_str());
  if (handle == nullptr) return files;
  struct dirent *fd;
  while ((fd = readdir(handle)) != nullptr) {
    std::string entry(fd->d_name);
    if ((entry == "..") || (entry == ".")) continue;
    files.push_back(entry);
  }
  closedir(handle);
#endif
  return files;
}

/* ----------------------------------------------------------------------
   Change current directory
------------------------------------------------------------------------- */

int platform::chdir(const std::string &path)
{
#if defined(_WIN32)
  return ::_chdir(path.c_str());
#else
  return ::chdir(path.c_str());
#endif
}

/* ----------------------------------------------------------------------
   Create a directory. Create entire path if necessary.
------------------------------------------------------------------------- */

int platform::mkdir(const std::string &path)
{
  std::deque<std::string> dirlist = {path};
  std::string dirname = path_dirname(path);

  while ((dirname != ".") && (dirname != "")) {
    dirlist.push_front(dirname);
    dirname = path_dirname(dirname);
  }

  int rv;
  for (const auto &dir : dirlist) {
    if (!path_is_directory(dir)) {
#if defined(_WIN32)
      rv = ::_mkdir(dir.c_str());
#else
      rv = ::mkdir(dir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP);
#endif
      if (rv != 0) return rv;
    }
  }
  return 0;
}

/* ----------------------------------------------------------------------
   Delete a directory and its contents recursively
------------------------------------------------------------------------- */

int platform::rmdir(const std::string &path)
{
  // recurse through directory tree deleting files and directories
  auto entries = list_directory(path);
  for (const auto &entry : entries) {
    const auto newpath = path_join(path, entry);
    if (path_is_directory(newpath))
      rmdir(newpath);
    else
      unlink(newpath);
  }
#if defined(_WIN32)
  return ::_rmdir(path.c_str());
#else
  return ::rmdir(path.c_str());
#endif
}

/* ----------------------------------------------------------------------
   Delete a file
------------------------------------------------------------------------- */

int platform::unlink(const std::string &path)
{
#if defined(_WIN32)
  return ::_unlink(path.c_str());
#else
  return ::unlink(path.c_str());
#endif
}

/* ----------------------------------------------------------------------
   Get current file stream position
------------------------------------------------------------------------- */

bigint platform::ftell(FILE *fp)
{
#if defined(_WIN32)
  return (bigint)::_ftelli64(fp);
#else
  return (bigint)::ftell(fp);
#endif
}

/* ----------------------------------------------------------------------
   Set current file stream position
------------------------------------------------------------------------- */

int platform::fseek(FILE *fp, bigint pos)
{
#if defined(_WIN32)
  if (pos == platform::END_OF_FILE)
    return ::_fseeki64(fp, 0, SEEK_END);
  else
    return ::_fseeki64(fp, (__int64) pos, SEEK_SET);
#else
  if (pos == platform::END_OF_FILE)
    return ::fseek(fp, 0, SEEK_END);
  else
    return ::fseek(fp, (long) pos, SEEK_SET);
#endif
}

/* ----------------------------------------------------------------------
   Truncate opened file to given length
------------------------------------------------------------------------- */

int platform::ftruncate(FILE *fp, bigint length)
{
#if defined(_WIN32)
  HANDLE h = (HANDLE) _get_osfhandle(_fileno(fp));
  LARGE_INTEGER li_start, li_length;
  li_start.QuadPart = (int64_t) 0;
  li_length.QuadPart = (int64_t) length;
  if (SetFilePointerEx(h, li_start, NULL, FILE_CURRENT) &&
      SetFilePointerEx(h, li_length, NULL, FILE_BEGIN) && SetEndOfFile(h)) {
    return 0;
  } else {
    return 1;
  }
#else
  platform::fseek(fp, length);
  return ::ftruncate(fileno(fp), (off_t) length);
#endif
}

/* ----------------------------------------------------------------------
   open pipe
------------------------------------------------------------------------- */

FILE *platform::popen(const std::string &cmd, const std::string &mode)
{
  FILE *fp = nullptr;
#if defined(_WIN32)
  if (mode == "r")
    fp = ::_popen(cmd.c_str(), "rb");
  else if (mode == "w")
    fp = ::_popen(cmd.c_str(), "wb");
#else
  if (mode == "r")
    fp = ::popen(cmd.c_str(), "r");
  else if (mode == "w")
    fp = ::popen(cmd.c_str(), "w");
#endif
  return fp;
}

/* ----------------------------------------------------------------------
   close pipe
------------------------------------------------------------------------- */

int platform::pclose(FILE *fp)
{
#if defined(_WIN32)
  return ::_pclose(fp);
#else
  return ::pclose(fp);
#endif
}

/* ----------------------------------------------------------------------
   strip off leading part of path, return just the filename
------------------------------------------------------------------------- */

std::string platform::path_basename(const std::string &path)
{
  size_t start = path.find_last_of(platform::filepathsep);

  if (start == std::string::npos) {
    start = 0;
  } else {
    start += 1;
  }

  return path.substr(start);
}

/* ----------------------------------------------------------------------
   Return only the leading part of a path, return just the directory
------------------------------------------------------------------------- */

std::string platform::path_dirname(const std::string &path)
{
  size_t start = path.find_last_of(platform::filepathsep);

  if (start == std::string::npos) return ".";

  return path.substr(0, start);
}

/* ----------------------------------------------------------------------
   join two paths.
   if one of the two is an empty string just return the other unmodified
   if the first string ends in the separator or the second begins with one, trim them
------------------------------------------------------------------------- */

std::string platform::path_join(const std::string &a, const std::string &b)
{
  if (a.empty()) return b;
  if (b.empty()) return a;

  // remove trailing separator(s) in first part
  std::string joined = a;
  while (joined.find_last_of(platform::filepathsep) == joined.size() - 1) {
    for (const auto &s : platform::filepathsep)
      if (joined.back() == s) joined.pop_back();
  }

  // skip over leading separator(s) in second part
  std::size_t skip = 0;
  while (b.find_first_of(platform::filepathsep, skip) == skip) ++skip;

  // combine and return
  joined += platform::filepathsep[0] + b.substr(skip);
  return joined;
}

/* ----------------------------------------------------------------------
   try to open file for reading to prove if it exists and is accessible
------------------------------------------------------------------------- */

bool platform::file_is_readable(const std::string &path)
{
  FILE *fp = fopen(path.c_str(), "r");
  if (fp) {
    fclose(fp);
    return true;
  }
  return false;
}
/* ----------------------------------------------------------------------
   determine available disk space, if supported. Return -1 if not.
------------------------------------------------------------------------- */

double platform::disk_free(const std::string &path)
{
  double bytes_free = -1.0;

#if defined(__linux__) || defined(__APPLE__) || defined(__FreeBSD__) || defined(__DragonFly__) || \
    defined(__OpenBSD__) || defined(__NetBSD__)
  struct statvfs fs;

  if (path.size()) {
    int rv = statvfs(path.c_str(), &fs);
    if (rv == 0) {
#if defined(__linux__)
      bytes_free = fs.f_bavail * fs.f_bsize;
#elif defined(__APPLE__) || defined(__FreeBSD__) || defined(__DragonFly__) || \
    defined(__OpenBSD__) || defined(__NetBSD__)
      bytes_free = fs.f_bavail * fs.f_frsize;
#endif
    }
  }
#elif defined(_WIN32)
  uint64_t is_free = 0;
  if (GetDiskFreeSpaceEx(path.c_str(), (PULARGE_INTEGER) &is_free, nullptr, nullptr))
    bytes_free = is_free;
#endif
  return bytes_free;
}

/* ----------------------------------------------------------------------
   check if filename has a known compression extension
------------------------------------------------------------------------- */

bool platform::has_compress_extension(const std::string &file)
{
  return find_compress_type(file).style != ::compress_info::NONE;
}

/* ----------------------------------------------------------------------
   open pipe to read a compressed file
------------------------------------------------------------------------- */

FILE *platform::compressed_read(const std::string &file)
{
  FILE *fp = nullptr;

#if defined(LAMMPS_GZIP)
  const auto &compress = find_compress_type(file);
  if (compress.style == ::compress_info::NONE) return nullptr;

  if (find_exe_path(compress.command).size())
    // put quotes around file name so that they may contain blanks
    fp = popen((compress.command + compress.uncompressflags + "\"" + file + "\""), "r");
#endif
  return fp;
}

/* ----------------------------------------------------------------------
   open pipe to write a compressed file
------------------------------------------------------------------------- */

FILE *platform::compressed_write(const std::string &file)
{
  FILE *fp = nullptr;

#if defined(LAMMPS_GZIP)
  const auto &compress = find_compress_type(file);
  if (compress.style == ::compress_info::NONE) return nullptr;

  if (find_exe_path(compress.command).size())
    // put quotes around file name so that they may contain blanks
    fp = popen((compress.command + compress.compressflags + "\"" + file + "\""), "w");
#endif
  return fp;
}

/* ---------------------------------------------------------------------- */
