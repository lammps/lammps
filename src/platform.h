/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_PLATFORM_H
#define LMP_PLATFORM_H

/*! \file platform.h */

#include "lmptype.h"

#include <cstdio>
#include <string>
#include <vector>

namespace LAMMPS_NS {
namespace platform {

  /*! Return the consumed CPU time for the current process in seconds
   *
   * This is a wrapper around the POSIX function getrusage() and its Windows equivalent.
   * It is to be used in a similar fashion than MPI_Wtime().  Its resolution may be rather
   * low so it can only be trusted when observing processes consuming CPU time of at least
   * a few seconds.
   *
   *  \return used CPU time in seconds */

  double cputime();

  /*! Return the wall clock state for the current process in seconds
   *
   * This this clock is counting continuous time and is initialized during
   * load of the executable/library.  Its absolute value must be considered
   * arbitrary and thus elapsed wall times are measured in taking differences.
   * It is therefore to be used in a similar fashion as MPI_Wtime() but
   * has a different offset, usually leading to better resolution.
   *
   *  \return wall clock time in seconds */

  double walltime();

  /*! Suspend execution for a microsecond interval
   *
   * This emulates the usleep(3) BSD function call also mentioned in POSIX.1-2001.
   * This is not a precise delay; it may be longer, but not shorter.
   *
   * \param  usec  length of delay in microseconds */

  void usleep(int usec);

  /*! Return string with the operating system version and architecture info
   *
   *  \return  string with info about the OS and the platform is is running on */

  std::string os_info();

  /*! Return string with C++ standard version used to compile LAMMPS.
   *
   * This function uses predefined compiler macros to identify
   * the C++ standard version used to compile LAMMPS with.
   *
   *  \return  string with the C++ standard version or "unknown" */

  std::string cxx_standard();

  /*! Return string with compiler version info
   *
   * This function uses predefined compiler macros to identify
   * Compilers and their version and configuration info.
   *
   *  \return  string with the compiler information text */

  std::string compiler_info();

  /*! Return string with OpenMP standard version info
   *
   * This function uses predefined compiler macros to identify
   * OpenMP support and the supported version of the standard.
   *
   *  \return  string with the openmp information text */

  std::string openmp_standard();

  /*! Return string with MPI vendor info
   *
   * This function uses predefined macros to identify
   * the vendor of the MPI library used.
   *
   *  \return  string with the MPI vendor information text */

  std::string mpi_vendor();

  /*! Return string with MPI version info
   *
   * This function uses predefined macros and MPI function
   * calls to identify the version of the MPI library used.
   *
   *  \param  major  major version of the MPI standard (set on exit)
   *  \param  minor  minor version of the MPI standard (set on exit)
   *  \return  string with the MPI version information text */

  std::string mpi_info(int &major, int &minor);

  /*! Return string with list of available compression types and executables
   *
   * This function tests which of the supported compression executables
   * are available for reading or writing compressed files where supported.
   *
   *  \return  string with list of available compression tools */

  std::string compress_info();

  /*! Add variable to the environment
   *
   * \param  vardef  variable name or variable definition (NAME=value)
   * \return -1 if failure otherwise 0 */

  int putenv(const std::string &vardef);

  /*! Delete variable from the environment
   *
   * \param  variable  variable name
   * \return -1 if failure otherwise 0 */

  int unsetenv(const std::string &variable);

  /*! Get list of entries in a path environment variable
   *
   * This provides a list of strings of the entries in an environment
   * variable that is containing a "path" like "PATH" or "LD_LIBRARY_PATH".
   *
   * \param   var  name of the environment variable
   * \return  vector with strings of all entries in that path variable */

  std::vector<std::string> list_pathenv(const std::string &var);

  /*! Open a shared object file or library
   *
   * \param  fname   name or path of the shared object
   * \return  handle to the shared object or null */

  void *dlopen(const std::string &fname);

  /*! Obtain error diagnostic info after dynamic linking function calls
   *
   * Return a human-readable string describing the most recent error that
   * occurred when using one of the functions for dynamic loading objects
   * the last call to this function. The string is empty, if there was no error.
   *
   * \return  string with error message or empty */

  std::string dlerror();

  /*! Close a shared object
   *
   * This releases the object corresponding to the provided handle.
   * Resolved symbols associated with this handle may not be used
   * after this call
   *
   * \param  handle  handle to an opened shared object
   * \return 0 if successful, non-zero of not */

  int dlclose(void *handle);

  /*! Resolve a symbol in shared object
   *
   * \param  handle  handle to an opened shared object
   * \param  symbol  name of the symbol to extract
   * \return  pointer to the resolved symbol or null */

  void *dlsym(void *handle, const std::string &symbol);

  /*! Platform specific file path component separator
   *
   * This is a string with the character that separates directories and filename in paths on
   * a platform. If multiple are characters are provided, the first is the preferred one. */

#if !defined(_WIN32)
  constexpr char filepathsep[] = "/";
#else
  constexpr char filepathsep[] = "\\/";
#endif

  /*! Platform specific path environment variable component separator
   *
   * This is the character that separates entries in "PATH"-style environment variables. */

#if !defined(_WIN32)
  constexpr char pathvarsep = ':';
#else
  constexpr char pathvarsep = ';';
#endif

  /*! Try to detect pathname from FILE pointer
   *
   * Currently only supported on Linux, MacOS, and Windows. Otherwise will report "(unknown)".
   *
   *  \param  fp   FILE pointer struct from STDIO library for which we want to detect the name
   *  \param  buf  storage buffer for pathname. output will be truncated if not large enough
   *  \param  len  size of storage buffer. output will be truncated to this length - 1
   *  \return  pointer to the storage buffer with path or a NULL pointer if buf is invalid
   *  or the buffer size is too small */

  const char *guesspath(FILE *fp, char *buf, int len);

  /*! Check if a file pointer may be connected to a console
   *
   * \param  fp file pointer
   * \return true if the file pointer is flagged as a TTY */

  bool is_console(FILE *fp);

  /*! Get string with path to the current directory
   *
   * \return path to the current directory or empty string */

  std::string current_directory();

  /*! Check if a path is a directory
   *
   * \param  path  directory path
   * \return true if the directory exists */

  bool path_is_directory(const std::string &path);

  /*! Get list of entries in a directory
   *
   * This provides a list of strings of the entries in the directory
   * without the leading path name while also skipping over ".." and ".".
   *
   * \param   dir  path to directory
   * \return  vector with strings of all directory entries */

  std::vector<std::string> list_directory(const std::string &dir);

  /*! Find pathname of an executable in the standard search path
   *
   * This function will traverse the list of directories in the PATH
   * environment variable and look for the executable *cmd*.  If the
   * file exists and is executable the full path is returned as string,
   * otherwise and empty string is returned.
   *
   * On Windows the *cmd* string must not include and extension as
   * this function will automatically append the extensions ".exe",
   * ".com" and ".bat" and look for those paths. On Windows also the
   * current directory is checked (and first), while otherwise not unless
   * "." exists in the PATH environment variable.
   *
   * Because of the nature of the check, this will not detect shell functions
   * built-in command or aliases.
   *
   * \param   cmd  name of command
   * \return  vector with strings of all directory entries */

  std::string find_exe_path(const std::string &cmd);

  /*! Change current directory
   *
   * \param  path  new current working directory path
   * \return -1 if unsuccessful, otherwise >= 0  */

  int chdir(const std::string &path);

  /*! Create a directory
   *
   * Unlike the the ``mkdir()`` or ``_mkdir()`` functions of the
   * C library, this function will also try to create non-existing sub-directories
   * in case they don't exist, and thus behave like the ``mkdir -p`` command rather
   * than plain ``mkdir`` or ``md`.
   *
   * \param  path  directory path
   * \return -1 if unsuccessful, otherwise >= 0  */

  int mkdir(const std::string &path);

  /*! Delete a directory
   *
   * Unlike the the ``rmdir()`` or ``_rmdir()`` functions of the
   * C library, this function will check for the contents of the
   * folder and recurse into any sub-folders, if necessary and
   * delete all contained folders and their contents before
   * deleting the folder *path*.
   *
   * \param  path  directory path
   * \return -1 if unsuccessful, otherwise >= 0  */

  int rmdir(const std::string &path);

  /*! Delete a file
   *
   *  \param   path    path to file to be deleted
   *  \return  0 on success, -1 on error */

  int unlink(const std::string &path);

  /*! Get current file position
   *
   *  \param   fp      FILE pointer of the given file
   *  \return  current FILE pointer position cast to a bigint */

  bigint ftell(FILE *fp);

  /*! constant to seek to the end of the file */
  constexpr bigint END_OF_FILE = -1;

  /*! Set absolute file position
   *
   * If the absolute position is END_OF_FILE, then position at the end of the file.
   *
   *  \param   fp      FILE pointer of the given file
   *  \param   pos     new position of the FILE pointer
   *  \return  0 if successful, otherwise -1 */

  int fseek(FILE *fp, bigint pos);

  /*! Truncate file to a given length and re-position file pointer
   *
   *  \param   fp      FILE pointer of the given file
   *  \param   length  length to which the file is being truncated to
   *  \return  0 if successful, otherwise -1 */

  int ftruncate(FILE *fp, bigint length);

  /*! Open a pipe to a command for reading or writing
   *
   *  \param  cmd   command for the pipe
   *  \param  mode  "r" for reading from *cmd* or "w" for writing to *cmd*
   *  \return  file pointer to the pipe if successful or null */

  FILE *popen(const std::string &cmd, const std::string &mode);

  /*! Close a previously opened pipe
   *
   *  \param  fp   FILE pointer for the pipe
   *  \return  exit status of the pipe command or -1 in case of errors */

  int pclose(FILE *fp);

  /*! Strip off leading part of path, return just the filename
   *
   * \param path file path
   * \return file name */

  std::string path_basename(const std::string &path);

  /*! Return the directory part of a path. Return "." if empty
   *
   * \param path file path
   * \return directory name */

  std::string path_dirname(const std::string &path);

  /*! Join two pathname segments
   *
   * This uses the forward slash '/' character unless LAMMPS is compiled
   * for Windows where it uses the backward slash '\\'
   *
   * \param   a  first path
   * \param   b  second path
   * \return     combined path */

  std::string path_join(const std::string &a, const std::string &b);

  /*! Check if file exists and is readable
   *
   * \param path file path
   * \return true if file exists and is readable */

  bool file_is_readable(const std::string &path);

  /*! Check if a file name ends in a known extension for a compressed file format
   *
   * Currently supported file extensions are: .gz, .bz2, .zst, .xz, .lzma, lz4
   *
   *  \param  file  name of the file to check
   *  \return  true if the file has a known extension, otherwise false  */

  bool has_compress_extension(const std::string &file);

  /*! Open pipe to compressed text file for reading
   *
   *  \param  file  name of the file to open
   *  \return  FILE pointer to pipe using for reading the compressed file. */

  FILE *compressed_read(const std::string &file);

  /*! Open pipe to compressed text file for writing
   *
   *  \param  file  name of the file to open
   *  \return  FILE pointer to pipe using for reading the compressed file. */

  FILE *compressed_write(const std::string &file);

}    // namespace platform
}    // namespace LAMMPS_NS
#endif
