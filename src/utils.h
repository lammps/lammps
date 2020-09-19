/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_UTILS_H
#define LMP_UTILS_H

/*! \file utils.h */

#include "lmptype.h"
#include <string>
#include <vector>
#include <cstdio>

namespace LAMMPS_NS {

  // forward declarations
  class Error;
  class LAMMPS;

  namespace utils {

    /** Match text against a simplified regex pattern
     *
     *  \param text the text to be matched against the pattern
     *  \param pattern the search pattern, which may contain regexp markers
     *  \return true if the pattern matches, false if not */

    bool strmatch(const std::string &text, const std::string &pattern);

    /** Send message to screen and logfile, if available
     *
     *  \param lmp   pointer to LAMMPS class instance
     *  \param mesg  message to be printed */

    void logmesg(LAMMPS *lmp, const std::string &mesg);

    /** return a string representing the current system error status
     *
     *  This is a wrapper around calling strerror(errno).
     *
     *  \return  error string */

    std::string getsyserror();

    /** safe wrapper around fgets() which aborts on errors
     *  or EOF and prints a suitable error message to help debugging
     *
     *  \param srcname  name of the calling source file (from FLERR macro)
     *  \param srcline  line in the calling source file (from FLERR macro)
     *  \param s        buffer for storing the result of fgets()
     *  \param size     size of buffer s (max number of bytes read by fgets())
     *  \param fp       file pointer used by fgets()
     *  \param filename file name associated with fp (may be a null pointer; then LAMMPS will try to detect)
     *  \param error    pointer to Error class instance (for abort) */

    void sfgets(const char *srcname, int srcline, char *s, int size,
                FILE *fp, const char *filename, Error *error);

    /** safe wrapper around fread() which aborts on errors
     *  or EOF and prints a suitable error message to help debugging
     *
     *  \param srcname  name of the calling source file (from FLERR macro)
     *  \param srcline  line in the calling source file (from FLERR macro)
     *  \param s        buffer for storing the result of fread()
     *  \param size     size of data elements read by fread()
     *  \param num      number of data elements read by fread()
     *  \param fp       file pointer used by fread()
     *  \param filename file name associated with fp (may be a null pointer; then LAMMPS will try to detect)
     *  \param error    pointer to Error class instance (for abort) */

    void sfread(const char *srcname, int srcline, void *s, size_t size,
                size_t num, FILE *fp, const char *filename, Error *error);

    /** Report if a requested style is in a package or may have a typo
     *
     *  \param style type of style that is to be checked for
     *  \param name  name of style that was not found
     *  \param lmp   pointer to top-level LAMMPS class instance
     *  \return string usable for error messages */

    std::string check_packages_for_style(const std::string &style,
                                         const std::string &name, LAMMPS *lmp);

    /** Convert a string to a floating point number while checking
     *  if it is a valid floating point or integer number
     *
     *  \param file     name of source file for error message
     *  \param line     line number in source file for error message
     *  \param str      string to be converted to number
     *  \param do_abort determines whether to call Error::one() or Error::all()
     *  \param lmp      pointer to top-level LAMMPS class instance
     *  \return         double precision floating point number
     */
    double numeric(const char *file, int line, const char *str,
                   bool do_abort, LAMMPS *lmp);

    /** Convert a string to an integer number while checking
     *  if it is a valid integer number (regular int)
     *
     *  \param file     name of source file for error message
     *  \param line     line number in source file for error message
     *  \param str      string to be converted to number
     *  \param do_abort determines whether to call Error::one() or Error::all()
     *  \param lmp      pointer to top-level LAMMPS class instance
     *  \return         integer number (regular int)  */

    int inumeric(const char *file, int line, const char *str,
                 bool do_abort, LAMMPS *lmp);

    /** Convert a string to an integer number while checking
     *  if it is a valid integer number (bigint)
     *
     *  \param file     name of source file for error message
     *  \param line     line number in source file for error message
     *  \param str      string to be converted to number
     *  \param do_abort determines whether to call Error::one() or Error::all()
     *  \param lmp      pointer to top-level LAMMPS class instance
     *  \return         integer number (bigint) */

    bigint bnumeric(const char *file, int line, const char *str,
                    bool do_abort, LAMMPS *lmp);

    /** Convert a string to an integer number while checking
     *  if it is a valid integer number (tagint)
     *
     * \param file     name of source file for error message
     * \param line     line number in source file for error message
     * \param str      string to be converted to number
     * \param do_abort determines whether to call Error::one() or Error::all()
     * \param lmp      pointer to top-level LAMMPS class instance
     * \return         integer number (tagint) */

    tagint tnumeric(const char *file, int line, const char *str,
                    bool do_abort, LAMMPS *lmp);

    /** Compute index bounds derived from a string with a possible wildcard
     *
     * This functions processes the string in *str* and set the values of *nlo*
     * and *nhi* according to the following five cases:
     *
     * - a single number, i: nlo = i; nhi = i;
     * - a single asterisk, \*: nlo = nmin; nhi = nmax;
     * - a single number followed by an asterisk, i\*: nlo = i; nhi = nmax;
     * - a single asterisk followed by a number, \*i: nlo = nmin; nhi = i;
     * - two numbers with an asterisk in between. i\*j: nlo = i; nhi = j;
     *
     * \param file     name of source file for error message
     * \param line     line number in source file for error message
     * \param str      string to be processed
     * \param nmin     smallest possible lower bound
     * \param nmax     largest allowed upper bound
     * \param nlo      lower bound
     * \param nhi      upper bound
     * \param error    pointer to Error class for out-of-bounds messages */

    template <typename TYPE>
    void bounds(const char *file, int line, const std::string &str,
                bigint nmin, bigint nmax, TYPE &nlo, TYPE &nhi, Error *error);

    /** Expand list of arguments when containing fix/compute wildcards
     *
     *  This function searches the list of arguments in *arg* for strings
     *  of the kind c_ID[*] or f_ID[*] referring to computes or fixes.
     *  Any such strings are replaced by one or more strings with the
     *  '*' character replaced by the corresponding possible numbers as
     *  determined from the fix or compute instance.  Other strings are
     *  just copied. If the *mode* parameter is set to 0, expand global
     *  vectors, but not global arrays; if it is set to 1, expand global
     *  arrays (by column) but not global vectors.
     *
     *  If any expansion happens, the earg list and all its
     *  strings are new allocations and must be freed explicitly by the
     *  caller. Otherwise arg and earg will point to the same address
     *  and no explicit de-allocation is needed by the caller.
     *
     * \param file  name of source file for error message
     * \param line  line number in source file for error message
     * \param narg  number of arguments in current list
     * \param arg   argument list, possibly containing wildcards
     * \param mode  select between global vectors(=0) and arrays (=1)
     * \param earg  new argument list with wildcards expanded
     * \param lmp   pointer to top-level LAMMPS class instance
     * \return      number of arguments in expanded list */

    int expand_args(const char *file, int line, int narg, char **arg,
                    int mode, char **&earg, LAMMPS *lmp);

    /** Trim leading and trailing whitespace. Like TRIM() in Fortran.
     *
     * \param line string that should be trimmed
     * \return new string without whitespace (string) */

    std::string trim(const std::string &line);

    /** Return string with anything from '#' onward removed
     *
     * \param line string that should be trimmed
     * \return new string without comment (string) */

    std::string trim_comment(const std::string &line);

    /** Count words in string with custom choice of separating characters
     *
     * \param text string that should be searched
     * \param separators string containing characters that will be treated as whitespace
     * \return number of words found */

    size_t count_words(const std::string &text, const std::string &separators);

    /** Count words in string, ignore any whitespace matching " \t\r\n\f"
     *
     * \param text string that should be searched
     * \return number of words found */

    size_t count_words(const std::string &text);

    /** Count words in C-string, ignore any whitespace matching " \t\r\n\f"
     *
     * \param text string that should be searched
     * \return number of words found */

    size_t count_words(const char *text);

    /** Count words in a single line, trim anything from '#' onward
     *
     * \param text string that should be trimmed and searched
     * \param separators string containing characters that will be treated as whitespace
     * \return number of words found */

    size_t trim_and_count_words(const std::string &text, const std::string &separators = " \t\r\n\f");

    /** Take text and split into non-whitespace words.
     *
     * This can handle strings with single and double quotes, escaped quotes,
     * and escaped codes within quotes, but due to using an STL container and
     * STL strings is rather slow because of making copies. Designed for parsing
     * command lines and similar text and not for time critical processing.
     * Use a tokenizer class for that.
     *
\verbatim embed:rst

*See also*
   :cpp:class:`Tokenizer`, :cpp:class:`ValueTokenizer`

\endverbatim
     *
     * \param text string that should be split
     * \return STL vector with the words */

    std::vector<std::string> split_words(const std::string &text);

    /** Check if string can be converted to valid integer
     *
     * \param str string that should be checked
     * \return true, if string contains valid a integer, false otherwise */

    bool is_integer(const std::string &str);

    /** Check if string can be converted to valid floating-point number
     *
     * \param str string that should be checked
     * \return true, if string contains valid number, false otherwise */

    bool is_double(const std::string &str);

    /** Try to detect pathname from FILE pointer.
     *
     * Currently only supported on Linux, otherwise will report "(unknown)".
     *
     *  \param buf  storage buffer for pathname. output will be truncated if not large enough
     *  \param len  size of storage buffer. output will be truncated to this length - 1
     *  \param fp   FILE pointer struct from STDIO library for which we want to detect the name
     *  \return pointer to the storage buffer, i.e. buf */

    const char *guesspath(char *buf, int len, FILE *fp);

    /** Strip off leading part of path, return just the filename
     *
     * \param path file path
     * \return file name */

    std::string path_basename(const std::string &path);

    /** Join two pathname segments
     *
     * This uses the forward slash '/' character unless LAMMPS is compiled
     * for Windows where it used the equivalent backward slash '\\'.
     *
     * \param   a  first path
     * \param   b  second path
     * \return     combined path */

    std::string path_join(const std::string &a, const std::string &b);

    /** Check if file exists and is readable
     *
     * \param path file path
     * \return true if file exists and is readable */

    bool file_is_readable(const std::string &path);

    /** Determine full path of potential file. If file is not found in current directory,
     *  search directories listed in LAMMPS_POTENTIALS environment variable
     *
     * \param path file path
     * \return full path to potential file */

    std::string get_potential_file_path(const std::string &path);

    /** Read potential file and return DATE field if it is present
     *
     * \param path file path
     * \param potential_name name of potential that is being read
     * \return DATE field if present */

    std::string get_potential_date(const std::string &path,
                                   const std::string &potential_name);

    /** Read potential file and return UNITS field if it is present
     *
     * \param path file path
     * \param potential_name name of potential that is being read
     * \return UNITS field if present */

    std::string get_potential_units(const std::string &path,
                                    const std::string &potential_name);

    enum { NOCONVERT = 0, METAL2REAL = 1, REAL2METAL = 1<<1 };
    enum { UNKNOWN = 0, ENERGY };

    /** Return bitmask of available conversion factors for a given property
     *
     * \param property property to be converted
     * \return bitmask indicating available conversions */
    int get_supported_conversions(const int property);

    /** Return unit conversion factor for given property and selected from/to units
     *
     * \param property property to be converted
     * \param conversion constant indicating the conversion
     * \return conversion factor */

    double get_conversion_factor(const int property, const int conversion);

    /** Open a potential file as specified by *name*
     *
     * If opening the file directly fails, the function will search for
     * it in the list of folder pointed to by the environment variable
     * ``LAMMPS_POTENTIALS`` (if it is set).
     *
     * If the potential file has a ``UNITS`` tag in the first line, the
     * tag's value is compared to the current unit style setting.
     * The behavior of the function then depends on the value of the
     * *auto_convert* parameter.  If it is a null pointer, then the unit
     * values must match or else the open will fail with an error.  Otherwise
     * the bitmask that *auto_convert* points to is used check for
     * compatibility with possible automatic conversions by the calling
     * function.  If compatible, the bitmask is set to the required
     * conversion or ``utils::NOCONVERT``.
     *
     * \param name          file- or pathname of the potential file
     * \param lmp           pointer to top-level LAMMPS class instance
     * \param auto_convert  pointer to unit conversion bitmask or ``nullptr``
     * \return              FILE pointer of the opened potential file or ``nullptr`` */

    FILE *open_potential(const std::string &name, LAMMPS *lmp, int *auto_convert);

    /** Convert a time string to seconds
     *
     * The strings "off" and "unlimited" result in -1
     *
     * \param timespec a string in the following format: ([[HH:]MM:]SS)
     * \return total in seconds */

    double timespec2seconds(const std::string &timespec);

    /** Convert a LAMMPS version date to a number
     *
     * This will generate a number YYYYMMDD from a date string
     * (with or without blanks) that is suitable for numerical
     * comparisons, i.e. later dates will generate a larger number.
     *
     * The day may or may not have a leading zero, the month
     * is identified by the first 3 letters (so there may be more)
     * and the year may be 2 or 4 digits (the missing 2 digits will
     * be assumed as 20. That is 04 corresponds to 2004).
     *
     * No check is made whether the date is valid.
     *
     * \param  date  string in the format (Day Month Year)
     * \return       date code */
    int date2num(const std::string &date);

    /** Custom merge sort implementation
     *
     * This function provides a custom upward hybrid merge sort
     * implementation with support to pass an opaque pointer to
     * the comparison function, e.g. for access to class members.
     * This avoids having to use global variables.  For improved
     * performance, it uses an in-place insertion sort on initial
     * chunks of up to 64 elements and switches to merge sort from
     * then on.
     *
     * \param  index  Array with indices to be sorted
     * \param  num    Length of the index array
     * \param  ptr    Pointer to opaque object passed to comparison function
     * \param  comp   Pointer to comparison function */

    void merge_sort(int *index, int num, void *ptr,
                    int (*comp)(int, int, void *));
  }
}

#endif

/* ERROR/WARNING messages:

*/
