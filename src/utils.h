/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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

#include "fmt/format.h"
#include "lmptype.h"

#include <mpi.h>

#include <vector>    // IWYU pragma: export

namespace LAMMPS_NS {

// forward declarations
class Error;
class LAMMPS;

namespace utils {

  /*! Match text against a simplified regex pattern
   *
   *  \param text the text to be matched against the pattern
   *  \param pattern the search pattern, which may contain regexp markers
   *  \return true if the pattern matches, false if not */

  bool strmatch(const std::string &text, const std::string &pattern);

  /*! Find sub-string that matches a simplified regex pattern
   *
   *  \param text the text to be matched against the pattern
   *  \param pattern the search pattern, which may contain regexp markers
   *  \return the string that matches the pattern or an empty one */

  std::string strfind(const std::string &text, const std::string &pattern);

  /*! Print error message about missing arguments for command
   *
   * This function simplifies the repetitive reporting missing arguments to a command.
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param cmd      name of the failing command
   *  \param error    pointer to Error class instance (for abort) or nullptr */

  void missing_cmd_args(const std::string &file, int line, const std::string &cmd, Error *error);

  /* Internal function handling the argument list for logmesg(). */

  void fmtargs_logmesg(LAMMPS *lmp, fmt::string_view format, fmt::format_args args);

  /*! Send formatted message to screen and logfile, if available
   *
   * This function simplifies the repetitive task of outputting some
   * message to both the screen and/or the log file. The template
   * wrapper with fmtlib format and argument processing allows
   * this function to work similar to ``fmt::print()``.
   *
   *  \param lmp    pointer to LAMMPS class instance
   *  \param format format string of message to be printed
   *  \param args   arguments to format string */

  template <typename... Args> void logmesg(LAMMPS *lmp, const std::string &format, Args &&...args)
  {
    fmtargs_logmesg(lmp, format, fmt::make_format_args(args...));
  }

  /*! \overload
   *
   *  \param lmp    pointer to LAMMPS class instance
   *  \param mesg   string with message to be printed */

  void logmesg(LAMMPS *lmp, const std::string &mesg);

  /*! Return text redirecting the user to a specific paragraph in the manual
   *
   * The LAMMPS manual contains detailed detailed explanations for errors and
   * warnings where a simple error message may not be sufficient.  These can
   * be reached through URLs with a numeric code.  This function creates the
   * corresponding text to be included into the error message that redirects
   * the user to that URL.
   *
   *  \param errorcode   number pointing to a paragraph in the manual */

  std::string errorurl(int errorcode);

  /*! Flush output buffers
   *
   *  This function calls fflush() on screen and logfile FILE pointers
   *  if available and thus tells the operating system to output all
   *  currently buffered data. This is local operation and independent
   *  from buffering by a file system or an MPI library.
   */

  void flush_buffers(LAMMPS *lmp);

  /*! Return a string representing the current system error status
   *
   *  This is a wrapper around calling strerror(errno).
   *
   *  \return  error string */

  std::string getsyserror();

  /*! Wrapper around fgets() which reads whole lines but truncates the
   *  data to the buffer size and ensures a newline char at the end.
   *
   *  This function is useful for reading line based text files with
   *  possible comments that should be parsed later. This applies to
   *  data files, potential files, atomfile variable files and so on.
   *  It is used instead of fgets() by utils::read_lines_from_file().
   *
   *  \param s        buffer for storing the result of fgets()
   *  \param size     size of buffer s (max number of bytes returned)
   *  \param fp       file pointer used by fgets() */

  char *fgets_trunc(char *s, int size, FILE *fp);

  /*! Safe wrapper around fgets() which aborts on errors
   *  or EOF and prints a suitable error message to help debugging.
   *
   *  Use nullptr as the error parameter to avoid the abort on EOF or error.
   *
   *  \param srcname  name of the calling source file (from FLERR macro)
   *  \param srcline  line in the calling source file (from FLERR macro)
   *  \param s        buffer for storing the result of fgets()
   *  \param size     size of buffer s (max number of bytes read by fgets())
   *  \param fp       file pointer used by fgets()
   *  \param filename file name associated with fp (may be a null pointer; then LAMMPS will try to detect)
   *  \param error    pointer to Error class instance (for abort) or nullptr */

  void sfgets(const char *srcname, int srcline, char *s, int size, FILE *fp, const char *filename,
              Error *error);

  /*! Safe wrapper around fread() which aborts on errors
   *  or EOF and prints a suitable error message to help debugging.
   *
   *  Use nullptr as the error parameter to avoid the abort on EOF or error.
   *
   *  \param srcname  name of the calling source file (from FLERR macro)
   *  \param srcline  line in the calling source file (from FLERR macro)
   *  \param s        buffer for storing the result of fread()
   *  \param size     size of data elements read by fread()
   *  \param num      number of data elements read by fread()
   *  \param fp       file pointer used by fread()
   *  \param filename file name associated with fp (may be a null pointer; then LAMMPS will try to detect)
   *  \param error    pointer to Error class instance (for abort) or nullptr */

  void sfread(const char *srcname, int srcline, void *s, size_t size, size_t num, FILE *fp,
              const char *filename, Error *error);

  /*! Read N lines of text from file into buffer and broadcast them
   *
   * This function uses repeated calls to fread() to fill a buffer with
   * newline terminated text.  If a line does not end in a newline (e.g.
   * at the end of a file), it is added.  The caller has to allocate an
   * nlines by nmax sized buffer for storing the text data.
   * Reading is done by MPI rank 0 of the given communicator only, and
   * thus only MPI rank 0 needs to provide a valid file pointer.
   *
   *  \param fp       file pointer used by fread
   *  \param nlines   number of lines to be read
   *  \param nmax     maximum length of a single line
   *  \param buffer   buffer for storing the data.
   *  \param me       MPI rank of calling process in MPI communicator
   *  \param comm     MPI communicator for broadcast
   *  \return         1 if the read was short, 0 if read was successful */

  int read_lines_from_file(FILE *fp, int nlines, int nmax, char *buffer, int me, MPI_Comm comm);

  /*! Report if a requested style is in a package or may have a typo
   *
   *  \param style type of style that is to be checked for
   *  \param name  name of style that was not found
   *  \param lmp   pointer to top-level LAMMPS class instance
   *  \return string usable for error messages */

  std::string check_packages_for_style(const std::string &style, const std::string &name,
                                       LAMMPS *lmp);

  /*! Convert a string to a boolean while checking whether it is a valid boolean term.
   *  Valid terms are 'yes', 'no', 'true', 'false', 'on', 'off', and '1', '0'. Only
   *  lower case is accepted.
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to logical
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param lmp      pointer to top-level LAMMPS class instance
   *  \return         1 if string resolves to "true", otherwise 0 */

  int logical(const char *file, int line, const std::string &str, bool do_abort, LAMMPS *lmp);

  /*! \overload
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to logical
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param lmp      pointer to top-level LAMMPS class instance
   *  \return         1 if string resolves to "true", otherwise 0 */

  int logical(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp);

  /*! Convert a string to a floating point number while checking
   *  if it is a valid floating point or integer number
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to number
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param lmp      pointer to top-level LAMMPS class instance
   *  \return         double precision floating point number */

  double numeric(const char *file, int line, const std::string &str, bool do_abort, LAMMPS *lmp);

  /*! \overload
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to number
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param lmp      pointer to top-level LAMMPS class instance
   *  \return         double precision floating point number */

  double numeric(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp);

  /*! Convert a string to an integer number while checking
   *  if it is a valid integer number (regular int)
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to number
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param lmp      pointer to top-level LAMMPS class instance
   *  \return         integer number (regular int)  */

  int inumeric(const char *file, int line, const std::string &str, bool do_abort, LAMMPS *lmp);

  /*! \overload
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to number
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param lmp      pointer to top-level LAMMPS class instance
   *  \return         double precision floating point number */

  int inumeric(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp);

  /*! Convert a string to an integer number while checking
   *  if it is a valid integer number (bigint)
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to number
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param lmp      pointer to top-level LAMMPS class instance
   *  \return         integer number (bigint) */

  bigint bnumeric(const char *file, int line, const std::string &str, bool do_abort, LAMMPS *lmp);

  /*! \overload
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to number
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param lmp      pointer to top-level LAMMPS class instance
   *  \return         double precision floating point number */

  bigint bnumeric(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp);

  /*! Convert a string to an integer number while checking
   *  if it is a valid integer number (tagint)
   *
   * \param file     name of source file for error message
   * \param line     line number in source file for error message
   * \param str      string to be converted to number
   * \param do_abort determines whether to call Error::one() or Error::all()
   * \param lmp      pointer to top-level LAMMPS class instance
   * \return         integer number (tagint) */

  tagint tnumeric(const char *file, int line, const std::string &str, bool do_abort, LAMMPS *lmp);

  /*! \overload
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to number
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param lmp      pointer to top-level LAMMPS class instance
   *  \return         double precision floating point number */

  tagint tnumeric(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp);

  /*! Compute index bounds derived from a string with a possible wildcard
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
  void bounds(const char *file, int line, const std::string &str, bigint nmin, bigint nmax,
              TYPE &nlo, TYPE &nhi, Error *error);

  /*! Expand list of arguments when containing fix/compute wildcards
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

  int expand_args(const char *file, int line, int narg, char **arg, int mode, char **&earg,
                  LAMMPS *lmp);

  /*! Expand type label string into its equivalent numeric type
   *
   *  This function checks if a given string may be a type label and
   *  then searches the labelmap type indicated by the *mode* argument
   *  for the corresponding numeric type.  If this is found a copy of
   *  the numeric type string is made and returned. Otherwise a null
   *  pointer is returned.
   *  If a string is returned, the calling code must free it with delete[].
   *
   * \param file  name of source file for error message
   * \param line  line number in source file for error message
   * \param str   type string to be expanded
   * \param mode  select labelmap using constants from Atom class
   * \param lmp   pointer to top-level LAMMPS class instance
   * \return      pointer to expanded string or null pointer */

  char *expand_type(const char *file, int line, const std::string &str, int mode, LAMMPS *lmp);

  /*! Make C-style copy of string in new storage
   *
   * This allocates a storage buffer and copies the C-style or
   * C++ style string into it.  The buffer is allocated with "new"
   * and thus needs to be deallocated with "delete[]".
   *
   * \param text  string that should be copied
   * \return new buffer with copy of string */

  char *strdup(const std::string &text);

  /*! Convert string to lowercase
   *
   * \param line  string that should be converted
   * \return new string with all lowercase characters */

  std::string lowercase(const std::string &line);

  /*! Convert string to uppercase
   *
   * \param line  string that should be converted
   * \return new string with all uppercase characters */

  std::string uppercase(const std::string &line);

  /*! Trim leading and trailing whitespace. Like TRIM() in Fortran.
   *
   * \param line  string that should be trimmed
   * \return new string without whitespace (string) */

  std::string trim(const std::string &line);

  /*! Return string with anything from the first '#' character onward removed
   *
   * \param line  string that should be trimmed
   * \return new string without comment (string) */

  std::string trim_comment(const std::string &line);

  /*! Replace first '*' character in a string with a number, optionally zero-padded
   *
   * If there is no '*' character in the string, return the original string.
   * If the number requires more characters than the value of the *pad*
   * argument, do not add zeros; otherwise add as many zeroes as needed to
   * the left to make the the number representation *pad* characters wide.
   *
   * \param name  string with file containing a '*' (or not)
   * \param step  step number to replace the (first) '*'
   * \param pad   zero-padding (may be zero)
   * \return  processed string */

  std::string star_subst(const std::string &name, bigint step, int pad);

  /*! Check if a string will likely have UTF-8 encoded characters
   *
   * UTF-8 uses the 7-bit standard ASCII table for the first 127 characters and
   * all other characters are encoded as multiple bytes.  For the multi-byte
   * characters the first byte has either the highest two, three, or four bits
   * set followed by a zero bit and followed by one, two, or three more bytes,
   * respectively, where the highest bit is set and the second highest bit set
   * to 0.  The remaining bits combined are the character code, which is thus
   * limited to 21-bits.
   *
   * For the sake of efficiency this test only checks if a character in the string
   * has the highest bit set and thus is very likely an UTF-8 character.  It will
   * not be able to tell this this is a valid UTF-8 character or whether it is a
   * 2-byte, 3-byte, or 4-byte character.
   *
\verbatim embed:rst

*See also*
   :cpp:func:`utils::utf8_subst`

\endverbatim
   * \param line  string that should be checked
   * \return true if string contains UTF-8 encoded characters (bool) */

  inline bool has_utf8(const std::string &line)
  {
    for (auto c : line)
      if (c & 0x80U) return true;
    return false;
  }

  /*! Replace known UTF-8 characters with ASCII equivalents
   *
\verbatim embed:rst

*See also*
   :cpp:func:`utils::has_utf8`

\endverbatim
   * \param line  string that should be converted
   * \return new string with ascii replacements (string) */

  std::string utf8_subst(const std::string &line);

  /*! Count words in string with custom choice of separating characters
   *
   * \param text string that should be searched
   * \param separators string containing characters that will be treated as whitespace
   * \return number of words found */

  size_t count_words(const std::string &text, const std::string &separators);

  /*! Count words in string, ignore any whitespace matching " \t\r\n\f"
   *
   * \param text string that should be searched
   * \return number of words found */

  size_t count_words(const std::string &text);

  /*! Count words in C-string, ignore any whitespace matching " \t\r\n\f"
   *
   * \param text string that should be searched
   * \return number of words found */

  size_t count_words(const char *text);

  /*! Count words in a single line, trim anything from '#' onward
   *
   * \param text string that should be trimmed and searched
   * \param separators string containing characters that will be treated as whitespace
   * \return number of words found */

  size_t trim_and_count_words(const std::string &text, const std::string &separators = " \t\r\n\f");

  /*! Take list of words and join them with a given separator text.
   *
   * This is the inverse operation of what the split_words() function
   * Tokenizer classes do.
   *
   * \param words  STL vector with strings
   * \param sep    separator string (may be empty)
   * \return  string with the concatenated words and separators */

  std::string join_words(const std::vector<std::string> &words, const std::string &sep);

  /*! Take text and split into non-whitespace words.
   *
   * This can handle strings with single and double quotes, escaped quotes,
   * and escaped codes within quotes, but due to using an STL container and
   * STL strings is rather slow because of making copies. Designed for
   * parsing command lines and similar text and not for time critical
   * processing.  Use a tokenizer class if performance matters.
   *
\verbatim embed:rst

*See also*
   :cpp:class:`Tokenizer`, :cpp:class:`ValueTokenizer`

\endverbatim
   * \param text string that should be split
   * \return STL vector with the words */

  std::vector<std::string> split_words(const std::string &text);

  /*! Take multi-line text and split into lines
   *
   * \param text string that should be split
   * \return STL vector with the lines */
  std::vector<std::string> split_lines(const std::string &text);

  /*! Check if string can be converted to valid integer
   *
   * \param str string that should be checked
   * \return true, if string contains valid a integer, false otherwise */

  bool is_integer(const std::string &str);

  /*! Check if string can be converted to valid floating-point number
   *
   * \param str string that should be checked
   * \return true, if string contains valid number, false otherwise */

  bool is_double(const std::string &str);

  /*! Check if string is a valid ID
   * ID strings may contain only letters, numbers, and underscores.
   *
   * \param str string that should be checked
   * \return true, if string contains valid id, false otherwise */

  bool is_id(const std::string &str);

  /*! Check if string is a valid type label, or numeric type, or numeric type range.
   * Numeric type or type range may only contain digits or the '*' character.
   * Type label strings may not contain a digit, or a '*', or a '#' character as the
   * first character to distinguish them from comments and numeric types or type ranges.
   * They also may not contain any whitespace. If the string is a valid numeric type
   * or type range the function returns 0, if it is a valid type label the function
   * returns 1, otherwise it returns -1.
   *
   * \param str string that should be checked
   * \return 0, 1, or -1, depending on whether the string is valid numeric type, valid type label or neither, respectively */

  int is_type(const std::string &str);

  /*! Determine full path of potential file. If file is not found in current directory,
   *  search directories listed in LAMMPS_POTENTIALS environment variable
   *
   * \param path file path
   * \return full path to potential file */

  std::string get_potential_file_path(const std::string &path);

  /*! Read potential file and return DATE field if it is present
   *
   * \param path file path
   * \param potential_name name of potential that is being read
   * \return DATE field if present */

  std::string get_potential_date(const std::string &path, const std::string &potential_name);

  /*! Read potential file and return UNITS field if it is present
   *
   * \param path file path
   * \param potential_name name of potential that is being read
   * \return UNITS field if present */

  std::string get_potential_units(const std::string &path, const std::string &potential_name);

  enum { NOCONVERT = 0, METAL2REAL = 1, REAL2METAL = 1 << 1 };
  enum { UNKNOWN = 0, ENERGY };

  /*! Return bitmask of available conversion factors for a given property
   *
   * \param property property to be converted
   * \return bitmask indicating available conversions */

  int get_supported_conversions(const int property);

  /*! Return unit conversion factor for given property and selected from/to units
   *
   * \param property property to be converted
   * \param conversion constant indicating the conversion
   * \return conversion factor */

  double get_conversion_factor(const int property, const int conversion);

  /*! Open a potential file as specified by *name*
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

  /*! Convert a time string to seconds
   *
   * The strings "off" and "unlimited" result in -1
   *
   * \param timespec a string in the following format: ([[HH:]MM:]SS)
   * \return total in seconds */

  double timespec2seconds(const std::string &timespec);

  /*! Convert a LAMMPS version date to a number
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

  /*! Return current date as string
   *
   * This will generate a string containing the current date in YYYY-MM-DD format.
   *
   * \return       string with current date */

  std::string current_date();

  /*! Binary search in a vector of ascending doubles of length N
   *
   * If the value is smaller than the smallest value in the vector, 0 is returned.
   * If the value is larger or equal than the largest value in the vector, N-1 is returned.
   * Otherwise the index that satisfies the condition
   *
   * haystack[index] <= value < haystack[index+1]
   *
   * is returned, i.e. a value from 1 to N-2. Note that if there are tied values in the
   * haystack, always the larger index is returned as only that satisfied the condition.
   *
   * \param  needle    search value for which are are looking for the closest index
   * \param  n         size of the haystack array
   * \param  haystack  array with data in ascending order.
   * \return           index of value in the haystack array smaller or equal to needle */
  int binary_search(const double needle, const int n, const double *haystack);

  /*! Custom merge sort implementation
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

  void merge_sort(int *index, int num, void *ptr, int (*comp)(int, int, void *));
}    // namespace utils
}    // namespace LAMMPS_NS

#endif
