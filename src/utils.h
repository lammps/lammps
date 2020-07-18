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

    /** \brief Match text against a simplified regex pattern
     *
     *  \param text the text to be matched against the pattern
     *  \param pattern the search pattern, which may contain regexp markers
     *  \return true if the pattern matches, false if not
     */
    bool strmatch(std::string text, std::string pattern);

    /** \brief Send message to screen and logfile, if available
     *
     *  \param lmp   pointer to LAMMPS class instance
     *  \param mesg  message to be printed
     */
    void logmesg(LAMMPS *lmp, const std::string &mesg);

    /** \brief return a string representing the current system error status
     *
     *  This is a wrapper around calling strerror(errno).
     *
     *  \return  error string
     */
    std::string getsyserror();

    /** \brief safe wrapper around fgets() which aborts on errors
     *  or EOF and prints a suitable error message to help debugging
     *
     *  \param srcname  name of the calling source file (from FLERR macro)
     *  \param srcline  line in the calling source file (from FLERR macro)
     *  \param s        buffer for storing the result of fgets()
     *  \param size     size of buffer s (max number of bytes read by fgets())
     *  \param fp       file pointer used by fgets()
     *  \param filename file name associated with fp (may be NULL; then LAMMPS will try to detect)
     *  \param error    pointer to Error class instance (for abort)
     */
    void sfgets(const char *srcname, int srcline, char *s, int size,
                FILE *fp, const char *filename, Error *error);

    /** \brief safe wrapper around fread() which aborts on errors
     *  or EOF and prints a suitable error message to help debugging
     *
     *  \param srcname  name of the calling source file (from FLERR macro)
     *  \param srcline  line in the calling source file (from FLERR macro)
     *  \param s        buffer for storing the result of fread()
     *  \param size     size of data elements read by fread()
     *  \param num      number of data elements read by fread()
     *  \param fp       file pointer used by fread()
     *  \param filename file name associated with fp (may be NULL; then LAMMPS will try to detect)
     *  \param error    pointer to Error class instance (for abort)
     */
    void sfread(const char *srcname, int srcline, void *s, size_t size,
                size_t num, FILE *fp, const char *filename, Error *error);

    /** \brief Report if a requested style is in a package or may have a typo
     *
     *  \param style type of style that is to be checked for
     *  \param name  name of style that was not found
     *  \param lmp   pointer to top-level LAMMPS class instance
     *  \return string usable for error messages
     */
    std::string check_packages_for_style(const std::string &style,
                                         const std::string &name, LAMMPS *lmp);

    /** \brief Convert a string to a floating point number while checking
        if it is a valid floating point or integer number
     *
     *  \param file name of source file for error message
     *  \param line in source file for error message
     *  \param str  string to be converted to number
     *  \param do_abort determines whether to call Error::one() or Error::all()
     *  \param lmp   pointer to top-level LAMMPS class instance
     *  \return double precision floating point number
     */
    double numeric(const char *file, int line, const char *str,
                   bool do_abort, LAMMPS *lmp);

    /** \brief Convert a string to an integer number while checking
        if it is a valid integer number (regular int)
     *
     *  \param file name of source file for error message
     *  \param line in source file for error message
     *  \param str  string to be converted to number
     *  \param do_abort determines whether to call Error::one() or Error::all()
     *  \param lmp   pointer to top-level LAMMPS class instance
     *  \return integer number (regular int)
     */
    int inumeric(const char *file, int line, const char *str,
                 bool do_abort, LAMMPS *lmp);

    /** \brief Convert a string to an integer number while checking
        if it is a valid integer number (bigint)
     *
     *  \param file name of source file for error message
     *  \param line in source file for error message
     *  \param str  string to be converted to number
     *  \param do_abort determines whether to call Error::one() or Error::all()
     *  \param lmp   pointer to top-level LAMMPS class instance
     *  \return integer number (bigint)
     */
    bigint bnumeric(const char *file, int line, const char *str,
                    bool do_abort, LAMMPS *lmp);

    /** \brief Convert a string to an integer number while checking
        if it is a valid integer number (tagint)
     *
     *  \param file name of source file for error message
     *  \param line in source file for error message
     *  \param str  string to be converted to number
     *  \param do_abort determines whether to call Error::one() or Error::all()
     *  \param lmp   pointer to top-level LAMMPS class instance
     *  \return integer number (tagint)
     */
    tagint tnumeric(const char *file, int line, const char *str,
                    bool do_abort, LAMMPS *lmp);

    /**
     * \brief Trim anything from '#' onward
     * \param line string that should be trimmed
     * \return new string without comment (string)
     */
    std::string trim_comment(const std::string & line);

    /**
     * \brief Count words in string
     * \param text string that should be searched
     * \param separators string containing characters that will be treated as whitespace
     * \return number of words found
     */
    size_t count_words(const std::string & text, const std::string & separators);

    /**
     * \brief Count words in string, ignore any whitespace matching " \t\r\n\f"
     * \param text string that should be searched
     * \param separators string containing characters that will be treated as whitespace
     * \return number of words found
     */
    size_t count_words(const std::string & text);

    /**
     * \brief Count words in C-string, ignore any whitespace matching " \t\r\n\f"
     * \param text string that should be searched
     * \param separators string containing characters that will be treated as whitespace
     * \return number of words found
     */
    size_t count_words(const char * text);

    /**
     * \brief Count words in a single line, trim anything from '#' onward
     * \param text string that should be trimmed and searched
     * \param separators string containing characters that will be treated as whitespace
     * \return number of words found
     */
    size_t trim_and_count_words(const std::string & text, const std::string & separators = " \t\r\n\f");

    /**
     * \brief Take text and split into non-whitespace words.
     *
     * This can handle single and double quotes, escaped quotes,
     * and escaped codes within quotes, but due to using an STL
     * container and STL strings is rather slow because of making
     * copies. Designed for parsing command lines and similar text
     * and not for time critical processing. Use a tokenizer for that.
     *
     * \param text string that should be split
     * \return STL vector with the words
     */
    std::vector<std::string> split_words(const std::string &text);

    /**
     * \brief Check if string can be converted to valid integer
     * \param text string that should be checked
     * \return true, if string contains valid integer, false otherwise
     */
    bool is_integer(const std::string & str);

    /**
     * \brief Check if string can be converted to valid floating-point number
     * \param text string that should be checked
     * \return true, if string contains valid floating-point number, false otherwise
     */
    bool is_double(const std::string & str);

    /** \brief try to detect pathname from FILE pointer.
     *
     * Currently only supported on Linux, otherwise will report "(unknown)".
     *
     *  \param buf  storage buffer for pathname. output will be truncated if not large enough
     *  \param len  size of storage buffer. output will be truncated to this length - 1
     *  \param fp   FILE pointer structe from STDIO library for which we want to detect the name
     *  \return pointer to the storage buffer, i.e. buf
     */
    const char *guesspath(char *buf, int len, FILE *fp);

    /**
     * \brief Strip off leading part of path, return just the filename
     * \param path file path
     * \return file name
     */
    std::string path_basename(const std::string & path);

    /**
     * \brief Join two paths
     * \param a first path
     * \param b second path
     * \return combined path
     */
    std::string path_join(const std::string & a, const std::string & b);

    /**
     * \brief Check if file exists and is readable
     * \param path file path
     * \return true if file exists and is readable
     */
    bool file_is_readable(const std::string & path);

    /**
     * \brief Determine full path of potential file
     *        If file is not found in current directory, search LAMMPS_POTENTIALS folder
     * \param path file path
     * \return full path to potential file
     */
    std::string get_potential_file_path(const std::string& path);

    /**
     * \brief Read potential file and return DATE field if it is present
     * \param path file path
     * \param potential_name name of potential that is being read
     * \return DATE field if present
     */
    std::string get_potential_date(const std::string & path, const std::string & potential_name);

    /**
     * \brief Read potential file and return UNITS field if it is present
     * \param path file path
     * \param potential_name name of potential that is being read
     * \return UNITS field if present
     */
    std::string get_potential_units(const std::string & path, const std::string & potential_name);

    enum { NOCONVERT = 0, METAL2REAL = 1, REAL2METAL = 1<<1 };
    enum { UNKNOWN = 0, ENERGY };

    /**
     * \brief Return bitmask of available conversion factors for a given propert
     * \param property property to be converted
     * \return bitmask indicating available conversions
     */
    int get_supported_conversions(const int property);

    /**
     * \brief Return unit conversion factor for given property and selected from/to units
     * \param property property to be converted
     * \param conversion constant indicating the conversion
     * \return conversion factor
     */
    double get_conversion_factor(const int property, const int conversion);
  }
}

#endif

/* ERROR/WARNING messages:

*/
