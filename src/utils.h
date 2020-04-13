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
    std::string check_packages_for_style(std::string style,
                                         std::string name, LAMMPS *lmp);

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
  }
}

#endif

/* ERROR/WARNING messages:

*/
