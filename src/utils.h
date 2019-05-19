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

#include <string>
#include <cstdio>

namespace LAMMPS_NS {

  // forward declarations
  class Error;

  namespace utils {

    /** \brief Match text against a simplified regex pattern
     *
     *  \param text the text to be matched against the pattern
     *  \param pattern the search pattern, which may contain regexp markers
     *  \return true if the pattern matches, false if not
     */
    bool strmatch(std::string text, std::string pattern);

    /** Categories of special arguments for cfvarg() function
     *
     * Enum starts from 100 to avoid conflicts with other local define flags
     */
    enum {NONE=100,              /// does not match any category
          COMPUTE,               /// processed a compute
          FIX,                   /// processed a fix
          VARIABLE               /// processed a variable
    };

    /** \brief Convenience function to process 'c_', 'f_', and 'v_' arguments
     *
     *  \param mode types to search for. 1-3 char string from 'c', 'f', or 'v'
     *  \param arg  argument string to test against the prefixes
     *  \param cfv_id name or ID of the compute, fix, or variable
     *  \return utils::COMPUTE, utils::FIX, utils::VARIABLE or utils::NONE
     */
    int cfvarg(std::string mode, const char *arg, char *&cfv_id);

    /** \brief safe wrapper around fgets() which aborts on errors
     *  or EOF and prints a suitable error message to help debugging
     *
     *  \param srcname  name of the calling source file (from FLERR macro)
     *  \param srcline  line in the calling source file (from FLERR macro)
     *  \param s        buffer for storing the result of fgets()
     *  \param size     size of buffer s (max number of bytes read by fgets())
     *  \param fp       file pointer used by fgets()
     *  \param filename file name associated with fp (for error message)
     *  \param error    pointer to Error class instance (for abort)
     */
    void sfgets(const char *srcname, int srcline, char *s, int size,
                FILE *fp, const char *filename, Error *error);

    int kim_simulator_json_parse(int argc, char **argv);
  }
}

#endif

/* ERROR/WARNING messages:

*/
