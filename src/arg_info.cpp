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

#include "arg_info.h"

#include <cstring>
#include <stdexcept>

using namespace LAMMPS_NS;

/** Class for processing references to fixes, computes and variables
 *
 * This class provides an abstraction for the repetitive task of
 * parsing arguments that may contain references to fixes, computes,
 * variables, or custom per-atom properties. It will identify the name
 * and the index value in the first and second dimension, if present.
 *
 * \param arg      string with possible reference
 * \param allowed  integer with bitmap of allowed types of references */

ArgInfo::ArgInfo(const std::string &arg, int allowed) : type(NONE), dim(0), index1(-1), index2(-1)
{
  if (((arg.size() > 3) && (arg[1] == '2') && (arg[2] == '_')) ||
      ((arg.size() > 2) && (arg[1] == '_'))) {
    if ((arg[0] == 'c') && (allowed & COMPUTE))
      type = COMPUTE;
    else if ((arg[0] == 'f') && (allowed & FIX))
      type = FIX;
    else if ((arg[0] == 'v') && (allowed & VARIABLE))
      type = VARIABLE;
    else if ((arg[0] == 'd') && (allowed & DNAME))
      type = DNAME;
    else if ((arg[0] == 'i') && (allowed & INAME))
      type = INAME;
    else {
      index1 = 0;
      name = arg;
      return;
    }
    const int offset = (arg[1] == '_') ? 2 : 3;

    std::size_t has_idx1 = arg.find('[', offset);
    if (has_idx1 != std::string::npos) {
      name = arg.substr(offset, has_idx1 - offset);
      dim = 1;

      std::size_t has_idx2 = arg.find('[', has_idx1 + 1);
      if (has_idx2 != std::string::npos) {
        dim = 2;

        if (arg[arg.size() - 1] != ']') {
          type = UNKNOWN;
        } else {
          try {
            index2 = std::stoi(arg.substr(has_idx2 + 1, arg.size() - (has_idx2 + 2)));
          } catch (std::invalid_argument &) {
            type = UNKNOWN;
          }
        }
      } else
        has_idx2 = arg.size();

      if ((arg[has_idx2 - 1] != ']') || ((dim == 1) && (arg.find(']') != has_idx2 - 1))) {
        type = UNKNOWN;
      } else {
        try {
          index1 = std::stoi(arg.substr(has_idx1 + 1, arg.size() - (has_idx1 + 2)));
        } catch (std::invalid_argument &) {
          type = UNKNOWN;
        }
      }
    } else {
      index1 = 0;
      name = arg.substr(offset);
    }
  } else {
    index1 = 0;
    name = arg;
  }
}

/* ---------------------------------------------------------------------- */

/*! make copy of the ID of the reference as C-style string
 *
 * The ID is copied into a buffer allocated with "new" and thus
 * must be later deleted with "delete []" to avoid a memory leak.
 * Because it is a full copy in a newly allocated buffer, the
 * lifetime of this string extends beyond the the time the ArgInfo
 * class is in scope.
 *
 * \return copy of string as char * */

char *ArgInfo::copy_name()
{
  char *dest = new char[name.size() + 1];
  strcpy(dest, name.c_str());    // NOLINT
  return dest;
}
