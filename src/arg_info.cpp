/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "arg_info.h"

#include <stdexcept>
#include <cstring>

using namespace LAMMPS_NS;

ArgInfo::ArgInfo(const std::string &arg, int allowed)
  : type(NONE), dim(0), index1(-1), index2(-1)
{
  if ((arg.size() > 2) && (arg[1] == '_')) {
    if ((arg[0] == 'c') && (allowed & COMPUTE))  type = COMPUTE;
    else if ((arg[0] == 'f') && (allowed & FIX)) type = FIX;
    else if ((arg[0] == 'v') && (allowed & VARIABLE)) type = VARIABLE;
    else {
      index1 = 0;
      name = arg;
      return;
    }

    std::size_t has_idx1 = arg.find('[',2);
    if (has_idx1 != std::string::npos) {
      name = arg.substr(2,has_idx1-2);
      dim = 1;

      std::size_t has_idx2 = arg.find('[',has_idx1+1);
      if (has_idx2 != std::string::npos) {
        dim = 2;

        if (arg[arg.size()-1] != ']') {
          type = UNKNOWN;
        } else {
          try {
            index2 = std::stoi(arg.substr(has_idx2+1,arg.size()-(has_idx2+2)));
          } catch (std::invalid_argument &) {
            type = UNKNOWN;
          }
        }
      } else has_idx2 = arg.size();

      if (arg[has_idx2-1] != ']') {
        type = UNKNOWN;
      } else {
        try {
          index1 = std::stoi(arg.substr(has_idx1+1,arg.size()-(has_idx2+2)));
        } catch (std::invalid_argument &) {
          type = UNKNOWN;
        }
      }
    } else {
      index1 = 0;
      name = arg.substr(2);
    }
  } else {
    index1 = 0;
    name = arg;
  }
}

/* ---------------------------------------------------------------------- */

char *ArgInfo::copy_name()
{
  char *dest = new char[name.size()+1];
  strcpy(dest,name.c_str());
  return dest;
}

