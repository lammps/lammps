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

#ifndef LMP_PYTHON_WRAPPER_H
#define LMP_PYTHON_WRAPPER_H

// true interface to embedded Python
// used when PYTHON package is installed

#ifdef LMP_PYTHON

#include "python.h"

#else

// dummy interface to PYTHON
// needed for compiling when PYTHON is not installed

namespace LAMMPS_NS {

class Python {
 public:
  int python_exists;

  Python(class LAMMPS *) {python_exists = 0;}
  ~Python() {}
  void command(int, char **) {}
  void invoke_function(int, char *) {}
  int find(char *) {return -1;}
  int variable_match(char *, char *, int) {return -1;}

};

}

#endif
#endif
