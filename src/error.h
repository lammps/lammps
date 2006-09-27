/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef ERROR_H
#define ERROR_H

#include "stdio.h"
#include "lammps.h"

class Error : public LAMMPS {
 public:
  Error() {}
  ~Error() {}

  void universe_all(char *);
  void universe_one(char *);

  void all(char *);
  void one(char *);
  void warning(char *);
};

#endif
