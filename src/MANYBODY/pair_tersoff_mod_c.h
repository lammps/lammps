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

#ifdef PAIR_CLASS

PairStyle(tersoff/mod/c,PairTersoffMODC)

#else

#ifndef LMP_PAIR_TERSOFF_MOD_C_H
#define LMP_PAIR_TERSOFF_MOD_C_H

#include "pair_tersoff_mod.h"

namespace LAMMPS_NS {

class PairTersoffMODC : public PairTersoffMOD {
 public:
  PairTersoffMODC(class LAMMPS *lmp) : PairTersoffMOD(lmp) {};
  ~PairTersoffMODC() {}

 protected:
  void read_file(char *);
  void repulsive(Param *, double, double &, int, double &);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot open Tersoff potential file %s

The specified potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in Tersoff potential file

Incorrect number of words per line in the potential file.

E: Illegal Tersoff parameter

One or more of the coefficients defined in the potential file is
invalid.

U: Potential file has duplicate entry

The potential file has more than one entry for the same element.

U: Potential file is missing an entry

The potential file does not have a needed entry.

*/
