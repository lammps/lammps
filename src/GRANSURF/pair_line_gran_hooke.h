/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(line/gran/hooke,PairLineGranHooke)

#else

#ifndef LMP_PAIR_LINE_GRAN_HOOKE_H
#define LMP_PAIR_LINE_GRAN_HOOKE_H

#include "pair_line_gran_hooke_history.h"

namespace LAMMPS_NS {

class PairLineGranHooke : public PairLineGranHookeHistory {
 public:
  PairLineGranHooke(class LAMMPS *);
  ~PairLineGranHooke() {}
  void compute(int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Pair gran/line has inconsitent internal size

UNDOCUMENTED

E: Pair gran/line iteraction between 2 lines

UNDOCUMENTED

E: Pair gran/line iteraction between 2 spheres

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair gran/line requires atom style line

UNDOCUMENTED

E: Pair gran/line requires ghost atoms store velocity

UNDOCUMENTED

E: Pair granular with shear history requires newton pair off

This is a current restriction of the implementation of pair
granular styles with history.

E: Could not find pair fix ID

UNDOCUMENTED

E: Pair gran/line found %g line end pts with more than 2 lines

UNDOCUMENTED

U: Pair granular requires atom style sphere

Self-explanatory.

U: Pair granular requires ghost atoms store velocity

Use the comm_modify vel yes command to enable this.

*/
