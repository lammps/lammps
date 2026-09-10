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

#ifdef FIX_CLASS
// clang-format off
FixStyle(graphics/lines,FixGraphicsLines);
// clang-format on
#else

#ifndef LMP_FIX_GRAPHICS_LINES_H
#define LMP_FIX_GRAPHICS_LINES_H

#include "fix.h"

namespace LAMMPS_NS {
class Compute;
class FixStoreAtom;

class FixGraphicsLines : public Fix {
 public:
  FixGraphicsLines(class LAMMPS *, int, char **);
  ~FixGraphicsLines() override;
  void post_constructor() override;
  int setmask() override;
  void setup(int) override;
  void end_of_step() override;

  void write_restart(FILE *) override;
  void restart(char *) override;

  double compute_scalar() override;
  int image(int *&, double **&) override;

 protected:
  int nrepeat, nfreq, nlength;
  int nvalues, ivalue, firstflag;

  // for storing graphics objects
  int numobjs;
  int *imgobjs;
  double **imgparms;

  char *id_cprop;          // internal compute property/atom ID
  Compute *cprop;          // pointer to internal compute property/atom instance
  char *id_fave;           // internal fix ave/atom ID
  Fix *fave;               // pointer to internal fix ave/atom instance
  char *id_fstore;         // internal fix STORE/ATOM ID
  FixStoreAtom *fstore;    // pointer to internal fix STORE/ATOM instance
};
}    // namespace LAMMPS_NS
#endif
#endif
