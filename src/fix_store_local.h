/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(store/local,FixStoreLocal);
// clang-format on
#else

#ifndef LMP_FIX_STORE_LOCAL_H
#define LMP_FIX_STORE_LOCAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStoreLocal : public Fix {
 public:
  FixStoreLocal(class LAMMPS *, int, char **);
  ~FixStoreLocal();
  int setmask();
  void post_force(int);
  double memory_usage();
  void add_data(double *, int, int);
  int nvalues;

 private:
  int nmax;

  double *vector;
  double **array;

  int ncount;
  int nreset;

  void reallocate(int);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid keyword in fix store/local command

Self-explanatory.

E: Unused instance of fix store/local

Instance of fix store/local is not associated with any other LAMMPS 
class such as a bond style, pair style, etc.

*/
