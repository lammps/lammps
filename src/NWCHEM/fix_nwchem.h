/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nwchem,FixNWChem);
// clang-format on
#else

#ifndef LMP_FIX_NWCHEM_H
#define LMP_FIX_NWCHEM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNWChem : public Fix {
 public:
  FixNWChem(class LAMMPS *, int, char **);
  virtual ~FixNWChem();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void pre_force(int);
  void post_force(int);
  void min_post_force(int);
  double compute_scalar();

 protected:
  char *id_pe;
  int pbcflag;

  int nqm;
  tagint *qmIDs;
  double **xqm,**fqm;
  double *qpotential,*qqm;
  int *qm2lmp;
  double qmenergy;
  
  double lmp2qm_distance,lmp2qm_energy,qm2lmp_force;

  class Compute *c_pe;
  class Pair *pair_coul;

  int pspw_minimizer(MPI_Comm, int, double *, double *, double *, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Must use units metal with fix nwchem command

UNDOCUMENTED

E: Fix nwchem currently runs only in serial

UNDOCUMENTED

E: LAMMPS is linked against incompatible NWChem library

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix nwchem does not yet support a LAMMPS calculation of a Coulomb potential

UNDOCUMENTED

E: Could not find fix nwchem compute ID

UNDOCUMENTED

E: Fix nwchem compute ID does not compute pe/atom

UNDOCUMENTED

E: Fix nwchem requires 3d problem

UNDOCUMENTED

E: Fix nwchem cannot compute Coulomb potential

UNDOCUMENTED

E: Fix nwchem requires 3d simulation

UNDOCUMENTED

W: Fix nwchem should come after all other integration fixes

UNDOCUMENTED

E: Internal NWChem problem

UNDOCUMENTED

*/
