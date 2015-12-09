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

#ifdef COMPUTE_CLASS

ComputeStyle(dilatation/atom,ComputeDilatationAtom)

#else

#ifndef LMP_COMPUTE_DILATATION_ATOM_H
#define LMP_COMPUTE_DILATATION_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeDilatationAtom : public Compute {
  friend class PairPeriPMB;
  friend class PairPeriLPS;
  friend class PairPeriVES;
  friend class PairPeriEPS;
 public:
  ComputeDilatationAtom(class LAMMPS *, int, char **);
  ~ComputeDilatationAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double *dilatation;
  int isPMB,isLPS,isVES,isEPS;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: More than one compute dilatation/atom

Self-explanatory.

E: Compute dilatation/atom cannot be used with this pair style

Self-explanatory.

E: Compute dilatation/atom requires Peridynamic pair style

Self-explanatory.

*/
