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

#ifdef IMPROPER_CLASS

ImproperStyle(harmonic,ImproperHarmonic)

#else

#ifndef LMP_IMPROPER_HARMONIC_H
#define LMP_IMPROPER_HARMONIC_H

#include "stdio.h"
#include "improper.h"

namespace LAMMPS_NS {

class ImproperHarmonic : public Improper {
 public:
  ImproperHarmonic(class LAMMPS *);
  virtual ~ImproperHarmonic();
  virtual void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

 protected:
  double *k,*chi;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Improper problem: %d %ld %d %d %d %d

Conformation of the 4 listed improper atoms is extreme; you may want
to check your simulation geometry.

E: Incorrect args for improper coefficients

Self-explanatory.  Check the input script or data file.

*/
