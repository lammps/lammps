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

ImproperStyle(cossq,ImproperCossq)

#else

#ifndef LMP_IMPROPER_COSSQ_H
#define LMP_IMPROPER_COSSQ_H

#include "improper.h"

namespace LAMMPS_NS {

class ImproperCossq : public Improper {
 public:
  ImproperCossq(class LAMMPS *);
  virtual ~ImproperCossq();
  virtual void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

 protected:
  double *k, *chi;
  void allocate();
};

}

#endif
#endif
