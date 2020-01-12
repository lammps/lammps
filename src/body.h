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

#ifndef LMP_BODY_H
#define LMP_BODY_H

#include "pointers.h"
#include "atom_vec_body.h"

namespace LAMMPS_NS {

class Body : protected Pointers {
 public:
  MyPoolChunk<int> *icp;
  MyPoolChunk<double> *dcp;

  char *style;
  int size_forward;           // max extra values packed for comm
  int size_border;            // max extra values packed for border comm
  int maxexchange;            // max size of exchanged atom

  AtomVecBody *avec;          // ptr to class that stores body bonus info

  Body(class LAMMPS *, int, char **);
  virtual ~Body();

  // methods implemented by child classes

  virtual int pack_comm_body(struct AtomVecBody::Bonus *, double *) {return 0;}
  virtual int unpack_comm_body(struct AtomVecBody::Bonus *, double *) {return 0;}
  virtual int pack_border_body(struct AtomVecBody::Bonus *, double *) {return 0;}
  virtual int unpack_border_body(struct AtomVecBody::Bonus *,
                                 double *) {return 0;}

  virtual void data_body(int, int, int, int*, double *) = 0;
  virtual int noutrow(int) = 0;
  virtual int noutcol() = 0;
  virtual void output(int, int, double *) = 0;
  virtual int image(int, double, double, int *&, double **&) = 0;

  virtual double radius_body(int, int, int *, double *) {return 0.0;}
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
