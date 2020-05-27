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

#ifdef BODY_CLASS

BodyStyle(rounded/polyhedron,BodyRoundedPolyhedron)

#else

#ifndef LMP_BODY_ROUNDED_POLYHEDRON_H
#define LMP_BODY_ROUNDED_POLYHEDRON_H

#include "body.h"
#include "atom_vec_body.h"

namespace LAMMPS_NS {

class BodyRoundedPolyhedron : public Body {
 public:
  BodyRoundedPolyhedron(class LAMMPS *, int, char **);
  ~BodyRoundedPolyhedron();
  int nsub(struct AtomVecBody::Bonus *);
  double *coords(struct AtomVecBody::Bonus *);
  int nedges(struct AtomVecBody::Bonus *);
  double *edges(struct AtomVecBody::Bonus *);
  int nfaces(struct AtomVecBody::Bonus *);
  double *faces(struct AtomVecBody::Bonus *);
  double enclosing_radius(struct AtomVecBody::Bonus *);
  double rounded_radius(struct AtomVecBody::Bonus *);

  int pack_border_body(struct AtomVecBody::Bonus *, double *);
  int unpack_border_body(struct AtomVecBody::Bonus *, double *);
  void data_body(int, int, int, int *, double *);
  double radius_body(int, int, int *, double *);

  int noutrow(int);
  int noutcol();
  void output(int, int, double *);
  int image(int, double, double, int *&, double **&);

 private:
  int *imflag;
  double **imdata;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid body rounded/polyhedron command

Arguments in atom-style command are not correct.

E: Invalid format in Bodies section of data file

The specified number of integer or floating point values does not
appear.

E: Incorrect # of integer values in Bodies section of data file

See doc page for body style.

E: Incorrect integer value in Bodies section of data file

See doc page for body style.

E: Incorrect # of floating-point values in Bodies section of data file

See doc page for body style.

E: Insufficient Jacobi rotations for body rounded/polyhedron

Eigensolve for rigid body was not sufficiently accurate.

*/
