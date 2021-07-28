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

#ifdef BODY_CLASS
// clang-format off
BodyStyle(rounded/polygon,BodyRoundedPolygon);
// clang-format on
#else

#ifndef LMP_BODY_ROUNDED_POLYGON_H
#define LMP_BODY_ROUNDED_POLYGON_H

#include "atom_vec_body.h"
#include "body.h"

namespace LAMMPS_NS {

class BodyRoundedPolygon : public Body {
 public:
  BodyRoundedPolygon(class LAMMPS *, int, char **);
  ~BodyRoundedPolygon();
  int nsub(struct AtomVecBody::Bonus *);
  double *coords(struct AtomVecBody::Bonus *);
  int nedges(struct AtomVecBody::Bonus *);
  double *edges(struct AtomVecBody::Bonus *);
  double enclosing_radius(struct AtomVecBody::Bonus *);
  double rounded_radius(struct AtomVecBody::Bonus *);

  int pack_border_body(struct AtomVecBody::Bonus *, double *);
  int unpack_border_body(struct AtomVecBody::Bonus *, double *);
  void data_body(int, int, int, int *, double *);
  int pack_data_body(tagint, int, double *);
  int write_data_body(FILE *, double *);
  double radius_body(int, int, int *, double *);

  int noutrow(int);
  int noutcol();
  void output(int, int, double *);
  int image(int, double, double, int *&, double **&);

 private:
  int *imflag;
  double **imdata;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid body rounded/polygon command

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

E: Insufficient Jacobi rotations for body nparticle

Eigensolve for rigid body was not sufficiently accurate.

*/
