/* ----------------------------------------------------------------------
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

BodyStyle(nparticle,BodyNparticle)

#else

#ifndef LMP_BODY_NPARTICLE_H
#define LMP_BODY_NPARTICLE_H

#include "body.h"
#include "atom_vec_body.h"

namespace LAMMPS_NS {

class BodyNparticle : public Body {
 public:
  BodyNparticle(class LAMMPS *, int, char **);
  ~BodyNparticle();
  int nsub(class AtomVecBody::Bonus *);
  double *coords(class AtomVecBody::Bonus *);

  int pack_border_body(class AtomVecBody::Bonus *, double *);
  int unpack_border_body(class AtomVecBody::Bonus *, double *);
  void data_body(int, int, int, char **, char **);

  int noutrow(int);
  int noutcol();
  void output(int, int, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid body nparticle command

Arguments in atom-style command are not correct.

E: Insufficient Jacobi rotations for body nparticle

Eigensolve for rigid body was not sufficiently accurate.

*/
