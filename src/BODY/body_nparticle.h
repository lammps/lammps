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

#ifdef BODY_CLASS
// clang-format off
BodyStyle(nparticle,BodyNparticle);
// clang-format on
#else

#ifndef LMP_BODY_NPARTICLE_H
#define LMP_BODY_NPARTICLE_H

#include "atom_vec_body.h"
#include "body.h"

namespace LAMMPS_NS {

class BodyNparticle : public Body {
 public:
  BodyNparticle(class LAMMPS *, int, char **);
  ~BodyNparticle() override;
  int nsub(struct AtomVecBody::Bonus *);
  double *coords(struct AtomVecBody::Bonus *);

  int pack_border_body(struct AtomVecBody::Bonus *, double *) override;
  int unpack_border_body(struct AtomVecBody::Bonus *, double *) override;
  void data_body(int, int, int, int *, double *) override;
  int pack_data_body(tagint, int, double *) override;
  int write_data_body(FILE *, double *) override;
  double radius_body(int, int, int *, double *) override;

  int noutrow(int) override;
  int noutcol() override;
  void output(int, int, double *) override;
  int image(int, double, double, int *&, double **&) override;

 private:
  int *imflag;
  double **imdata;
};

}    // namespace LAMMPS_NS

#endif
#endif
