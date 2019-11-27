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

#ifdef ATOM_CLASS

AtomStyle(sphere,AtomVecSphere)

#else

#ifndef LMP_ATOM_VEC_SPHERE_H
#define LMP_ATOM_VEC_SPHERE_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecSphere : public AtomVec {
 public:
  AtomVecSphere(class LAMMPS *);
  void process_args(int, char **);
  void init();
  void create_atom_post(int);
  void data_atom_post(int);
  void pack_data_pre(int);
  void pack_data_post(int);

 private:
  int radvary;
  double radius,rmass;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid radius in Atoms section of data file

Radius must be >= 0.0.

E: Invalid density in Atoms section of data file

Density value cannot be <= 0.0.

*/
