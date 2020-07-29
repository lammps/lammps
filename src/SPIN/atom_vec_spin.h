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

AtomStyle(spin,AtomVecSpin)

#else

#ifndef LMP_ATOM_VEC_SPIN_H
#define LMP_ATOM_VEC_SPIN_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecSpin : public AtomVec {
 public:
  AtomVecSpin(class LAMMPS *);

  void grow_pointers();
  void force_clear(int, size_t);
  void data_atom_post(int);

 private:
  double **sp,**fm,**fm_long;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
