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

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(spin,AtomVecSpin);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_SPIN_H
#define LMP_ATOM_VEC_SPIN_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecSpin : virtual public AtomVec {
 public:
  AtomVecSpin(class LAMMPS *);

  void grow_pointers() override;
  void force_clear(int, size_t) override;
  void data_atom_post(int) override;

 protected:
  double **sp, **fm, **fm_long;
};

}    // namespace LAMMPS_NS

#endif
#endif
