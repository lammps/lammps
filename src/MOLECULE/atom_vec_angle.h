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

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(angle,AtomVecAngle);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_ANGLE_H
#define LMP_ATOM_VEC_ANGLE_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecAngle : public AtomVec {
 public:
  AtomVecAngle(class LAMMPS *);
  ~AtomVecAngle() override;

  void grow_pointers() override;
  void pack_restart_pre(int) override;
  void pack_restart_post(int) override;
  void unpack_restart_init(int) override;
  void data_atom_post(int) override;

 private:
  int *num_bond, *num_angle;
  int **bond_type, **angle_type;
  int **nspecial;

  int any_bond_negative, any_angle_negative;
  int bond_per_atom, angle_per_atom;
  int *bond_negative, *angle_negative;
};

}    // namespace LAMMPS_NS

#endif
#endif
