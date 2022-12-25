/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(dielectric,AtomVecDielectric);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_DIELECTRIC_H
#define LMP_ATOM_VEC_DIELECTRIC_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecDielectric : virtual public AtomVec {
  friend class PairLJCutCoulDebyeDielectric;
  friend class PairLJLongCoulLongDielectric;

 public:
  AtomVecDielectric(class LAMMPS *);

  void init() override;
  void grow_pointers() override;
  void create_atom_post(int) override;
  void data_atom_post(int) override;
  void unpack_restart_init(int) override;
  int property_atom(const std::string &) override;
  void pack_property_atom(int, double *, int, int) override;

 protected:
  int *num_bond, *num_angle, *num_dihedral, *num_improper;
  int **bond_type, **angle_type, **dihedral_type, **improper_type;
  int **nspecial;

  int bond_per_atom, angle_per_atom, dihedral_per_atom, improper_per_atom;

  double **mu;
  double *area, *ed, *em, *epsilon, *curvature, *q_scaled;
};

}    // namespace LAMMPS_NS

#endif
#endif
