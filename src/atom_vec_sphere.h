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
AtomStyle(sphere,AtomVecSphere);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_SPHERE_H
#define LMP_ATOM_VEC_SPHERE_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecSphere : virtual public AtomVec {
 public:
  AtomVecSphere(class LAMMPS *);
  void process_args(int, char **) override;
  void init() override;

  void grow_pointers() override;
  void create_atom_post(int) override;
  void data_atom_post(int) override;
  void pack_data_pre(int) override;
  void pack_data_post(int) override;

 protected:
  double *radius, *rmass;
  double **omega;

  int radvary;
  double radius_one, rmass_one;
};

}    // namespace LAMMPS_NS

#endif
#endif
