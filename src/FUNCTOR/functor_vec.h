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

#ifndef LMP_FUNCTOR_VEC_H
#define LMP_FUNCTOR_VEC_H

// Helper used by the FUNCTOR pair driver: the per-atom position and force
// arrays (atom->x, atom->f) are allocated as one contiguous block of doubles
// with a separate array of row pointers.  Casting the first row pointer to a
// vec3_t* lets the compiler treat the data as a contiguous array of structs
// (stride 3 doubles), which improves vectorization and cache behavior.  This
// is the same trick used in src/OPT/pair_lj_cut_opt.cpp.

namespace LAMMPS_NS::functor {

struct vec3_t {
  double x, y, z;
};

}    // namespace LAMMPS_NS::functor

#endif
