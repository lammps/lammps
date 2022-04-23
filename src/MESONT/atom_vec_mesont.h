/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(mesont,AtomVecMesoNT);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_MESONT_H
#define LMP_ATOM_VEC_MESONT_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecMesoNT : public AtomVec {
 public:
  AtomVecMesoNT(class LAMMPS *);
};

}    // namespace LAMMPS_NS

#endif
#endif
