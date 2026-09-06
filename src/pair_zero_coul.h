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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(zero/coul,PairZeroCoul);
// clang-format on
#else

#ifndef LMP_PAIR_ZERO_COUL_H
#define LMP_PAIR_ZERO_COUL_H

#include "pair_zero.h"

namespace LAMMPS_NS {

/* pair zero/coul is a variant of pair style zero that presents itself as
   a Coulombic pair style: it declares compatibility with kspace styles
   and provides the real-space Coulomb cutoff, but computes no
   interactions, just like pair style zero.  This allows computing only
   the k-space contribution to forces and energies, or satisfying styles
   that require the presence of a Coulombic pair style, e.g. for testing
   and debugging purposes. */

class PairZeroCoul : public PairZero {
 public:
  PairZeroCoul(class LAMMPS *);
  void *extract(const char *, int &) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
