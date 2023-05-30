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
AtomStyle(wavepacket,AtomVecWavepacket);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_WAVEPACKET_H
#define LMP_ATOM_VEC_WAVEPACKET_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecWavepacket : virtual public AtomVec {
 public:
  AtomVecWavepacket(class LAMMPS *);

  void grow_pointers() override;
  void force_clear(int, size_t) override;
  void create_atom_post(int) override;
  void data_atom_post(int) override;
  int property_atom(const std::string &) override;
  void pack_property_atom(int, double *, int, int) override;

 private:
  int *spin;
  double *q, *eradius, *ervel, *erforce;
};

}    // namespace LAMMPS_NS

#endif
#endif
