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
PairStyle(hybrid/kk,PairHybridKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_HYBRID_KOKKOS_H
#define LMP_PAIR_HYBRID_KOKKOS_H

#include "pair_hybrid.h"
#include "pair_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class PairHybridKokkos : public PairHybrid {
  friend class FixGPU;
  friend class FixIntel;
  friend class FixOMP;
  friend class Force;
  friend class Respa;
  friend class Info;
 public:
  typedef LMPDeviceType device_type;

  PairHybridKokkos(class LAMMPS *);

  void compute(int, int) override;
  void init_style() override;

 private:
  DAT::t_x_array_randomread x;
  DAT::t_f_array f;
  friend void pair_virial_fdotr_compute<PairHybridKokkos>(PairHybridKokkos*);
};

}

#endif
#endif

