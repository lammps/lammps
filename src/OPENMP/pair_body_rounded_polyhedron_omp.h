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

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(body/rounded/polyhedron/omp,PairBodyRoundedPolyhedronOMP);
// clang-format on
#else

#ifndef LMP_PAIR_BODY_ROUNDED_POLYHEDRON_OMP_H
#define LMP_PAIR_BODY_ROUNDED_POLYHEDRON_OMP_H

#include "pair_body_rounded_polyhedron.h"
#include "thr_omp.h"

#include <vector>

namespace LAMMPS_NS {

// threaded compute of pair style body/rounded/polyhedron and of the pair
// styles derived from it, given as BASE, which provide the contact models.
// The member functions are defined in pair_body_rounded_polyhedron_omp_impl.h
// and instantiated in the source file of each threaded style.

template <class BASE> class PairBodyRoundedPolyhedronOMPT : public BASE, public ThrOMP {

 public:
  PairBodyRoundedPolyhedronOMPT(class LAMMPS *);
  ~PairBodyRoundedPolyhedronOMPT() override;

  void compute(int, int) override;
  double memory_usage() override;

 protected:
  // set up the per-thread data of the contact models before the threaded loop
  virtual void init_thr(int /*nthreads*/) {}

 private:
  double **fnc_thr;    // per-thread fnc, see PairBodyRoundedPolyhedron
  int nmax_thr;        // allocated number of atoms per thread in fnc_thr
  std::vector<PairBodyRoundedPolyhedron::Scratch> scratch_thr;    // per-thread scratch space

  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData *const thr, double **fnc_t,
            PairBodyRoundedPolyhedron::Scratch &s);
};

class PairBodyRoundedPolyhedronOMP : public PairBodyRoundedPolyhedronOMPT<PairBodyRoundedPolyhedron> {
 public:
  PairBodyRoundedPolyhedronOMP(class LAMMPS *);
};

extern template class PairBodyRoundedPolyhedronOMPT<PairBodyRoundedPolyhedron>;

}    // namespace LAMMPS_NS

#endif
#endif
