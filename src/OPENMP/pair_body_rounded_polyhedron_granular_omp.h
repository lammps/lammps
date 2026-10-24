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
PairStyle(body/rounded/polyhedron/granular/omp,PairBodyRoundedPolyhedronGranularOMP);
// clang-format on
#else

#ifndef LMP_PAIR_BODY_ROUNDED_POLYHEDRON_GRANULAR_OMP_H
#define LMP_PAIR_BODY_ROUNDED_POLYHEDRON_GRANULAR_OMP_H

#include "pair_body_rounded_polyhedron_granular.h"
#include "pair_body_rounded_polyhedron_omp.h"

namespace LAMMPS_NS {

class PairBodyRoundedPolyhedronGranularOMP :
    public PairBodyRoundedPolyhedronOMPT<PairBodyRoundedPolyhedronGranular> {

 public:
  PairBodyRoundedPolyhedronGranularOMP(class LAMMPS *);
  ~PairBodyRoundedPolyhedronGranularOMP() override;

  void init_style() override;
  double memory_usage() override;

 protected:
  void init_thr(int nthreads) override;
  Granular_NS::GranularModel *contact_model(int itype, int jtype) override;

 private:
  // a copy of all models for each thread, since the models hold the state
  // of the contact whose forces they compute
  std::vector<std::vector<Granular_NS::GranularModel *>> models_thr;
  int models_stale;    // 1 if the copies must be made again, after init_style()

  void destroy_models_thr();
};

extern template class PairBodyRoundedPolyhedronOMPT<PairBodyRoundedPolyhedronGranular>;

}    // namespace LAMMPS_NS

#endif
#endif
