/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(airebo/omp,PairAIREBOOMP);
// clang-format on
#else

#ifndef LMP_PAIR_AIREBO_OMP_H
#define LMP_PAIR_AIREBO_OMP_H

#include "pair_airebo.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairAIREBOOMP : public PairAIREBO, public ThrOMP {

 public:
  PairAIREBOOMP(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

 protected:
  double bondorder_thr(int i, int j, double rij[3], double rijmag, double VA, ThrData *const thr);
  double bondorderLJ_thr(int i, int j, double rij[3], double rijmag, double VA, double rij0[3],
                         double rijmag0, ThrData *const thr);

  // threaded mirror of PairAIREBO::bondorder_Pij_cross (bond-centric P cross
  // force); no-op for atom-centric P, so airebo/omp is unchanged
  void bondorder_Pij_cross_thr(int i, int j, int itype, int jtype, double VA, double tmppij,
                               double tmppji, const double dN2PIJ[2], const double dN2PJI[2],
                               double *const *const f, ThrData *const thr);

  void FREBO_thr(int ifrom, int ito, int eflag, double *pv0, ThrData *const thr);
  void FLJ_thr(int ifrom, int ito, int eflag, double *pv1, ThrData *const thr);
  void TORSION_thr(int ifrom, int ito, int eflag, double *pv2, ThrData *const thr);
  void REBO_neigh_thr();
};

}    // namespace LAMMPS_NS

#endif
#endif
