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
PairStyle(body/rounded/polygon/omp,PairBodyRoundedPolygonOMP);
// clang-format on
#else

#ifndef LMP_PAIR_BODY_ROUNDED_POLYGON_OMP_H
#define LMP_PAIR_BODY_ROUNDED_POLYGON_OMP_H

#include "pair_body_rounded_polygon.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairBodyRoundedPolygonOMP : public PairBodyRoundedPolygon, public ThrOMP {

 public:
  PairBodyRoundedPolygonOMP(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData *const thr);

  void sphere_against_sphere_thr(int i, int j, double delx, double dely, double delz, double rsq,
                                 double k_n, double k_na, const double *const *v, dbl3_t *f,
                                 ThrData *thr);

  int vertex_against_edge_thr(int i, int j, double k_n, double k_na, const double *const *x,
                              dbl3_t *f, dbl3_t *torque, const tagint *tag, Contact *contact_list,
                              int &num_contacts, double &evdwl, double *facc, double *vforce_thr,
                              double *eforce_thr, ThrData *thr);

  void contact_forces_thr(Contact &contact, double j_a, const double *const *x,
                          const double *const *v, const double *const *angmom, dbl3_t *f,
                          dbl3_t *torque, double &evdwl, double *facc, double *vforce_thr,
                          double *eforce_thr, ThrData *thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
