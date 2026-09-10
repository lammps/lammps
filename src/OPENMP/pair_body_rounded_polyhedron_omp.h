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

namespace LAMMPS_NS {

class PairBodyRoundedPolyhedronOMP : public PairBodyRoundedPolyhedron, public ThrOMP {

 public:
  PairBodyRoundedPolyhedronOMP(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

 private:
  template <int EVFLAG, int EFLAG>
  void eval(int ifrom, int ito, ThrData *const thr);

  void sphere_against_sphere_thr(int ibody, int jbody, int itype, int jtype,
                                 double delx, double dely, double delz, double rsq,
                                 double **v, dbl3_t *f, int evflag, ThrData *thr);

  void sphere_against_edge_thr(int ibody, int jbody, int itype, int jtype, double **x,
                               double **v, dbl3_t *f, dbl3_t *torque, double **angmom,
                               double *vflag_thr, int evflag, ThrData *thr);

  void sphere_against_face_thr(int ibody, int jbody, int itype, int jtype, double **x,
                               double **v, dbl3_t *f, dbl3_t *torque, double **angmom,
                               int evflag, ThrData *thr);

  int edge_against_edge_thr(int ibody, int jbody, int itype, int jtype, double **x,
                            Contact *contact_list, int &num_contacts, double &evdwl,
                            double *facc, dbl3_t *f, dbl3_t *torque, double **v,
                            double **angmom);

  int edge_against_face_thr(int ibody, int jbody, int itype, int jtype, double **x,
                            Contact *contact_list, int &num_contacts, double &evdwl,
                            double *facc, dbl3_t *f, dbl3_t *torque, double **v,
                            double **angmom, double *vflag_thr);

  int interaction_face_to_edge_thr(int ibody, int face_index, double *xmi,
                                   double rounded_radius_i, int jbody, int edge_index,
                                   double *xmj, double rounded_radius_j, int itype,
                                   int jtype, double cut_inner, Contact *contact_list,
                                   int &num_contacts, double &energy, double *facc,
                                   dbl3_t *f, dbl3_t *torque, double **x, double **v,
                                   double **angmom, double *vflag_thr);

  int interaction_edge_to_edge_thr(int ibody, int edge_index_i, double *xmi,
                                   double rounded_radius_i, int jbody, int edge_index_j,
                                   double *xmj, double rounded_radius_j, int itype,
                                   int jtype, double cut_inner, Contact *contact_list,
                                   int &num_contacts, double &energy, double *facc,
                                   dbl3_t *f, dbl3_t *torque, double **x, double **v,
                                   double **angmom);

  void pair_force_and_torque_thr(int ibody, int jbody, double *pi, double *pj, double r,
                                 double contact_dist, int itype, int jtype, double **x,
                                 double **v, dbl3_t *f, dbl3_t *torque, double **angmom,
                                 int jflag, double &energy, double *facc);

  void contact_forces_thr(int ibody, int jbody, double *xi, double *xj, double delx,
                          double dely, double delz, double fx, double fy, double fz,
                          double **x, double **v, double **angmom, dbl3_t *f,
                          dbl3_t *torque, double *facc);

  void rescale_cohesive_forces_thr(double **x, dbl3_t *f, dbl3_t *torque,
                                   Contact *contact_list, int &num_contacts, int itype,
                                   int jtype, double *facc);
};

}    // namespace LAMMPS_NS

#endif
#endif
