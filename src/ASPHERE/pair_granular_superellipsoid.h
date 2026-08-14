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
   Contributing author: Jacopo Bilotto (EPFL), Jibril B. Coulibaly
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(granular/superellipsoid,PairGranularSuperellipsoid);
// clang-format on
#else

#ifndef LMP_PAIR_GRANULAR_SUPERELLIPSOID_H
#define LMP_PAIR_GRANULAR_SUPERELLIPSOID_H

#include "pair.h"

#include "atom_vec_ellipsoid.h"

namespace LAMMPS_NS {

class PairGranularSuperellipsoid : public Pair {
 public:
  PairGranularSuperellipsoid(class LAMMPS *);
  ~PairGranularSuperellipsoid() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void reset_dt() override;
  double single(int, int, int, int, double, double, double, double &) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double memory_usage() override;
  void transfer_history(double *, double *, int, int) override;

 protected:
  int freeze_group_bit;

  int neighprev;
  double *onerad_dynamic, *onerad_frozen;
  double *maxrad_dynamic, *maxrad_frozen;

  class FixDummy *fix_dummy;
  class FixNeighHistory *fix_history;

  // storage of rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, null pointer if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  // Model variables
  double dt;
  int **normal_model;
  int **damping_model;
  int **tangential_model;
  int **limit_damping;
  int default_hist_size;
  int contact_radius_flag;

  // Normal coefficients
  double **kn, **gamman;    // Hooke + Hertz

  // Tangential coefficients
  double **kt, **xt, **xmu;    // linear_history

  // Intermediate values for contact model
  int history_update, touchjj, itype, jtype;
  double Fnormal, forces[3], torquesi[3], torquesj[3];
  double radi, radj, meff, Fntot, contact_radius;
  double *xi, *xj, *vi, *vj;
  double fs[3], ft[3];
  double dx[3], nx[3], r, rsq, rinv, Reff, radsum, delta, dR;
  double vr[3], vn[3], vnnr, vt[3], wr[3], vtr[3], vrel;

  double *quati, *quatj, *angmomi, *angmomj, *inertiai, *inertiaj;
  double X0[4], nij[3], Ri[3][3], Rj[3][3];
  double shapei0[3], blocki0[3], shapej0[3], blockj0[3], shapei[3], blocki[3], shapej[3], blockj[3];
  double *history_data, *xref;
  AtomVecEllipsoid::BlockType flagi, flagj;
  tagint tagi, tagj;

  void allocate();
  double mix_geom(double, double);
  double mix_mean(double, double);
  bool check_contact();
  void calculate_forces();

 private:
  int size_history;
  int heat_flag;

  // optional user-specified global cutoff, per-type user-specified cutoffs
  double **cutoff_type;
  double cutoff_global;
  int contact_formulation;
  int bounding_box;
  int curvature_model;

  int extra_svector;

  void rotate_rescale_vec(double *hislocal, double *n);

  // Below not implemented. Placeholder if we decide not to compute local hessian in line search
  static double
  shape_and_gradient_local(const double *, const double *, const double *,
                           double *);    // would return a vector of temporary variables
  static double hessian_local(
      const double *, const double *, const double *,
      double *);    // would use the above vector of temporary variables to compute local hessian
};

}    // namespace LAMMPS_NS

#endif
#endif
