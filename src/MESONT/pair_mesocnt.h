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
   Contributing author: Philipp Kloza (University of Cambridge)
                        pak37@cam.ac.uk
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(mesocnt, PairMesoCNT);
// clang-format on
#else

#ifndef LMP_PAIR_MESOCNT_H
#define LMP_PAIR_MESOCNT_H

#include "pair.h"

namespace LAMMPS_NS {
class PotentialFileReader;
class PairMesoCNT : public Pair {
 public:
  PairMesoCNT(class LAMMPS *);
  ~PairMesoCNT() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 protected:
  int nend_types;
  int uinf_points, gamma_points, phi_points, usemi_points;
  int *end_types, *reduced_nlist, *numchainlist, *selfid;
  int **reduced_neighlist, **nchainlist, **endlist, **selfpos;
  int ***chainlist;

  bool segment_flag, neigh_flag;

  double ang, ang_inv, eunit, funit;
  double delta1, delta2;
  double r, rsq, d, rc, rcsq, rc0, cutoff, cutoffsq, neigh_cutoff;
  double r_ang, rsq_ang, d_ang, rc_ang, rcsq_ang, cutoff_ang, cutoffsq_ang;
  double sig, sig_ang, comega, ctheta;
  double hstart_uinf, hstart_gamma, hstart_phi, psistart_phi, hstart_usemi, xistart_usemi;
  double delh_uinf, delh_gamma, delh_phi, delpsi_phi, delh_usemi, delxi_usemi;

  double p1[3], p2[3], p[3], m[3];
  double *param, *w, *wnode;
  double **dq_w;
  double ***q1_dq_w, ***q2_dq_w;
  double *gl_nodes_finf, *gl_nodes_fsemi;
  double *gl_weights_finf, *gl_weights_fsemi;
  double *uinf_data, *gamma_data, **phi_data, **usemi_data;
  double **uinf_coeff, **gamma_coeff, ****phi_coeff, ****usemi_coeff;
  double **flocal, **fglobal, **basis;

  bool match_end(int);

  int count_chains(int *, int);

  void allocate();
  void bond_neigh_id();
  void bond_neigh_topo();
  void neigh_common(int, int, int &, int *);
  void chain_lengths(int *, int, int *);
  void chain_split(int *, int, int *, int **, int *);
  void sort(int *, int);

  void read_file(const char *);
  void read_data(PotentialFileReader &, double *, double &, double &, int);
  void read_data(PotentialFileReader &, double **, double &, double &, double &, double &, int);

  void spline_coeff(double *, double **, double, int);
  void spline_coeff(double **, double ****, double, double, int);

  void geometry(const double *, const double *, const double *, const double *, const double *,
                double *, double *, double *, double **);

  void finf(const double *, double &, double **);
  void fsemi(const double *, double &, double &, double **);

  // Legendre-Gauss integration

  double legendre(int, double);
  void gl_init_nodes(int, double *);
  void gl_init_weights(int, double *, double *);

  // inlined functions for efficiency

  inline void weight(const double *, const double *, const double *, const double *, double &,
                     double *, double *, double *, double *);

  inline double spline(double, double, double, double **, int);
  inline double dspline(double, double, double, double **, int);
  inline double spline(double, double, double, double, double, double, double ****, int);
  inline double dxspline(double, double, double, double, double, double, double ****, int);
  inline double dyspline(double, double, double, double, double, double, double ****, int);

  inline double heaviside(double x)
  {
    if (x > 0)
      return 1.0;
    else
      return 0.0;
  }

  inline double s(double x)
  {
    return heaviside(-x) + heaviside(x) * heaviside(1 - x) * (1 - x * x * (3 - 2 * x));
  }

  inline double ds(double x) { return 6 * heaviside(x) * heaviside(1 - x) * x * (x - 1); }

  inline double s5(double x)
  {
    double x2 = x * x;
    return heaviside(-x) + heaviside(x) * heaviside(1 - x) * (1 - x2 * x * (6 * x2 - 15 * x + 10));
  }

  inline double ds5(double x)
  {
    double x2 = x * x;
    return -30 * heaviside(x) * heaviside(1 - x) * x2 * (x2 - 2 * x + 1);
  }
};

}    // namespace LAMMPS_NS

#endif
#endif
