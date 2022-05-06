/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(adp,PairADP);
// clang-format on
#else

#ifndef LMP_PAIR_ADP_H
#define LMP_PAIR_ADP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairADP : public Pair {
 public:
  PairADP(class LAMMPS *);
  ~PairADP() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;

 protected:
  int nmax;    // allocated size of per-atom arrays
  double cutforcesq, cutmax;

  // per-atom arrays

  double *rho, *fp;
  double **mu, **lambda;

  // potentials as array data

  int nrho, nr;
  int nfrho, nrhor, nz2r;
  int nu2r, nw2r;
  double **frho, **rhor, **z2r;
  double **u2r, **w2r;
  int *type2frho, **type2rhor, **type2z2r;
  int **type2u2r, **type2w2r;

  // potentials in spline form used for force computation

  double dr, rdr, drho, rdrho;
  double ***rhor_spline, ***frho_spline, ***z2r_spline;
  double ***u2r_spline, ***w2r_spline;

  // potentials as file data

  struct Setfl {
    char **elements;
    int nelements, nrho, nr;
    double drho, dr, cut;
    double *mass;
    double **frho, **rhor, ***z2r;
    double ***u2r, ***w2r;
  };
  Setfl *setfl;

  void allocate();
  virtual void array2spline();
  void interpolate(int, double, double *, double **);

  void read_file(char *);
  virtual void file2array();
};

}    // namespace LAMMPS_NS

#endif
#endif
