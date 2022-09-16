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
PairStyle(eam,PairEAM);
// clang-format on
#else

#ifndef LMP_PAIR_EAM_H
#define LMP_PAIR_EAM_H

#include "pair.h"

namespace LAMMPS_NS {

class PairEAM : public Pair {
 public:
  friend class FixSemiGrandCanonicalMC;    // Alex Stukowski option

  // public variables so ATC package can access them

  double cutmax;

  // potentials as array data

  int nrho, nr;
  int nfrho, nrhor, nz2r;
  double **frho, **rhor, **z2r;
  int *type2frho, **type2rhor, **type2z2r;

  // potentials in spline form used for force computation

  double dr, rdr, drho, rdrho, rhomax, rhomin;
  double ***rhor_spline, ***frho_spline, ***z2r_spline;

  PairEAM(class LAMMPS *);
  ~PairEAM() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;
  void *extract_peratom(const char *, int &) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;
  void swap_eam(double *, double **) override;

 protected:
  int nmax;    // allocated size of per-atom arrays
  double cutforcesq;
  double **scale;
  bigint embedstep;    // timestep, the embedding term was computed

  // per-atom arrays

  double *rho, *fp;
  int *numforce;

  // potentials as file data

  struct Funcfl {
    char *file;
    int nrho, nr;
    double drho, dr, cut, mass;
    double *frho, *rhor, *zr;
  };
  Funcfl *funcfl;
  int nfuncfl;

  struct Setfl {
    char **elements;
    int nelements, nrho, nr;
    double drho, dr, cut;
    double *mass;
    double **frho, **rhor, ***z2r;
  };
  Setfl *setfl;

  struct Fs {
    char **elements;
    int nelements, nrho, nr;
    double drho, dr, cut;
    double *mass;
    double **frho, ***rhor, ***z2r;
  };
  Fs *fs;

  virtual void allocate();
  virtual void array2spline();
  void interpolate(int, double, double *, double **);

  virtual void read_file(char *);
  virtual void file2array();
};

}    // namespace LAMMPS_NS

#endif
#endif
