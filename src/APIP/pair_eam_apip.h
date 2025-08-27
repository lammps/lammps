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
   Contributing authors: Stephen Foiles (SNL), Murray Daw (SNL) (EAM)
     David Immel (d.immel@fz-juelich.de, FZJ, Germany) (APIP)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(eam/apip,PairEAMAPIP);
// clang-format on
#else

#ifndef LMP_PAIR_EAM_APIP_H
#define LMP_PAIR_EAM_APIP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairEAMAPIP : public Pair {
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

  PairEAMAPIP(class LAMMPS *);
  ~PairEAMAPIP() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void setup() override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void *extract(const char *, int &) override;

  double memory_usage() override;
  void swap_eam(double *, double **) override;

 protected:
  int nmax;    // allocated size of per-atom arrays
  double cutforcesq;
  double **scale;

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

  // stats required for load balancing
  int n_non_complex_accumulated;    // number of calculated atoms
  double time_wall_accumulated;     // time required for the atom calculation
  double time_per_atom;             // time for one calculation

  void calculate_time_per_atom();

  bool lambda_thermostat;    // true/false there is one/no fix lambda_thermostat
};

}    // namespace LAMMPS_NS

#endif
#endif
