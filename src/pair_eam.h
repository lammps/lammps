/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef PAIR_EAM_H
#define PAIR_EAM_H

#include "pair.h"

class PairEAM : public Pair {
  friend class FixEnergy;

 public:
  double *rho,*fp;
  int nmax;

  PairEAM();
  virtual ~PairEAM();
  void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void single(int, int, int, int, double, double, double, int, One &);

  int pack_comm(int, int *, double *, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int memory_usage();

 protected:
  double cutforcesq,cutmax;

  // potentials as file data

  int *map;                   // which element each atom type maps to

  struct Funcfl {
    char *file;
    int nrho,nr;
    double drho,dr,cut,mass;
    double *frho,*rhor,*zr;
  };
  Funcfl *funcfl;
  int nfuncfl;

  struct Setfl {
    char **elements;
    int nelements,nrho,nr;
    double drho,dr,cut;
    double *mass;
    double **frho,**rhor,***z2r;
  };
  Setfl *setfl;

  struct Fs {
    char **elements;
    int nelements,nrho,nr;
    double drho,dr,cut;
    double *mass;
    double **frho,***rhor,***z2r;
  };
  Fs *fs;

  // potentials as array data

  int nrho,nr;
  int nfrho,nrhor,nz2r;
  double **frho,**rhor,**z2r;
  int *type2frho,**type2rhor,**type2z2r;
  
  // potentials in spline form used for force computation

  double dr,rdr,drho,rdrho;
  double ***rhor_spline,***frho_spline,***z2r_spline;

  void allocate();
  void array2spline();
  void interpolate(int, double, double *, double **);
  void grab(FILE *, int, double *);
  void single_rho(int, int, double, double &, double &);
  void single_embed(int, int, double &, int, double &);

  virtual void read_file(char *);
  virtual void file2array();
};

#endif
