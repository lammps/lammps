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
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  virtual double init_one(int, int);
  virtual void init_style();
  virtual void single(int, int, int, int, double, double, double, int, One &);

  int pack_comm(int, int *, double *, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int memory_usage();

 protected:
  double cutforcesq,cutmax;
  int **tabindex;

  // potential as read in from file

  struct Table {
    char *filename;           // file it was read from
    int ith,jth;              // for setfl, which i,j entry in file
    int nrho,nr;              // array lengths
    double drho,dr,cut,mass;  // array spacings and cutoff, mass
    double *frho,*rhor;       // array values
    double *zr,*z2r;          // zr set for funcfl, z2r set for setfl
  };
  int ntables;
  Table *tables;

  // potential stored in multi-type setfl format
  // created from read-in tables

  int nrho,nr;
  double drho,dr;
  double **frho,**rhor,**zrtmp;
  double ***z2r;

  // potential in spline form used for force computation
  // created from multi-type setfl format by interpolate()

  double rdr,rdrho;
  double **rhor_0,**rhor_1,**rhor_2,**rhor_3,**rhor_4,**rhor_5,**rhor_6;
  double **frho_0,**frho_1,**frho_2,**frho_3,**frho_4,**frho_5,**frho_6;
  double ***z2r_0,***z2r_1,***z2r_2,***z2r_3,***z2r_4,***z2r_5,***z2r_6;

  void allocate();
  int read_funcfl(char *);
  void convert_funcfl();
  virtual void interpolate();
  virtual void interpolate_deallocate();
  void grab(FILE *, int, double *);
  void skip(FILE *, int);
  void single_embed(int, int, double &, int, double &);
};

#endif
