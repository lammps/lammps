/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef BOND_TABLE_H
#define BOND_TABLE_H

#include "stdio.h"
#include "bond.h"

namespace LAMMPS_NS {

class BondTable : public Bond {
 public:
  BondTable(class LAMMPS *);
  ~BondTable();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, double, int, int);

 private:
  int tabstyle,n,nm1;
  double *r0;

  struct Table {
    int ninput,fpflag;
    double fplo,fphi,r0;
    double *rfile,*efile,*ffile;
    double *e2file,*f2file;
    double delta,invdelta,deltasq6;
    double *r,*e,*de,*f,*df,*e2,*f2;
  };

  int ntables;
  Table *tables;
  int *tabindex;
  
  void allocate();
  void null_table(Table *);
  void free_table(Table *);
  void read_table(Table *, char *, char *);
  void bcast_table(Table *);
  void spline_table(Table *);
  void compute_table(Table *);

  void param_extract(Table *, char *);
  void spline(double *, double *, int, double, double, double *);
  double splint(double *, double *, double *, int, double);

  void uf_lookup(int, double, double &, double &);
  void u_lookup(int, double, double &);
};

}

#endif
