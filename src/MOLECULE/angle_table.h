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

#ifndef ANGLE_TABLE_H
#define ANGLE_TABLE_H

#include "stdio.h"
#include "angle.h"

namespace LAMMPS_NS {

class AngleTable : public Angle {
 public:
  AngleTable(class LAMMPS *);
  ~AngleTable();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int);

 private:
  int tabstyle,n,nm1;
  double *theta0;

  struct Table {
    int ninput,fpflag;
    double fplo,fphi,theta0;
    double *afile,*efile,*ffile;
    double *e2file,*f2file;
    double delta,invdelta,deltasq6;
    double *ang,*e,*de,*f,*df,*e2,*f2;
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
