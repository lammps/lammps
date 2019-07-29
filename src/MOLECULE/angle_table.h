/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS

AngleStyle(table,AngleTable)

#else

#ifndef LMP_ANGLE_TABLE_H
#define LMP_ANGLE_TABLE_H

#include <cstdio>
#include "angle.h"

namespace LAMMPS_NS {

class AngleTable : public Angle {
 public:
  AngleTable(class LAMMPS *);
  virtual ~AngleTable();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int);

 protected:
  int tabstyle,tablength;
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
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Unknown table style in angle style table

Self-explanatory.

E: Illegal number of angle table entries

There must be at least 2 table entries.

E: Invalid angle table length

Length must be 2 or greater.

E: Angle table must range from 0 to 180 degrees

Self-explanatory.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

E: Did not find keyword in table file

Keyword used in pair_coeff command was not found in table file.

E: Invalid keyword in angle table parameters

Self-explanatory.

E: Angle table parameters did not set N

List of angle table parameters must include N setting.

E: Illegal angle in angle style table

UNDOCUMENTED

*/
