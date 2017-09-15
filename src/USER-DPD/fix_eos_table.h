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

#ifdef FIX_CLASS

FixStyle(eos/table,FixEOStable)

#else

#ifndef LMP_FIX_EOS_TABLE_H
#define LMP_FIX_EOS_TABLE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEOStable : public Fix {
 public:
  FixEOStable(class LAMMPS *, int, char **);
  virtual ~FixEOStable();
  int setmask();
  virtual void init();
  virtual void post_integrate();
  virtual void end_of_step();
  void energy_lookup(double, double &);
  void temperature_lookup(double, double &);

 protected:
  enum{LINEAR};

  int tabstyle,tablength;
  struct Table {
    int ninput;
    double lo,hi;
    double *rfile,*efile;
    double *e2file;
    double delta,invdelta,deltasq6;
    double *r,*e,*de,*e2;
  };
  int ntables;
  Table *tables;

  void allocate();
  void null_table(Table *);
  void free_table(Table *);
  void read_table(Table *, Table *, char *, char *);
  void bcast_table(Table *);
  void spline_table(Table *);
  void compute_table(Table *);

  void param_extract(Table *, Table *, char *);
  void spline(double *, double *, int, double, double, double *);
  double splint(double *, double *, double *, int, double);

  };
}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Unknown table style in fix eos/table

Style of table is invalid for use with fix eos/table command.

E: Illegal number of eos/table entries

There must be at least 2 table entries.

E: Invalid eos/table length

Length of read-in fix eos/table is invalid

E: eos/table values are not increasing

The EOS must be a monotonically, increasing function

E:  FixEOStable requires atom_style with internal temperature and energies (e.g. dpd)

Self-explanatory.

E: Internal temperature < zero

Self-explanatory.  EOS may not be valid under current simulation conditions.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Did not find keyword in table file

Keyword used in fix eos/table command was not found in table file.

E: Invalid keyword in fix eos/table parameters

Keyword used in list of table parameters is not recognized.

E: fix eos/table parameters did not set N

List of fix eos/table parameters must include N setting.

E: Temperature is not within table cutoffs

The internal temperature does not lie with the minimum
and maximum temperature cutoffs of the table

E: Energy is not within table cutoffs

The internal energy does not lie with the minimum
and maximum energy cutoffs of the table

*/
