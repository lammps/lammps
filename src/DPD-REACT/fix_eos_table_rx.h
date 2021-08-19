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

#ifdef FIX_CLASS
// clang-format off
FixStyle(eos/table/rx,FixEOStableRX);
// clang-format on
#else

#ifndef LMP_FIX_EOS_TABLE_RX_H
#define LMP_FIX_EOS_TABLE_RX_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEOStableRX : public Fix {
 public:
  FixEOStableRX(class LAMMPS *, int, char **);
  virtual ~FixEOStableRX();
  int setmask();
  void setup(int);
  virtual void init();
  virtual void post_integrate();
  virtual void end_of_step();
  void energy_lookup(int, double, double &);
  void temperature_lookup(int, double, double &);

 protected:
  enum { LINEAR };

  int tabstyle, tablength;
  struct Table {
    int ninput;
    double lo, hi;
    double *rfile, *efile;
    double *e2file;
    double delta, invdelta, deltasq6;
    double *r, *e, *de, *e2;
  };
  int ntables;
  Table *tables, *tables2;

  void allocate();
  void null_table(Table *);
  void free_table(Table *);
  void read_table(Table *, Table *, char *, char *);
  void bcast_table(Table *);
  void spline_table(Table *);
  void compute_table(Table *);

  void param_extract(Table *, char *);
  void spline(double *, double *, int, double, double, double *);
  double splint(double *, double *, double *, int, double);

  int nspecies;

  void read_file(char *);

  double *dHf, *energyCorr, *tempCorrCoeff, *moleculeCorrCoeff;

  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);

  int *eosSpecies;
  int ncolumn;
  bool rx_flag;
};
}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: FixEOStableRX requires a fix rx command.

The fix rx command must come before the pair style command in the input file

E:  There are no rx species specified

There must be at least one species specified through the fix rx command

E:  Invalid eos/table/rx length

The eos/table/rx table must have more than one entry.

E:  eos/table/rx values are not increasing

The equation-of-state must an increasing function

E:  FixEOStableRX requires atom_style with internal temperature and energies (e.g. dpd)

Self-explanatory.

E:  Internal temperature <= zero.

Self-explanatory.

E:  Cannot open eos table/rx potential file %s

Self-explanatory.

E:  Incorrect format in eos table/rx file

Self-explanatory.

E:  Cannot open file %s

Self-explanatory.

E:  Did not find keyword in table file

Self-explanatory.

E:  Illegal fix eos/table/rx command

Incorrect number of arguments specified for the fix eos/table/rx command.

E:  Invalid keyword in fix eos/table/rx parameters

Self-explanatory.

E:  The number of columns in fix eos/table/rx does not match the number of species.

Self-explanatory.  Check format for fix eos/table/rx file.

E:  fix eos/table/rx parameters did not set N

The number of table entries was not set in the eos/table/rx file

W:  Secant solver did not converge because table bounds were exceeded

The secant solver failed to converge, resulting in the lower or upper table bound temperature to be returned

E: NaN detected in secant solver.

Self-explanatory.

E: Maxit exceeded in secant solver

The maximum number of iterations was exceeded in the secant solver

*/
