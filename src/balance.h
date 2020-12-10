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

#ifdef COMMAND_CLASS

CommandStyle(balance,Balance)

#else

#ifndef LMP_BALANCE_H
#define LMP_BALANCE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Balance : protected Pointers {
 public:
  class RCB *rcb;
  class FixStore *fixstore;       // per-atom weights stored in FixStore
  int wtflag;                     // 1 if particle weighting is used
  int varflag;                    // 1 if weight style var(iable) is used
  int outflag;                    // 1 for output of balance results to file

  Balance(class LAMMPS *);
  ~Balance();
  void command(int, char **);
  void options(int, int, char **);
  void weight_storage(char *);
  void init_imbalance(int);
  void set_weights();
  double imbalance_factor(double &);
  void shift_setup(char *, int, double);
  int shift();
  int *bisection(int sortflag = 0);
  void dumpout(bigint);

 private:
  int me,nprocs;

  double thresh;                                    // threshold to perform LB
  int style;                                        // style of LB
  int xflag,yflag,zflag;                            // xyz LB flags
  double *user_xsplit,*user_ysplit,*user_zsplit;    // params for xyz LB
  int oldrcb;                                    // use old-style RCB compute

  int nitermax;              // params for shift LB
  double stopthresh;
  char bstr[4];

  int shift_allocate;        // 1 if SHIFT vectors have been allocated
  int ndim;                  // length of balance string bstr
  int *bdim;                 // XYZ for each character in bstr
  double *onecost;           // work vector of counts in one dim
  double *allcost;           // counts for slices in one dim
  double *sum;               // cumulative count for slices in one dim
  double *target;            // target sum for slices in one dim
  double *lo,*hi;            // lo/hi split coords that bound each target
  double *losum,*hisum;      // cumulative counts at lo/hi coords
  int rho;                   // 0 for geometric recursion
                             // 1 for density weighted recursion

  double *proccost;          // particle cost per processor
  double *allproccost;       // proccost summed across procs

  int nimbalance;                 // number of user-specified weight styles
  class Imbalance **imbalances;   // list of Imb classes, one per weight style
  double *weight;                 // ptr to FixStore weight vector

  FILE *fp;                  // balance output file
  int firststep;

  double imbalance_splits();
  void shift_setup_static(char *);
  void tally(int, int, double *);
  int adjust(int, double *);
  int binary(double, int, double *);
#ifdef BALANCE_DEBUG
  void debug_shift_output(int, int, int, double *);
#endif
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Balance command before simulation box is defined

The balance command cannot be used before a read_data, read_restart,
or create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot balance in z dimension for 2d simulation

Self-explanatory.

E: Balance shift string is invalid

The string can only contain the characters "x", "y", or "z".

E: Balance rcb cannot be used with comm_style brick

Comm_style tiled must be used instead.

E: Lost atoms via balance: original %ld current %ld

This should not occur.  Report the problem to the developers.

E: Unknown (fix) balance weight method

UNDOCUMENTED

E: Cannot open (fix) balance output file

UNDOCUMENTED

E: Balance produced bad splits

This should not occur.  It means two or more cutting plane locations
are on top of each other or out of order.  Report the problem to the
developers.

U: Cannot open balance output file

Self-explanatory.

*/
