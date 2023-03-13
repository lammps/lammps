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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(balance,Balance);
// clang-format on
#else

#ifndef LMP_BALANCE_H
#define LMP_BALANCE_H

#include "command.h"

namespace LAMMPS_NS {

class Balance : public Command {
 public:
  class RCB *rcb;
  class FixStorePeratom *fixstore;    // per-atom weights stored in FixStorePeratom
  int wtflag;                         // 1 if particle weighting is used
  int varflag;                        // 1 if weight style var(iable) is used
  int sortflag;                       // 1 if sorting of comm messages is done
  int outflag;                        // 1 for output of balance results to file

  Balance(class LAMMPS *);
  ~Balance() override;
  void command(int, char **) override;
  void options(int, int, char **, int);
  void weight_storage(char *);
  void init_imbalance(int);
  void set_weights();
  double imbalance_factor(double &);
  void shift_setup(char *, int, double);
  int shift();
  int *bisection();
  void dumpout(bigint);

  static constexpr int BSTR_SIZE = 3;

 private:
  int me, nprocs;

  double thresh;                                      // threshold to perform LB
  int style;                                          // style of LB
  int xflag, yflag, zflag;                            // xyz LB flags
  double *user_xsplit, *user_ysplit, *user_zsplit;    // params for xyz LB
  int oldrcb;                                         // use old-style RCB compute

  int nitermax;    // params for shift LB
  double stopthresh;
  char bstr[BSTR_SIZE + 1];

  int shift_allocate;       // 1 if SHIFT vectors have been allocated
  int ndim;                 // length of balance string bstr
  int *bdim;                // XYZ for each character in bstr
  double *onecost;          // work vector of counts in one dim
  double *allcost;          // counts for slices in one dim
  double *sum;              // cumulative count for slices in one dim
  double *target;           // target sum for slices in one dim
  double *lo, *hi;          // lo/hi split coords that bound each target
  double *losum, *hisum;    // cumulative counts at lo/hi coords
  int rho;                  // 0 for geometric recursion
                            // 1 for density weighted recursion

  double *proccost;       // particle cost per processor
  double *allproccost;    // proccost summed across procs

  int nimbalance;                  // number of user-specified weight styles
  class Imbalance **imbalances;    // list of Imb classes, one per weight style
  double *weight;                  // ptr to FixStore weight vector

  FILE *fp;    // balance output file
  int firststep;

  double imbalance_splits();
  void shift_setup_static(char *);
  void tally(int, int, double *);
  int adjust(int, double *);
#ifdef BALANCE_DEBUG
  void debug_shift_output(int, int, int, double *);
#endif
};
}    // namespace LAMMPS_NS
#endif
#endif
