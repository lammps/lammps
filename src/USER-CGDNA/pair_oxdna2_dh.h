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

#ifdef PAIR_CLASS

PairStyle(oxdna2/dh,PairOxdna2Dh)

#else

#ifndef LMP_PAIR_OXDNA2_DH_H
#define LMP_PAIR_OXDNA2_DH_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOxdna2Dh : public Pair {
 public:
  PairOxdna2Dh(class LAMMPS *);
  virtual ~PairOxdna2Dh();
  virtual void compute_interaction_sites(double *, double *, double *);
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  void *extract(const char *, int &);

 protected:

  double **qeff_dh_pf,**kappa_dh;
  double **b_dh,**cut_dh_ast,**cutsq_dh_ast,**cut_dh_c,**cutsq_dh_c;

  virtual void allocate();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
