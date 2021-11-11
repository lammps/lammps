/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This pair style is written by Daniele Rapetti (iximiel@gmail.com)
   ------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(smatb/single,PairSMATBSingle)
// clang-format on
#else

#ifndef LMP_PAIR_SMATBSINGLE_H
#define LMP_PAIR_SMATBSINGLE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSMATBSingle : public Pair {
 public:
  PairSMATBSingle(class LAMMPS *);
  virtual ~PairSMATBSingle();
  void compute(int, int);    //workhorse routine that computes pairwise interactions
  /*
    void compute_inner();
    void compute_middle();
    void compute_outer(int, int);
    */
  void settings(int, char **);    //reads the input script line with arguments you define
  void coeff(int, char **);       //set coefficients for one i,j type pair
  //void init_style();//initialization specific to this pair style
  double init_one(int, int);    //perform initialization for one i,j type pair

  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual void write_data(FILE *);
  virtual void write_data_all(FILE *);

  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  virtual int pack_reverse_comm(int, int, double *);
  virtual void unpack_reverse_comm(int, int *, double *);

 protected:
  virtual void allocate();
  int nmax;                         // allocated size of per-atom arrays
  double *on_eb;                    //allocated to store up caclulation values
  double r0;                        // interaction radius, user-given
  double p, A, q, QSI;              // parameters user-given
  double cutOffStart, cutOffEnd;    //cut offs, user given
  double cutOffEnd2;                //squared cut off end, calculated
  double a3, a4, a5;                //polynomial for cutoff linking to zero:   Ae^p substitution
  double x3, x4, x5;                //polynomial for cutoff linking to zero: QSIe^q substitution
};

}    // namespace LAMMPS_NS

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
