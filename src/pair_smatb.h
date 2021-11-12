/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This style is written by Daniele Rapetti (iximiel@gmail.com)
   ------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(smatb,PairSMATB)
// clang-format on
#else

#ifndef LMP_PAIR_SMATB_H
#define LMP_PAIR_SMATB_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSMATB : public Pair {
 public:
  PairSMATB(class LAMMPS *);
  ~PairSMATB() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  virtual void allocate();
  // allocated size of per-atom arrays
  int nmax;
  //allocated to store up calculation values
  double *on_eb{nullptr};
  // interaction radius, user-given
  double **r0{nullptr};
  // parameters user-given
  double **p{nullptr};
  double **A{nullptr};
  double **q{nullptr};
  double **QSI{nullptr};
  //extremes of the cut off, user given
  double **cutOffStart{nullptr};
  double **cutOffEnd{nullptr};
  //squared cut off end, calculated
  double **cutOffEnd2{nullptr};
  //polynomial for cutoff linking to zero:   Ae^p substitution
  double **a3{nullptr};
  double **a4{nullptr};
  double **a5{nullptr};
  //polynomial for cutoff linking to zero: QSIe^q substitution
  double **x3{nullptr};
  double **x4{nullptr};
  double **x5{nullptr};
  /* latex form of the potential (R_c is cutOffEnd, \Xi is QSI):

       E_i =
       \sum_{j,R_{ij}\leq R_c} A  e^{-p \left(\frac{R_{ij}}{R_{0}}-1\right)}
       -\sqrt{\sum_{j,R_{ij}\leq R_c}\Xi^2 e^{-2q\left(\frac{R_{ij}}{R_{0}}-1\right)}}.

       NB::this form does not have the polynomial link to 0 for the cut off
     */
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
