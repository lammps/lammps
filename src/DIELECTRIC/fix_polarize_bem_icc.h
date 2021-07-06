/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(polarize/bem/icc,FixPolarizeBEMICC);
// clang-format on
#else

#ifndef LMP_FIX_POLARIZE_BEM_ICC_H
#define LMP_FIX_POLARIZE_BEM_ICC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPolarizeBEMICC : public Fix {
 public:
  FixPolarizeBEMICC(class LAMMPS *, int, char **);
  ~FixPolarizeBEMICC();
  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void pre_force(int);
  virtual double compute_vector(int);
  int modify_param(int, char **);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);

  virtual void compute_induced_charges();
  void set_dielectric_params(double, double, double, double, int, double);

  class AtomVecDielectric *avec;

 protected:
  int nevery;                // to be invoked every time steps
  double **efield_pair;      // electrical field at position of atom i due to pair contribution
  double **efield_kspace;    // electrical field at position of atom i due to kspace contribution
  int kspaceflag;            // 1 if kspace is used for the induced charge computation
  int torqueflag, extraflag;

  void force_clear();

 private:
  int iterations;             // actual number of iterations
  int itr_max;                // maximum number of outer iterations
  double tol_abs, tol_rel;    // tolerance for convergence
  double rho;                 // current error
  double omega;               // iterative weight
  int randomized;             // 1 if generating random induced charges, 0 otherwise
  double ave_charge;          // average random charge
  int seed_charge;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
