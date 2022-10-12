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
FixStyle(polarize/bem/gmres,FixPolarizeBEMGMRES);
// clang-format on
#else

#ifndef LMP_FIX_POLARIZE_BEM_GMRES_H
#define LMP_FIX_POLARIZE_BEM_GMRES_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPolarizeBEMGMRES : public Fix {
 public:
  FixPolarizeBEMGMRES(class LAMMPS *, int, char **);
  ~FixPolarizeBEMGMRES() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void pre_force(int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  double compute_vector(int) override;

  int modify_param(int, char **) override;
  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  virtual void allocate();
  virtual void deallocate();

  virtual void compute_induced_charges();
  void set_dielectric_params(double, double, double, double, int, double);

  class AtomVecDielectric *avec;

 protected:
  int nmax;
  int nevery;    // to be invoked every time steps
  int *
      tag2mat;    // tag2mat[atom->tag[i]] = the index of the atoms in the induced charge arrays from atom tags
  int *
      mat2tag;    // mat2tag[idx] = the atom tag of the induced charge idx in the induced charge arrays
  int *induced_charge_idx;    // return the index of the atoms in the induced charge arrays
  int num_induced_charges;    // total number of induced charges
  double *induced_charges;    // values of induced charges
  double *buffer;             // buffer of size num_induced_charges
  double *q_backup;           // backup for the real charges
  int allocated;
  double **efield_pair;      // electrical field at position of atom i due to pair contribution
  double **efield_kspace;    // electrical field at position of atom i due to kspace contribution
  int kspaceflag;            // 1 if kspace is used for the induced charge computation
  int torqueflag, extraflag;

  void force_clear();
  double vec_dot(const double *, const double *,
                 int);    // dot product between two vectors of length n

 private:
  int mat_dim;          // matrix dimension = total number of induced charges
  int mr;               // number of vectors used to span the Krylov space
  int iterations;       // actual number of iterations
  int itr_max;          // maximum number of outer iterations
  int randomized;       // 1 if generating random induced charges, 0 otherwise
  double ave_charge;    // average random charge
  int seed_charge;
  double epsilon0e2q;    // convert epsilon0 times efield to unit of charge per area

  double *c, *g, *h, *r, *s, *v, *y;    // vectors used by the solver
  double *rhs;                          // right-hand side vector of the equation Ax = b
  double tol_abs, tol_rel;              // tolerance for convergence
  double normb;                         // norm of the rhs vector b
  double rho;                           // norm of (b - Ax)
  int first;    // 1 if first time invoked (initializing induced charges with zero)

  void gmres_solve(double *, double *);            // GMRES workhorse
  void apply_operator(double *, double *, int);    // compute Ax without explicitly storing A
  void update_residual(double *, double *,
                       int);    // compute (b - Ax) directly (without computing b and Ax explcitly)

  // Givens rotations
  inline void mult_givens(double c, double s, int k, double *g)
  {
    double g1 = c * g[k] - s * g[k + 1];
    double g2 = s * g[k] + c * g[k + 1];
    g[k] = g1;
    g[k + 1] = g2;
  }
};

}    // namespace LAMMPS_NS

#endif
#endif
