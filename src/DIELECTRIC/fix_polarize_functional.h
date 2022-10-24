/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(polarize/functional,FixPolarizeFunctional);
// clang-format on
#else

#ifndef LMP_FIX_POLARIZE_FUNCTIONAL_H
#define LMP_FIX_POLARIZE_FUNCTIONAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPolarizeFunctional : public Fix {
 public:
  FixPolarizeFunctional(class LAMMPS *, int, char **);
  ~FixPolarizeFunctional() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup(int) override;
  void setup_pre_force(int vflag) override;
  void pre_force(int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  int modify_param(int, char **) override;
  double memory_usage() override;
  void allocate();
  void deallocate();
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;

 protected:
  int nmax;
  class AtomVecDielectric *avec;
  class NeighList *list;

  void set_dielectric_params(double, double, double, double, int, double);
  void charge_rescaled(int);
  void update_induced_charges();

  double **inverse_matrix;
  double **G1ww, **ndotGww, **G2ww, **G3ww, **Rww;

  int *
      tag2mat;    // tag2mat[atom->tag[i]] = the index of the atoms in the induced charge arrays from atom tags
  int *
      mat2tag;    // mat2tag[idx] = the atom tag of the induced charge idx in the induced charge arrays
  int *induced_charge_idx;    // return the index of the atoms in the induced charge arrays
  int num_induced_charges;    // total number of induced charges
  double *induced_charges;    // values of induced charges
  int *
      tag2mat_ions;    // tag2mat_ions[atom->tag[i]] returns the index of the atoms in the ion arrays from atom tags
  int *mat2tag_ions;    // mat2tag_ions[idx] returns the atom tag of the ion idx in the ion arrays
  int *ion_idx;         // return the index of the atoms in the ion arrays
  int num_ions;         // total number of ions
  double *rhs1;
  double *rhs2;
  double **buffer1;
  double **buffer2;

  int allocated;
  int kspaceflag;            // 1 if kspace is used for the induced charge computation
  double **efield_pair;      // electrical field at position of atom i due to pair contribution
  double **efield_kspace;    // electrical field at position of atom i due to kspace contribution
  int torqueflag, extraflag;
  double g_ewald;
  int includingG3ww;

  void calculate_Rww_cutoff();
  void calculate_qiRqw_cutoff();

  // qw, qq ion-interface terms

  double *qiRqwVector;
  double **G1qw_real;
  double *sum2G2wq;
  double *sum1G2qw;
  double *sum1G1qw_epsilon;
  double *sum2ndotGwq_epsilon;

  // conjugate gradient solver

  double *cg_r;
  double *cg_p;
  double *cg_Ap;
  double **cg_A;
  double tolerance;

  void calculate_matrix_multiply_vector(double **, double *, double *, int);
  double inner_product(double *, double *, int);
  void cg_solver(double **, double *, double *, int);

  inline double greens_real(double);
  inline double grad_greens_real_factor(double);
  inline void calculate_grad_greens_real(double *, double, double, double);
  inline double calculate_greens_ewald(double, double, double);
  inline void calculate_grad_greens_ewald(double *, double, double, double);
  inline double calculate_greens_ewald_self_vertex(double);
  inline double calculate_ndotgreens_ewald_self_vertex(double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif
