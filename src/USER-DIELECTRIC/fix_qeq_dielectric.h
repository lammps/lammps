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

/* ----------------------------------------------------------------------
   Contributing author:
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(qeq/dielectric,FixQEqDielectric)

#else

#ifndef LMP_FIX_QEQ_DIELECTRIC_H
#define LMP_FIX_QEQ_DIELECTRIC_H

#include <vector>
#include "atom.h"
#include "fix.h"
//#include "ewald_dielectric.h"

namespace LAMMPS_NS {

class FixQEqDielectric : public Fix {
 public:
  FixQEqDielectric(class LAMMPS *, int, char **);
  ~FixQEqDielectric();
  int setmask();
  void init();
  void init_list(int,class NeighList *);
  void setup_pre_force(int);
  void pre_force(int);
  int modify_param(int, char**);

  void setup_pre_force_respa(int, int);
  void pre_force_respa(int, int, int);
  void min_setup_pre_force(int);
  void min_pre_force(int);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  double memory_usage();

  void calculate_Rww_full();
  void calculate_Rww_cutoff();

  void calculate_qiRqw_full();
  void calculate_qiRqw_cutoff();

  void compute_induced_charges();
  void set_dielectric_params(double, double, double, double);

 protected:
  class AtomVecDielectric* avec;

  int nevery;
  int nlevels_respa;
  class NeighList *list;

  bigint ngroup;
  int full;
  int includingG3ww;

  // Ewald-related functions and data structures
//  class EwaldDielectric *ewaldDielectric;

  // interface term, ww
  // matrix M*M
  double **inverse_matrix;
  double **G1ww, **ndotGww, **G2ww, **G3ww, **Rww;
  double ***gradG1ww; //matrix of vectors M*M*3
  double *q_induced_charges;

  // qw, qq ion-interface terms
  double *qiRqwVector;
  double  **G1qq_real, **G1qw_real;
  double  ***gradG1wq_real; //matrix of vectors M*M*3
  double *sum2G2wq;
  // temp data for each steps
  double *sum1G2qw;
  double *sum1G3qw;
  double *sum1G1qw_epsilon;
  double *sum2ndotGwq_epsilon;

  double greens_real(double);
  double grad_greens_real_factor(double);
  void calculate_grad_greens_real(double *, double, double, double);
  double calculate_greens_ewald(double, double, double);
  void calculate_grad_greens_ewald(double *, double, double, double);
  void calculate_inverse_matrix(double **, double **, int);
  void calculate_matrix_multiply_vector(double **, double *, double *, int);
  double calculate_greens_ewald_self_vertex(double);
  double calculate_ndotgreens_ewald_self_vertex(double, double);
  double g_ewald;

 private:
  // to store the global id (tag) of different particles.
  int n_induced_charges, n_ions;
  int *tags_interface, *tags_ions;
  void setup_tags();
  inline int get_index_interface(int num) {
    int id = tags_interface[num]; // id = fixed global id
    return atom->map(id);
  };
  inline int get_index_ions(int num) {
    int id = tags_ions[num]; // id = fixed global id
    return atom->map(id);
  };
  void change_all_ions_q_to_real_q();
  void change_all_ions_q_to_scaled_q();

  int n_induced_charges_local, n_ions_local;
  std::vector<int> tags_interface_local, tags_ions_local;
  std::vector<int> matrix_index_from_local_index;
  std::vector<int> matrix_index_from_global_id;
  void setup_tags_local();
  inline int get_index_interface_local(int num) {
    int id = tags_interface_local[num]; // id = fixed global id
    return atom->map(id);
  };
  inline int get_index_ions_local(int num) {
    int id = tags_ions_local[num]; // id = fixed global id
    return atom->map(id);
  };

  inline int get_matrix_index_from_local_index(int k) {
    // for atom with local index k, what is the index in matrix (Rww + RwwT)^-1 for the whole system.
    // notice this local index include ghost index too!
    int id = atom->tag[k]; // id = fixed global id
    return matrix_index_from_global_id[id];
  }

  // read all em, ed, norm of interface, and epsilon of ions
  void print_all_properties();

};

}

#endif
#endif
