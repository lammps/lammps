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
FixStyle(orient/eco,FixOrientECO);
// clang-format on
#else

#ifndef LMP_FIX_ORIENT_ECO_H
#define LMP_FIX_ORIENT_ECO_H

#include "fix.h"

namespace LAMMPS_NS {

class FixOrientECO : public Fix {
 public:
  FixOrientECO(class LAMMPS *, int, char **);
  ~FixOrientECO();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  double compute_scalar();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double memory_usage();

 private:
  struct Nbr;    // forward declaration. private struct for managing precomputed terms

  int me;              // this processors rank
  int nmax;            // maximal # of owned + ghost atoms on this processor
  int ilevel_respa;    // used for RESPA integrator only

  int sign;                              // from sign of u
  double u_0;                            // synthetic potential energy
  double half_u;                         // half synthetic potential energy
  double eta;                            // threshold for thermal effects
  double inv_eta;                        // inverse threshold for thermal effects
  double r_cut;                          // cutoff radius
  double squared_cutoff;                 // squared cutoff radius
  double inv_squared_cutoff;             // inverse squared cutoff radius
  char *dir_filename;                    // filename of reference grain input
  double dir_vec[6][3];                  // direct lattice vectors
  double reciprocal_vectors[2][3][3];    // reciprocal lattice vectors

  double added_energy;    // energy added by fix

  double **order;    // order parameter and normalized order
                     // parameter per atom

  double norm_fac;        // normalization constant
  double inv_norm_fac;    // inverse normalization constant

  Nbr *nbr;                 // pointer on array of precomputed terms
  class NeighList *list;    // LAMMPS' neighbor list

  void get_reciprocal();    // calculate reciprocal lattice vectors
  int get_norm();           // compute normalization factor
};

}    // namespace LAMMPS_NS

#endif    // LMP_FIX_ORIENT_ECO_H
#endif    // FIX_CLASS
