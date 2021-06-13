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

/* ----------------------------------------------------------------------
    Contributing author: Lukas Fath (KIT)
    some subroutines are from fix_shake.cpp
  ------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(filter/corotate,FixFilterCorotate);
// clang-format on
#else

#ifndef LMP_FIX_FILTER_COROTATE_H
#define LMP_FIX_FILTER_COROTATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFilterCorotate : public Fix {
 public:
  FixFilterCorotate(class LAMMPS *, int, char **);
  ~FixFilterCorotate();
  void setup(int);
  void setup_pre_neighbor();
  void pre_neighbor();
  void setup_pre_force_respa(int, int);
  //    void setup_post_force_respa(int,int);
  void pre_force_respa(int, int, int);
  void post_force_respa(int, int, int);

  void init();
  int setmask();

  double compute_array(int, int);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  void grow_arrays(int);
  double memory_usage();

  void copy_arrays(int, int, int);
  void set_arrays(int);
  void update_arrays(int, int);

  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

 protected:
  int me, nprocs;

  int flevel;    //filtered respa level

  double **help2;      //temp derivative
  double **x_store;    //temp for atom->x
  double *g;           //temp for derivative

  double *n1, *n2, *n3, *del1, *del2, *del3;

  double **dn1dx, **dn2dx, **dn3dx;

  int *bond_flag, *angle_flag;    // bond/angle types to constrain
  int *type_flag;                 // constrain bonds to these types
  double *mass_list;              // constrain bonds to these masses
  int nmass;                      // # of masses in mass_list

  int molecular;                             // copy of atom->molecular
  double *bond_distance, *angle_distance;    // constraint distances

  int nlevels_respa;    // copies of needed rRESPA variables

  double **x, **v, **f;    // local ptrs to atom class quantities
  double *mass, *rmass;
  int *type;
  int nlocal;
  // atom-based arrays
  int *shake_flag;    // 0 if atom not in SHAKE cluster
  // 1 = size 3 angle cluster
  // 2,3,4 = size of bond-only cluster
  tagint **shake_atom;    // global IDs of atoms in cluster
  // central atom is 1st
  // lowest global ID is 1st for size 2
  int **shake_type;    // bondtype of each bond in cluster
  // for angle cluster, 3rd value
  //   is angletype
  int *nshake;    // count

  int *list;                               // list of clusters to SHAKE
  int nlist, maxlist;                      // size and max-size of list
  double ***clist_derv;                    //stores derivative
  double **clist_q0;                       //stores reference config
  int *clist_nselect1, *clist_nselect2;    //stores length of each selec. list
  int **clist_select1, **clist_select2;    //stores selection lists

  void find_clusters();
  int masscheck(double);

  void filter_inner();
  void filter_outer();

  void general_cluster(int, int);

  void stats();
  int bondtype_findset(int, tagint, tagint, int);
  int angletype_findset(int, tagint, tagint, int);

  // callback functions for ring communication

  static void ring_bonds(int, char *, void *);
  static void ring_nshake(int, char *, void *);
  static void ring_shake(int, char *, void *);

  int sgn(double val) { return (0 < val) - (val < 0); };
};

}    // namespace LAMMPS_NS

#endif
#endif
