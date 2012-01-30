/* ----------------------------------------------------------------------
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

PairStyle(kim,PairKIM)

#else

#ifndef LMP_PAIR_KIM_H
#define LMP_PAIR_KIM_H

class KIM_API_model;
#include "pair.h"

namespace LAMMPS_NS {

class PairKIM : public Pair {
 public:
  PairKIM(class LAMMPS *);
  ~PairKIM();

  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

 private:
  double cut_global;

  char **elements;              // names of unique elements
  int *map;                     // mapping from atom types to elements
  int nelements;                // # of unique elements

  void allocate();

  // KIM data

  KIM_API_model *pkim;
  int coordinates_ind,numberOfAtoms_ind,numberAtomTypes_ind;
  int atomTypes_ind,compute_ind;
  int get_half_neigh_ind,get_full_neigh_ind,neighObject_ind;
  int cutoff_ind,energy_ind;
  int energyPerAtom_ind,force_ind,forces_ind,virialGlobal_ind;
  int virialPerAtom_ind,process_d1Edr_ind;
  int localnall;
  bool support_atypes;
  bool support_Rij;

  char *testname;
  char *modelname;
  char testfile[160];
  char modelfile[160];
  char *test_descriptor_string;

  int *atypeMapKIM;
  int *atomTypes;      // atom types converted to KIM indices

  void my_init();
  void my_init2();
  void my_free();
  void create_atypeMapKIM();
  void atypeMap2KIM();
  void set_statics();
  void set_volatiles();
  void init2zero(KIM_API_model *, int *);

  // static methods used as callbacks from KIM

  static void neigh_iterator(void *,int **,int *, int *);
  static void neigh_iterator2(void*, int **, int *, int *,
			      double **, double **);
  static int get_neigh(void **,int *, int *, int *, int *, int **, double **);

  static void process_d1Edr_init(KIM_API_model **,
				 double *, double *, double **,
				 int *, int *, int *);
  static void process_d1Edr_fast();
  static void process_d1Edr(KIM_API_model **, double *, double *,
			    double **, int *, int *, int *);
  static void err_treat(const char *, int);
  static void err_treat(int, const char *, int);
};

}

#endif
#endif
