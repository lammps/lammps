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

#ifndef PAIR_REAX_H
#define PAIR_REAX_H

#include "pair.h"

namespace LAMMPS_NS {

class PairREAX : public Pair {
 public:
  PairREAX(class LAMMPS *);
  ~PairREAX();

  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  
 private:
 
  double cutmax;                // max cutoff for all elements

  void allocate();
  void read_files(char *, char *);
  void neigh_f2c(int, int *, int *, int **);
  void neigh_c2f(int, int *, int *, int **);

  // For ReaxFF Verlet list

  double rcutvsq, rcutbsq; 
  int iprune,ihb;
  
  // For the midpoint method in REAX
  bool Lmidpoint; 

  // Cutoff for hydrogen bond, Taper value for Coulomb interactions
  double hbcut,swb;

  void write_reax_positions(); // particle positions into reax fortran
  void write_reax_vlist();  // Verlet neighbor lists into reax fortran
  void read_reax_forces();  // forces in reax fortran into LAMMPS
  void read_reax_atom_virial(); // virial stress per atom in reax fortran into LAMMPS

  // For charge equilibration
  double swa;
  double swc0, swc1, swc2, swc3, swc4, swc5, swc6, swc7;
  double precision;
  struct ff_params {double rcutsq; int np; double* params;};
  ff_params* param_list;
  int nentries;
  double chpot;
  // For communication of w[i] in CG solve
  double *w;

  void taper_setup();
  double taper_E(const double &, const double &);
  double taper_F(const double &, const double &);

  void compute_charge(double &);
  void sparse_product(const int &, const int &, const int &, double[],
		      int[], int[], double[], double[]);
  void cg_solve(const int &, const int &, double[], int[], 
		       int[], double[], double[]);

  void charge_reax(const int &, const int &, double[],
		   double[], int[], int[], double[]);
};

}

#endif
