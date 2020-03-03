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

#ifdef PAIR_CLASS

PairStyle(cac,PairCAC)

#else

#ifndef LMP_PAIR_CAC_H
#define LMP_PAIR_CAC_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCAC : public Pair {
 public:
	double cutmax;                // max cutoff for all elements
  int pre_force_flag;           // set to 1 if computing something before force   
   
  PairCAC(class LAMMPS *);
  virtual ~PairCAC();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **){}
  virtual void init_style();
  virtual double init_one(int, int){ return 0.0; }

  //functions for Asa_Data and class to obtain
  double shape_function(double, double, double,int,int);
  void interpolation(int);
  double shape_function_derivative(double, double, double,int,int,int);
 
  //set of shape functions
  double quad_shape_one(double s, double t, double w){ return (1-s)*(1-t)*(1-w)/8;}
  double quad_shape_two(double s, double t, double w){ return (1+s)*(1-t)*(1-w)/8;}
  double quad_shape_three(double s, double t, double w){ return (1+s)*(1+t)*(1-w)/8;}
  double quad_shape_four(double s, double t, double w){ return (1-s)*(1+t)*(1-w)/8;}
  double quad_shape_five(double s, double t, double w){ return (1-s)*(1-t)*(1+w)/8;}
  double quad_shape_six(double s, double t, double w){ return (1+s)*(1-t)*(1+w)/8;}
  double quad_shape_seven(double s, double t, double w){ return (1+s)*(1+t)*(1+w)/8;}
  double quad_shape_eight(double s, double t, double w){ return (1-s)*(1+t)*(1+w)/8;}
 
  virtual double memory_usage();

 protected:
  int outer_neighflag;
  int sector_flag;
  int ghost_quad;
  double cutforcesq;
  double **scale;
  int *current_element_scale;
  double mapped_volume;
  int reneighbor_time;
  int max_nodes_per_element, neigh_poly_count;
  double virial_density[6];

  //stores quadrature point coordinates and calculation coefficients
  double **quadrature_point_data;
  int quadrature_point_max;
  int *quadrature_counts, *e2quad_index;
  int max_quad_per_element;
  double *quadrature_weights;
  double *quadrature_abcissae;
  int quadrature_node_count;
	
  double cut_global_s;
  double cutoff_skin;
  int quad_allocated;
  int one_layer_flag;
  //used to call set of shape functions
  typedef double(PairCAC::*Shape_Functions)(double s, double t, double w);
  Shape_Functions *shape_functions;

  double element_energy;
  int quad_eflag;
  double quadrature_energy;
  double **mass_matrix;
  double **mass_copy;
  double **force_column;
  double *current_nodal_forces;
  double *current_force_column;
  double *current_x;
  int   *pivot;
  double **shape_quad_result;
  double **current_nodal_positions;
  double **neighbor_copy_ucell;
  int **neighbor_copy_index;
  int neighbor_element_type;
  int old_atom_count,old_all_atom_count, old_quad_count;
  int *old_atom_etype, *old_all_atom_etype;
  int ***inner_quad_lists_index;
  double ***inner_quad_lists_ucell;
  int *inner_quad_lists_counts;
  int ***outer_quad_lists_index;
  double ***outer_quad_lists_ucell;
  int *outer_quad_lists_counts;
  int atomic_flag;
  int nmax, nmax_surf;
  int neighrefresh;
  int maxneigh;
  int maxneigh_quad_inner, maxneigh_quad_outer;
  int maxneigh2;
  int current_element_type, current_poly_count;
  int natomic, atomic_counter;
  int *type_array;
  int poly_counter;
  int qi, pqi;
  int local_inner_max;
  int local_outer_max;
  int densemax;
  double **inner_neighbor_coords;
  double **outer_neighbor_coords;
  int *inner_neighbor_types;
  int *outer_neighbor_types;
  double *inner_neighbor_charges;

	virtual void allocate();
  void copy_vectors(int copymode);

  //further CAC functions 
  void check_existence_flags();
  void compute_mass_matrix();
  void compute_forcev(int);
  void set_shape_functions();
  void LUPSolve(double **A, int *P, double *b, int N, double *x);
  int LUPDecompose(double **A, int N, double Tol, int *P);
  double shape_product(int,int);
  void quadrature_init(int degree);
  virtual void pre_force_densities() {}
  virtual void force_densities(int, double, double, double, double, double
    &fx, double &fy, double &fz) {}
  int mldivide3(const double mat[3][3], const double *vec, double *ans);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: LU matrix is degenerate

For some reason the mass matrix for one of the finite elements defined
in your CAC model is ill conditioned.

E: Unexpected argument in CAC pair style invocation

Self-explanatory.  Check the input script or data file. See the documention
for appropriate syntax.

E: CAC Pair style requires a CAC atom style

Self-explanatory.

E: CAC Pair styles cannot be invoked with hybrid; also 
don't use the word hybrid in your style name to avoid this error

Self-explanatory.

E: CAC Pair style requires a CAC comm style

Self-explanatory.

E: minimum points exceed element domain

The minimization algorithm used to locate adjacent finite element 
geometry and complete the virtual atom neighborhood is encountering
an invalid solution. Request help for this supposedly impossible scenario.
See the documentation for contact info.

*/
