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
  class Asa_Data *asa_pointer; 
  //variable for Asa_Data to obtain
  int poly_min; 
  int surf_select[2];
  int dof_surf_list[4];
  double quad_r[3];
  int neigh_nodes_per_element; 
  int *neighbor_element_scale;
  double ***neighbor_element_positions;    
   
  PairCAC(class LAMMPS *);
  virtual ~PairCAC();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **){}
  virtual void init_style();
  virtual double init_one(int, int){ return 0.0; }

  //functions for Asa_Data and class to obtain
  double shape_function(double, double, double,int,int);
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
  double cutforcesq;
  double **scale;
  int *current_element_scale;
  double mapped_volume;
  int reneighbor_time;
  int max_nodes_per_element, neigh_poly_count;

  //stores quadrature point coordinates and calculation coefficients for all nlocal
  double ** quadrature_point_data;
  double quadrature_point_max;
  double *quadrature_counts;
	
  double cut_global_s;
  int   quadrature_node_count;
  double cutoff_skin;
  double cell_vectors[3][3];
  double interior_scale[3];
  int **surf_set;
  int **dof_set;
  int **sort_surf_set;
  int **sort_dof_set;
  double shape_args[3];
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
  double *quad_node,quad_weight;  
  double *quadrature_weights;
  double *quadrature_abcissae;
  double *quadrature_result;
  double **shape_quad_result;
  double ***current_nodal_positions;
  double ***current_nodal_gradients;
  double **neighbor_copy_ucell;
  int **neighbor_copy_index;
  int neighbor_element_type;
  int old_atom_count, old_quad_count;
  int *old_atom_etype;
  int ****inner_quad_lists_index;
  double ****inner_quad_lists_ucell;
  int **inner_quad_lists_counts;
  int ****outer_quad_lists_index;
  double ****outer_quad_lists_ucell;
  int **outer_quad_lists_counts;
  double **old_quad_minima;
  double *old_minima_neighbors;
  double **interior_scales;
  int **surface_counts;
  int atomic_flag;
  int nmax;
  int expansion_count_inner, expansion_count_outer, max_expansion_count_inner, max_expansion_count_outer;
  int neighrefresh;
  int maxneigh;
  int maxneigh_quad_inner, maxneigh_quad_outer;
  int maxneigh2;
  int surface_counts_max[3];
  int surface_counts_max_old[3];
  int current_element_type, current_poly_count;
  int natomic, atomic_counter;
  int *type_array;
  int poly_counter;
  int current_list_index;
  int interior_flag;
  int neigh_quad_counter;
  int quad_list_counter;
  int local_inner_max;
  int local_outer_max;
  int densemax;

	virtual void allocate();
	virtual void read_file(char *) {}
  virtual void array2spline(){}
  virtual void file2array(){}
  void copy_vectors(int copymode);

  //further CAC functions 
  void quadrature_init(int degree);
  void check_existence_flags();
  //void init_asa_cg();
  void allocate_quad_neigh_list(int,int,int,int);
  void allocate_surface_counts();
  void compute_mass_matrix();
  void compute_quad_neighbors(int);
  void compute_forcev(int);
  void grow_quad_data();
  void neigh_list_cord(double& coordx, double& coordy, double& coordz, int, int, double, double, double);
  void set_shape_functions();
  void compute_surface_depths(double &x, double &y, double &z, 
    int &xb, int &yb, int &zb, int flag);
  void LUPSolve(double **A, int *P, double *b, int N, double *x);
  void neighbor_accumulate(double,double,double,int, int,int);
  int LUPDecompose(double **A, int N, double Tol, int *P);
  double shape_product(int,int);
  void quad_list_build(int, double, double, double);
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
