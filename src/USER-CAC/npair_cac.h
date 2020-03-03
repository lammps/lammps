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

#ifdef NPAIR_CLASS

NPairStyle(cac,
           NPairCAC,
	 NP_BIN | NP_ATOMONLY |
	NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_CAC)

#else

#ifndef LMP_NPAIR_CAC_H
#define LMP_NPAIR_CAC_H

#include "npair.h"

namespace LAMMPS_NS {

class NPairCAC : public NPair {
 public:

  NPairCAC(class LAMMPS *);
  ~NPairCAC();
  void build(class NeighList *);
  void compute_surface_depths(double &x, double &y, double &z,
	  int &xb, int &yb, int &zb, int eindex, int etype);

  int CAC_decide_quad2element(double *, int);
  int CAC_decide_element2element(int, int);

  //functions for Asa_Data and NPair to use
  double shape_function(double, double, double,int,int);
  void interpolation(int);
  double shape_function_derivative(double, double, double,int,int,int);
  int quad_sector_select(double, double, double,int,int);
  void allocate_quad_int(int **&x, int);
  void allocate_quad_double(double **&x, int);
  
  class Asa_Data *asa_pointer; 
  //variables for Asa_Data and others to obtain
  int poly_min; 
  int surf_select[2];
  int dof_surf_list[4];
  double quad_r[3];
  int neigh_nodes_per_element; 
  int *neighbor_element_scale, neighbor_element_type;
  double ***neighbor_element_positions; 
  double **current_nodal_positions;
  double cutneighmax;
  int quadrature_node_count;
  int poly_counter;
  
  double *quadrature_weights;
  double *quadrature_abcissae;
  virtual bigint memory_usage();
  

protected:
  int max_atom_count;
  int neigh_allocated, quad_allocated;
  int maxneigh, quad_count_max;
  int nmax;
  int sector_flag, outer_neigh_flag, ghost_quad_flag;
  double ***inner_quad_lists_ucell, ***outer_quad_lists_ucell;
  int **surface_counts, ***inner_quad_lists_index, *inner_quad_lists_counts,
    ***outer_quad_lists_index, *outer_quad_lists_counts;
  int *inner_quad_neigh_maxes, *outer_quad_neigh_maxes;
  int surface_counts_max[3];
  double **interior_scales;
	int surface_counts_max_old[3];
  int pqi, qi, neigh_count;
  int **list_container, *list_maxes, *interface_flags;
  int *scan_flags,*bin_scan_flags,*bins_searched, nbins_searched, maxbins_searched;
  double **quadrature_point_data, cut_global;
  int quad_scan_list_max, quadrature_point_max, quadrature_poly_max, quadrature_poly_count;
  int *quadrature_counts, *e2quad_index;
  int *current_quad_list;
  int max_quad_per_element;
  int **surf_set;
  int **dof_set;
  int **sort_surf_set;
  int **sort_dof_set;
  int **neighbor_weights;

  void quadrature_init(int degree);
  void allocate_neigh_list();
  int compute_quad_points(int);
  void allocate_local_arrays();
  void allocate_quad_neigh_list();
  void quad_list_build(int, double, double, double);
  void neighbor_accumulate(double,double,double,int, int,int);
  void compute_quad_neighbors(int);
  void grow_quad_data();
  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: excessive/negative bin index

Internal debug error message that may be removed in the near future.

*/
