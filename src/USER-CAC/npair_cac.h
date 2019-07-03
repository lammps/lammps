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
  double shape_function(double, double, double, int, int);
  void compute_surface_depths(double &x, double &y, double &z,
	  int &xb, int &yb, int &zb, int flag);

  int CAC_decide_quad2element(int);

  double ***current_nodal_positions;
  double cutneighmax;
  int   quadrature_node_count;
  int   quadrature_point_count;
  int **element_scale;
  int current_element_scale[3];
  int current_poly_counter;
  double *quadrature_weights;
  double *quadrature_abcissae;
  void quadrature_init(int degree);
  virtual bigint memory_usage();
  

protected:
  int old_atom_count, old_quad_count;
  int expansion_count, max_expansion_count;
	int *old_atom_etype;
  int quad_allocated;
  int maxneigh_quad;
  int nmax;
  double **current_element_quad_points;
  double current_quad_point[3];
  int **surface_counts;
  int surface_counts_max[3];
  double **interior_scales;
	int surface_counts_max_old[3];
	int *neighbor_copy_index;
  int quadrature_counter;
  int bad_bin_flag;
  int ***quad_list_container;
  int *scan_flags;
  int max_quad_alloc;
  void allocate_quad_neigh_list(int,int,int,int);
  int compute_quad_points(int);
  void allocate_surface_counts();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: excessive/negative bin index

Internal debug error message that may be removed in the near future.

*/
