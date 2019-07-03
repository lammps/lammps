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

#ifdef NBIN_CLASS

NBinStyle(cac,
          NBinCAC,
          NB_CAC)

#else

#ifndef LMP_NBIN_CAC_H
#define LMP_NBIN_CAC_H

#include "nbin.h"

namespace LAMMPS_NS {

class NBinCAC : public NBin {
 public:
  NBinCAC(class LAMMPS *);
  ~NBinCAC();
  void setup_bins(int);
  void bin_atoms();
  virtual void bin_atoms_setup(int);
  
  protected:
  int *bin_overlap_limits;
  int first_alloc;
  int max_bin_expansion_count;
  int *bin_expansion_counts;
  int nmax;
  double **current_element_quad_points;
  double current_quad_point[3];
  int **surface_counts;
  int surface_counts_max[3];
  double **interior_scales;
	int surface_counts_max_old[3];
	int *neighbor_copy_index;
  int quadrature_counter;
  int   quadrature_node_count;
  double ***current_nodal_positions;
  int **element_scale;
  int current_element_scale[3];
  int current_poly_counter;
  double *quadrature_weights;
  double *quadrature_abcissae;
  int quad_rule_initialized;
  double bsubbox[3],bsubboxlo[3],bsubboxhi[3];
  int max_nall;                  //upper bound on number of local + ghost atoms that have been encountered
  int *max_nbin_overlap;         //upper bound on how many bins an element has overlapped
  int foreign_boxes;
  virtual int quad2bins(double *quad_position);
  virtual int element2bins(int element_index);
  virtual void CAC_bin_atoms_setup(int);
  int compute_quad_points(int);
  void CAC_setup_bins(int);
  void allocate_surface_counts();
  void quadrature_init(int degree);
  void expand_overlap_arrays(int size);
  void compute_surface_depths(double &x, double &y, double &z,
	  int &xb, int &yb, int &zb, int flag);
  double shape_function(double, double, double, int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot use the CAC method without the CAC comm style

Self-explanatory

E: Domain too large for neighbor bins

The domain is deemed excessively demanding in terms of bin allocation due to memory resources etc.

E: excessive/negatic bin index

Internal debugging check in case of logical failure in the algorithm. Contact author for support

E: Cannot use neighbor bins - box size << cutoff

UNDOCUMENTED

E: Too many neighbor bins

UNDOCUMENTED

*/
