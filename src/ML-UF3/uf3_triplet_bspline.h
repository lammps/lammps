/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pointers.h"

#include "uf3_bspline_basis2.h"
#include "uf3_bspline_basis3.h"
#include "uf3_pair_bspline.h"

#include <vector>

#ifndef UF3_TRIPLET_BSPLINE_H
#define UF3_TRIPLET_BSPLINE_H

namespace LAMMPS_NS {
class uf3_triplet_bspline {
 private:
  LAMMPS *lmp;
  int knot_vect_size_ij, knot_vect_size_ik, knot_vect_size_jk;
  std::vector<std::vector<std::vector<double>>> coeff_matrix, dncoeff_matrix_ij, dncoeff_matrix_ik,
      dncoeff_matrix_jk;
  std::vector<std::vector<double>> knot_matrix;
  std::vector<uf3_bspline_basis3> bsplines_ij, bsplines_ik, bsplines_jk;
  std::vector<uf3_bspline_basis2> dnbsplines_ij, dnbsplines_ik, dnbsplines_jk;
  int get_starting_index_uniform(const std::vector<double>, int, double, double);
  int get_starting_index_nonuniform(const std::vector<double>, int, double, double);
  int (uf3_triplet_bspline::*get_starting_index)(const std::vector<double>, int, double, double);
  //double knot_spacing_ij=0,knot_spacing_ik=0,knot_spacing_jk=0;
  //double _alignvar(, 8) ret_val[4]; // Force memory alignment on 8 byte boundaries
  double ret_val[4];

  int starting_knot(const std::vector<double>, int, double);

 public:
  //Dummy Constructor
  uf3_triplet_bspline();
  uf3_triplet_bspline(LAMMPS *ulmp, const std::vector<std::vector<double>> &uknot_matrix,
                      const std::vector<std::vector<std::vector<double>>> &ucoeff_matrix,
                      const int &uknot_spacing_type);
  
  uf3_triplet_bspline(LAMMPS *ulmp, double **uknot_array, const int *uknot_array_size,
                      double ***ucoeff_array, const int *ucoeff_array_size,
                      const int &uknot_spacing_type);

  ~uf3_triplet_bspline();
  int knot_spacing_type;
  double knot_spacing_ij=0,knot_spacing_ik=0,knot_spacing_jk=0;
  double *eval(double value_rij, double value_rik, double value_rjk);

  double memory_usage();
};
}    // namespace LAMMPS_NS
#endif
