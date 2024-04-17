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

#include <vector>

#ifndef UF3_PAIR_BSPLINE_H
#define UF3_PAIR_BSPLINE_H

namespace LAMMPS_NS {

class uf3_pair_bspline {
 private:
  int knot_vect_size, coeff_vect_size;
  std::vector<double> knot_vect, dnknot_vect;
  std::vector<double> coeff_vect, dncoeff_vect;
  std::vector<uf3_bspline_basis3> bspline_bases;
  std::vector<uf3_bspline_basis2> dnbspline_bases;
  int get_starting_index_uniform(double), get_starting_index_nonuniform(double);
  int (uf3_pair_bspline::*get_starting_index)(double);
  //double knot_spacing=0;
  LAMMPS *lmp;

 public:
  // dummy constructor
  uf3_pair_bspline();
  uf3_pair_bspline(LAMMPS *ulmp, const std::vector<double> &uknot_vect,
                   const std::vector<double> &ucoeff_vect,
                   const int &uknot_spacing_type);
  
  uf3_pair_bspline(LAMMPS *ulmp, const double* uknot_array,
                   const int uknot_array_size,
                   const double* ucoeff_array,
                   const int ucoeff_array_size,
                   const int uknot_spacing_type);

  ~uf3_pair_bspline();
  int knot_spacing_type;
  double knot_spacing=0;
  double ret_val[2];
  double *eval(double value_rij);
  double memory_usage();
};
}    // namespace LAMMPS_NS
#endif
