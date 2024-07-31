/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Ajinkya Hire (Univ. of Florida),
                        Hendrik Kraß (Univ. of Constance),
                        Matthias Rupp (Luxembourg Institute of Science and Technology),
                        Richard Hennig (Univ of Florida)
---------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(uf3,PairUF3);
// clang-format on
#else

#ifndef LMP_PAIR_UF3_H
#define LMP_PAIR_UF3_H

#include "pair.h"

namespace LAMMPS_NS {

class PairUF3 : public Pair {
 public:
  PairUF3(class LAMMPS *);
  ~PairUF3() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  void init_list(int, class NeighList *) override;    // needed for ptr to full neigh list
  double init_one(int, int) override;                 // needed for cutoff radius for neighbour list
  double single(int, int, int, int, double, double, double, double &) override;
  double memory_usage() override;

 protected:
  int ***setflag_3b, **knot_spacing_type_2b, ***knot_spacing_type_3b;
  double **cut, ***cut_3b, **cut_3b_list, ****min_cut_3b;
  double **knot_spacing_2b, ****knot_spacing_3b;

  double ***n2b_knots_array, ***n2b_coeff_array;
  int **n2b_knots_array_size, **n2b_coeff_array_size;
  double ****cached_constants_2b, ****cached_constants_2b_deri;

  int ***map_3b;
  double ***n3b_knots_array, ****n3b_coeff_array;
  int **n3b_knots_array_size, **n3b_coeff_array_size;
  double ****coeff_for_der_jk, ****coeff_for_der_ik, ****coeff_for_der_ij;
  double ****cached_constants_3b, ****cached_constants_3b_deri;

  int *neighshort, maxshort;    // short neighbor list array for 3body interaction

  void uf3_read_unified_pot_file(char *potf_name);
  void communicate();
  int bsplines_created;
  bool pot_3b;

  virtual void allocate();
  void create_bsplines();
  void create_cached_constants_2b();
  void create_cached_constants_3b();

  int get_starting_index_uniform_2b(int i, int j, double r);
  int get_starting_index_uniform_3b(int i, int j, int k, double r, int knot_dim);

  int get_starting_index_nonuniform_2b(int i, int j, double r);
  int get_starting_index_nonuniform_3b(int i, int j, int k, double r, int knot_dim);

  int (PairUF3::*get_starting_index_2b)(int i, int j, double r);
  int (PairUF3::*get_starting_index_3b)(int i, int j, int k, double r, int knot_dim);

  int nbody_flag;
  int max_num_knots_2b;
  int max_num_coeff_2b;
  int max_num_knots_3b;
  int max_num_coeff_3b;
  int tot_interaction_count_3b;
};

}    // namespace LAMMPS_NS

#endif
#endif
