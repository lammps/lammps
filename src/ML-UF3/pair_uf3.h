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
  void uf3_read_pot_file(char *potf_name);
  void uf3_read_pot_file(int i, int j, char *potf_name);
  void uf3_read_pot_file(int i, int j, int k, char *potf_name);
  void uf3_read_unified_pot_file(char *potf_name);
  void communicate();
  int nbody_flag, n2body_pot_files, n3body_pot_files, tot_pot_files;
  int bsplines_created;
  bool pot_3b;
  int ***setflag_3b, **knot_spacing_type_2b, ***knot_spacing_type_3b;
  double **cut, ***cut_3b, **cut_3b_list, ****min_cut_3b;
  virtual void allocate();
  void create_bsplines();
  struct UF3Impl *uf3_impl; //PIMPLE (pointer-to-implementation)
  UF3Impl *get_UF3Impl();

  int max_num_knots_2b = 0;
  int max_num_coeff_2b = 0;
  int max_num_knots_3b = 0;
  int max_num_coeff_3b = 0;
  double ***n2b_knots_array, ***n2b_coeff_array;
  int **n2b_knots_array_size, **n2b_coeff_array_size;

  int ***map_3b, tot_interaction_count_3b;
  double ***n3b_knots_array, ****n3b_coeff_array;
  int **n3b_knots_array_size, **n3b_coeff_array_size;

  /*void uf3_read_2b_pot_block(int itype, int jtype, std::string iele,
                             std::string jele,
                             TextFileReader &txtfilereader);
  
  void uf3_read_3b_pot_block(int itype, int jtype, int ktype,
                                      std::string iele, std::string jele,
                                      std::string kele,
                                      TextFileReader &txtfilereader);*/

  //Accessor function called by pair_uf3_kokkos.cpp
  //Will probably be removed once std::vector are converted to arrays
  std::vector<std::vector<std::vector<double>>>& get_n2b_knot();
  std::vector<std::vector<std::vector<double>>>& get_n2b_coeff();
  std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& get_n3b_knot_matrix();
  std::vector<std::vector<std::vector<double>>>& get_n3b_coeff_matrix_key(std::string key);
  double get_knot_spacing_2b(int i, int j);
  double get_knot_spacing_3b_ij(int i, int j, int k);
  double get_knot_spacing_3b_ik(int i, int j, int k);
  double get_knot_spacing_3b_jk(int i, int j, int k);
  int *neighshort, maxshort;    // short neighbor list array for 3body interaction
};

}    // namespace LAMMPS_NS

#endif
#endif

