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

#ifdef FIX_CLASS
// clang-format off
FixStyle(neighbor/swap,FixNeighborSwap);
// clang-format on
#else

#ifndef LMP_FIX_NEIGH_SWAP_H
#define LMP_FIX_NEIGH_SWAP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNeighborSwap : public Fix {
 public:
  FixNeighborSwap(class LAMMPS *, int, char **);
  ~FixNeighborSwap();
  int setmask();
  void init();
  void pre_exchange();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double compute_vector(int);
  double memory_usage();
  void write_restart(FILE *);
  void restart(char *);

 private:
  int nevery, seed;
  int ke_flag;       // yes = conserve ke, no = do not conserve ke
  int diff_flag;     // yes = simulate diffusion of central atom, no = swap only to certain types
  int rates_flag;    // yes = use modified type rates, no = swap rates are equivilent across types
  int ncycles;
  int niswap, njswap;                  // # of i,j swap atoms on all procs
  int niswap_local, njswap_local;      // # of swap atoms on this proc
  int niswap_before, njswap_before;    // # of swap atoms on procs < this proc
  // int global_i_ID;                     // global id of selected i atom
  class Region *region;    // swap region
  char *idregion;          // swap region id

  int nswaptypes;
  int jtype_selected;
  int id_center;
  double x_center;
  double y_center;
  double z_center;
  int *type_list;
  double *rate_list;

  double nswap_attempts;
  double nswap_successes;

  bool unequal_cutoffs;

  int atom_swap_nmax;
  double beta, r_0;
  double local_probability;     // Total swap probability stored on this proc
  double global_probability;    // Total swap probability across all proc
  double prev_probability;      // Swap probability on proc < this proc
  double *qtype;
  double energy_stored;
  double **sqrt_mass_ratio;
  double **voro_neighbor_list;
  int *local_swap_iatom_list;
  int *local_swap_neighbor_list;
  int *local_swap_type_list;    // Type list index of atoms stored on this proc
  double *local_swap_probability;

  class RanPark *random_equal;

  class Compute *c_voro;
  class Compute *c_pe;

  void options(int, char **);
  int attempt_swap();
  double energy_full();
  int pick_i_swap_atom();
  int pick_j_swap_neighbor(int);
  double get_distance(double[3], double[3]);
  void build_i_neighbor_list(int);
  void update_iswap_atoms_list();
};

}    // namespace LAMMPS_NS

#endif
#endif
