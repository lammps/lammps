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

#ifdef FIX_CLASS

FixStyle(atom/swap,FixAtomSwap)

#else

#ifndef LMP_FIX_MCSWAP_H
#define LMP_FIX_MCSWAP_H

#include <cstdio>
#include "fix.h"

namespace LAMMPS_NS {

class FixAtomSwap : public Fix {
 public:
  FixAtomSwap(class LAMMPS *, int, char **);
  ~FixAtomSwap();
  int setmask();
  void init();
  void pre_exchange();
  int attempt_semi_grand();
  int attempt_swap();
  double energy_full();
  int pick_semi_grand_atom();
  int pick_i_swap_atom();
  int pick_j_swap_atom();
  void update_semi_grand_atoms_list();
  void update_swap_atoms_list();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double compute_vector(int);
  double memory_usage();
  void write_restart(FILE *);
  void restart(char *);

 private:
  int nevery,seed;
  int conserve_ke_flag;                   // yes = conserve ke, no = do not conserve ke
  int semi_grand_flag;                    // yes = semi-grand canonical, no = constant composition
  int ncycles;
  int niswap,njswap;                      // # of i,j swap atoms on all procs
  int niswap_local,njswap_local;          // # of swap atoms on this proc
  int niswap_before,njswap_before;        // # of swap atoms on procs < this proc
  int nswap;                              // # of swap atoms on all procs
  int nswap_local;                        // # of swap atoms on this proc
  int nswap_before;                       // # of swap atoms on procs < this proc
  int regionflag;                         // 0 = anywhere in box, 1 = specific region
  int iregion;                            // swap region
  char *idregion;                         // swap region id

  int nswaptypes,nmutypes;
  int *type_list;
  double *mu;

  double nswap_attempts;
  double nswap_successes;

  bool unequal_cutoffs;

  int atom_swap_nmax;
  double beta;
  double *qtype;
  double energy_stored;
  double **sqrt_mass_ratio;
  int *local_swap_iatom_list;
  int *local_swap_jatom_list;
  int *local_swap_atom_list;

  class RanPark *random_equal;
  class RanPark *random_unequal;

  class Compute *c_pe;

  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix atom/swap does not exist

Self-explanatory.

E: Must specify at least 2 types in fix atom/swap command

Self-explanatory.

E: Need nswaptypes mu values in fix atom/swap command

Self-explanatory.

E: Only 2 types allowed when not using semi-grand in fix atom/swap command

Self-explanatory.

E: Mu not allowed when not using semi-grand in fix atom/swap command

Self-explanatory.

E: Invalid atom type in fix atom/swap command

The atom type specified in the atom/swap command does not exist.

E: All atoms of a swapped type must have the same charge.

Self-explanatory.

E: At least one atom of each swapped type must be present to define charges.

Self-explanatory.

E: All atoms of a swapped type must have same charge.

Self-explanatory.

E: Cannot do atom/swap on atoms in atom_modify first group

This is a restriction due to the way atoms are organized in a list to
enable the atom_modify first command.

*/
