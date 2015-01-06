/* ----------------------------------------------------------------------
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

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixAtomSwap : public Fix {
 public:
  FixAtomSwap(class LAMMPS *, int, char **);
  ~FixAtomSwap();
  int setmask();
  void init();
  void pre_exchange();
  int attempt_swap();
  double energy_full();
  int pick_i_swap_atom();
  int pick_j_swap_atom();
  void update_swap_atoms_list();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double compute_vector(int);
  double memory_usage();
  void write_restart(FILE *);
  void restart(char *);

 private:
  int atom_swap_itype,atom_swap_jtype,nevery,seed;
  int ncycles;
  int niswap,njswap;                      // # of swap atoms on all procs
  int niswap_local,njswap_local;          // # of swap atoms on this proc
  int niswap_before,njswap_before;        // # of swap atoms on procs < this proc
  int regionflag;                         // 0 = anywhere in box, 1 = specific region
  int iregion;                            // swap region
  char *idregion;                         // swap region id

  double nswap_attempts;
  double nswap_successes;

  bool unequal_cutoffs;
  
  int atom_swap_nmax;
  double beta;
  double qitype,qjtype;
  double energy_stored;
  int *local_swap_iatom_list;
  int *local_swap_jatom_list;

  class RanPark *random_equal;
  
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

E: Invalid atom type in fix atom/swap command

The atom type specified in the atom/swap command does not exist.

E: Illegal fix atom/swap gas mass <= 0

The computed mass of the designated atom type was less 
than or equal to zero.

E: Cannot do atom/swap on atoms in atom_modify first group

This is a restriction due to the way atoms are organized in a list to
enable the atom_modify first command.


*/
