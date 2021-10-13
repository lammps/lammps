/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(mol/swap,FixMolSwap);
// clang-format on
#else

#ifndef LMP_FIX_MOL_SWAP_H
#define LMP_FIX_MOL_SWAP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMolSwap : public Fix {
 public:
  FixMolSwap(class LAMMPS *, int, char **);
  ~FixMolSwap();
  int setmask();
  void init();
  void pre_exchange();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double compute_vector(int);
  void write_restart(FILE *);
  void restart(char *);

 private:
  int nevery, ncycles, seed;
  int itype, jtype;
  double temperature;

  int ke_flag;            // 1 if kinetic energy is also swapped
  double i2j_vscale;      // scale factors for velocity to keep KE constant
  double j2i_vscale;

  int qflag;              // 1 if charge is also swapped
  double iq,jq;           // charge values for all itype,jtype atoms

  bool unequal_cutoffs;   // 1 if itype and jtype have any different cutoffs
  tagint minmol,maxmol;   // range of mol IDs selected for swaps

  double nswap_attempt;   // cummulative stats on MC attempts and accepts
  double nswap_accept;

  double beta;            // 1/kT
  double energy_stored;   // energy of current state as swaps are accepted

  class RanPark *random;
  class Compute *c_pe;

  int attempt_swap();
  double energy_full();
};

}    // namespace LAMMPS_NS

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
