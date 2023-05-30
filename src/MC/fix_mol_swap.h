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
  ~FixMolSwap() override;
  int setmask() override;
  void init() override;
  void pre_exchange() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double compute_vector(int) override;
  void write_restart(FILE *) override;
  void restart(char *) override;

 private:
  int nevery, ncycles, seed;
  int itype, jtype;
  double temperature;

  int ke_flag;          // 1 if kinetic energy is also swapped
  double i2j_vscale;    // scale factors for velocity to keep KE constant
  double j2i_vscale;

  int qflag;        // 1 if charge is also swapped
  double iq, jq;    // charge values for all itype,jtype atoms

  bool unequal_cutoffs;     // 1 if itype and jtype have any different cutoffs
  tagint minmol, maxmol;    // range of mol IDs selected for swaps

  double nswap_attempt;    // cummulative stats on MC attempts and accepts
  double nswap_accept;

  double beta;             // 1/kT
  double energy_stored;    // energy of current state as swaps are accepted

  class RanPark *random;
  class Compute *c_pe;

  int attempt_swap();
  double energy_full();
};

}    // namespace LAMMPS_NS

#endif
#endif
