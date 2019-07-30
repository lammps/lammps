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

#ifdef PAIR_CLASS

PairStyle(hybrid/kk,PairHybridKokkos)

#else

#ifndef LMP_PAIR_HYBRID_KOKKOS_H
#define LMP_PAIR_HYBRID_KOKKOS_H

#include "pair_hybrid.h"
#include "pair_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class PairHybridKokkos : public PairHybrid {
  friend class FixGPU;
  friend class FixIntel;
  friend class FixOMP;
  friend class Force;
  friend class Respa;
  friend class Info;
 public:
  typedef LMPDeviceType device_type;

  PairHybridKokkos(class LAMMPS *);
  virtual ~PairHybridKokkos();
  void compute(int, int);

 private:
  DAT::t_x_array_randomread x;
  DAT::t_f_array f;
  friend void pair_virial_fdotr_compute<PairHybridKokkos>(PairHybridKokkos*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Pair style hybrid cannot have hybrid as an argument

Self-explanatory.

E: Pair style hybrid cannot have none as an argument

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair coeff for hybrid has invalid style

Style in pair coeff must have been listed in pair_style command.

E: Pair hybrid sub-style is not used

No pair_coeff command used a sub-style specified in the pair_style
command.

E: Pair_modify special setting for pair hybrid incompatible with global special_bonds setting

Cannot override a setting of 0.0 or 1.0 or change a setting between
0.0 and 1.0.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Invoked pair single on pair style none

A command (e.g. a dump) attempted to invoke the single() function on a
pair style none, which is illegal.  You are probably attempting to
compute per-atom quantities with an undefined pair style.

E: Pair hybrid sub-style does not support single call

You are attempting to invoke a single() call on a pair style
that doesn't support it.

E: Pair hybrid single calls do not support per sub-style special bond values

Self-explanatory.

E: Unknown pair_modify hybrid sub-style

The choice of sub-style is unknown.

E: Coulomb cutoffs of pair hybrid sub-styles do not match

If using a Kspace solver, all Coulomb cutoffs of long pair styles must
be the same.

*/
