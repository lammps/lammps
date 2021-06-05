/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(efield/atom,ComputeEfieldAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_EFIELD_ATOM_H
#define LMP_COMPUTE_EFIELD_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeEfieldAtom : public Compute {
 public:
  ComputeEfieldAtom(class LAMMPS *, int, char **);
  ~ComputeEfieldAtom();
  void init();
  void setup();
  void compute_peratom();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int pairflag;
  int kspaceflag;
  double **efield_pair, **efield_kspace;

  int nmax;
  double **efield;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute stress/atom temperature ID

Self-explanatory.

E: Compute stress/atom temperature ID does not compute temperature

The specified compute must compute temperature.

E: Per-atom virial was not tallied on needed timestep

You are using a thermo keyword that requires potentials to have
tallied the virial, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
