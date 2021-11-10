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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(fragment/atom,ComputeFragmentAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_FRAGMENT_ATOM_H
#define LMP_COMPUTE_FRAGMENT_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeFragmentAtom : public Compute {
 public:
  ComputeFragmentAtom(class LAMMPS *, int, char **);
  ~ComputeFragmentAtom();
  void init();
  void compute_peratom();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double memory_usage();

 private:
  int nmax, commflag, singleflag;
  int *stack, *clist, *markflag;
  double *fragmentID;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute fragment/atom used when bonds are not allowed

UNDOCUMENTED

E: Cannot use compute fragment/atom unless atoms have IDs

Atom IDs are used to identify fragments.

E: Compute fragment/atom requires a bond style to be defined

This is so that a bond list is generated which is used to find fragments.

W: More than one compute fragment/atom

It is not efficient to use compute fragment/atom  more than once.

*/
