/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef ATOM_ANGLE_H
#define ATOM_ANGLE_H

#include "atom.h"

class AtomAngle : public Atom {
 public:
  AtomAngle(int, char **);
  ~AtomAngle() {}
  void copy(int, int);
  void pack_comm(int, int *, double *, int *);
  void unpack_comm(int, int, double *);
  void pack_reverse(int, int, double *);
  void unpack_reverse(int, int *, double *);
  void pack_border(int, int *, double *, int *);
  void unpack_border(int, int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
};

#endif
