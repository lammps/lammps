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

#ifndef FIX_PERI_NEIGH_H
#define FIX_PERI_NEIGH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPeriNeigh : public Fix {
  friend class PairPeriPMB;
  friend class ComputeDamageAtom;

 public:
  FixPeriNeigh(class LAMMPS *,int, char **);
  ~FixPeriNeigh();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

 private:
  int first;                 // flag for first time initialization
  int maxpartner;            // max # of peridynamic neighs for any atom
  int *npartner;             // # of neighbors for each atom
  int **partner;             // neighs for each atom, stored as global IDs
  double **r0;               // initial distance to partners
  double *vinter;            // sum of vfrac for bonded neighbors

  class NeighList *list;
};

}

#endif
