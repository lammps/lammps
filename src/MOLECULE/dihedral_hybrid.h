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

#ifndef DIHEDRAL_HYBRID_H
#define DIHEDRAL_HYBRID_H

#include "stdio.h"
#include "dihedral.h"

class DihedralHybrid : public Dihedral {
 public:
  DihedralHybrid();
  ~DihedralHybrid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, int, char **);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  int memory_usage();

 private:
  int nstyles;                  // # of different dihedral styles
  Dihedral **styles;            // class list for each Dihedral style
  char **keywords;              // keyword for each dihedral style
  int *map;                     // which style each dihedral type points to

  int *ndihedrallist;           // # of dihedrals in sub-style dihedrallists
  int *maxdihedral;             // max # of dihedrals sub-style lists can store
  int ***dihedrallist;          // dihedrallist for each sub-style
  
  void allocate();
};

#endif
