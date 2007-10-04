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

#ifndef DIHEDRAL_HYBRID_H
#define DIHEDRAL_HYBRID_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralHybrid : public Dihedral {
 public:
  DihedralHybrid(class LAMMPS *);
  ~DihedralHybrid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, int, char **);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  double memory_usage();

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

}

#endif
