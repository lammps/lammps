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

#ifdef DIHEDRAL_CLASS
// clang-format off
DihedralStyle(hybrid,DihedralHybrid);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_HYBRID_H
#define LMP_DIHEDRAL_HYBRID_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralHybrid : public Dihedral {
 public:
  int nstyles;          // # of different dihedral styles
  Dihedral **styles;    // class list for each Dihedral style
  char **keywords;      // keyword for each dihedral style

  DihedralHybrid(class LAMMPS *);
  ~DihedralHybrid() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double memory_usage() override;

 private:
  int *map;    // which style each dihedral type points to

  int *ndihedrallist;     // # of dihedrals in sub-style dihedrallists
  int *maxdihedral;       // max # of dihedrals sub-style lists can store
  int ***dihedrallist;    // dihedrallist for each sub-style

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
