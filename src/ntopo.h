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

#ifndef LMP_NTOPO_H
#define LMP_NTOPO_H

#include "pointers.h"    // IWYU pragma: keep

namespace LAMMPS_NS {

class NTopo : protected Pointers {
 public:
  int nbondlist, nanglelist, ndihedrallist, nimproperlist;
  int **bondlist, **anglelist, **dihedrallist, **improperlist;

  NTopo(class LAMMPS *);
  ~NTopo() override;

  virtual void build() = 0;

  void add_temporary_bond(int, int, int);
  double memory_usage();

 protected:
  int me, nprocs;
  int maxbond, maxangle, maxdihedral, maximproper;
  int cluster_check;    // copy from Neighbor

  void allocate_bond();
  void allocate_angle();
  void allocate_dihedral();
  void allocate_improper();

  void bond_check();
  void angle_check();
  void dihedral_check(int, int **);
};

}    // namespace LAMMPS_NS

#endif
