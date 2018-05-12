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

#ifndef LMP_NTOPO_H
#define LMP_NTOPO_H

#include "pointers.h"

namespace LAMMPS_NS {

class NTopo : protected Pointers {
 public:
  int nbondlist,nanglelist,ndihedrallist,nimproperlist;
  int **bondlist,**anglelist,**dihedrallist,**improperlist;

  NTopo(class LAMMPS *);
  virtual ~NTopo();

  virtual void build() = 0;

  bigint memory_usage();

 protected:
  int me,nprocs;
  int maxbond,maxangle,maxdihedral,maximproper;
  int cluster_check;             // copy from Neighbor

  void allocate_bond();
  void allocate_angle();
  void allocate_dihedral();
  void allocate_improper();

  void bond_check();
  void angle_check();
  void dihedral_check(int, int **);
};

}

#endif

/* ERROR/WARNING messages:

E: Bond extent > half of periodic box length

UNDOCUMENTED

E: Angle extent > half of periodic box length

UNDOCUMENTED

E: Dihedral/improper extent > half of periodic box length

UNDOCUMENTED

*/
