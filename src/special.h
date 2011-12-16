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

#ifndef LMP_SPECIAL_H
#define LMP_SPECIAL_H

#include "pointers.h"

namespace LAMMPS_NS {

class Special : protected Pointers {
 public:
  Special(class LAMMPS *);
  ~Special();
  void build();

 private:
  int me,nprocs;
  int **onetwo,**onethree,**onefour;
  int dihedral_flag;

  void combine();
  void angle_trim();
  void dihedral_trim();
};

}

#endif

/* ERROR/WARNING messages:

E: 1-3 bond count is inconsistent

An inconsistency was detected when computing the number of 1-3
neighbors for each atom.  This likely means something is wrong with
the bond topologies you have defined.

E: 1-4 bond count is inconsistent

An inconsistency was detected when computing the number of 1-4
neighbors for each atom.  This likely means something is wrong with
the bond topologies you have defined.

*/
