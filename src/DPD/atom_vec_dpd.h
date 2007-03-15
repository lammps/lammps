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

#ifndef ATOM_VEC_DPD_H
#define ATOM_VEC_DPD_H

#include "atom_vec_atomic.h"

namespace LAMMPS_NS {

class AtomVecDPD : public AtomVecAtomic {
 public:
  AtomVecDPD(class LAMMPS *, int, char **);
  void zero_ghost(int, int);
  int pack_comm(int, int *, double *, int, int *);
  int pack_comm_one(int, double *);
  void unpack_comm(int, int, double *);
  int unpack_comm_one(int, double *);
  int pack_border(int, int *, double *, int, int *);
  int pack_border_one(int, double *);
  void unpack_border(int, int, double *);
  int unpack_border_one(int, double *);
};

}

#endif
