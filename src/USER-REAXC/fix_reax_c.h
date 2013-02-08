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

/* ----------------------------------------------------------------------
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

   Please cite the related publication:
   H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
   "Parallel Reactive Molecular Dynamics: Numerical Methods and
   Algorithmic Techniques", Parallel Computing, in press.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(REAXC,FixReaxC)

#else

#ifndef LMP_FIX_REAXC_H
#define LMP_FIX_REAXC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixReaxC : public Fix {
  friend class PairReaxC;

 public:
  FixReaxC(class LAMMPS *,int, char **);
  ~FixReaxC();
  int setmask();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);

 private:
  int maxbonds;              // max # of bonds for any atom
  int maxhbonds;             // max # of Hbonds for any atom
  int *num_bonds;            // # of bonds for each atom
  int *num_hbonds;           // # of Hbonds for each atom
};

}

#endif
#endif
