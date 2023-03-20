// clang-format off
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

/* ----------------------------------------------------------------------
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

   Please cite the related publication:
   H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
   "Parallel Reactive Molecular Dynamics: Numerical Methods and
   Algorithmic Techniques", Parallel Computing, in press.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(REAXFF,FixReaxFF);
// clang-format on
#else

#ifndef LMP_FIX_REAXFF_H
#define LMP_FIX_REAXFF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixReaxFF : public Fix {
  friend class PairReaxFF;
  friend class PairReaxFFOMP;

 public:
  FixReaxFF(class LAMMPS *, int, char **);
  ~FixReaxFF() override;
  int setmask() override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:
  int maxbonds;       // max # of bonds for any atom
  int maxhbonds;      // max # of Hbonds for any atom
  int *num_bonds;     // # of bonds for each atom
  int *num_hbonds;    // # of Hbonds for each atom
  int oldnmax;        // arrays' size before growing
};

}    // namespace LAMMPS_NS

#endif
#endif
