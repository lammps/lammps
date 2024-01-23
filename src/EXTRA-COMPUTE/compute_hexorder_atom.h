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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(hexorder/atom,ComputeHexOrderAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_HEXORDER_ATOM_H
#define LMP_COMPUTE_HEXORDER_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeHexOrderAtom : public Compute {
 public:
  ComputeHexOrderAtom(class LAMMPS *, int, char **);
  ~ComputeHexOrderAtom() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int nmax, maxneigh, ncol, nnn, ndegree;
  double cutsq;
  class NeighList *list;
  double *distsq;
  int *nearest;

  double **qnarray;

  void calc_qn_complex(double, double, double &, double &);
  void calc_qn_trig(double, double, double &, double &);
  void select2(int, int, double *, int *);
};

}    // namespace LAMMPS_NS

#endif
#endif
