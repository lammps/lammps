/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(orientorder/atom,ComputeOrientOrderAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_ORIENTORDER_ATOM_H
#define LMP_COMPUTE_ORIENTORDER_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeOrientOrderAtom : public Compute {
 public:
  ComputeOrientOrderAtom(class LAMMPS *, int, char **);
  ~ComputeOrientOrderAtom() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  double memory_usage() override;
  double cutsq;
  int iqlcomp, qlcomp, qlcompflag, wlflag, wlhatflag;
  int *qlist;
  int nqlist;
  double *qnormfac, *qnormfac2;

 protected:
  int nmax, maxneigh, ncol, nnn;
  class NeighList *list;
  double *distsq;
  int *nearest;
  double **rlist;
  int qmax;
  double **qnarray;
  double **qnm_r;
  double **qnm_i;

  void select3(int, int, double *, int *, double **);
  void calc_boop(double **rlist, int numNeighbors, double qn[], int nlist[], int nnlist);

  double polar_prefactor(int, int, double);
  double associated_legendre(int, int, double);

  virtual void init_wigner3j();
  double triangle_coeff(const int a, const int b, const int c);
  double w3j(const int L, const int j1, const int j2, const int j3);
  double *w3jlist;    // Wigner coeffs
  int widx_max;
  int chunksize;
};

}    // namespace LAMMPS_NS

#endif
#endif
