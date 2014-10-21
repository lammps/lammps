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

#ifdef COMPUTE_CLASS

ComputeStyle(snav/atom,ComputeSNAVAtom)

#else

#ifndef LMP_COMPUTE_SNAV_ATOM_H
#define LMP_COMPUTE_SNAV_ATOM_H

#include "compute.h"
 
namespace LAMMPS_NS {

class ComputeSNAVAtom : public Compute {
 public:
  ComputeSNAVAtom(class LAMMPS *, int, char **);
  ~ComputeSNAVAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int nmax, njmax, diagonalstyle;
  int ncoeff,nvirial;
  double **cutsq;
  class NeighList *list;
  double **snav;
  double rcutfac;
  double *radelem;
  double *wjelem;

  class SNA** snaptr;

};

}

#endif
#endif
