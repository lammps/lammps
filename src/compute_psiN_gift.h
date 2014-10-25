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

#ifdef COMPUTE_CLASS

ComputeStyle(PsiNGift/atom,ComputePsiNGift)

#else

#ifndef LMP_COMPUTE_PSIN_GIFT_H
#define LMP_COMPUTE_PSIN_GIFT_H

#include "compute.h"
#include "domain.h"
#include "complex"
#include "vector"


namespace LAMMPS_NS {

class ComputePsiNGift : public Compute {
 public:
  ComputePsiNGift(class LAMMPS *, int, char **);
  ~ComputePsiNGift();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double Npsi;
  
  double sphere_step;
  double cutsq;
  class NeighList *list;
  double **PsiN;
  
  double calculatePhi(double delx,double dely);
  std::vector<int> createHull(int i,std::vector<int> keys,std::vector<double> distance);
  std::vector<double> sort(std::vector<double> array);
  void sort2arraysbyvalue(std::vector<int> &keys,std::vector<double> &distance);
    /*
  std::complex<> psi;
  std::complex j1;
  */ 
};

}

#endif
#endif
