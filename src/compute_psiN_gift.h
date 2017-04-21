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

  typedef double coord_t;         // coordinate type
  typedef double coord2_t;  // must be big enough to hold 2*max(|coordinate|)^2  

  int nmax;
  double Npsi;
    
  double sphere_step;
  double cutsq;
  class NeighList *list;
  double **PsiN;
 
  struct Point{
  int key;
  ComputePsiNGift::coord_t x, y;
	bool operator <(const Point &p) const {
		return x < p.x || (x == p.x && y < p.y);
	}
  };

  ComputePsiNGift::coord2_t cross(const Point &O, const Point &A, const Point &B);
  std::vector<ComputePsiNGift::Point> convex_hull(std::vector<ComputePsiNGift::Point> P);

  double calculatePhi(double delx,double dely);
  std::vector<int> createHull(int i,std::vector<int> keys,std::vector<double> distance, std::vector<double> angle);
  std::vector<int> createHull(int i,std::vector<int> keys,std::vector<double> distance, double** &x);
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
