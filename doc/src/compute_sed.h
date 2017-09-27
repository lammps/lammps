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
   Contributing authors: Tianli Feng, Divya Chalise, Xiulin Ruan
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(sed,ComputeSed)

#else

#ifndef LMP_COMPUTE_Sed_H
#define LMP_COMPUTE_Sed_H

#ifndef __EIGEN__
#define __EIGEN__
typedef struct eigen{
  double re[3];
  double im[3];
} EIGEN;
#endif

#include "compute.h"

namespace LAMMPS_NS {
	
class ComputeSed : public Compute {
 public:
  ComputeSed(class LAMMPS *, int, char **);
  ~ComputeSed();
  void init();
  void compute_vector();
 private:
  int nk, nkpp, nkst, nked;   // Number definition for k-points 
  int nlines, nlinemax, nlinemin;
  int natoms, nbasis;
  double delt;
  int gap;
  int *idmid20, *at2bs, *indat2cell;
  double **kpt, **ratom;
  double *mass;
  EIGEN *eigenv;
};

}

#endif
#endif