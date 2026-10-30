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

#ifdef FIX_CLASS

FixStyle(spin/heat/flux,FixSpinheatFlux)

#else

#ifndef FIX_SPINHEAT_FLUX_H //Define a  macro
#define FIX_SPINHEAT_FLUX_H

#include "stdio.h"
#include "stdlib.h"
#include "fix.h"
#include <fftw3.h>

using namespace std;

namespace LAMMPS_NS {

class FixSpinheatFlux : public Fix {
 public:
  FixSpinheatFlux(class LAMMPS *, int, char **);
  ~FixSpinheatFlux();
  int setmask();
  void setup(int);
  void end_of_step();
  void init();
  virtual double compute_scalar();
 // void write_restart(FILE *);
 // void restart(char *);
  //int size_restart(int);
  //int  pack_restart(int, double *);
  //void unpack_restart(int, int);  

  private:
  int nevery, nsample;
  int ntimestep;
  bigint nvalid;           //Calculate nmodes if ts >= nvalid
  int lastindex;       // index in values ring of latest time 
  //int *itag_old;       //Global atoms IDs
  int    *mask;
  double *mass;
  double **x, **sp, **fm;
  char   *group2;
  double *Jin, *Jout;
  double dx, dy, dz, Jij, r_delta, jv, J_mag, alpha, gamma, sigma, r_cut, temp, hbar, kb;
  int i, j, hevi_step, atompgroupi, atompgroupj, jgroup,jgroupbit, me, nump, n, nlocal;
  int *Jd,*mpi_Jd,*Jd0;
  FILE *fp_W;
  bigint nextvalid();
};

}

#endif
#endif
