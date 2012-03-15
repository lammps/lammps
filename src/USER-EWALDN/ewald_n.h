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

#ifdef KSPACE_CLASS

KSpaceStyle(ewald/n,EwaldN)

#else

#ifndef LMP_EWALD_N_H
#define LMP_EWALD_N_H

#include "kspace.h"
#include "math_complex.h"

namespace LAMMPS_NS {

#define EWALD_NORDER	6
#define EWALD_NFUNCS	4
#define EWALD_MAX_NSUMS	10
#define EWALD_NSUMS	{1, 1, 7, 1}

typedef struct cvector { complex x, y, z; } cvector;
typedef struct hvector { double x, y, z; } hvector;
typedef struct kvector { long x, y, z; } kvector;

class EwaldN : public KSpace {
 public:
  EwaldN(class LAMMPS *, int, char **);
  ~EwaldN();
  void init();
  void setup();
  void compute(int, int);
  double memory_usage() {return bytes;}

 private:
  double unit[6];
  int function[EWALD_NFUNCS], first_output;

  int nkvec, nkvec_max, nevec, nevec_max,
      nbox, nfunctions, nsums, sums;
  int peratom_allocate_flag;
  int nmax;
  double bytes;
  double gsqmx,q2,b2;
  double *kenergy, energy_self[EWALD_NFUNCS];
  double *kvirial, virial_self[EWALD_NFUNCS];
  double **energy_self_peratom;
  double **virial_self_peratom;
  cvector *ekr_local;
  hvector *hvec;
  kvector *kvec;

  double mumurd2e, dielectric, *B, volume;
  struct Sum { double x, x2; } sum[EWALD_MAX_NSUMS];
  complex *cek_local, *cek_global;
 
  double rms(int, double, bigint, double, double);
  void reallocate();
  void allocate_peratom();
  void reallocate_atoms();
  void deallocate();
  void deallocate_peratom();
  void coefficients();
  void init_coeffs();
  void init_coeff_sums();
  void init_self();
  void init_self_peratom();
  void compute_ek();
  void compute_force();
  void compute_surface();
  void compute_energy();
  void compute_energy_peratom();
  void compute_virial();
  void compute_virial_peratom();
  void compute_slabcorr();
  double NewtonSolve(double, double, bigint, double, double);
  double f(double, double, bigint, double, double);
  double derivf(double, double, bigint, double, double);
};

}

#endif
#endif
