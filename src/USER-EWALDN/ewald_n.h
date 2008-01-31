/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef EWALD_N_H
#define EWALD_N_H

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
  double bytes;
  double precision, g2_max;
  double *kenergy, energy_self[EWALD_NFUNCS];
  double *kvirial, virial_self[EWALD_NFUNCS];
  cvector *ekr_local;
  hvector *hvec;
  kvector *kvec;

  double qqrd2e, mumurd2e, dielectric, *B, volume;
  struct Sum { double x, x2; } sum[EWALD_MAX_NSUMS];
  complex *cek_local, *cek_global;
 
  void reallocate();
  void reallocate_atoms();
  void deallocate();
  void coefficients();
  void init_coeffs();
  void init_coeff_sums();
  void init_self();
  void compute_ek();
  void compute_force();
  void compute_surface();
  void compute_energy(int);
  void compute_virial(int);
  void compute_slabcorr(int);
};

}

#endif

