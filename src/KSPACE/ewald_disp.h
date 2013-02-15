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

KSpaceStyle(ewald/disp,EwaldDisp)

#else

#ifndef LMP_EWALD_DISP_H
#define LMP_EWALD_DISP_H

#include "kspace.h"
#include "math_complex.h"

namespace LAMMPS_NS {

#define EWALD_NORDER        6
#define EWALD_NFUNCS        4
#define EWALD_MAX_NSUMS     10
#define EWALD_NSUMS        {1, 1, 7, 1}

typedef struct cvector { complex x, y, z; } cvector;
typedef struct hvector { double x, y, z; } hvector;
typedef struct kvector { long x, y, z; } kvector;

class EwaldDisp : public KSpace {
 public:
  EwaldDisp(class LAMMPS *, int, char **);
  ~EwaldDisp();
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

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use EwaldDisp with 2d simulation

This is a current restriction of this command.

E: Cannot use nonperiodic boundaries with EwaldDisp

UNDOCUMENTED

E: Incorrect boundaries with slab EwaldDisp

UNDOCUMENTED

E: KSpace style is incompatible with Pair style

UNDOCUMENTED

E: Unsupported mixing rule in kspace_style ewald/disp

UNDOCUMENTED

E: Unsupported order in kspace_style ewald/disp

UNDOCUMENTED

E: Cannot use Ewald/disp solver on system with no charge or LJ particles

UNDOCUMENTED

W: System is not charge neutral, net charge = %g

The total charge on all atoms on the system is not 0.0, which
is not valid for Ewald or PPPM.

E: KSpace accuracy too large to estimate G vector

UNDOCUMENTED

W: Ewald/disp Newton solver failed, using old method to estimate g_ewald

UNDOCUMENTED

E: KSpace accuracy too low

UNDOCUMENTED

E: epsilon or sigma reference not set by pair style in ewald/n

UNDOCUMENTED

*/
