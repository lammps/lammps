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
  double gsqmx,q2,b2,M2;
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

  double rms(int, double, bigint, double, double, double);
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
  void compute_virial_dipole();
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

For kspace style ewald/disp, all 3 dimensions must have periodic
boundaries unless you use the kspace_modify command to define a 2d
slab with a non-periodic z dimension.

E: Incorrect boundaries with slab EwaldDisp

Must have periodic x,y dimensions and non-periodic z dimension to use
2d slab option with Ewald.

E: KSpace style is incompatible with Pair style

Setting a kspace style requires that a pair style with a long-range
Coulombic or dispersion component be used.

E: Unsupported mixing rule in kspace_style ewald/disp

Only geometric mixing is supported.

E: Unsupported order in kspace_style ewald/disp

Only 1/r^6 dispersion or dipole terms are supported.

E: Cannot (yet) use 'electron' units with dipoles

This feature is not yet supported.

E: Cannot use Ewald/disp solver on system with no charge, dipole, or LJ particles

No atoms in system have a non-zero charge or dipole, or are LJ particles.  Change
charges/dipoles or change options of the kspace solver/pair style.

W: System is not charge neutral, net charge = %g

The total charge on all atoms on the system is not 0.0, which
is not valid for the long-range Coulombic solvers.

W: Ewald/disp Newton solver failed, using old method to estimate g_ewald

Self-explanatory. Choosing a different cutoff value may help.

E: KSpace accuracy too low

Requested accuracy must be less than 1.0.

E: Epsilon or sigma reference not set by pair style in ewald/n

The pair style is not providing the needed epsilon or sigma values.

E: Cannot (yet) use kspace slab correction with long-range dipoles and non-neutral systems or per-atom energy

This feature is not yet supported.

*/
