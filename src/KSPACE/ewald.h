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

KSpaceStyle(ewald,Ewald)

#else

#ifndef LMP_EWALD_H
#define LMP_EWALD_H

#include "kspace.h"

namespace LAMMPS_NS {

class Ewald : public KSpace {
 public:
  Ewald(class LAMMPS *, int, char **);
  virtual ~Ewald();
  void init();
  void setup();
  virtual void compute(int, int);
  double memory_usage();

  void compute_group_group(int, int, int);

 protected:
  int kxmax,kymax,kzmax;
  int kcount,kmax,kmax3d,kmax_created;
  double gsqmx,qsum,qsqsum,q2,volume;
  int nmax;

  double unitk[3];
  int *kxvecs,*kyvecs,*kzvecs;
  double *ug;
  double **eg,**vg;
  double **ek;
  double *sfacrl,*sfacim,*sfacrl_all,*sfacim_all;
  double ***cs,***sn;

  // group-group interactions

  int group_allocate_flag;
  double *sfacrl_A,*sfacim_A,*sfacrl_A_all,*sfacim_A_all;
  double *sfacrl_B,*sfacim_B,*sfacrl_B_all,*sfacim_B_all;

  double rms(int, double, bigint, double);
  virtual void eik_dot_r();
  void coeffs();
  virtual void allocate();
  void deallocate();
  void slabcorr();

  // group-group interactions

  void allocate_groups();
  void deallocate_groups();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use Ewald with triclinic box

This feature is not yet supported.

E: Cannot use Ewald with 2d simulation

The kspace style ewald cannot be used in 2d simulations.  You can use
2d Ewald in a 3d simulation; see the kspace_modify command.

E: Kspace style requires atom attribute q

The atom style defined does not have these attributes.

E: Cannot use nonperiodic boundaries with Ewald

For kspace style ewald, all 3 dimensions must have periodic boundaries
unless you use the kspace_modify command to define a 2d slab with a
non-periodic z dimension.

E: Incorrect boundaries with slab Ewald

Must have periodic x,y dimensions and non-periodic z dimension to use
2d slab option with Ewald.

E: KSpace style is incompatible with Pair style

Setting a kspace style requires that a pair style with a long-range
Coulombic component be selected.

E: Cannot use kspace solver on system with no charge

No atoms in system have a non-zero charge.

W: System is not charge neutral, net charge = %g

The total charge on all atoms on the system is not 0.0, which
is not valid for Ewald or PPPM.

*/
