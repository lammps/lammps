/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/dipole,PPPMDipole);
// clang-format on
#else

#ifndef LMP_PPPM_DIPOLE_H
#define LMP_PPPM_DIPOLE_H

#include "pppm.h"

namespace LAMMPS_NS {

class PPPMDipole : public PPPM {
 public:
  PPPMDipole(class LAMMPS *);
  virtual ~PPPMDipole();
  void init();
  void setup();
  void setup_grid();
  void compute(int, int);
  int timing_1d(int, double &);
  int timing_3d(int, double &);
  double memory_usage();

 protected:
  void set_grid_global();
  double newton_raphson_f();

  void allocate();
  void allocate_peratom();
  void deallocate();
  void deallocate_peratom();
  void compute_gf_denom();

  void slabcorr();

  // grid communication

  void pack_forward_grid(int, void *, int, int *);
  void unpack_forward_grid(int, void *, int, int *);
  void pack_reverse_grid(int, void *, int, int *);
  void unpack_reverse_grid(int, void *, int, int *);

  // dipole

  FFT_SCALAR ***densityx_brick_dipole, ***densityy_brick_dipole, ***densityz_brick_dipole;
  FFT_SCALAR ***vdxx_brick_dipole, ***vdyy_brick_dipole, ***vdzz_brick_dipole;
  FFT_SCALAR ***vdxy_brick_dipole, ***vdxz_brick_dipole, ***vdyz_brick_dipole;
  FFT_SCALAR ***ux_brick_dipole, ***uy_brick_dipole, ***uz_brick_dipole;
  FFT_SCALAR ***v0x_brick_dipole, ***v1x_brick_dipole, ***v2x_brick_dipole;
  FFT_SCALAR ***v3x_brick_dipole, ***v4x_brick_dipole, ***v5x_brick_dipole;
  FFT_SCALAR ***v0y_brick_dipole, ***v1y_brick_dipole, ***v2y_brick_dipole;
  FFT_SCALAR ***v3y_brick_dipole, ***v4y_brick_dipole, ***v5y_brick_dipole;
  FFT_SCALAR ***v0z_brick_dipole, ***v1z_brick_dipole, ***v2z_brick_dipole;
  FFT_SCALAR ***v3z_brick_dipole, ***v4z_brick_dipole, ***v5z_brick_dipole;
  FFT_SCALAR *work3, *work4;
  FFT_SCALAR *densityx_fft_dipole, *densityy_fft_dipole, *densityz_fft_dipole;

  class GridComm *gc_dipole;

  int only_dipole_flag;
  double musum, musqsum, mu2;

  double find_gewald_dipole(double, double, bigint, double, double);
  double newton_raphson_f_dipole(double, double, bigint, double, double);
  double derivf_dipole(double, double, bigint, double, double);
  double compute_df_kspace_dipole();
  double compute_qopt_dipole();
  void compute_gf_dipole();
  void make_rho_dipole();
  void brick2fft_dipole();
  void poisson_ik_dipole();
  void poisson_peratom_dipole();
  void fieldforce_ik_dipole();
  void fieldforce_peratom_dipole();
  double final_accuracy_dipole();
  void musum_musq();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot (yet) use charges with Kspace style PPPMDipole

Charge-dipole interactions are not yet implemented in PPPMDipole so this
feature is not yet supported.

E: Must redefine kspace_style after changing to triclinic box

Self-explanatory.

E: Kspace style requires atom attribute mu

The atom style defined does not have this attribute.

E: Cannot (yet) use kspace_modify diff ad with dipoles

This feature is not yet supported.

E: Cannot (yet) use 'electron' units with dipoles

This feature is not yet supported.

E: Cannot yet use triclinic cells with PPPMDipole

This feature is not yet supported.

E: Cannot yet use TIP4P with PPPMDipole

This feature is not yet supported.

E: Cannot use nonperiodic boundaries with PPPM

For kspace style pppm, all 3 dimensions must have periodic boundaries
unless you use the kspace_modify command to define a 2d slab with a
non-periodic z dimension.

E: Incorrect boundaries with slab PPPM

Must have periodic x,y dimensions and non-periodic z dimension to use
2d slab option with PPPM.

E: PPPM order cannot be < 2 or > than %d

This is a limitation of the PPPM implementation in LAMMPS.

E: KSpace style is incompatible with Pair style

Setting a kspace style requires that a pair style with matching
long-range dipole components be used.

W: Reducing PPPM order b/c stencil extends beyond nearest neighbor processor

This may lead to a larger grid than desired. See the kspace_modify overlap
command to prevent changing of the PPPM order.

E: PPPM order < minimum allowed order

The default minimum order is 2. This can be reset by the
kspace_modify minorder command.

E: PPPM grid stencil extends beyond nearest neighbor processor

This is not allowed if the kspace_modify overlap setting is no.

E: Cannot (yet) compute per-atom virial with kspace style pppm/dipole

This feature is not yet supported.

E: KSpace accuracy must be > 0

The kspace accuracy designated in the input must be greater than zero.

E: Could not compute grid size

The code is unable to compute a grid size consistent with the desired
accuracy. This error should not occur for typical problems. Please
send an email to the developers.

E: PPPM grid is too large

The global PPPM grid is larger than OFFSET in one or more dimensions.
OFFSET is currently set to 4096. You likely need to decrease the
requested accuracy.

E: Could not compute g_ewald

The Newton-Raphson solver failed to converge to a good value for
g_ewald. This error should not occur for typical problems. Please
send an email to the developers.

E: Non-numeric box dimensions - simulation unstable

The box size has apparently blown up.

E: Out of range atoms - cannot compute PPPM

One or more atoms are attempting to map their charge to a PPPM grid
point that is not owned by a processor. This is likely for one of two
reasons, both of them bad. First, it may mean that an atom near the
boundary of a processor's sub-domain has moved more than 1/2 the
"neighbor skin distance"_neighbor.html without neighbor lists being
rebuilt and atoms being migrated to new processors. This also means
you may be missing pairwise interactions that need to be computed.
The solution is to change the re-neighboring criteria via the
"neigh_modify"_neigh_modify command. The safest settings are "delay 0
every 1 check yes". Second, it may mean that an atom has moved far
outside a processor's sub-domain or even the entire simulation box.
This indicates bad physics, e.g. due to highly overlapping atoms, too
large a timestep, etc.

E: Using kspace solver PPPMDipole on system with no dipoles

Must have non-zero dipoles with PPPMDipole.

E: Must use kspace_modify gewald for system with no dipoles

Self-explanatory.

*/
