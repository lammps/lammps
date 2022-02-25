/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS

// clang-format off
KSpaceStyle(pppm/electrode, PPPMElectrode)
// clang-format on

#else

#ifndef LMP_PPPM_ELECTRODE_H
#define LMP_PPPM_ELECTRODE_H

#include "boundary_correction.h"
#include "electrode_kspace.h"
#include "pppm.h"
#include <algorithm>

#if defined(FFT_FFTW3)
#define LMP_FFT_LIB "FFTW3"
#elif defined(FFT_MKL)
#define LMP_FFT_LIB "MKL FFT"
#elif defined(FFT_CUFFT)
#define LMP_FFT_LIB "cuFFT"
#else
#define LMP_FFT_LIB "KISS FFT"
#endif

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define LMP_FFT_PREC "single"
#define MPI_FFT_SCALAR MPI_FLOAT
#else

typedef double FFT_SCALAR;
#define LMP_FFT_PREC "double"
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

namespace LAMMPS_NS {

class PPPMElectrode : public PPPM, public ElectrodeKSpace {
 public:
  PPPMElectrode(class LAMMPS *);
  virtual ~PPPMElectrode();
  virtual void init();
  virtual void setup();
  virtual void setup_grid();
  virtual void compute(int, int);

  void compute_vector(bigint *, double *);
  void compute_vector_corr(bigint *, double *);
  void compute_matrix(bigint *, double **);
  void compute_matrix_corr(bigint *, double **);

  virtual void compute_group_group(int, int, int);

 protected:
  FFT_SCALAR ***electrolyte_density_brick;
  FFT_SCALAR *electrolyte_density_fft;
  class BoundaryCorrection *boundcorr;

  virtual void set_grid_global();
  void set_grid_local();

  virtual void allocate();
  virtual void deallocate();
  virtual void allocate_peratom();
  double compute_df_kspace();
  // double estimate_ik_error(double, double, bigint);
  virtual double compute_qopt();
  virtual void compute_gf_ik();
  virtual void compute_gf_ad();

  /* ----------------------------------------------------------------------
     denominator for Hockney-Eastwood Green's function
       of x,y,z = sin(kx*deltax/2), etc

              inf                 n-1
     S(n,k) = Sum  W(k+pi*j)**2 = Sum b(l)*(z*z)**l
             j=-inf               l=0

            = -(z*z)**n /(2n-1)! * (d/dx)**(2n-1) cot(x)  at z = sin(x)
     gf_b = denominator expansion coeffs
  ------------------------------------------------------------------------- */

  inline double gf_denom(const double &x, const double &y, const double &z) const
  {
    double sx, sy, sz;
    sz = sy = sx = 0.0;
    for (int l = order - 1; l >= 0; l--) {
      sx = gf_b[l] + sx * x;
      sy = gf_b[l] + sy * y;
      sz = gf_b[l] + sz * z;
    }
    double s = sx * sy * sz;
    return s * s;
  };

 private:
  int compute_step;
  void start_compute();
  void make_rho_in_brick(bigint *, FFT_SCALAR ***, bool);
  void project_psi(bigint *, double *vec);
  void one_step_multiplication(bigint *, std::vector<double>, double **, double **, int const);
  void two_step_multiplication(bigint *, std::vector<double>, double **, double **, int const);
  bool compute_vector_called;
  bigint *imat_cached;
};

}    // namespace LAMMPS_NS

#endif
#endif

    /* ERROR/WARNING messages:

    E: Illegal ... command

    Self-explanatory.  Check the input script syntax and compare to the
    documentation for the command.  You can use -echo screen as a
    command-line option when running LAMMPS to see the offending line.

    E: Must redefine kspace_style after changing to triclinic box

    UNDOCUMENTED

    E: Cannot (yet) use PPPM with triclinic box and kspace_modify diff ad

    This feature is not yet supported.

    E: Cannot (yet) use PPPM with triclinic box and slab correction

    This feature is not yet supported.

    E: Cannot use PPPM with 2d simulation

    The kspace style pppm cannot be used in 2d simulations.  You can use
    2d PPPM in a 3d simulation; see the kspace_modify command.

    E: PPPM can only currently be used with comm_style brick

    This is a current restriction in LAMMPS.

    E: KSpace style requires atom attribute q

    The atom style defined does not have these attributes.

    E: Cannot use non-periodic boundaries with PPPM

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
    long-range Coulombic or dispersion components be used.

    E: Pair style is incompatible with TIP4P KSpace style

    The pair style does not have the requires TIP4P settings.

    E: Bond and angle potentials must be defined for TIP4P

    Cannot use TIP4P pair potential unless bond and angle potentials
    are defined.

    E: Bad TIP4P angle type for PPPM/TIP4P

    Specified angle type is not valid.

    E: Bad TIP4P bond type for PPPM/TIP4P

    Specified bond type is not valid.

    W: Reducing PPPM order b/c stencil extends beyond nearest neighbor processor

    This may lead to a larger grid than desired.  See the kspace_modify overlap
    command to prevent changing of the PPPM order.

    E: PPPM order < minimum allowed order

    The default minimum order is 2.  This can be reset by the
    kspace_modify minorder command.

    E: PPPM grid stencil extends beyond nearest neighbor processor

    This is not allowed if the kspace_modify overlap setting is no.

    E: KSpace accuracy must be > 0

    The kspace accuracy designated in the input must be greater than zero.

    E: Must use kspace_modify gewald for uncharged system

    UNDOCUMENTED

    E: Could not compute grid size

    The code is unable to compute a grid size consistent with the desired
    accuracy.  This error should not occur for typical problems.  Please
    send an email to the developers.

    E: PPPM grid is too large

    The global PPPM grid is larger than OFFSET in one or more dimensions.
    OFFSET is currently set to 4096.  You likely need to decrease the
    requested accuracy.

    E: Could not compute g_ewald

    The Newton-Raphson solver failed to converge to a good value for
    g_ewald.  This error should not occur for typical problems.  Please
    send an email to the developers.

    E: Non-numeric box dimensions - simulation unstable

    The box size has apparently blown up.

    E: Out of range atoms - cannot compute PPPM

    One or more atoms are attempting to map their charge to a PPPM grid
    point that is not owned by a processor.  This is likely for one of two
    reasons, both of them bad.  First, it may mean that an atom near the
    boundary of a processor's sub-domain has moved more than 1/2 the
    "neighbor skin distance"_neighbor.html without neighbor lists being
    rebuilt and atoms being migrated to new processors.  This also means
    you may be missing pairwise interactions that need to be computed.
    The solution is to change the re-neighboring criteria via the
    "neigh_modify"_neigh_modify command.  The safest settings are "delay 0
    every 1 check yes".  Second, it may mean that an atom has moved far
    outside a processor's sub-domain or even the entire simulation box.
    This indicates bad physics, e.g. due to highly overlapping atoms, too
    large a timestep, etc.

    E: Cannot (yet) use K-space slab correction with compute group/group for
    triclinic systems

    This option is not yet supported.

    E: Cannot (yet) use kspace_modify diff ad with compute group/group

    This option is not yet supported.

    U: Cannot (yet) use PPPM with triclinic box and TIP4P

    This feature is not yet supported.

    */
