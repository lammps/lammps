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

#ifndef LMP_KSPACE_H
#define LMP_KSPACE_H

#include "pointers.h"    // IWYU pragma: export

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_FLOAT
#else
typedef double FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

namespace LAMMPS_NS {

class KSpace : protected Pointers {
  friend class ThrOMP;
  friend class FixOMP;

 public:
  double energy;    // accumulated energies
  double energy_1, energy_6;
  double virial[6];          // accumulated virial: xx,yy,zz,xy,xz,yz
  double *eatom, **vatom;    // accumulated per-atom energy/virial
  double e2group;            // accumulated group-group energy
  double f2group[3];         // accumulated group-group force
  int triclinic_support;     // 1 if supports triclinic geometries

  int ewaldflag;         // 1 if a Ewald solver
  int pppmflag;          // 1 if a PPPM solver
  int msmflag;           // 1 if a MSM solver
  int dispersionflag;    // 1 if a LJ/dispersion solver
  int tip4pflag;         // 1 if a TIP4P solver
  int dipoleflag;        // 1 if a dipole solver
  int spinflag;          // 1 if a spin solver
  int differentiation_flag;
  int neighrequest_flag;    // used to avoid obsolete construction
                            // of neighbor lists
  int mixflag;              // 1 if geometric mixing rules are enforced
                            // for LJ coefficients
  bool conp_one_step;       // calculate A matrix in one step with pppm
  int slabflag, wireflag;
  int scalar_pressure_flag;    // 1 if using MSM fast scalar pressure
  double slab_volfactor, wire_volfactor;

  int warn_nonneutral;    // 0 = error if non-neutral system
                          // 1 = warn once if non-neutral system
                          // 2 = warn, but already warned
  int warn_nocharge;      // 0 = already warned
                          // 1 = warn if zero charge

  int order, order_6, order_allocated;
  double accuracy;             // accuracy of KSpace solver (force units)
  double accuracy_absolute;    // user-specified accuracy in force units
  double accuracy_relative;    // user-specified dimensionless accuracy
                               // accuracy = acc_rel * two_charge_force
  double accuracy_real_6;      // real space accuracy for
                               // dispersion solver (force units)
  double accuracy_kspace_6;    // reciprocal space accuracy for
                               // dispersion solver (force units)
  int auto_disp_flag;          // use automatic parameter generation for pppm/disp
  double two_charge_force;     // force in user units of two point
                               // charges separated by 1 Angstrom

  double g_ewald, g_ewald_6;
  int nx_pppm, ny_pppm, nz_pppm;          // global FFT grid for Coulombics
  int nx_pppm_6, ny_pppm_6, nz_pppm_6;    // global FFT grid for dispersion
  int nx_msm_max, ny_msm_max, nz_msm_max;

  int group_group_enable;    // 1 if style supports group/group calculation

  int centroidstressflag;    // centroid stress compared to two-body stress
                             // CENTROID_SAME = same as two-body stress
                             // CENTROID_AVAIL = different and implemented
                             // CENTROID_NOTAVAIL = different, not yet implemented

  // KOKKOS host/device flag and data masks

  ExecutionSpace execution_space;
  unsigned int datamask_read, datamask_modify;
  int copymode;

  int compute_flag;       // 0 if skip compute()
  int fftbench;           // 0 if skip FFT timing
  int collective_flag;    // 1 if use MPI collectives for FFT/remap
  int stagger_flag;       // 1 if using staggered PPPM grids

  double splittol;    // tolerance for when to truncate splitting

  KSpace(class LAMMPS *);
  ~KSpace() override;
  void two_charge();
  void triclinic_check();
  void modify_params(int, char **);
  void *extract(const char *);
  void compute_dummy(int, int);

  // triclinic

  void x2lamdaT(double *, double *);
  void lamda2xT(double *, double *);
  void lamda2xvector(double *, double *);
  void kspacebbox(double, double *);

  // public so can be called by commands that change charge

  void qsum_qsq(int warning_flag = 1);

  // general child-class methods

  virtual void settings(int, char **){};
  virtual void init() = 0;
  virtual void setup() = 0;
  virtual void setup_grid(){};
  virtual void compute(int, int) = 0;
  virtual void compute_group_group(int, int, int){};

  virtual void pack_forward_grid(int, void *, int, int *){};
  virtual void unpack_forward_grid(int, void *, int, int *){};
  virtual void pack_reverse_grid(int, void *, int, int *){};
  virtual void unpack_reverse_grid(int, void *, int, int *){};

  virtual int timing(int, double &, double &) { return 0; }
  virtual int timing_1d(int, double &) { return 0; }
  virtual int timing_3d(int, double &) { return 0; }

  virtual int modify_param(int, char **) { return 0; }
  virtual double memory_usage() { return 0.0; }

  /* ----------------------------------------------------------------------
   compute gamma for MSM and pair styles
   see Eq 4 from Parallel Computing 35 (2009) 164-177
------------------------------------------------------------------------- */

  double gamma(const double &rho) const
  {
    if (rho <= 1.0) {
      const int split_order = order / 2;
      const double rho2 = rho * rho;
      double g = gcons[split_order][0];
      double rho_n = rho2;
      for (int n = 1; n <= split_order; n++) {
        g += gcons[split_order][n] * rho_n;
        rho_n *= rho2;
      }
      return g;
    } else
      return (1.0 / rho);
  }

  /* ----------------------------------------------------------------------
   compute the derivative of gamma for MSM and pair styles
   see Eq 4 from Parallel Computing 35 (2009) 164-177
------------------------------------------------------------------------- */

  double dgamma(const double &rho) const
  {
    if (rho <= 1.0) {
      const int split_order = order / 2;
      const double rho2 = rho * rho;
      double dg = dgcons[split_order][0] * rho;
      double rho_n = rho * rho2;
      for (int n = 1; n < split_order; n++) {
        dg += dgcons[split_order][n] * rho_n;
        rho_n *= rho2;
      }
      return dg;
    } else
      return (-1.0 / rho / rho);
  }

  double **get_gcons() { return gcons; }
  double **get_dgcons() { return dgcons; }

 protected:
  int gridflag, gridflag_6;
  int gewaldflag, gewaldflag_6;
  int minorder, overlap_allowed;
  int adjust_cutoff_flag;
  int suffix_flag;    // suffix compatibility flag
  bigint natoms_original;
  double scale, qqrd2e;
  double qsum, qsqsum, q2;
  double **gcons, **dgcons;    // accumulated per-atom energy/virial

  int evflag, evflag_atom;
  int eflag_either, eflag_global, eflag_atom;
  int vflag_either, vflag_global, vflag_atom;
  int maxeatom, maxvatom;

  int kewaldflag;                      // 1 if kspace range set for Ewald sum
  int kx_ewald, ky_ewald, kz_ewald;    // kspace settings for Ewald sum

  void pair_check();
  void ev_init(int eflag, int vflag, int alloc = 1)
  {
    if (eflag || vflag)
      ev_setup(eflag, vflag, alloc);
    else
      evflag = evflag_atom = eflag_either = eflag_global = eflag_atom = vflag_either =
          vflag_global = vflag_atom = 0;
  }
  void ev_setup(int, int, int alloc = 1);
  double estimate_table_accuracy(double, double);
};

}    // namespace LAMMPS_NS

#endif
