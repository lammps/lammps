/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
KSpaceStyle(pppm/electrode, PPPMElectrode);
// clang-format on
#else

#ifndef LMP_PPPM_ELECTRODE_H
#define LMP_PPPM_ELECTRODE_H

#include "electrode_kspace.h"
#include "pppm.h"

namespace LAMMPS_NS {

class PPPMElectrode : public PPPM, public ElectrodeKSpace {
 public:
  PPPMElectrode(class LAMMPS *);
  ~PPPMElectrode() override;
  void init() override;
  void setup() override;
  void reset_grid() override;
  void compute(int, int) override;

  void compute_vector(double *, int, int, bool) override;
  void compute_vector_corr(double *, int, int, bool) override;
  void compute_matrix(bigint *, double **, bool) override;
  void compute_matrix_corr(bigint *, double **) override;

  void compute_group_group(int, int, int) override;

 protected:
  FFT_SCALAR ***electrolyte_density_brick;
  FFT_SCALAR *electrolyte_density_fft;
  class BoundaryCorrection *boundcorr;

  void set_grid_global() override;
  void set_grid_local() override;

  void allocate() override;
  void deallocate() override;
  void allocate_peratom() override;
  double compute_df_kspace() override;
  double compute_qopt() override;
  void compute_gf_ik() override;
  void compute_gf_ad() override;

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
  int last_source_grpbit;
  bool last_invert_source;
  void start_compute();
  void make_rho_in_brick(int, FFT_SCALAR ***, bool);
  void project_psi(double *, int);
  void one_step_multiplication(bigint *, double *, double **, double **, int const, bool);
  void two_step_multiplication(bigint *, double *, double **, double **, int const, bool);
  void build_amesh(int, int, int, double *, double *);
  bool compute_vector_called;
};

}    // namespace LAMMPS_NS

#endif
#endif
