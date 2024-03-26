/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(box/relax,FixBoxRelax);
// clang-format on
#else

#ifndef LMP_FIX_BOX_RELAX_H
#define LMP_FIX_BOX_RELAX_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBoxRelax : public Fix {
 public:
  FixBoxRelax(class LAMMPS *, int, char **);
  ~FixBoxRelax() override;
  int setmask() override;
  void init() override;

  double min_energy(double *) override;
  void min_store() override;
  void min_clearstore() override;
  void min_pushstore() override;
  void min_popstore() override;
  int min_reset_ref() override;
  void min_step(double, double *) override;
  double max_alpha(double *) override;
  int min_dof() override;

  int modify_param(int, char **) override;

 private:
  int p_flag[6];
  int pstyle, pcouple, allremap;
  int dimension;
  double p_target[6], p_current[6];
  double vol0, xprdinit, yprdinit, zprdinit;
  double vmax, pv2e, pflagsum;
  int kspace_flag;

  static constexpr int MAX_LIFO_DEPTH = 2;
  int current_lifo;                      // LIFO stack pointer
  double boxlo0[MAX_LIFO_DEPTH][3];      // low box bounds at start of line search
  double boxhi0[MAX_LIFO_DEPTH][3];      // high box bounds at start of line search
  double boxtilt0[MAX_LIFO_DEPTH][3];    // xy,xz,yz tilts at start of line search
  double ds[6];                          // increment in scale matrix

  int scaleyz;    // 1 if yz scaled with lz
  int scalexz;    // 1 if xz scaled with lz
  int scalexy;    // 1 if xy scaled with ly

  double fixedpoint[3];    // Location of dilation fixed-point

  char *id_temp, *id_press;
  class Compute *temperature, *pressure;
  int tflag, pflag;

  std::vector<Fix *> rfix;

  double sigma[6];        // scaled target stress
  double utsigma[3];      // weighting for upper-tri elements
                          // of modified sigma
  int sigmamod_flag;      // 1 if modified sigma to be used
  double fdev[6];         // Deviatoric force on cell
  int deviatoric_flag;    // 0 if target stress tensor is hydrostatic
  double h0[6];           // h_inv of reference (zero strain) box
  double h0_inv[6];       // h_inv of reference (zero strain) box
  int nreset_h0;          // interval for resetting h0
  double p_hydro;         // hydrostatic component of target stress

  void remap();
  void couple();

  void compute_sigma();
  void compute_deviatoric();
  double compute_strain_energy();
  void compute_press_target();
  double compute_scalar() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
