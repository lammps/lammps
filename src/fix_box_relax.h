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

#ifdef FIX_CLASS

FixStyle(box/relax,FixBoxRelax)

#else

#ifndef LMP_FIX_BOX_RELAX_H
#define LMP_FIX_BOX_RELAX_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBoxRelax : public Fix {
 public:
  FixBoxRelax(class LAMMPS *, int, char **);
  ~FixBoxRelax();
  int setmask();
  void init();

  double min_energy(double *);
  void min_store();
  void min_clearstore();
  void min_pushstore();
  void min_popstore();
  int min_reset_ref();
  void min_step(double, double *);
  double max_alpha(double *);
  int min_dof();

  int modify_param(int, char **);

 private:
  int p_flag[6];
  int pstyle,pcouple,allremap;
  int dimension;
  double p_target[6],p_current[6];
  double vol0,xprdinit,yprdinit,zprdinit;
  double vmax,pv2e,pflagsum;

  int current_lifo;              // LIFO stack pointer
  double boxlo0[2][3];           // box bounds at start of line search
  double boxhi0[2][3];
  double boxtilt0[2][3];         // xy,xz,yz tilts at start of line search
  double s0[3];                  // scale matrix at start of line search
  double ds[6];                  // increment in scale matrix

  char *id_temp,*id_press;
  class Compute *temperature,*pressure;
  int tflag,pflag;

  int nrigid;
  int *rfix;

  double sigma[6];                 // scaled target stress
  double utsigma[3];               // weighting for upper-tri elements 
                                   // of modified sigma
  int sigmamod_flag;               // 1 if modified sigma to be used
  double fdev[6];                  // Deviatoric force on cell
  int deviatoric_flag;             // 0 if target stress tensor is hydrostatic
  double h0[6];                    // h_inv of reference (zero strain) box
  double h0_inv[6];                // h_inv of reference (zero strain) box
  int nreset_h0;                   // interval for resetting h0
  double p_hydro;                  // hydrostatic component of target stress

  void remap();
  void couple();

  void compute_sigma();
  void compute_deviatoric();
  double compute_strain_energy();
  void compute_press_target();
  double compute_scalar();
};

}

#endif
#endif
