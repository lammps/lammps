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

#ifdef FIX_CLASS
// clang-format off
FixStyle(msst,FixMSST);
// clang-format on
#else

#ifndef FIX_MSST_H
#define FIX_MSST_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMSST : public Fix {
 public:
  FixMSST(class LAMMPS *, int, char **);
  ~FixMSST() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  double compute_scalar() override;
  double compute_vector(int) override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  int modify_param(int, char **) override;
  double memory_usage() override;

 private:
  double dtv, dtf, dthalf;        // full and half step sizes
  double boltz, nktv2p, mvv2e;    // Boltzmann factor and unit conversions
  double total_mass;              // mass of the computational cell

  double omega[3];    // time derivative of the volume
  double p_current[3], dilation[3];
  double qmass;     // effective cell mass
  double mu;        // effective cell viscosity
  double tscale;    // converts thermal energy to compressive
                    // strain ke at simulation start
  int dftb;         // flag for use with DFTB+

  double velocity_sum;                  // sum of the velocities squared
  double damping;                       // damping function for TS force term at
                                        //   small volume difference (v0 - vol)
  double T0S0;                          // initial TS term for DFTB+ simulations
  double S_elec, S_elec_1, S_elec_2;    // time history of electron entropy
                                        //   for DFTB+ simulaitons
  double TS_dot;                        // time derivative of TS term for
                                        //   DFTB+ simulations

  double **old_velocity;    // saved velocities

  int kspace_flag;    // 1 if KSpace invoked, 0 if not
  int nrigid;         // number of rigid fixes
  int *rfix;          // indices of rigid fixes

  char *id_temp, *id_press;    // strings with identifiers of
  char *id_pe;                 // created computes

  class Compute *temperature;    // computes created to evaluate
  class Compute *pressure;       // thermodynamic quantities
  class Compute *pe;
  int tflag, pflag, vsflag, peflag;    // flags to keep track of computes that
                                       // were created

  // shock initial conditions

  double e0;                     // initial energy
  double v0;                     // initial volume
  double p0;                     // initial pressure
  double velocity;               // velocity of the shock
  double lagrangian_position;    // Lagrangian location of computational cell
  int direction;                 // direction of shock
  int p0_set;                    // is pressure set
  int v0_set;                    // is volume set
  int e0_set;                    // is energy set
  double TS_int;                 // needed for conserved quantity
                                 //   with thermal electronic excitations
  double beta;                   // energy conservation scaling factor

  int maxold;                         // allocated size of old_velocity
  class FixExternal *fix_external;    // ptr to fix external

  // functions

  void couple();
  void remap(int);
  double compute_etotal();
  double compute_vol();
  double compute_hugoniot();
  double compute_rayleigh();
  double compute_lagrangian_speed();
  double compute_lagrangian_position();
  double compute_vsum();
};

}    // namespace LAMMPS_NS

#endif
#endif
