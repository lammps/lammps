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
/* Implementation of the Multi-Scale Shock Method.
   See Reed, Fried, Joannopoulos, Phys. Rev. Lett., 90, 235503(2003).
   Implementation by Laurence Fried, LLNL, 4/2007.
*/

#ifdef FIX_CLASS

FixStyle(msst,FixMSST)

#else

#ifndef FIX_MSST_H
#define FIX_MSST_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMSST : public Fix {
 public:
  FixMSST(class LAMMPS *, int, char **);
  ~FixMSST();
  int setmask();
  void init();
  void setup(int);
  void initial_integrate(int);
  void final_integrate();
  double compute_scalar();
  double compute_vector(int);
  void write_restart(FILE *);
  void restart(char *);
  int modify_param(int, char **);

 private:
  double dtv,dtf,dthalf;           // Full and half step sizes.
  double boltz,nktv2p, mvv2e;      // Boltzmann factor and unit conversions.
  double total_mass;               // Mass of the computational cell.

  double omega[3];                 // Time derivative of the volume.
  double p_current[3],dilation[3];
  double qmass;                    // Effective cell mass.
  double mu;                       // Effective cell viscosity.
  double tscale;                   // Converts thermal energy to compressive 
                                   // strain ke at simulation start
  
  double velocity_sum;             // Sum of the velocities squared.

  double **old_velocity;           // Saved velocities.

  int kspace_flag;                 // 1 if KSpace invoked, 0 if not
  int nrigid;                      // number of rigid fixes
  int *rfix;                       // indices of rigid fixes

  char *id_temp,*id_press;         // Strings with identifiers of
  char *id_pe;                     // created computes.

  class Compute *temperature;      // Computes created to evaluate 
  class Compute *pressure;         // thermodynamic quantities.
  class Compute *pe;
  int tflag,pflag,vsflag,peflag;   // Flags to keep track of computes that
                                   // were created.

  // shock initial conditions.

  double e0;                       // Initial energy
  double v0;                       // Initial volume
  double p0;                       // Initial pressure
  double velocity;                 // Velocity of the shock.
  double lagrangian_position;      // Lagrangian location of computational cell  
  int direction;                   // Direction of shock
  int p0_set;                      // Is pressure set.
  int v0_set;                      // Is volume set.
  int e0_set;                      // Is energy set.
    
  int atoms_allocated;             // The number of allocated atoms in old_velocity.

  // functions 

  void couple();
  void remap(int);
  void check_alloc(int n);
  double compute_etotal();
  double compute_vol();
  double compute_hugoniot();
  double compute_rayleigh();
  double compute_lagrangian_speed();
  double compute_lagrangian_position();
  double compute_vsum();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix msst tscale must satisfy 0 <= tscale < 1

Self-explanatory.

E: Fix msst requires a periodic box

Self-explanatory.

E: Cannot use fix msst without per-type mass defined

Self-explanatory.

E: Could not find fix msst compute ID

Self-explanatory.

E: Fix msst compute ID does not compute temperature

Self-explanatory.

E: Fix msst compute ID does not compute pressure

Self-explanatory.

E: Fix msst compute ID does not compute potential energy

Self-explanatory.

E: Could not find fix_modify temperature ID

The compute ID for computing temperature does not exist.

E: Fix_modify temperature ID does not compute temperature

The compute ID assigned to the fix must compute temperature.

W: Temperature for MSST is not for group all

User-assigned temperature to MSST fix does not compute temperature for
all atoms.  Since MSST computes a global pressure, the kinetic energy
contribution from the temperature is assumed to also be for all atoms.
Thus the pressure used by MSST could be inaccurate.

E: Could not find fix_modify pressure ID

The compute ID for computing pressure does not exist.

E: Fix_modify pressure ID does not compute pressure

The compute ID assigned to the fix must compute pressure.

*/
