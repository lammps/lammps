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
  int kspace_flag;

  int current_lifo;              // LIFO stack pointer
  double boxlo0[2][3];           // box bounds at start of line search
  double boxhi0[2][3];
  double boxtilt0[2][3];         // xy,xz,yz tilts at start of line search
  double ds[6];                  // increment in scale matrix

  int scaleyz;                   // 1 if yz scaled with lz
  int scalexz;                   // 1 if xz scaled with lz
  int scalexy;                   // 1 if xy scaled with ly

  double fixedpoint[3];          // Location of dilation fixed-point

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

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid fix box/relax command for a 2d simulation

Fix box/relax styles involving the z dimension cannot be used in
a 2d simulation.

E: Invalid fix box/relax command pressure settings

If multiple dimensions are coupled, those dimensions must be specified.

E: Cannot use fix box/relax on a non-periodic dimension

When specifying a diagonal pressure component, the dimension must be
periodic.

E: Cannot use fix box/relax on a 2nd non-periodic dimension

When specifying an off-diagonal pressure component, the 2nd of the two
dimensions must be periodic.  E.g. if the xy component is specified,
then the y dimension must be periodic.

E: Cannot use fix box/relax with tilt factor scaling on a 2nd non-periodic dimension

When specifying scaling on a tilt factor component, the 2nd of the two
dimensions must be periodic.  E.g. if the xy component is specified,
then the y dimension must be periodic.

E: Cannot use fix box/relax with both relaxation and scaling on a tilt factor

When specifying scaling on a tilt factor component, that component can not
also be controlled by the barostat. E.g. if scalexy yes is specified and
also keyword tri or xy, this is wrong.

E: Can not specify Pxy/Pxz/Pyz in fix box/relax with non-triclinic box

Only triclinic boxes can be used with off-diagonal pressure components.
See the region prism command for details.

E: Invalid fix box/relax pressure settings

Settings for coupled dimensions must be the same.

E: Temperature ID for fix box/relax does not exist

Self-explanatory.

E: Pressure ID for fix box/relax does not exist

The compute ID needed to compute pressure for the fix does not
exist.

E: Attempt to push beyond stack limit in fix box/relax

Internal LAMMPS error.  Please report it to the developers.

E: Attempt to pop empty stack in fix box/relax

Internal LAMMPS error.  Please report it to the developers.

E: Fix box/relax generated negative box length

The pressure being applied is likely too large.  Try applying
it incrementally, to build to the high pressure.

E: Could not find fix_modify temperature ID

The compute ID for computing temperature does not exist.

E: Fix_modify temperature ID does not compute temperature

The compute ID assigned to the fix must compute temperature.

W: Temperature for fix modify is not for group all

The temperature compute is being used with a pressure calculation
which does operate on group all, so this may be inconsistent.

E: Pressure ID for fix modify does not exist

Self-explanatory.

E: Could not find fix_modify pressure ID

The compute ID for computing pressure does not exist.

E: Fix_modify pressure ID does not compute pressure

The compute ID assigned to the fix must compute pressure.

*/
