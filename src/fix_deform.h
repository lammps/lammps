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

FixStyle(deform,FixDeform)

#else

#ifndef LMP_FIX_DEFORM_H
#define LMP_FIX_DEFORM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDeform : public Fix {
 public:
  int remapflag;                   // whether x,v are remapped across PBC
  int dimflag[6];                  // which dims are deformed

  FixDeform(class LAMMPS *, int, char **);
  ~FixDeform();
  int setmask();
  void init();
  void pre_exchange();
  void end_of_step();

 private:
  int triclinic,scaleflag;
  int flip,flipx,flipy;
  double *h_rate,*h_ratelo;
  int varflag;                     // 1 if VARIABLE option is used, 0 if not
  int kspace_flag;                 // 1 if KSpace invoked, 0 if not
  int nrigid;                      // number of rigid fixes
  int *rfix;                       // indices of rigid fixes
  class Irregular *irregular;      // for migrating atoms after box flips

  double TWOPI;

  struct Set {
    int style,substyle;
    double flo,fhi,ftilt;
    double dlo,dhi,dtilt;
    double scale,vel,rate;
    double amplitude,tperiod;
    double lo_initial,hi_initial;
    double lo_start,hi_start,lo_stop,hi_stop,lo_target,hi_target;
    double tilt_initial,tilt_start,tilt_stop,tilt_target,tilt_flip;
    double tilt_min,tilt_max;
    double vol_initial,vol_start;
    int fixed,dynamic1,dynamic2;
    char *hstr,*hratestr;
    int hvar,hratevar;
  };
  Set *set;

  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix deform tilt factors require triclinic box

Cannot deform the tilt factors of a simulation box unless it
is a triclinic (non-orthogonal) box.

E: Cannot use fix deform on a non-periodic boundary

When specifying a change is a box dimension, the dimension must be
periodic.

E: Cannot use fix deform on a 2nd non-periodic boundary

When specifying a tilt factor change, the 2nd of the two dimensions
must be periodic.  E.g. if the xy tilt is specified, then the y
dimension must be periodic.

E: Use of fix deform with undefined lattice

A lattice must be defined to use fix deform with units = lattice.

E: Fix deform volume setting is invalid

Cannot use volume style unless other dimensions are being controlled.

E: More than one fix deform

Only one fix deform can be defined at a time.

E: Final box dimension due to fix deform is < 0.0

Self-explanatory.

E: Cannot use fix deform trate on a box with zero tilt

The trate style alters the current strain.

E: Fix deform is changing yz by too much with changing xy

When both yz and xy are changing, it induces changes in xz if the
box must flip from one tilt extreme to another.  Thus it is not
allowed for yz to grow so much that a flip is induced.

*/
