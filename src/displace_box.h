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

#ifdef COMMAND_CLASS

CommandStyle(displace_box,DisplaceBox)

#else

#ifndef LMP_DISPLACE_BOX_H
#define LMP_DISPLACE_BOX_H

#include "pointers.h"

namespace LAMMPS_NS {

class DisplaceBox : protected Pointers {
 public:
  DisplaceBox(class LAMMPS *);
  void command(int, char **);

 private:
  int remapflag,scaleflag;

  struct Set {
    int style,substyle;
    double flo,fhi,ftilt;
    double dlo,dhi,dtilt;
    double scale;
    double lo_start,hi_start;
    double lo_stop,hi_stop;
    double tilt_start,tilt_stop;
    double vol_start;
    int fixed,dynamic1,dynamic2;
  };
  Set *set;

  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Displace_box command before simulation box is defined

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot displace_box after reading restart file with per-atom info

This is because the restart file info cannot be migrated with the
atoms.  You can get around this by performing a 0-timestep run which
will assign the restart file info to actual atoms.

E: Could not find displace_box group ID

Group ID used in the displace_box command does not exist.

E: Displace_box tilt factors require triclinic box

Cannot use tilt factors unless the simulation box is
non-orthogonal.

E: Cannot displace_box on a non-periodic boundary

Self-explanatory.

E: Use of displace_box with undefined lattice

Must use lattice command with displace_box command if units option is
set to lattice.

E: Fix deform volume setting is invalid

Cannot use volume style unless other dimensions are being controlled.

E: Induced tilt by displace_box is too large

The final tilt value must be between -1/2 and 1/2 of the perpendicular
box length.

E: Lost atoms via displace_box: original %ld current %ld

UNDOCUMENTED

*/
