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

CommandStyle(change_box,ChangeBox)

#else

#ifndef LMP_CHANGE_BOX_H
#define LMP_CHANGE_BOX_H

#include "pointers.h"

namespace LAMMPS_NS {

class ChangeBox : protected Pointers {
 public:
  ChangeBox(class LAMMPS *);
  void command(int, char **);

 private:
  int scaleflag;
  double scale[3];

  struct Operation {
    int style,flavor;
    int dim,boundindex;
    int vdim1,vdim2;
    double flo,fhi,ftilt;
    double dlo,dhi,dtilt;
    double scale;
  };

  Operation *ops;
  int nops;

  double boxlo[3],h_inv[6];

  void options(int, char **);
  void save_box_state();
  void volume_preserve(int, int, double);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Change_box command before simulation box is defined

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot change_box after reading restart file with per-atom info

UNDOCUMENTED

E: Could not find change_box group ID

UNDOCUMENTED

E: Cannot change_box in z dimension for 2d simulation

UNDOCUMENTED

E: Change_box volume used incorrectly

UNDOCUMENTED

E: Cannot change_box in xz or yz for 2d simulation

UNDOCUMENTED

E: Use of change_box with undefined lattice

UNDOCUMENTED

E: Cannot change box tilt factors for orthogonal box

UNDOCUMENTED

E: Cannot run 2d simulation with nonperiodic Z dimension

UNDOCUMENTED

E: Cannot change box to orthogonal when tilt is non-zero

UNDOCUMENTED

E: Cannot change box ortho/triclinic with dumps defined

UNDOCUMENTED

E: Cannot change box ortho/triclinic with certain fixes defined

UNDOCUMENTED

W: Lost atoms via change_box: original %ld current %ld

UNDOCUMENTED

U: Displace_box command before simulation box is defined

Self-explanatory.

U: Cannot displace_box after reading restart file with per-atom info

This is because the restart file info cannot be migrated with the
atoms.  You can get around this by performing a 0-timestep run which
will assign the restart file info to actual atoms.

U: Could not find displace_box group ID

Group ID used in the displace_box command does not exist.

U: Displace_box tilt factors require triclinic box

Cannot use tilt factors unless the simulation box is
non-orthogonal.

U: Cannot displace_box on a non-periodic boundary

Self-explanatory.

U: Use of displace_box with undefined lattice

Must use lattice command with displace_box command if units option is
set to lattice.

U: Fix deform volume setting is invalid

Cannot use volume style unless other dimensions are being controlled.

U: Induced tilt by displace_box is too large

The final tilt value must be between -1/2 and 1/2 of the perpendicular
box length.

U: Lost atoms via displace_box: original %ld current %ld

UNDOCUMENTED

*/
