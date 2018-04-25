/* -*- c++ -*- ----------------------------------------------------------
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

FixStyle(mscg,FixMSCG)

#else

#ifndef LMP_FIX_MSCG_H
#define LMP_FIX_MSCG_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMSCG : public Fix {
 public:
  FixMSCG(class LAMMPS *, int, char **);
  ~FixMSCG();
  int setmask();
  void post_constructor();
  void init();
  void end_of_step();
  void post_run();

 private:
  int range_flag,name_flag,me,nprocs;
  int nframes,n_frames,block_size,n_cg_sites,n_cg_types,*cg_site_types;
  int max_partners_bond,max_partners_angle,max_partners_dihedral;
  unsigned *n_partners_bond,*n_partners_angle,*n_partners_dihedral;
  unsigned **partners_bond,**partners_angle,**partners_dihedral;
  double *x1d,*f1d,**f;
  double box_half_lengths[3];
  char **type_names;
  void *mscg_struct;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix mscg does not yet support parallel use via MPI

UNDOCUMENTED

E: Fix mscg must be used with 32-bit atom IDs

UNDOCUMENTED

E: Fix mscg does not yet support triclinic geometries

Self-explanatory.

E: Bond list overflow, boost fix_mscg max

UNDOCUMENTED

E: Angle list overflow, boost fix_mscg max

UNDOCUMENTED

E: Dihedral list overflow, boost fix_mscg max

UNDOCUMENTED

W: Fix mscg n_frames is inconsistent with control.in

The control.in file read by the MSCG lib has a parameter n_frames
that should be equal to the number of frames processed by the
fix mscg command. If not equal, the fix will still run, but the
calculated residuals may be normalized incorrectly.

W: Fix mscg n_frames is not divisible by block_size in control.in

The control.in file read by the MSCG lib has a parameter block_size
that should be a divisor of the number of frames processed by the
fix mscg command. If not, the fix will still run, but some frames may
not be included in the MSCG calculations.

U: Fix mscg does not yet support mpi

Self-explanatory.

U: Bond/Angle/Dihedral list overflow, boost fix_mscg max

A site has more bond/angle/dihedral partners that the maximum and
has overflowed the bond/angle/dihedral partners list. Increase the
corresponding fix_mscg max arg.

*/
