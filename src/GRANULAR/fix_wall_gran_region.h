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

FixStyle(wall/gran/region,FixWallGranRegion)

#else

#ifndef LMP_FIX_WALL_GRAN_REGION_H
#define LMP_FIX_WALL_GRAN_REGION_H

#include "fix_wall_gran.h"

namespace LAMMPS_NS {

class FixWallGranRegion : public FixWallGran {
 public:
  FixWallGranRegion(class LAMMPS *, int, char **);
  ~FixWallGranRegion();
  void post_force(int);
  void write_restart(FILE *);
  void restart(char* );
  void init();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

 private:
  class Region *region;
  char *region_style;
  int nregion;

  // shear history for multiple contacts per particle

  int tmax;              // max # of region walls one particle can touch
  int *ncontact;         // # of shear contacts per particle
  int **walls;           // which wall each contact is with
  double ***shearmany;   // shear history per particle per contact
  int *c2r;              // contact to region mapping
                         // c2r[i] = index of Ith contact in
                         //   region-contact[] list of contacts
  int motion_resetflag;  // used by restart to indicate that region
                         //    vel info is to be reset

  void update_contacts(int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Region ID for fix wall/gran/region does not exist

UNDOCUMENTED

W: Region properties for region %s changed between runs, resetting its motion

UNDOCUMENTED

W: Region properties for region %s are inconsistent with restart file, resetting its motion

UNDOCUMENTED

E: Too many wall/gran/region contacts for one particle

UNDOCUMENTED

U: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

U: Fix wall/gran requires atom style sphere

Self-explanatory.

U: Cannot use wall in periodic dimension

Self-explanatory.

U: Cannot wiggle and shear fix wall/gran

Cannot specify both options at the same time.

U: Invalid wiggle direction for fix wall/gran

Self-explanatory.

U: Invalid shear direction for fix wall/gran

Self-explanatory.

U: Fix wall/gran is incompatible with Pair style

Must use a granular pair style to define the parameters needed for
this fix.

*/
