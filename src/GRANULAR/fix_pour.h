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

FixStyle(pour,FixPour)

#else

#ifndef LMP_FIX_POUR_H
#define LMP_FIX_POUR_H

#include "fix.h"
#include "near_list.h"

namespace LAMMPS_NS {

class FixPour : public Fix {
 public:
  FixPour(class LAMMPS *, int, char **);
  ~FixPour();
  int setmask();
  void init();
  void pre_exchange();
  void reset_dt();
  void *extract(const char *, int &);

 protected:
  virtual void init_near_lists(int nnew, INearList*&, IDistributedNearList*&);
  virtual void cleanup_near_lists(INearList*&, IDistributedNearList*&);

  int ninsert,ntype,seed;
  int iregion,mode,idnext,dstyle,npoly,rigidflag,shakeflag;
  int ignoreflag,ignoreline,ignoretri;
  double radius_one,radius_max;
  double radius_lo,radius_hi;
  double *radius_poly,*frac_poly;
  double density_lo,density_hi;
  double volfrac;
  int maxattempt;
  int region_style;
  double rate;
  double vxlo,vxhi,vylo,vyhi,vy,vz;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double xc,yc,rc;
  double grav;
  char *idrigid,*idshake;

  class Molecule **onemols;
  int nmol,natom_max;
  double molradius_max;
  double *molfrac;
  double **coords;
  imageint *imageflags;
  class Fix *fixrigid,*fixshake;
  double oneradius;

  int nfreq,nfirst,ninserted,nper;
  double lo_current,hi_current;
  tagint maxtag_all,maxmol_all;
  class RanPark *random,*random2;

  void find_maxid();
  int overlap(int);
  bool outside(int, double, double, double);
  void xyz_random(double, double *);
  double radius_sample();
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

E: Cannot yet use fix pour with the KOKKOS package

This feature is not yet supported.

E: Fix pour requires atom attributes radius, rmass

The atom style defined does not have these attributes.

E: Invalid atom type in fix pour command

Self-explanatory.

E: Must specify a region in fix pour

Self-explanatory.

E: Fix pour region does not support a bounding box

Not all regions represent bounded volumes.  You cannot use
such a region with the fix pour command.

E: Fix pour region cannot be dynamic

Only static regions can be used with fix pour.

E: Insertion region extends outside simulation box

Self-explanatory.

E: Must use a z-axis cylinder region with fix pour

Self-explanatory.

E: Must use a block or cylinder region with fix pour

Self-explanatory.

E: Must use a block region with fix pour for 2d simulations

Self-explanatory.

E: Cannot use fix_pour unless atoms have IDs

Self-explanatory.

E: Fix pour molecule must have coordinates

The defined molecule does not specify coordinates.

E: Fix pour molecule must have atom types

The defined molecule does not specify atom types.

E: Invalid atom type in fix pour mol command

The atom types in the defined molecule are added to the value
specified in the create_atoms command, as an offset.  The final value
for each atom must be between 1 to N, where N is the number of atom
types.

E: Fix pour molecule template ID must be same as atom style template ID

When using atom_style template, you cannot pour molecules that are
not in that template.

E: Cannot use fix pour rigid and not molecule

Self-explanatory.

E: Cannot use fix pour shake and not molecule

Self-explanatory.

E: Cannot use fix pour rigid and shake

These two attributes are conflicting.

E: No fix gravity defined for fix pour

Gravity is required to use fix pour.

E: Fix pour insertion count per timestep is 0

Self-explanatory.

E: Cannot use fix pour with triclinic box

This option is not yet supported.

E: Gravity must point in -z to use with fix pour in 3d

Self-explanatory.

E: Gravity must point in -y to use with fix pour in 2d

Self-explanatory.

E: Gravity changed since fix pour was created

The gravity vector defined by fix gravity must be static.

E: Fix pour rigid fix does not exist

Self-explanatory.

E: Fix pour and fix rigid/small not using same molecule template ID

Self-explanatory.

E: Fix pour shake fix does not exist

Self-explanatory.

E: Fix pour and fix shake not using same molecule template ID

Self-explanatory.

W: Less insertions than requested

The fix pour command was unsuccessful at finding open space
for as many particles as it tried to insert.

E: Too many total atoms

See the setting for bigint in the src/lmptype.h file.

E: New atom IDs exceed maximum allowed ID

See the setting for tagint in the src/lmptype.h file.

E: Fix pour region ID does not exist

Self-explanatory.

E: Molecule template ID for fix pour does not exist

Self-explanatory.

E: Fix pour polydisperse fractions do not sum to 1.0

Self-explanatory.

E: Cannot change timestep with fix pour

This is because fix pour pre-computes the time delay for particles to
fall out of the insertion volume due to gravity.

*/
