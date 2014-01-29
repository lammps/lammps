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

#ifndef LMP_DOMAIN_H
#define LMP_DOMAIN_H

#include "math.h"
#include "pointers.h"

namespace LAMMPS_NS {

class Domain : protected Pointers {
 public:
  int box_exist;                         // 0 = not yet created, 1 = exists
  int dimension;                         // 2 = 2d, 3 = 3d
  int nonperiodic;                       // 0 = periodic in all 3 dims
                                         // 1 = periodic or fixed in all 6
                                         // 2 = shrink-wrap in any of 6
  int xperiodic,yperiodic,zperiodic;     // 0 = non-periodic, 1 = periodic
  int periodicity[3];                    // xyz periodicity as array

  int boundary[3][2];                    // settings for 6 boundaries
                                         // 0 = periodic
                                         // 1 = fixed non-periodic
                                         // 2 = shrink-wrap non-periodic
                                         // 3 = shrink-wrap non-per w/ min

  int triclinic;                         // 0 = orthog box, 1 = triclinic
  int tiltsmall;                         // 1 if limit tilt, else 0

                                         // orthogonal box
  double xprd,yprd,zprd;                 // global box dimensions
  double xprd_half,yprd_half,zprd_half;  // half dimensions
  double prd[3];                         // array form of dimensions
  double prd_half[3];                    // array form of half dimensions

                                         // triclinic box
                                         // xprd,xprd_half,prd,prd_half =
                                         // same as if untilted
  double prd_lamda[3];                   // lamda box = (1,1,1)
  double prd_half_lamda[3];              // lamda half box = (0.5,0.5,0.5)

  double boxlo[3],boxhi[3];              // orthogonal box global bounds

                                         // triclinic box
                                         // boxlo/hi = same as if untilted
  double boxlo_lamda[3],boxhi_lamda[3];  // lamda box = (0,1)
  double boxlo_bound[3],boxhi_bound[3];  // bounding box of tilted domain
  double corners[8][3];                  // 8 corner points

                                         // orthogonal box & triclinic box
  double minxlo,minxhi;                  // minimum size of global box
  double minylo,minyhi;                  // when shrink-wrapping
  double minzlo,minzhi;                  // tri only possible for non-skew dims

                                         // orthogonal box
  double sublo[3],subhi[3];              // sub-box bounds on this proc

                                         // triclinic box
                                         // sublo/hi = undefined
  double sublo_lamda[3],subhi_lamda[3];  // bounds of subbox in lamda

                                         // triclinic box
  double xy,xz,yz;                       // 3 tilt factors
  double h[6],h_inv[6];                  // shape matrix in Voigt notation
  double h_rate[6],h_ratelo[3];          // rate of box size/shape change

  int box_change;                // 1 if any of next 3 flags are set, else 0
  int box_change_size;           // 1 if box size changes, 0 if not
  int box_change_shape;          // 1 if box shape changes, 0 if not
  int box_change_domain;         // 1 if proc sub-domains change, 0 if not

  int deform_flag;                // 1 if fix deform exist, else 0
  int deform_vremap;              // 1 if fix deform remaps v, else 0
  int deform_groupbit;            // atom group to perform v remap for

  class Lattice *lattice;                  // user-defined lattice

  int nregion;                             // # of defined Regions
  int maxregion;                           // max # list can hold
  class Region **regions;                  // list of defined Regions

  Domain(class LAMMPS *);
  virtual ~Domain();
  virtual void init();
  void set_initial_box();
  virtual void set_global_box();
  virtual void set_lamda_box();
  virtual void set_local_box();
  virtual void reset_box();
  virtual void pbc();
  void image_check();
  void box_too_small_check();
  void minimum_image(double &, double &, double &);
  void minimum_image(double *);
  int closest_image(int, int);
  void closest_image(const double * const, const double * const,
                     double * const);
  void remap(double *, imageint &);
  void remap(double *);
  void remap_near(double *, double *);
  void unmap(double *, imageint);
  void unmap(double *, imageint, double *);
  void image_flip(int, int, int);
  void set_lattice(int, char **);
  void add_region(int, char **);
  void delete_region(int, char **);
  int find_region(char *);
  void set_boundary(int, char **, int);
  void set_box(int, char **);
  void print_box(const char *);
  void boundary_string(char *);

  virtual void lamda2x(int);
  virtual void x2lamda(int);
  virtual void lamda2x(double *, double *);
  virtual void x2lamda(double *, double *);
  void x2lamda(double *, double *, double *, double *);
  void bbox(double *, double *, double *, double *);
  void box_corners();

  // minimum image convention check
  // return 1 if any distance > 1/2 of box size
  // indicates a special neighbor is actually not in a bond, 
  //   but is a far-away image that should be treated as an unbonded neighbor
  // inline since called from neighbor build inner loop
  // 
  inline int minimum_image_check(double dx, double dy, double dz) {
    if (xperiodic && fabs(dx) > xprd_half) return 1;
    if (yperiodic && fabs(dy) > yprd_half) return 1;
    if (zperiodic && fabs(dz) > zprd_half) return 1;
    return 0;
  }

 protected:
  double small[3];                  // fractions of box lengths
};

}

#endif

/* ERROR/WARNING messages:

E: Box bounds are invalid

The box boundaries specified in the read_data file are invalid.  The
lo value must be less than the hi value for all 3 dimensions.

E: Cannot skew triclinic box in z for 2d simulation

Self-explanatory.

E: Triclinic box skew is too large

The displacement in a skewed direction must be less than half the box
length in that dimension.  E.g. the xy tilt must be between -half and
+half of the x box length.  This constraint can be relaxed by using
the box tilt command.

W: Triclinic box skew is large

The displacement in a skewed direction is normally required to be less
than half the box length in that dimension.  E.g. the xy tilt must be
between -half and +half of the x box length.  You have relaxed the
constraint using the box tilt command, but the warning means that a
LAMMPS simulation may be inefficient as a result.

E: Illegal simulation box

The lower bound of the simulation box is greater than the upper bound.

E: Bond atom missing in image check

The 2nd atom in a particular bond is missing on this processor.
Typically this is because the pairwise cutoff is set too short or the
bond has blown apart and an atom is too far away.

W: Inconsistent image flags

The image flags for a pair on bonded atoms appear to be inconsistent.
Inconsistent means that when the coordinates of the two atoms are
unwrapped using the image flags, the two atoms are far apart.
Specifically they are further apart than half a periodic box length.
Or they are more than a box length apart in a non-periodic dimension.
This is usually due to the initial data file not having correct image
flags for the 2 atoms in a bond that straddles a periodic boundary.
They should be different by 1 in that case.  This is a warning because
inconsistent image flags will not cause problems for dynamics or most
LAMMPS simulations.  However they can cause problems when such atoms
are used with the fix rigid or replicate commands.

W: Bond atom missing in image check

The 2nd atom in a particular bond is missing on this processor.
Typically this is because the pairwise cutoff is set too short or the
bond has blown apart and an atom is too far away.

E: Bond atom missing in box size check

The 2nd atoms needed to compute a particular bond is missing on this
processor.  Typically this is because the pairwise cutoff is set too
short or the bond has blown apart and an atom is too far away.

W: Bond atoms missing in box size check

The 2nd atoms needed to compute a particular bond is missing on this
processor.  Typically this is because the pairwise cutoff is set too
short or the bond has blown apart and an atom is too far away.

W: Bond/angle/dihedral extent > half of periodic box length

This is a restriction because LAMMPS can be confused about which image
of an atom in the bonded interaction is the correct one to use.
"Extent" in this context means the maximum end-to-end length of the
bond/angle/dihedral.  LAMMPS computes this by taking the maximum bond
length, multiplying by the number of bonds in the interaction (e.g. 3
for a dihedral) and adding a small amount of stretch.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Reuse of region ID

A region ID cannot be used twice.

E: Invalid region style

The choice of region style is unknown.

E: Delete region ID does not exist

Self-explanatory.

E: Both sides of boundary must be periodic

Cannot specify a boundary as periodic only on the lo or hi side.  Must
be periodic on both sides.

*/
