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

#ifndef LMP_LATTICE_H
#define LMP_LATTICE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Lattice : protected Pointers {
 public:
  enum { NONE, SC, BCC, FCC, HCP, DIAMOND, SQ, SQ2, HEX, CUSTOM };

  int style;                              // NONE,SC,FCC,etc
  double xlattice, ylattice, zlattice;    // lattice scale factors in 3 dims
  double a1[3], a2[3], a3[3];             // edge vectors of unit cell
  int nbasis;                             // # of basis atoms in unit cell
  double **basis;                         // fractional coords of each basis atom
                                          // within unit cell (0 <= coord < 1)

  Lattice(class LAMMPS *, int, char **);
  ~Lattice();
  void lattice2box(double &, double &, double &);
  void box2lattice(double &, double &, double &);
  void bbox(int, double, double, double, double &, double &, double &, double &, double &,
            double &);

 private:
  double scale;
  double origin[3];    // lattice origin
  int orientx[3];      // lattice orientation vecs
  int orienty[3];      // orientx = what lattice dir lies
  int orientz[3];      //           along x dim in box

  double primitive[3][3];    // lattice <-> box transform matrices
  double priminv[3][3];
  double rotaterow[3][3];
  double rotatecol[3][3];

  int orthogonal();
  int right_handed();
  int collinear();
  void setup_transform();
  void add_basis(double, double, double);
  double dot(double *, double *);
  void cross(double *, double *, double *);
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Lattice style incompatible with simulation dimension

2d simulation can use sq, sq2, or hex lattice.  3d simulation can use
sc, bcc, or fcc lattice.

E: Invalid option in lattice command for non-custom style

Certain lattice keywords are not supported unless the
lattice style is "custom".

E: No basis atoms in lattice

Basis atoms must be defined for lattice style user.

E: Lattice orient vectors are not orthogonal

The three specified lattice orientation vectors must be mutually
orthogonal.

E: Lattice orient vectors are not right-handed

The three specified lattice orientation vectors must create a
right-handed coordinate system such that a1 cross a2 = a3.

E: Lattice primitive vectors are collinear

The specified lattice primitive vectors do not for a unit cell with
non-zero volume.

E: Lattice settings are not compatible with 2d simulation

One or more of the specified lattice vectors has a non-zero z
component.

E: Lattice spacings are invalid

Each x,y,z spacing must be > 0.

E: Degenerate lattice primitive vectors

Invalid set of 3 lattice vectors for lattice command.

E: Zero-length lattice orient vector

Self-explanatory.

*/
