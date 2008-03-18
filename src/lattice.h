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

#ifndef LATTICE_H
#define LATTICE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Lattice : protected Pointers {
 public:
  int style;                           // enum list of NONE,SC,FCC,etc
  double xlattice,ylattice,zlattice;   // lattice scale factors in 3 dims
  double a1[3],a2[3],a3[3];            // vectors that bound unit cell
  int nbasis;                          // # of atoms in basis of unit cell
  double **basis;                      // fractional coords of each basis atom
                                       // within unit cell (0 <= coord < 1)

  Lattice(class LAMMPS *, int, char **);
  ~Lattice();
  void lattice2box(double &, double &, double &);
  void box2lattice(double &, double &, double &);
  void bbox(int, double, double, double,
	    double &, double &, double &, double &, double &, double &);

private:
  double scale;
  double origin[3];                    // lattice origin
  int orientx[3];                      // lattice orientation vecs
  int orienty[3];                      // orientx = what lattice dir lies
  int orientz[3];                      //           along x dim in box

  double primitive[3][3];              // lattice <-> box transform matrices
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

}

#endif
