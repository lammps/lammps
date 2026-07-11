/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
  ~Lattice() override;
  void lattice2box(double &, double &, double &);
  void box2lattice(double &, double &, double &);
  void bbox(int, double, double, double, double &, double &, double &, double &, double &,
            double &);
  int is_general_triclinic();

 private:
  int triclinic_general;    // 1 if general triclinic, else 0
  int oriented;             // 1 if non-default orient xyz, else 0
  double scale;
  double origin[3];    // lattice origin
  int orientx[3];      // lattice orientation vecs
  int orienty[3];      // orientx = what lattice dir lies
  int orientz[3];      //           along x dim in box

  double primitive[3][3];    // lattice <-> box transformation matrices
  double priminv[3][3];
  double rotaterow[3][3];
  double rotatecol[3][3];

  double a1_prime[3];    // a123 rotated to restricted triclinic orientation
  double a2_prime[3];
  double a3_prime[3];

  int orthogonal();
  int right_handed_orientation();
  int right_handed_primitive();
  int collinear();
  void setup_transform(double *, double *, double *);
  void add_basis(double, double, double);
  double dot(double *, double *);
  void cross(double *, double *, double *);
};

}    // namespace LAMMPS_NS

#endif
