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

//
// Contributing author, Richard Meng, Queen's University at Kingston, 22.11.24, contact@richardzjm.com
//

#ifndef LMP_MTP_RADIAL_BASIS_H
#define LMP_MTP_RADIAL_BASIS_H

#include "pointers.h"

namespace LAMMPS_NS {

class TextFileReader;

class RadialMTPBasis : protected Pointers {
 public:
  RadialMTPBasis(TextFileReader &tfr, LAMMPS *lmp);
  RadialMTPBasis(int size, LAMMPS *lmp);
  virtual ~RadialMTPBasis();

  virtual void calc_radial_basis(double dist) = 0;
  virtual void calc_radial_basis_ders(double dist) = 0;

  int size;             // The size of the radial basis functions
  double min_cutoff;    //  Minimum radius value
  double max_cutoff;    // Cutoff radius
  double scaling;       // All radial functions are multiplied by scaling

  // Values and derivatives for radial basis functions
  double *radial_basis_vals;
  double *radial_basis_ders;

 protected:
 private:
  //Specifically reads the basis properties (ie. cutoffs and size) and not the radial parameters
  void read_basis_properties(TextFileReader &tfr);
};

class RBChebyshev : public RadialMTPBasis {
 public:
  RBChebyshev(int size, LAMMPS *lmp) : RadialMTPBasis(size, lmp) {}
  RBChebyshev(TextFileReader &tfr, LAMMPS *lmp) : RadialMTPBasis(tfr, lmp) {}
  void calc_radial_basis(double val) override;
  void calc_radial_basis_ders(double val) override;
};

}    // namespace LAMMPS_NS

#endif    // LMP_MTP_RADIAL_BASIS_H