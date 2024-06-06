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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(create_atoms,CreateAtoms);
// clang-format on
#else

#ifndef LMP_CREATE_ATOMS_H
#define LMP_CREATE_ATOMS_H

#include "command.h"

namespace LAMMPS_NS {

class CreateAtoms : public Command {
 public:
  CreateAtoms(class LAMMPS *);
  void command(int, char **) override;

 private:
  int ntype, style, mode, nbasis, nrandom, seed;
  int remapflag;
  int maxtry;
  int quat_user;
  int overlapflag;
  double overlap;
  int subsetflag;
  bigint nsubset;
  double subsetfrac;
  int *basistype;
  double xone[3], quatone[4], **xmol;
  double radthresh, radscale, mesh_density;

  int varflag, vvar, xvar, yvar, zvar;
  char *vstr, *xstr, *ystr, *zstr;
  char *xstr_copy, *ystr_copy, *zstr_copy;

  int ilo, ihi, jlo, jhi, klo, khi;

  int nlatt;             // number of owned lattice sites
  int nlatt_overflow;    // 1 if local nlatt exceeds a 32-bit int

  int *flag;    // flag subset of particles to insert on lattice
  int *next;
  int mesh_style;

  class Region *region;
  class Molecule *onemol;
  class RanMars *ranmol;
  class RanMars *ranlatt;

  int triclinic;
  double sublo[3], subhi[3];    // epsilon-extended proc sub-box for adding atoms

  void add_single();
  void add_random();
  void add_mesh(const char *);
  int add_bisection(const double[3][3], tagint);
  int add_quasirandom(const double[3][3], tagint);
  void add_lattice();
  void loop_lattice(int);
  void get_xmol(double *);
  void add_molecule();
  int vartest(double *);    // evaluate a variable with new atom position
};

}    // namespace LAMMPS_NS

#endif
#endif
