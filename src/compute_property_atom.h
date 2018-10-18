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

#ifdef COMPUTE_CLASS

ComputeStyle(property/atom,ComputePropertyAtom)

#else

#ifndef LMP_COMPUTE_PROPERTY_ATOM_H
#define LMP_COMPUTE_PROPERTY_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePropertyAtom : public Compute {
 public:
  ComputePropertyAtom(class LAMMPS *, int, char **);
  ~ComputePropertyAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nvalues;
  int nmax;
  int *index;
  double *buf;
  class AtomVecEllipsoid *avec_ellipsoid;
  class AtomVecLine *avec_line;
  class AtomVecTri *avec_tri;
  class AtomVecBody *avec_body;

  typedef void (ComputePropertyAtom::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_id(int);
  void pack_molecule(int);
  void pack_proc(int);
  void pack_type(int);
  void pack_mass(int);

  void pack_x(int);
  void pack_y(int);
  void pack_z(int);
  void pack_xs(int);
  void pack_ys(int);
  void pack_zs(int);
  void pack_xs_triclinic(int);
  void pack_ys_triclinic(int);
  void pack_zs_triclinic(int);
  void pack_xu(int);
  void pack_yu(int);
  void pack_zu(int);
  void pack_xu_triclinic(int);
  void pack_yu_triclinic(int);
  void pack_zu_triclinic(int);
  void pack_ix(int);
  void pack_iy(int);
  void pack_iz(int);

  void pack_vx(int);
  void pack_vy(int);
  void pack_vz(int);
  void pack_fx(int);
  void pack_fy(int);
  void pack_fz(int);
  void pack_q(int);
  void pack_mux(int);
  void pack_muy(int);
  void pack_muz(int);
  void pack_mu(int);
  void pack_radius(int);
  void pack_diameter(int);

  void pack_spx(int);		// pack magnetic variables
  void pack_spy(int);
  void pack_spz(int);
  void pack_sp(int);
  void pack_fmx(int);
  void pack_fmy(int);
  void pack_fmz(int);

  void pack_omegax(int);
  void pack_omegay(int);
  void pack_omegaz(int);
  void pack_angmomx(int);
  void pack_angmomy(int);
  void pack_angmomz(int);
  void pack_shapex(int);
  void pack_shapey(int);
  void pack_shapez(int);
  void pack_quatw(int);
  void pack_quati(int);
  void pack_quatj(int);
  void pack_quatk(int);
  void pack_tqx(int);
  void pack_tqy(int);
  void pack_tqz(int);
  void pack_end1x(int);
  void pack_end1y(int);
  void pack_end1z(int);
  void pack_end2x(int);
  void pack_end2y(int);
  void pack_end2z(int);
  void pack_corner1x(int);
  void pack_corner1y(int);
  void pack_corner1z(int);
  void pack_corner2x(int);
  void pack_corner2y(int);
  void pack_corner2z(int);
  void pack_corner3x(int);
  void pack_corner3y(int);
  void pack_corner3z(int);

  void pack_nbonds(int);

  void pack_iname(int);
  void pack_dname(int);

  void pack_property_atom(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute property/atom for atom property that isn't allocated

Self-explanatory.

E: Compute property/atom integer vector does not exist

The command is accessing a vector added by the fix property/atom
command, that does not exist.

E: Compute property/atom floating point vector does not exist

The command is accessing a vector added by the fix property/atom
command, that does not exist.

E: Invalid keyword in compute property/atom command

Self-explanatory.

*/
