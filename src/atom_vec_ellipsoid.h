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

#ifdef ATOM_CLASS

AtomStyle(ellipsoid,AtomVecEllipsoid)

#else

#ifndef LMP_ATOM_VEC_ELLIPSOID_H
#define LMP_ATOM_VEC_ELLIPSOID_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecEllipsoid : public AtomVec {
 public:
  struct Bonus {
    double shape[3];
    double quat[4];
    int ilocal;
  };
  struct Bonus *bonus;

  AtomVecEllipsoid(class LAMMPS *);
  ~AtomVecEllipsoid();

  void copy_bonus(int, int, int);
  void clear_bonus();
  int pack_comm_bonus(int, int *, double *);
  void unpack_comm_bonus(int, int, double *);
  int pack_border_bonus(int, int *, double *);
  int unpack_border_bonus(int, int, double *);
  int pack_exchange_bonus(int, double *);
  int unpack_exchange_bonus(int, double *);
  int size_restart_bonus();
  int pack_restart_bonus(int, double *);
  int unpack_restart_bonus(int, double *);
  void data_atom_bonus(int, char **);
  bigint memory_usage_bonus();

  void create_atom_post(int);
  void data_atom_post(int);
  void pack_data_pre(int);
  void pack_data_post(int);

  // unique to AtomVecEllipsoid

  void set_shape(int, double, double, double);

  int nlocal_bonus;

 private:
  int nghost_bonus,nmax_bonus;
  int ellipsoid_flag;
  double rmass;

  void grow_bonus();
  void copy_bonus_all(int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

E: Invalid density in Atoms section of data file

Density value cannot be <= 0.0.

E: Assigning ellipsoid parameters to non-ellipsoid atom

Self-explanatory.

E: Invalid shape in Ellipsoids section of data file

Self-explanatory.

*/
