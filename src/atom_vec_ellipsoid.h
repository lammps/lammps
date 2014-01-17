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
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  int pack_comm(int, int *, double *, int, int *);
  int pack_comm_vel(int, int *, double *, int, int *);
  int pack_comm_hybrid(int, int *, double *);
  void unpack_comm(int, int, double *);
  void unpack_comm_vel(int, int, double *);
  int unpack_comm_hybrid(int, int, double *);
  int pack_reverse(int, int, double *);
  int pack_reverse_hybrid(int, int, double *);
  void unpack_reverse(int, int *, double *);
  int unpack_reverse_hybrid(int, int *, double *);
  int pack_border(int, int *, double *, int, int *);
  int pack_border_vel(int, int *, double *, int, int *);
  int pack_border_hybrid(int, int *, double *);
  void unpack_border(int, int, double *);
  void unpack_border_vel(int, int, double *);
  int unpack_border_hybrid(int, int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  void data_atom(double *, imageint, char **);
  int data_atom_hybrid(int, char **);
  void data_vel(int, char **);
  int data_vel_hybrid(int, char **);
  void pack_data(double **);
  int pack_data_hybrid(int, double *);
  void write_data(FILE *, int, double **);
  int write_data_hybrid(FILE *, double *);
  void pack_vel(double **);
  int pack_vel_hybrid(int, double *);
  void write_vel(FILE *, int, double **);
  int write_vel_hybrid(FILE *, double *);
  bigint memory_usage();

  // manipulate Bonus data structure for extra atom info

  void clear_bonus();
  void data_atom_bonus(int, char **);

  // unique to AtomVecEllipsoid

  void set_shape(int, double, double, double);

 private:
  tagint *tag;
  int *type,*mask;
  imageint *image;
  double **x,**v,**f;
  double *rmass;
  double **angmom,**torque;
  int *ellipsoid;

  int nlocal_bonus,nghost_bonus,nmax_bonus;

  void grow_bonus();
  void copy_bonus(int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom ID in Atoms section of data file

Atom IDs must be positive integers.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

E: Invalid density in Atoms section of data file

Density value cannot be <= 0.0.

E: Assigning ellipsoid parameters to non-ellipsoid atom

Self-explanatory.

E: Invalid shape in Ellipsoids section of data file

Self-explanatory.

*/
