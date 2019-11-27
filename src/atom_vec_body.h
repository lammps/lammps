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

AtomStyle(body,AtomVecBody)

#else

#ifndef LMP_ATOM_VEC_BODY_H
#define LMP_ATOM_VEC_BODY_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecBody : public AtomVec {
 public:
  class Body *bptr;

  struct Bonus {
    double quat[4];
    double inertia[3];
    int ninteger,ndouble;
    int iindex,dindex;
    int *ivalue;
    double *dvalue;
    int ilocal;
  };
  struct Bonus *bonus;

  AtomVecBody(class LAMMPS *);
  ~AtomVecBody();
  void process_args(int, char **);

  void copy_bonus(int, int, int);
  void clear_bonus();
  int pack_comm_bonus(int, int *, double *);
  void unpack_comm_bonus(int, int, double *);
  int pack_reverse_bonus(int, int, double *);
  int pack_border_bonus(int, int *, double *);
  int unpack_border_bonus(int, int, double *);
  int pack_exchange_bonus(int, double *);
  int unpack_exchange_bonus(int, double *);
  int size_restart_bonus();
  int pack_restart_bonus(int, double *);
  int unpack_restart_bonus(int, double *);
  void data_body(int, int, int, int *, double *);
  bigint memory_usage_bonus();

  void create_atom_post(int);
  void data_atom_post(int);
  void pack_data_pre(int);
  void pack_data_post(int);

  // methods used by other classes to query/set body info

  double radius_body(int, int, int *, double *);
  void set_quat(int, double *);

  int nlocal_bonus;

 private:
  int nghost_bonus,nmax_bonus;
  int intdoubleratio;            // sizeof(double) / sizeof(int)
  int body_flag;

  MyPoolChunk<int> *icp;
  MyPoolChunk<double> *dcp;

  void grow_bonus();
  void copy_bonus_all(int, int);
  // check(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Internal error in atom_style body

This error should not occur.  Contact the developers.

E: Invalid atom_style body command

No body style argument was provided.

E: Unrecognized body style

The choice of body style is unknown.

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

E: Invalid density in Atoms section of data file

Density value cannot be <= 0.0.

E: Assigning body parameters to non-body atom

Self-explanatory.

E: Assigning quat to non-body atom

Self-explanatory.

*/
