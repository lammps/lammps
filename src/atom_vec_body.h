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

AtomStyle(body,AtomVecBody)

#else

#ifndef LMP_ATOM_VEC_BODY_H
#define LMP_ATOM_VEC_BODY_H

#include "atom_vec.h"
#include "my_pool_chunk.h"

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
  void data_body(int, int, int, char **, char **);

 private:
  tagint *tag;
  int *type,*mask;
  imageint *image;
  double **x,**v,**f;
  double *rmass;
  double **angmom,**torque;
  int *body;

  int nlocal_bonus,nghost_bonus,nmax_bonus;
  int intdoubleratio;       // sizeof(double) / sizeof(int)

  MyPoolChunk<int> *icp;
  MyPoolChunk<double> *dcp;

  void grow_bonus();
  void copy_bonus(int, int);
  //void check(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Internal error in atom_style body

This error should not occur.  Contact the developers.

E: Invalid atom_style body command

No body style argument was provided.

E: Invalid body style

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

U: Invalid atom ID in Atoms section of data file

Atom IDs must be positive integers.

*/
