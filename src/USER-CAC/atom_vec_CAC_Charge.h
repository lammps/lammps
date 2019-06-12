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

AtomStyle(CAC/charge,AtomVecCAC_Charge)

#else

#ifndef LMP_ATOM_VEC_CAC_CHARGE_H
#define LMP_ATOM_VEC_CAC_CHARGE_H

#include "atom_vec_CAC.h"

namespace LAMMPS_NS {

class AtomVecCAC_Charge : public AtomVecCAC {
 public:
  AtomVecCAC_Charge(class LAMMPS *);
  virtual ~AtomVecCAC_Charge() {}
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  void process_args(int, char **);
  virtual int pack_comm(int, int *, double *, int, int *);
  virtual int pack_comm_vel(int, int *, double *, int, int *);
  virtual void unpack_comm(int, int, double *);
  virtual void unpack_comm_vel(int, int, double *);
  int pack_reverse(int, int, double *);
  void unpack_reverse(int, int *, double *);
  virtual int pack_border(int, int *, double *, int, int *);
  virtual int pack_border_vel(int, int *, double *, int, int *);
  virtual void unpack_border(int, int, double *);
  virtual void unpack_border_vel(int, int, double *);
  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  void data_atom(double *, imageint, char **);
  void pack_data(double **);
  void write_data(FILE *, int, double **);
  bigint memory_usage();
  virtual void force_clear(int, size_t);
  
 protected:
  tagint *tag;
  int *type,*mask;
  imageint *image;
   
  double **x,**v,**f;
  double ****nodal_positions,****nodal_velocities,****nodal_forces,
  ****nodal_gradients, **node_charges, ****initial_nodal_positions, *scale_search_range;
  int *poly_count, **node_types, *element_type,
	  **element_scale, scale_count, oneflag, *scale_list;
  int element_type_count;
  int search_range_max;
  int initial_size;
  char **element_names;
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

*/
