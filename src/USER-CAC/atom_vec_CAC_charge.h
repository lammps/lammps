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

AtomStyle(cac/charge,AtomVecCAC_Charge)

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

E: Invalid atom_style cac/charge command

Check the documentation for the correct arguments to the cac/charge atom style.

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: cac atom styles require a CAC comm style

Self-explanatory

E: poly count declared in data file was greater than maxpoly in input file

One of the input poly counts for an element in the data file was greater than
the maximum declared poly count that was passed as an arg to the atom style
cac invocation. Increase the argument value to be greater than any poly
count in your input file.

E: element type not yet defined, add definition in process_args function of atom_vec_CAC.cpp style

Self-explanatory. Contact author for detailed advice if encountering issues defining a new element type.

E: element type requires a greater number of nodes than the specified maximum nodes per element passed to atom style cac

One of the input element types in your data file requires more nodes than the maximum you specified as an
arg to the atom style cac command.

E: Invalid node index in CAC_Elements section of data file

A node index supplied as part of a CAC element's definition is not within the range of possible 
node indices (1-element_node_count) for this element type

E: Invalid poly index in CAC_Elements section of data file

A poly index supplied as part of a CAC element's definition is not within the range of possible 
poly indices (1-poly_count) for the poly_count declared in the element's header line

E: Invalid atom type in CAC_Elements section of data file

Atom types must range from 1 to specified # of types.

E: more than one type assigned to the same poly index in an element

The poly index represents an internal variable of the underlying crystal structure approximated 
by the given finite element. Since each such internal variable essentially represents a particle 
spanning a deforming lattice it is not correct to associate two different mass types to one poly index.

E: there are more nodes for one internal DOF than the element type admits

one of the poly indices has too many nodes associated with it in the declartion of 
element data. Make sure that the total number of line entries after your element
header is element_node_count*poly_count and that each poly index appears element_node_count
times.

E: cac atom style does not yet support writing data files

Self-explanatory

*/
