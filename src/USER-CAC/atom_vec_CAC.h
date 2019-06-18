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

AtomStyle(cac,AtomVecCAC)

#else

#ifndef LMP_ATOM_VEC_CAC_H
#define LMP_ATOM_VEC_CAC_H

#include "atom_vec.h"
#include "asa_user.h"

namespace LAMMPS_NS {

class AtomVecCAC : public AtomVec {
 public:
  //minimization algorithm structs for neighbor rebuild checks
  asacg_parm *cgParm;
  asa_parm *asaParm;
  asa_objective *Objective;


  AtomVecCAC(class LAMMPS *);
  virtual ~AtomVecCAC();
  virtual void grow(int);
  virtual void grow_reset();
  virtual void copy(int, int, int);
  virtual void process_args(int, char **);
  virtual int pack_comm(int, int *, double *, int, int *);
  virtual int pack_comm_vel(int, int *, double *, int, int *);
  virtual void unpack_comm(int, int, double *);
  virtual void unpack_comm_vel(int, int, double *);
  virtual int pack_reverse(int, int, double *);
  virtual void unpack_reverse(int, int *, double *);
  virtual int pack_border(int, int *, double *, int, int *);
  virtual int pack_border_vel(int, int *, double *, int, int *);
  virtual void unpack_border(int, int, double *);
  virtual void unpack_border_vel(int, int, double *);
  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(double *);
  virtual int size_restart();
  virtual int pack_restart(int, double *);
  virtual int unpack_restart(double *);
  virtual void create_atom(int, double *);
  virtual void data_atom(double *, imageint, char **);
  virtual void pack_data(double **);
  virtual void write_data(FILE *, int, double **);
  virtual bigint memory_usage();
  virtual void force_clear(int, size_t);
  virtual int check_distance_function(double deltasq); //specific neighbor rebuild check function 
  virtual void set_hold_properties(); //sets nodal positions at reneighboring step for comparison

  //CAC min objective functions for neighbor rebuild
  virtual double myvalue(asa_objective *asa);
  virtual void mygrad(asa_objective *asa);
  virtual double shape_function(double, double, double,int,int);
  virtual double shape_function_derivative(double, double, double,int,int,int);
  virtual void init();

 protected:
  tagint *tag;
  int *type,*mask;
  imageint *image;
  int nodes_per_element, maxpoly; //maximum number of nodes and atoms per unit cell per element in model
  double **x,**v,**f;
  double ****nodal_positions,****nodal_velocities,****nodal_forces,
    ****nodal_gradients, ****initial_nodal_positions, *scale_search_range;
  int *poly_count, **node_types, *element_type,
	  **element_scale, scale_count, oneflag, *scale_list;
  int element_type_count;
  int search_range_max;
  int initial_size;
  int poly_min;
  int min_nodes_per_element;
  int min_element_index;
  double deltasq_trigger;
  double ****hold_nodal_positions;
  int max_old;
  int *node_count_per_poly;
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
