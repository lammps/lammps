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

#ifndef LMP_COMM_H
#define LMP_COMM_H

#include "pointers.h"

namespace LAMMPS_NS {

class Comm : protected Pointers {
 public:
  int me,nprocs;                    // proc info
  int procgrid[3];                  // assigned # of procs in each dim
  int user_procgrid[3];             // user request for procs in each dim
  int myloc[3];                     // which proc I am in each dim
  int procneigh[3][2];              // my 6 neighboring procs
  int ghost_velocity;               // 1 if ghost atoms have velocity, 0 if not
  int maxforward_fix;               // comm sizes called from Fix,Pair
  int maxreverse_fix;
  int maxforward_pair;
  int maxreverse_pair;
  double cutghost[3];               // cutoffs used for acquiring ghost atoms
  double cutghostuser;              // user-specified ghost cutoff
  int ***grid2proc;                 // which proc owns i,j,k loc in 3d grid

  Comm(class LAMMPS *);
  virtual ~Comm();

  virtual void init();
  virtual void set_procs();                 // setup 3d grid of procs
  virtual void setup();                     // setup 3d communication pattern
  virtual void forward_comm(int dummy = 0); // forward communication of atom coords
  virtual void reverse_comm();              // reverse communication of forces
  virtual void exchange();                  // move atoms to new procs
  virtual void borders();                   // setup list of atoms to communicate

  virtual void forward_comm_pair(class Pair *);        // forward comm from a Pair
  virtual void reverse_comm_pair(class Pair *);        // reverse comm from a Pair
  virtual void forward_comm_fix(class Fix *);          // forward comm from a Fix
  virtual void reverse_comm_fix(class Fix *);          // reverse comm from a Fix
  virtual void forward_comm_compute(class Compute *);  // forward comm from a Compute
  virtual void reverse_comm_compute(class Compute *);  // reverse comm from a Compute

  virtual void set(int, char **);           // set communication style
  virtual bigint memory_usage();

 protected:
  int style;                        // single vs multi-type comm
  int nswap;                        // # of swaps to perform
  int need[3];                      // procs I need atoms from in each dim
  int triclinic;                    // 0 if domain is orthog, 1 if triclinic
  int maxswap;                      // max # of swaps memory is allocated for
  int size_forward;                 // # of per-atom datums in forward comm
  int size_reverse;                 // # of datums in reverse comm
  int size_border;                  // # of datums in forward border comm
  int *sendnum,*recvnum;            // # of atoms to send/recv in each swap
  int *sendproc,*recvproc;          // proc to send/recv to/from at each swap
  int *size_forward_recv;           // # of values to recv in each forward comm
  int *size_reverse_send;           // # to send in each reverse comm
  int *size_reverse_recv;           // # to recv in each reverse comm
  double *slablo,*slabhi;           // bounds of slab to send at each swap
  double **multilo,**multihi;       // bounds of slabs for multi-type swap
  double **cutghostmulti;           // cutghost on a per-type basis
  int *pbc_flag;                    // general flag for sending atoms thru PBC
  int **pbc;                        // dimension flags for PBC adjustments
  int comm_x_only,comm_f_only;      // 1 if only exchange x,f in for/rev comm
  int map_style;                    // non-0 if global->local mapping is done
  int bordergroup;                  // only communicate this group in borders

  int *firstrecv;                   // where to put 1st recv atom in each swap
  int **sendlist;                   // list of atoms to send in each swap
  int *maxsendlist;                 // max size of send list for each swap

  double *buf_send;                 // send buffer for all comm
  double *buf_recv;                 // recv buffer for all comm
  int maxsend,maxrecv;              // current size of send/recv buffer
  int maxforward,maxreverse;        // max # of datums in forward/reverse comm

  virtual void procs2box();                 // map procs to 3d box
  virtual void cross(double, double, double,
	     double, double, double,
	     double &, double &, double &);    // cross product
  virtual void grow_send(int,int);          // reallocate send buffer
  virtual void grow_recv(int);              // free/allocate recv buffer
  virtual void grow_list(int, int);         // reallocate one sendlist
  virtual void grow_swap(int);              // grow swap and multi arrays
  virtual void allocate_swap(int);          // allocate swap arrays
  virtual void allocate_multi(int);         // allocate multi arrays
  virtual void free_swap();                 // free swap arrays
  virtual void free_multi();                // free multi arrays
};

}

#endif
