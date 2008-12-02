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

#ifndef COMM_H
#define COMM_H

#include "pointers.h"

namespace LAMMPS_NS {

class Comm : protected Pointers {
 public:
  int me,nprocs;                    // proc info
  int style;                        // single vs multi-type comm
  int procgrid[3];                  // assigned # of procs in each dim
  int user_procgrid[3];             // user request for procs in each dim
  int myloc[3];                     // which proc I am in each dim
  int procneigh[3][2];              // my 6 neighboring procs
  int nswap;                        // # of swaps to perform
  int need[3];                      // procs I need atoms from in each dim
  int maxforward_fix;               // comm sizes called from Fix,Pair
  int maxreverse_fix;
  int maxforward_pair;
  int maxreverse_pair;
  double cutghost[3];               // cutoffs used for acquiring ghost atoms
  
  Comm(class LAMMPS *);
  ~Comm();

  void init();
  void set_procs();                 // setup 3d grid of procs
  void setup();                     // setup 3d communication pattern
  void communicate();               // communication of atom coords
  void reverse_communicate();       // reverse communication of forces
  void exchange();                  // move atoms to new procs
  void borders();                   // setup list of atoms to communicate

  void comm_pair(class Pair *);                // forward comm from a Pair
  void reverse_comm_pair(class Pair *);        // reverse comm from a Pair
  void comm_fix(class Fix *);                  // forward comm from a Fix
  void reverse_comm_fix(class Fix *);          // reverse comm from a Fix
  void comm_compute(class Compute *);          // forward comm from a Compute
  void reverse_comm_compute(class Compute *);  // reverse comm from a Compute

  void irregular();                 // irregular communication across all procs

  void set(int, char **);           // set communication style
  double memory_usage();

 private:
  int triclinic;                    // 0 if domain is orthog, 1 if triclinic
  int maxswap;                      // max # of swaps memory is allocated for
  int *sendnum,*recvnum;            // # of atoms to send/recv in each swap
  int *sendproc,*recvproc;          // proc to send/recv to/from at each swap
  int *size_comm_recv;              // # of values to recv in each forward comm
  int *size_reverse_send;           // # to send in each reverse comm
  int *size_reverse_recv;           // # to recv in each reverse comm
  double *slablo,*slabhi;           // bounds of slab to send at each swap
  double **multilo,**multihi;       // bounds of slabs for multi-type swap
  double **cutghostmulti;           // cutghost on a per-type basis
  int *pbc_flag;                    // general flag for sending atoms thru PBC
  int **pbc;                        // dimension flags for PBC adjustments
  int comm_x_only,comm_f_only;      // 1 if only exchange x,f in for/rev comm
  int map_style;                    // non-0 if global->local mapping is done
  int ***grid2proc;                 // which proc owns i,j,k loc in 3d grid
  int bordergroup;                  // only communicate this group in borders
  double cutghostuser;              // user-specified ghost cutoff

  int *firstrecv;                   // where to put 1st recv atom in each swap
  int **sendlist;                   // list of atoms to send in each swap
  int *maxsendlist;                 // max size of send list for each swap

  double *buf_send;                 // send buffer for all comm
  double *buf_recv;                 // recv buffer for all comm
  int maxsend,maxrecv;              // current size of send/recv buffer
  int maxforward,maxreverse;        // max # of datums in forward/reverse comm

  struct Plan {                // plan for irregular communication
    int nsend;                 // # of messages to send
    int nrecv;                 // # of messages to recv
    int sendmax;               // # of doubles in largest send message
    int *proc_send;            // procs to send to
    int *length_send;          // # of doubles to send to each proc
    int *num_send;             // # of datums to send to each proc
    int *index_send;           // list of which datums to send to each proc
    int *offset_send;          // where each datum starts in send buffer
    int *proc_recv;            // procs to recv from
    int *length_recv;          // # of doubles to recv from each proc
    MPI_Request *request;      // MPI requests for posted recvs
    MPI_Status *status;        // MPI statuses for WaitAll
  };

  void procs2box();                 // map procs to 3d box
  void cross(double, double, double,
	     double, double, double,
	     double &, double &, double &);    // cross product
  void grow_send(int,int);          // reallocate send buffer
  void grow_recv(int);              // free/allocate recv buffer
  void grow_list(int, int);         // reallocate one sendlist
  void grow_swap(int);              // grow swap and multi arrays
  void allocate_swap(int);          // allocate swap arrays
  void allocate_multi(int);         // allocate multi arrays
  void free_swap();                 // free swap arrays
  void free_multi();                // free multi arrays

  struct Plan *irregular_create(int, int *, int *, int *);
  void irregular_perform(Plan *, double *, int *, double *);
  void irregular_destroy(Plan *);
  int irregular_lookup(double *);
};

}

#endif
