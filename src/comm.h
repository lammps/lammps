/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// interprocessor communication

#ifndef COMM_H
#define COMM_H

#include "lammps.h"

class Fix;
class Pair;

class Comm : public LAMMPS {
 public:
  int me,nprocs;                    // proc info
  int procgrid[3];                  // assigned # of procs in each dim
  int user_procgrid[3];             // user request for procs in each dim
  int myloc[3];                     // which proc I am in each dim
  int procneigh[3][2];              // my 6 neighboring procs
  int nswap;                        // # of swaps to perform
  int need[3];                      // procs I need atoms from in each dim
  int maxforward_fix,maxreverse_fix; // comm sizes called from Fix,Pair
  int maxforward_pair,maxreverse_pair;
  
  Comm();
  ~Comm();

  void init();
  void set_procs();                 // setup 3d grid of procs
  void setup();                     // setup 3d communication pattern
  void communicate();               // communication of atom coords
  void reverse_communicate();       // reverse communication of forces
  void exchange();                  // move atoms to new procs
  void borders();                   // setup list of atoms to communicate
  int memory_usage();               // tally memory usage
  void comm_fix(Fix *);             // forward comm from a Fix
  void reverse_comm_fix(Fix *);     // reverse comm from a Fix
  void comm_pair(Pair *);           // forward comm from a Pair
  void reverse_comm_pair(Pair *);   // reverse comm from a Pair

 private:
  int maxswap;                      // max # of swaps memory is allocated for
  int *sendnum,*recvnum;            // # of atoms to send/recv in each swap
  int *sendproc,*recvproc;          // proc to send/recv to/from at each swap
  int *size_comm_send;              // # of values to send in each comm
  int *size_comm_recv;              // # of values to recv in each comm
  int *size_reverse_send;           // # to send in each reverse comm
  int *size_reverse_recv;           // # to recv in each reverse comm
  double *slablo,*slabhi;           // bounds of slab to send at each swap
  int **pbc_flags;                  // flags for sending atoms thru PBC
                                    // [0] = 1 if any dim is across PBC
                                    // [123] = 1 if dim 012 is across PBC
  int direct_flag;                  // 1 if only x,f are forward/reverse comm
  int map_style;                    // non-0 if global->local mapping is done

  int *firstrecv;                   // where to put 1st recv atom in each swap
  int **sendlist;                   // list of atoms to send in each swap
  int *maxsendlist;                 // max size of send list for each swap

  double *buf_send;                 // send buffer for all comm
  double *buf_recv;                 // recv buffer for all comm
  int maxsend,maxrecv;              // current size of send/recv buffer
  int maxforward,maxreverse;        // max # of datums in forward/reverse comm

  void procs2box();                 // map procs to 3d box
  void grow_send(int,int);          // reallocate send buffer
  void grow_recv(int);              // free/allocate recv buffer
  void grow_list(int, int);         // reallocate one sendlist
  void grow_swap();                 // grow swap arrays
  void allocate_swap(int);          // allocate swap arrays
  void free_swap();                 // free swap arrays
};

#endif
