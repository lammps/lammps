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
  int procgrid[3];                  // procs assigned in each dim of 3d grid
  int user_procgrid[3];             // user request for procs in each dim
  int myloc[3];                     // which proc I am in each dim
  int procneigh[3][2];              // my 6 neighboring procs, 0/1 = left/right
  int ghost_velocity;               // 1 if ghost atoms have velocity, 0 if not
  int uniform;                      // 1 = equal subdomains, 0 = load-balanced
  double *xsplit,*ysplit,*zsplit;   // fractional (0-1) sub-domain sizes
  double cutghost[3];               // cutoffs used for acquiring ghost atoms
  double cutghostuser;              // user-specified ghost cutoff
  int ***grid2proc;                 // which proc owns i,j,k loc in 3d grid
  int recv_from_partition;          // recv proc layout from this partition
  int send_to_partition;            // send my proc layout to this partition
                                    // -1 if no recv or send
  int other_partition_style;        // 0 = recv layout dims must be multiple of
                                    //     my layout dims
  int nthreads;                     // OpenMP threads per MPI process

  Comm(class LAMMPS *);
  virtual ~Comm();

  virtual void init();
  virtual void set_proc_grid(int outflag = 1); // setup 3d grid of procs
  virtual void setup();                        // setup 3d comm pattern
  virtual void forward_comm(int dummy = 0);    // forward comm of atom coords
  virtual void reverse_comm();                 // reverse comm of forces
  virtual void exchange();                     // move atoms to new procs
  virtual void borders();                      // setup list of atoms to comm

  virtual void forward_comm_pair(class Pair *);    // forward comm from a Pair
  virtual void reverse_comm_pair(class Pair *);    // reverse comm from a Pair
  virtual void forward_comm_fix(class Fix *);      // forward comm from a Fix
  virtual void reverse_comm_fix(class Fix *);      // reverse comm from a Fix
  virtual void forward_comm_variable_fix(class Fix *); // variable-size variant
  virtual void reverse_comm_variable_fix(class Fix *); // variable-size variant
  virtual void forward_comm_compute(class Compute *);  // forward from a Compute
  virtual void reverse_comm_compute(class Compute *);  // reverse from a Compute
  virtual void forward_comm_dump(class Dump *);    // forward comm from a Dump
  virtual void reverse_comm_dump(class Dump *);    // reverse comm from a Dump
  void forward_comm_array(int, double **);         // forward comm of array

  void ring(int, int, void *, int, void (*)(int, char *),   // ring comm
            void *, int self = 1);

  virtual void set(int, char **);         // set communication style
  void set_processors(int, char **);      // set 3d processor grid attributes

  virtual bigint memory_usage();

 protected:
  int style;                        // single vs multi-type comm
  int nswap;                        // # of swaps to perform = sum of maxneed
  int recvneed[3][2];               // # of procs away I recv atoms from
  int sendneed[3][2];               // # of procs away I send atoms to
  int maxneed[3];                   // max procs away any proc needs, per dim
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
  int gridflag;                     // option for creating 3d grid
  int mapflag;                      // option for mapping procs to 3d grid
  char xyz[4];                      // xyz mapping of procs to 3d grid
  char *customfile;                 // file with custom proc map
  char *outfile;                    // proc grid/map output file

  int otherflag;                    // 1 if this partition dependent on another
  int other_style;                  // style of dependency
  int other_procgrid[3];            // proc layout of another partition
  int other_coregrid[3];            // core layout of another partition
  int ncores;                       // # of cores per node
  int coregrid[3];                  // 3d grid of cores within a node
  int user_coregrid[3];             // user request for cores in each dim

  int *firstrecv;                   // where to put 1st recv atom in each swap
  int **sendlist;                   // list of atoms to send in each swap
  int *maxsendlist;                 // max size of send list for each swap

  double *buf_send;                 // send buffer for all comm
  double *buf_recv;                 // recv buffer for all comm
  int maxsend,maxrecv;              // current size of send/recv buffer
  int maxforward,maxreverse;        // max # of datums in forward/reverse comm

  int updown(int, int, int, double, int, double *);
                                            // compare cutoff to procs
  virtual void grow_send(int, int);         // reallocate send buffer
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

/* ERROR/WARNING messages:

W: OMP_NUM_THREADS environment is not set.

This environment variable must be set appropriately to use the
USER-OMP pacakge.

E: Bad grid of processors

The 3d grid of processors defined by the processors command does not
match the number of processors LAMMPS is being run on.

E: Processor count in z must be 1 for 2d simulation

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid group in communicate command

Self-explanatory.

E: Communicate group != atom_modify first group

Self-explanatory.

E: Invalid cutoff in communicate command

Specified cutoff must be >= 0.0.

E: Specified processors != physical processors

The 3d grid of processors defined by the processors command does not
match the number of processors LAMMPS is being run on.

E: Cannot use processors part command without using partitions

See the command-line -partition switch.

E: Invalid partitions in processors part command

Valid partitions are numbered 1 to N and the sender and receiver
cannot be the same partition.

E: Sending partition in processors part command is already a sender

Cannot specify a partition to be a sender twice.

E: Receiving partition in processors part command is already a receiver

Cannot specify a partition to be a receiver twice.

E: Processors grid numa and map style are incompatible

Using numa for gstyle in the processors command requires using
cart for the map option.

E: Processors part option and grid style are incompatible

Cannot use gstyle numa or custom with the part option.

*/
