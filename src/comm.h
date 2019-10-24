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

#ifndef LMP_COMM_H
#define LMP_COMM_H

#include "pointers.h"  // IWYU pragma: export

namespace LAMMPS_NS {

class Comm : protected Pointers {
 public:
  int style;     // comm pattern: 0 = 6-way stencil, 1 = irregular tiling
  int layout;    // LAYOUT_UNIFORM = equal-sized bricks
                 // LAYOUT_NONUNIFORM = logical bricks, but diff sizes via LB
                 // LAYOUT_TILED = general tiling, due to RCB LB
  enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};
  int mode;      // 0 = single cutoff, 1 = multi-type cutoff
  enum{SINGLE,MULTI};
  const char *comm_style;           //name of the comm style

  int me,nprocs;                    // proc info
  int ghost_velocity;               // 1 if ghost atoms have velocity, 0 if not
  double cutghost[3];               // cutoffs used for acquiring ghost atoms
  double cutghostuser;              // user-specified ghost cutoff (mode == 0)
  double *cutusermulti;            // per type user ghost cutoff (mode == 1)
  int recv_from_partition;          // recv proc layout from this partition
  int send_to_partition;            // send my proc layout to this partition
                                    // -1 if no recv or send
  int other_partition_style;        // 0 = recv layout dims must be multiple of
                                    //     my layout dims

  int nthreads;                // OpenMP threads per MPI process

  // public settings specific to layout = UNIFORM, NONUNIFORM

  int procgrid[3];                  // procs assigned in each dim of 3d grid
  int user_procgrid[3];             // user request for procs in each dim
  int myloc[3];                     // which proc I am in each dim
  int procneigh[3][2];              // my 6 neighboring procs, 0/1 = left/right
  double *xsplit,*ysplit,*zsplit;   // fractional (0-1) sub-domain sizes
  int ***grid2proc;                 // which proc owns i,j,k loc in 3d grid

  // public settings specific to layout = TILED

  int rcbnew;                       // 1 if just reset by rebalance, else 0
  double mysplit[3][2];             // fractional (0-1) bounds of my sub-domain
  double rcbcutfrac;                // fractional RCB cut by this proc
  int rcbcutdim;                    // dimension of RCB cut

  class NBin *bin_pointer;                //used to invoke possible binning for comm substyles
  class NStencil *stencil_pointer;        //used for possible binning in comm substyle

  // methods

  Comm(class LAMMPS *);
  virtual ~Comm();
  // NOTE: copy_arrays is called from a constructor and must not be made virtual
  void copy_arrays(class Comm *);
  virtual void init();
  void modify_params(int, char **);

  void set_processors(int, char **);      // set 3d processor grid attributes
  virtual void set_proc_grid(int outflag = 1); // setup 3d grid of procs

  double get_comm_cutoff();     // determine communication cutoff

  virtual void setup() = 0;                      // setup 3d comm pattern
  virtual void forward_comm(int dummy = 0) = 0;  // forward comm of atom coords
  virtual void reverse_comm() = 0;               // reverse comm of forces
  virtual void exchange() = 0;                   // move atoms to new procs
  virtual void borders() = 0;                    // setup list of atoms to comm

  // forward/reverse comm from a Pair, Fix, Compute, Dump

  virtual void forward_comm_pair(class Pair *) = 0;
  virtual void reverse_comm_pair(class Pair *) = 0;
  virtual void forward_comm_fix(class Fix *, int size=0) = 0;
  virtual void reverse_comm_fix(class Fix *, int size=0) = 0;
  virtual void reverse_comm_fix_variable(class Fix *) = 0;
  virtual void forward_comm_compute(class Compute *) = 0;
  virtual void reverse_comm_compute(class Compute *) = 0;
  virtual void forward_comm_dump(class Dump *) = 0;
  virtual void reverse_comm_dump(class Dump *) = 0;

  // forward comm of an array
  // exchange of info on neigh stencil
  // set processor mapping options

  virtual void forward_comm_array(int, double **) = 0;
  virtual int exchange_variable(int, double *, double *&) = 0;
  int binary(double, int, double *);

  // map a point to a processor, based on current decomposition

  virtual void coord2proc_setup() {}
  virtual int coord2proc(double *, int &, int &, int &);

  // memory usage

  virtual bigint memory_usage() = 0;

  // non-virtual functions common to all Comm styles

  void ring(int, int, void *, int, void (*)(int, char *, void *),
            void *, void *, int self = 1);
  int rendezvous(int, int, char *, int, int, int *,
                 int (*)(int, char *, int &, int *&, char *&, void *),
                 int, char *&, int, void *, int statflag=0);

  int read_lines_from_file(FILE *, int, int, char *);
  int read_lines_from_file_universe(FILE *, int, int, char *);

  // extract data useful to other classes
  virtual void *extract(const char *, int &) {return NULL;}
  
  //called to increase the size of bufextra by altering maxexchange variables
  void increase_max_atom(int);
  void increase_max_fix(int);

 protected:
  int bordergroup;           // only communicate this group in borders

  int triclinic;                    // 0 if domain is orthog, 1 if triclinic
  int map_style;                    // non-0 if global->local mapping is done
  int comm_x_only,comm_f_only;      // 1 if only exchange x,f in for/rev comm

  int size_forward;                 // # of per-atom datums in forward comm
  int size_reverse;                 // # of datums in reverse comm
  int size_border;                  // # of datums in forward border comm

  int maxforward,maxreverse;    // max # of datums in forward/reverse comm
  int maxexchange;              // max size of one exchanged atom
  int maxexchange_atom;         // contribution to maxexchange from AtomVec
  int maxexchange_fix;          // static contribution to maxexchange from Fixes
  int maxexchange_fix_dynamic;  // 1 if a fix has a dynamic contribution
  int bufextra;                 // augment send buf size for an exchange atom


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

  void init_exchange();
  int rendezvous_irregular(int, char *, int, int, int *,
                           int (*)(int, char *, int &, int *&, char *&, void *),
                           int, char *&, int, void *, int);
  int rendezvous_all2all(int, char *, int, int, int *,
                         int (*)(int, char *, int &, int *&, char *&, void *),
                         int, char *&, int, void *, int);
  void rendezvous_stats(int, int, int, int, int, int, bigint);

 public:
  enum{MULTIPLE};
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid group in comm_modify command

Self-explanatory.

E: Comm_modify group != atom_modify first group

Self-explanatory.

E: Use cutoff/multi keyword to set cutoff in multi mode

Mode is multi so cutoff keyword cannot be used.

E: Invalid cutoff in comm_modify command

Specified cutoff must be >= 0.0.

E: Use cutoff keyword to set cutoff in single mode

Mode is single so cutoff/multi keyword cannot be used.

E: Cannot set cutoff/multi before simulation box is defined

Self-explanatory.

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

E: Bad grid of processors

The 3d grid of processors defined by the processors command does not
match the number of processors LAMMPS is being run on.

E: Processor count in z must be 1 for 2d simulation

Self-explanatory.

E: Cannot put data on ring from NULL pointer

W: Communication cutoff is 0.0. No ghost atoms will be generated. Atoms may get lost.

The communication cutoff defaults to the maximum of what is inferred from pair and
bond styles (will be zero, if none are defined) and what is specified via
"comm_modify cutoff" (defaults to 0.0).  If this results to 0.0, no ghost atoms will
be generated and LAMMPS may lose atoms or use incorrect periodic images of atoms in
interaction lists.  To avoid, either define pair style zero with a suitable cutoff
or use comm_modify cutoff.

UNDOCUMENTED

U: OMP_NUM_THREADS environment is not set.

This environment variable must be set appropriately to use the
USER-OMP package.

*/
