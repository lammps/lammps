/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_COMM_H
#define LMP_COMM_H

#include "pointers.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class Comm : protected Pointers {
 public:
  enum { BRICK, TILED };
  int style;    // BRICK = 6-way stencil communication
                // TILED = irregular tiling communication

  enum { LAYOUT_UNIFORM, LAYOUT_NONUNIFORM, LAYOUT_TILED };
  int layout;    // LAYOUT_UNIFORM = equal-sized bricks
                 // LAYOUT_NONUNIFORM = logical bricks, but diff sizes via LB
                 // LAYOUT_TILED = general tiling, due to RCB LB
  enum { SINGLE, MULTI, MULTIOLD };
  int mode;    // SINGLE = single cutoff
               // MULTI = multi-collection cutoff
               // MULTIOLD = multiold-type cutoff

  int me, nprocs;               // proc info
  int ghost_velocity;           // 1 if ghost atoms have velocity, 0 if not
  double cutghost[3];           // cutoffs used for acquiring ghost atoms
  double cutghostuser;          // user-specified ghost cutoff (mode == SINGLE)
  double *cutusermulti;         // per collection user ghost cutoff (mode == MULTI)
  double *cutusermultiold;      // per type user ghost cutoff (mode == MULTIOLD)
  int ncollections;             // # of collections known by comm, used to test if # has changed
  int ncollections_cutoff;      // # of collections stored b cutoff/multi
  int recv_from_partition;      // recv proc layout from this partition
  int send_to_partition;        // send my proc layout to this partition
                                // -1 if no recv or send
  int other_partition_style;    // 0 = recv layout dims must be multiple of
                                //     my layout dims

  int nthreads;    // OpenMP threads per MPI process

  // public settings specific to layout = UNIFORM, NONUNIFORM

  int procgrid[3];                     // proc count assigned to each dim of 3d grid
  int user_procgrid[3];                // user request for proc counts in each dim
  int myloc[3];                        // which proc I am in each dim, 0 to N-1
  int procneigh[3][2];                 // my 6 neighboring procs, 0/1 = left/right
  double *xsplit, *ysplit, *zsplit;    // fractional (0-1) sub-domain sizes, includes 0/1
  int ***grid2proc;                    // which proc owns i,j,k loc in 3d grid

  // public settings specific to layout = TILED

  int rcbnew;              // 1 if just reset by rebalance, else 0
  double mysplit[3][2];    // fractional (0-1) bounds of my sub-domain
  double rcbcutfrac;       // fractional RCB cut by this proc
  int rcbcutdim;           // dimension of RCB cut

  // methods

  Comm(class LAMMPS *);
  ~Comm() override;
  // NOTE: copy_arrays is called from a constructor and must not be made virtual
  void copy_arrays(class Comm *);
  virtual void init();
  void modify_params(int, char **);

  void set_processors(int, char **);              // set 3d processor grid attributes
  virtual void set_proc_grid(int outflag = 1);    // setup 3d grid of procs

  double get_comm_cutoff();    // determine communication cutoff

  virtual void setup() = 0;                        // setup 3d comm pattern
  virtual void forward_comm(int dummy = 0) = 0;    // forward comm of atom coords
  virtual void reverse_comm() = 0;                 // reverse comm of forces
  virtual void exchange() = 0;                     // move atoms to new procs
  virtual void borders() = 0;                      // setup list of atoms to comm

  // forward/reverse comm from a Pair, Bond, Fix, Compute, Dump

  virtual void forward_comm(class Pair *) = 0;
  virtual void reverse_comm(class Pair *) = 0;
  virtual void forward_comm(class Bond *) = 0;
  virtual void reverse_comm(class Bond *) = 0;
  virtual void forward_comm(class Fix *, int size = 0) = 0;
  virtual void reverse_comm(class Fix *, int size = 0) = 0;
  virtual void reverse_comm_variable(class Fix *) = 0;
  virtual void forward_comm(class Compute *) = 0;
  virtual void reverse_comm(class Compute *) = 0;
  virtual void forward_comm(class Dump *) = 0;
  virtual void reverse_comm(class Dump *) = 0;

  // forward comm of an array

  virtual void forward_comm_array(int, double **) = 0;

  // map a point to a processor, based on current decomposition

  virtual void coord2proc_setup() {}
  virtual int coord2proc(double *, int &, int &, int &);

  // memory usage

  virtual double memory_usage() = 0;

  // non-virtual functions common to all Comm styles

  void ring(int, int, void *, int, void (*)(int, char *, void *), void *, void *, int self = 1);
  int rendezvous(int, int, char *, int, int, int *,
                 int (*)(int, char *, int &, int *&, char *&, void *), int, char *&, int, void *,
                 int statflag = 0);

  // extract data useful to other classes

  virtual void *extract(const char *, int &) { return nullptr; }

 protected:
  int bordergroup;    // only communicate this group in borders

  int triclinic;                   // 0 if domain is orthog, 1 if triclinic
  int map_style;                   // non-0 if global->local mapping is done
  int comm_x_only, comm_f_only;    // 1 if only exchange x,f in for/rev comm

  int size_forward;    // # of per-atom datums in forward comm
  int size_reverse;    // # of datums in reverse comm
  int size_border;     // # of datums in forward border comm

  int maxforward, maxreverse;     // max # of datums in forward/reverse comm
  int maxexchange;                // max size of one exchanged atom
  int maxexchange_atom;           // contribution to maxexchange from AtomVec
  int maxexchange_fix;            // static contribution to maxexchange from Fixes
  int maxexchange_fix_dynamic;    // 1 if a fix has a dynamic contribution
  int bufextra;                   // augment send buf size for an exchange atom

  int gridflag;        // option for creating 3d grid
  int mapflag;         // option for mapping procs to 3d grid
  char xyz[4];         // xyz mapping of procs to 3d grid
  char *customfile;    // file with custom proc map
  char *outfile;       // proc grid/map output file

  int otherflag;            // 1 if this partition dependent on another
  int other_style;          // style of dependency
  int other_procgrid[3];    // proc layout of another partition
  int other_coregrid[3];    // core layout of another partition
  int ncores;               // # of cores per node
  int coregrid[3];          // 3d grid of cores within a node
  int user_coregrid[3];     // user request for cores in each dim
  int multi_reduce;         // 1 if multi cutoff is intra-collection cutoff

  void init_exchange();
  int rendezvous_irregular(int, char *, int, int, int *,
                           int (*)(int, char *, int &, int *&, char *&, void *), int, char *&, int,
                           void *, int);
  int rendezvous_all2all(int, char *, int, int, int *,
                         int (*)(int, char *, int &, int *&, char *&, void *), int, char *&, int,
                         void *, int);
  void rendezvous_stats(int, int, int, int, int, int, bigint);

 public:
  enum { MULTIPLE };
};

}    // namespace LAMMPS_NS

#endif
