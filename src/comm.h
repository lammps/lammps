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
  int style;     // comm pattern: 0 = 6-way stencil, 1 = irregular tiling
  int layout;    // LAYOUT_UNIFORM = logical equal-sized bricks
                 // LAYOUT_NONUNIFORM = logical bricks, 
                 //                     but different sizes due to LB
                 // LAYOUT_TILED = general tiling, due to RCB LB

  int me,nprocs;                    // proc info
  int ghost_velocity;               // 1 if ghost atoms have velocity, 0 if not
  double cutghost[3];               // cutoffs used for acquiring ghost atoms
  double cutghostuser;              // user-specified ghost cutoff
  int recv_from_partition;          // recv proc layout from this partition
  int send_to_partition;            // send my proc layout to this partition
                                    // -1 if no recv or send
  int other_partition_style;        // 0 = recv layout dims must be multiple of
                                    //     my layout dims
  int maxexchange_atom;             // max contribution to exchange from AtomVec
  int maxexchange_fix;              // max contribution to exchange from Fixes
  int nthreads;                     // OpenMP threads per MPI process

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

  // methods

  Comm(class LAMMPS *);
  virtual ~Comm();
  void copy_arrays(class Comm *);
  virtual void init();
  void modify_params(int, char **);

  void set_processors(int, char **);      // set 3d processor grid attributes
  virtual void set_proc_grid(int outflag = 1); // setup 3d grid of procs

  virtual void setup() = 0;                      // setup 3d comm pattern
  virtual void forward_comm(int dummy = 0) = 0;  // forward comm of atom coords
  virtual void reverse_comm() = 0;               // reverse comm of forces
  virtual void exchange() = 0;                   // move atoms to new procs
  virtual void borders() = 0;                    // setup list of atoms to comm

  // forward/reverse comm from a Pair, Fix, Compute, Dump

  virtual void forward_comm_pair(class Pair *) = 0;
  virtual void reverse_comm_pair(class Pair *) = 0;
  virtual void forward_comm_fix(class Fix *) = 0;
  virtual void reverse_comm_fix(class Fix *) = 0;
  virtual void forward_comm_variable_fix(class Fix *) = 0;
  virtual void reverse_comm_variable_fix(class Fix *) = 0;
  virtual void forward_comm_compute(class Compute *) = 0;
  virtual void reverse_comm_compute(class Compute *) = 0;
  virtual void forward_comm_dump(class Dump *) = 0;
  virtual void reverse_comm_dump(class Dump *) = 0;

  // forward comm of an array
  // exchange of info on neigh stencil
  // set processor mapping options

  virtual void forward_comm_array(int, double **) = 0;  
  virtual int exchange_variable(int, double *, double *&) = 0;
  virtual bigint memory_usage() = 0;

  // non-virtual functions common to all Comm styles

  void ring(int, int, void *, int, void (*)(int, char *),
            void *, int self = 1);
  int read_lines_from_file(FILE *, int, int, char *);  
  int read_lines_from_file_universe(FILE *, int, int, char *);  

 protected:
  int mode;                  // 0 = single cutoff, 1 = multi-type cutoff
  int bordergroup;           // only communicate this group in borders

  int triclinic;                    // 0 if domain is orthog, 1 if triclinic
  int map_style;                    // non-0 if global->local mapping is done
  int comm_x_only,comm_f_only;      // 1 if only exchange x,f in for/rev comm

  int size_forward;                 // # of per-atom datums in forward comm
  int size_reverse;                 // # of datums in reverse comm
  int size_border;                  // # of datums in forward border comm

  int maxforward,maxreverse;        // max # of datums in forward/reverse comm
  int maxexchange;                  // max # of datums/atom in exchange comm 

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
};

}

#endif

/* ERROR/WARNING messages:

*/
