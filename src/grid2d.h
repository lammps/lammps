/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_GRID2D_H
#define LMP_GRID2D_H

#include "pointers.h"

namespace LAMMPS_NS {

class Grid2d : protected Pointers {
 public:
  enum { KSPACE = 0, PAIR = 1, FIX = 2 };    // calling classes

  Grid2d(class LAMMPS *, MPI_Comm, int, int, double, int, double,
         int &, int &, int &, int &,
         int &, int &, int &, int &);
  Grid2d(class LAMMPS *, MPI_Comm, int, int, int, int, int, int, int, int, int, int);
  Grid2d(class LAMMPS *, MPI_Comm, int, int, int, int, int, int, int, int, int, int,
         int, int, int, int, int);
  ~Grid2d() override;

  int identical(Grid2d *);
  void get_size(int &, int &);
  void get_bounds(int &, int &, int &, int &);
  void get_bounds_ghost(int &, int &, int &, int &);
  
  void setup(int &, int &);
  int ghost_adjacent();
  void forward_comm(int, void *, int, int, int, void *, void *, MPI_Datatype);
  void reverse_comm(int, void *, int, int, int, void *, void *, MPI_Datatype);

  void setup_remap(Grid2d *, int &, int &);
  void remap(int, void *, int, int, void *, void *, MPI_Datatype);

  void read_file(int, void *, FILE *, int, int);
  void write_file(int, void *, int, int, int, MPI_Datatype);

protected:
  int me, nprocs;
  int layout;           // REGULAR or TILED
  MPI_Comm gridcomm;    // communicator for this class
                        // usually world, but MSM calls with subset

  int ngrid[2];         // global grid size

  // inputs from caller via constructor

  int nx, ny;          // size of global grid in both dims
  int inxlo, inxhi;    // inclusive extent of my grid chunk
  int inylo, inyhi;    //   0 <= in <= N-1
  int outxlo, outxhi;      // inclusive extent of my grid chunk plus
  int outylo, outyhi;      //   ghost cells in all 4 directions
                           //   lo indices can be < 0, hi indices can be >= N
  int fullxlo, fullxhi;    // extent of grid chunk that caller stores
  int fullylo, fullyhi;    //   can be same as out indices or larger

  // -------------------------------------------
  // internal variables for BRICK layout
  // -------------------------------------------

  int procxlo, procxhi;    // 4 neighbor procs that adjoin me
  int procylo, procyhi;    // not used for comm_style = tiled

  int ghostxlo, ghostxhi;    // # of my owned grid planes needed
  int ghostylo, ghostyhi;    // by neighobr procs in each dir as their ghost planes

  // swap = exchange of owned and ghost grid cells between 2 procs, including self

  struct Swap {
    int sendproc;       // proc to send to for forward comm
    int recvproc;       // proc to recv from for forward comm
    int npack;          // # of datums to pack
    int nunpack;        // # of datums to unpack
    int *packlist;      // 3d array offsets to pack
    int *unpacklist;    // 3d array offsets to unpack
  };

  int nswap, maxswap;
  Swap *swap;

  // -------------------------------------------
  // internal variables for TILED layout
  // -------------------------------------------

  MPI_Request *requests;    // length of max messages this proc receives

  // request = sent to each proc whose owned cells overlap my ghost cells

  struct Request {
    int sender;    // sending proc
    int index;     // index of overlap on sender
    int box[4];    // box that overlaps receiver's owned cells
                   // wholly contained within global grid
  };

  Request *srequest, *rrequest;

  // response = reply from each proc whose owned cells overlap my ghost cells

  struct Response {
    int index;     // index of my overlap for the initial request
    int box[4];    // box that overlaps responder's owned cells
                   // wholly contained within global grid
                   // has to unwrapped by PBC to map to my ghost cells
  };

  Response *sresponse, *rresponse;

  // send = proc to send a subset of my owned cells to, for forward comm
  // for reverse comm, proc I receive ghost overlaps with my owned cells from
  // offset used in reverse comm to recv a message in middle of a large buffer

  struct Send {
    int proc;
    int npack;
    int *packlist;
    int offset;
  };

  // recv = proc to recv a subset of my ghost cells from, for forward comm
  // for reverse comm, proc I send a subset of my ghost cells to
  // offset used in forward comm to recv a message in middle of a large buffer

  struct Recv {
    int proc;
    int nunpack;
    int *unpacklist;
    int offset;
  };

  int adjacent;    // 0 on a proc who receives ghosts from a non-neighbor proc

  // copy = subset of my owned cells to copy into subset of my ghost cells
  // that describes forward comm, for reverse comm it is the opposite

  struct Copy {
    int npack;
    int nunpack;
    int *packlist;
    int *unpacklist;
  };

  int nsend, nrecv, ncopy;
  Send *send;
  Recv *recv;
  Copy *copy;
  
  // -------------------------------------------
  // internal variables for REMAP operation
  // -------------------------------------------

  MPI_Request *requests_remap;    // length of max messages this proc receives

  int nsend_remap, nrecv_remap, self_remap;
  Send *send_remap;
  Recv *recv_remap;
  Copy copy_remap;

  // -------------------------------------------
  // internal variables for OVERLAP operation
  // -------------------------------------------
  
  int *overlap_procs;       // length of Nprocs in communicator

  // BRICK decomposition

  double *xsplit,*ysplit,*zsplit;
  int ***grid2proc;

  // TILED decomposition

  // RCB tree of cut info
  // each proc contributes one value, except proc 0

  struct RCBinfo {
    int dim;    // 0,1 = which dim the cut is in
    int cut;    // grid index of lowest cell in upper half of cut
  };

  RCBinfo *rcbinfo;

  // overlap = a proc whose owned cells overlap with my owned or ghost box
  // includes overlaps across periodic boundaries, can also be self

  struct Overlap {
    int proc;      // proc whose cells overlap my cells
    int box[4];    // box of my cells which overlap proc's cells
                   // this box is wholly contained within global grid
    int pbc[2];    // PBC offsets to convert my box to a portion of my ghost box
                   // my ghost box may extend beyond global grid
  };

  int noverlap_list, maxoverlap_list;
  Overlap *overlap_list;

  // -------------------------------------------
  // internal methods
  // -------------------------------------------

  void store(int, int, int, int, int, int, int, int,
             int, int, int, int, int, int, int, int);
  
  virtual void setup_brick(int &, int &);
  virtual void setup_tiled(int &, int &);
  int ghost_adjacent_brick();
  int ghost_adjacent_tiled();
  
  template <class T> void forward_comm_brick(T *, int, int, int, void *, void *, MPI_Datatype);
  template <class T> void forward_comm_tiled(T *, int, int, int, void *, void *, MPI_Datatype);
  template <class T> void reverse_comm_brick(T *, int, int, int, void *, void *, MPI_Datatype);
  template <class T> void reverse_comm_tiled(T *, int, int, int, void *, void *, MPI_Datatype);

  template <class T> void remap_style(T *, int, int, void *, void *, MPI_Datatype);

  template <class T> void read_file_style(T *, FILE *, int, int);
  template <class T> void write_file_style(T *, int, int, int, MPI_Datatype);

  int compute_overlap(int *, int *, Overlap *&);
  void clean_overlap();
  void box_drop(int *, int *);
  void box_drop_grid(int *, int, int, int &, int *);

  virtual void grow_swap();
  void grow_overlap();
  void deallocate_remap();

  int indices(int *&, int, int, int, int);
  int find_proc_index(int, int, int, double *);
  void find_proc_box(int, int, int *);
};

}    // namespace LAMMPS_NS

#endif
