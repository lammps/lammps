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

#ifndef LMP_GRIDCOMM_H
#define LMP_GRIDCOMM_H

#include "pointers.h"

namespace LAMMPS_NS {

class GridComm : protected Pointers {
 public:
  enum { KSPACE = 0, PAIR = 1, FIX = 2 };    // calling classes

  GridComm(class LAMMPS *, MPI_Comm, int, int, int, int, int, int, int, int, int, int, int, int,
           int, int, int);
  GridComm(class LAMMPS *, MPI_Comm, int, int, int, int, int, int, int, int, int, int, int, int,
           int, int, int, int, int, int, int, int, int, int);
  ~GridComm() override;
  void setup(int &, int &);
  int ghost_adjacent();
  void forward_comm(int, void *, int, int, int, void *, void *, MPI_Datatype);
  void reverse_comm(int, void *, int, int, int, void *, void *, MPI_Datatype);
  void gather(int, void *, int, int, int, void *, MPI_Datatype);

 protected:
  int me, nprocs;
  int layout;           // REGULAR or TILED
  MPI_Comm gridcomm;    // communicator for this class
                        // usually world, but MSM calls with subset

  // inputs from caller via constructor

  int nx, ny, nz;      // size of global grid in all 3 dims
  int inxlo, inxhi;    // inclusive extent of my grid chunk
  int inylo, inyhi;    //   0 <= in <= N-1
  int inzlo, inzhi;
  int outxlo, outxhi;      // inclusive extent of my grid chunk plus
  int outylo, outyhi;      //   ghost cells in all 6 directions
  int outzlo, outzhi;      //   lo indices can be < 0, hi indices can be >= N
  int fullxlo, fullxhi;    // extent of grid chunk that caller stores
  int fullylo, fullyhi;    //   can be same as out indices or larger
  int fullzlo, fullzhi;

  // -------------------------------------------
  // internal variables for REGULAR layout
  // -------------------------------------------

  int procxlo, procxhi;    // 6 neighbor procs that adjoin me
  int procylo, procyhi;    // not used for comm_style = tiled
  int proczlo, proczhi;

  int ghostxlo, ghostxhi;    // # of my owned grid planes needed
  int ghostylo, ghostyhi;    // by neighobr procs in each dir as their ghost planes
  int ghostzlo, ghostzhi;

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

  int *overlap_procs;       // length of Nprocs in communicator
  MPI_Request *requests;    // length of max messages this proc receives

  // RCB tree of cut info
  // each proc contributes one value, except proc 0

  struct RCBinfo {
    int dim;    // 0,1,2 = which dim the cut is in
    int cut;    // grid index of lowest cell in upper half of cut
  };

  RCBinfo *rcbinfo;

  // overlap = a proc whose owned cells overlap with my extended ghost box
  // includes overlaps across periodic boundaries, can also be self

  struct Overlap {
    int proc;      // proc whose owned cells overlap my ghost cells
    int box[6];    // box that overlaps otherproc's owned cells
                   // this box is wholly contained within global grid
    int pbc[3];    // PBC offsets to convert box to a portion of my ghost box
                   // my ghost box may extend beyond global grid
  };

  int noverlap, maxoverlap;
  Overlap *overlap;

  // request = sent to each proc whose owned cells overlap my ghost cells

  struct Request {
    int sender;    // sending proc
    int index;     // index of overlap on sender
    int box[6];    // box that overlaps receiver's owned cells
                   // wholly contained within global grid
  };

  Request *srequest, *rrequest;

  // response = reply from each proc whose owned cells overlap my ghost cells

  struct Response {
    int index;     // index of my overlap for the initial request
    int box[6];    // box that overlaps responder's owned cells
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
  // internal methods
  // -------------------------------------------

  void initialize(MPI_Comm, int, int, int, int, int, int, int, int, int, int, int, int, int, int,
                  int, int, int, int, int, int, int, int, int, int, int, int, int);
  virtual void setup_regular(int &, int &);
  virtual void setup_tiled(int &, int &);
  void ghost_box_drop(int *, int *);
  void box_drop_grid(int *, int, int, int &, int *);

  int ghost_adjacent_regular();
  int ghost_adjacent_tiled();

  template <class T> void forward_comm_regular(T *, int, int, int, void *, void *, MPI_Datatype);
  template <class T> void forward_comm_tiled(T *, int, int, int, void *, void *, MPI_Datatype);
  template <class T> void reverse_comm_regular(T *, int, int, int, void *, void *, MPI_Datatype);
  template <class T> void reverse_comm_tiled(T *, int, int, int, void *, void *, MPI_Datatype);

  virtual void grow_swap();
  void grow_overlap();

  int indices(int *&, int, int, int, int, int, int);
};

}    // namespace LAMMPS_NS

#endif
