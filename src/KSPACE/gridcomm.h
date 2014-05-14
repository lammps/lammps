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

#ifndef LMP_GRIDCOMM_H
#define LMP_GRIDCOMM_H

#include "pointers.h"

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_FLOAT
#else
typedef double FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

namespace LAMMPS_NS {

class GridComm : protected Pointers {
 public:
  GridComm(class LAMMPS *, MPI_Comm, int, int,
           int, int, int, int, int, int,
           int, int, int, int, int, int,
           int, int, int, int, int, int);
  GridComm(class LAMMPS *, MPI_Comm, int, int,
           int, int, int, int, int, int,
           int, int, int, int, int, int,
           int, int, int, int, int, int,
           int, int, int, int, int, int);
  ~GridComm();
  void ghost_notify();
  int ghost_overlap();
  void setup();
  void forward_comm(class KSpace *, int);
  void reverse_comm(class KSpace *, int);
  double memory_usage();

 private:
  int me;
  int nforward,nreverse;
  MPI_Comm gridcomm;
  MPI_Request request;
  MPI_Status status;

  // in = inclusive indices of 3d grid chunk that I own
  // out = inclusive indices of 3d grid chunk I own plus ghosts I use
  // proc = 6 neighbor procs that surround me
  // ghost = # of my owned grid planes needed from me
  //         by each of 6 neighbor procs to become their ghost planes

  int inxlo,inxhi,inylo,inyhi,inzlo,inzhi;
  int outxlo,outxhi,outylo,outyhi,outzlo,outzhi;
  int outxlo_max,outxhi_max,outylo_max,outyhi_max,outzlo_max,outzhi_max;
  int procxlo,procxhi,procylo,procyhi,proczlo,proczhi;
  int ghostxlo,ghostxhi,ghostylo,ghostyhi,ghostzlo,ghostzhi;

  int nbuf;
  FFT_SCALAR *buf1,*buf2;

  struct Swap {
    int sendproc;       // proc to send to for forward comm
    int recvproc;       // proc to recv from for forward comm
    int npack;          // # of datums to pack
    int nunpack;        // # of datums to unpack
    int *packlist;      // 3d array offsets to pack
    int *unpacklist;    // 3d array offsets to unpack
  };

  int nswap;
  Swap *swap;

  int indices(int *&, int, int, int, int, int, int);
};

}

#endif
