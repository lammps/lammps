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

#include <mpi.h>

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_FLOAT
#else
typedef double FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

// details of how to do a 3d remap

struct remap_plan_3d {
  FFT_SCALAR *sendbuf;    // buffer for MPI sends
  FFT_SCALAR *scratch;    // scratch buffer for MPI recvs
  void (*pack)(FFT_SCALAR *, FFT_SCALAR *, struct pack_plan_3d *);
  // which pack function to use
  void (*unpack)(FFT_SCALAR *, FFT_SCALAR *, struct pack_plan_3d *);
  // which unpack function to use
  int *send_offset;                   // extraction loc for each send
  int *send_size;                     // size of each send message
  int *send_proc;                     // proc to send each message to
  struct pack_plan_3d *packplan;      // pack plan for each send message
  int *recv_offset;                   // insertion loc for each recv
  int *recv_size;                     // size of each recv message
  int *recv_proc;                     // proc to recv each message from
  int *recv_bufloc;                   // offset in scratch buf for each recv
  MPI_Request *request;               // MPI request for each posted recv
  struct pack_plan_3d *unpackplan;    // unpack plan for each recv message
  int nrecv;                          // # of recvs from other procs
  int nsend;                          // # of sends to other procs
  int self;                           // whether I send/recv with myself
  int memory;                         // user provides scratch space or not
  MPI_Comm comm;                      // group of procs performing remap
  int usecollective;                  // use collective or point-to-point MPI
  int commringlen;                    // length of commringlist
  int *commringlist;                  // ranks on communication ring of this plan
};

// collision between 2 regions

struct extent_3d {
  int ilo, ihi, isize;
  int jlo, jhi, jsize;
  int klo, khi, ksize;
};

// function prototypes

void remap_3d(FFT_SCALAR *, FFT_SCALAR *, FFT_SCALAR *, struct remap_plan_3d *);
struct remap_plan_3d *remap_3d_create_plan(MPI_Comm, int, int, int, int, int, int, int, int, int,
                                           int, int, int, int, int, int, int, int);
void remap_3d_destroy_plan(struct remap_plan_3d *);
int remap_3d_collide(struct extent_3d *, struct extent_3d *, struct extent_3d *);
