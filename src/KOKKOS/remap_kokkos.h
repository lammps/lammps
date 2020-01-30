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

#ifndef LMP_REMAP_KOKKOS_H
#define LMP_REMAP_KOKKOS_H

#include "pointers.h"
#include <mpi.h>
#include "fftdata_kokkos.h"
#include "remap.h"

namespace LAMMPS_NS {

// details of how to do a 3d remap

template<class DeviceType>
struct remap_plan_3d_kokkos {
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> AT;
  typename AT::t_FFT_SCALAR_1d d_sendbuf;                  // buffer for MPI sends
  typename AT::t_FFT_SCALAR_1d d_scratch;                  // scratch buffer for MPI recvs
  void (*pack)(typename AT::t_FFT_SCALAR_1d_um, int, typename AT::t_FFT_SCALAR_1d_um, int, struct pack_plan_3d *);
                                    // which pack function to use
  void (*unpack)(typename AT::t_FFT_SCALAR_1d_um, int, typename AT::t_FFT_SCALAR_1d_um, int, struct pack_plan_3d *);
                                    // which unpack function to use
  int *send_offset;                 // extraction loc for each send
  int *send_size;                   // size of each send message
  int *send_proc;                   // proc to send each message to
  struct pack_plan_3d *packplan;    // pack plan for each send message
  int *recv_offset;                 // insertion loc for each recv
  int *recv_size;                   // size of each recv message
  int *recv_proc;                   // proc to recv each message from
  int *recv_bufloc;                 // offset in scratch buf for each recv
  MPI_Request *request;             // MPI request for each posted recv
  struct pack_plan_3d *unpackplan;  // unpack plan for each recv message
  int nrecv;                        // # of recvs from other procs
  int nsend;                        // # of sends to other procs
  int self;                         // whether I send/recv with myself
  int memory;                       // user provides scratch space or not
  MPI_Comm comm;                    // group of procs performing remap
  int usecollective;                // use collective or point-to-point MPI
  int commringlen;                  // length of commringlist
  int *commringlist;                // ranks on communication ring of this plan
};

template<class DeviceType>
class RemapKokkos : protected Pointers {
 public:
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> AT;
  RemapKokkos(class LAMMPS *);
  RemapKokkos(class LAMMPS *, MPI_Comm,int,int,int,int,int,int,
        int,int,int,int,int,int,int,int,int,int,int);
  ~RemapKokkos();
  void perform(typename AT::t_FFT_SCALAR_1d, typename AT::t_FFT_SCALAR_1d, typename AT::t_FFT_SCALAR_1d);

  struct remap_plan_3d_kokkos<DeviceType> *plan;

  void remap_3d_kokkos(typename AT::t_FFT_SCALAR_1d, typename AT::t_FFT_SCALAR_1d, typename AT::t_FFT_SCALAR_1d, struct remap_plan_3d_kokkos<DeviceType> *);
  struct remap_plan_3d_kokkos<DeviceType> *remap_3d_create_plan_kokkos(MPI_Comm,
                                             int, int, int, int, int, int,
                                             int, int, int, int, int, int,
                                             int, int, int, int, int);
  void remap_3d_destroy_plan_kokkos(struct remap_plan_3d_kokkos<DeviceType> *);
};

}

#endif

/* ERROR/WARNING messages:

E: Could not create 3d remap plan

The FFT setup in pppm failed.

*/
