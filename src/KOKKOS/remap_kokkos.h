// clang-format off
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
  remap_plan_3d_kokkos() :
    pack(nullptr), unpack(nullptr), send_offset(nullptr), send_size(nullptr),
    send_proc(nullptr), send_bufloc(nullptr), packplan(nullptr), recv_offset(nullptr),
    recv_size(nullptr), recv_proc(nullptr), recv_bufloc(nullptr), isend_reqs(nullptr),
    request(nullptr), unpackplan(nullptr), nrecv(0), nsend(0), self(0), memory(0),
    comm(MPI_COMM_NULL), usecollective(0), usenonblocking(0), usegpu_aware(0),
    commringlen(0), commringlist(nullptr), sendcnts(nullptr), rcvcnts(nullptr),
    sdispls(nullptr), rdispls(nullptr), selfcommringloc(-1), selfnsendloc(-1),
    selfnrecvloc(-1) {}

  // frees every buffer managed with malloc()/free(), so that "delete plan" is
  // complete at any point during remap_3d_create_plan_kokkos(), however far it
  // got.  The Kokkos views clean up after themselves; the communicator does
  // not belong to the plan until create_plan() succeeds and is released by
  // remap_3d_destroy_plan_kokkos().

  ~remap_plan_3d_kokkos()
  {
#define SAFE_FREE(ptr) if (ptr) free(ptr)
    SAFE_FREE(send_offset);
    SAFE_FREE(send_size);
    SAFE_FREE(send_proc);
    SAFE_FREE(send_bufloc);
    SAFE_FREE(packplan);
    SAFE_FREE(recv_offset);
    SAFE_FREE(recv_size);
    SAFE_FREE(recv_proc);
    SAFE_FREE(recv_bufloc);
    SAFE_FREE(isend_reqs);
    SAFE_FREE(request);
    SAFE_FREE(unpackplan);
    SAFE_FREE(commringlist);
    SAFE_FREE(sendcnts);
    SAFE_FREE(rcvcnts);
    SAFE_FREE(sdispls);
    SAFE_FREE(rdispls);
#undef SAFE_FREE
  }

  remap_plan_3d_kokkos(const remap_plan_3d_kokkos &) = delete;
  remap_plan_3d_kokkos(remap_plan_3d_kokkos &&) = delete;
  remap_plan_3d_kokkos &operator=(const remap_plan_3d_kokkos &) = delete;
  remap_plan_3d_kokkos &operator=(remap_plan_3d_kokkos &&) = delete;

  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef FFTArrayTypes<DeviceType> FFT_AT;
  typename FFT_AT::t_FFT_SCALAR_1d d_sendbuf;                  // buffer for MPI sends
  FFT_HAT::t_FFT_SCALAR_1d h_sendbuf;                          // host buffer for MPI sends
  typename FFT_AT::t_FFT_SCALAR_1d d_scratch;                  // scratch buffer for MPI recvs
  FFT_HAT::t_FFT_SCALAR_1d h_scratch;                          // host scratch buffer for MPI recvs
  void (*pack)(typename FFT_AT::t_FFT_SCALAR_1d_um, int, typename FFT_AT::t_FFT_SCALAR_1d_um, int, struct pack_plan_3d *);
                                    // which pack function to use
  void (*unpack)(typename FFT_AT::t_FFT_SCALAR_1d_um, int, typename FFT_AT::t_FFT_SCALAR_1d_um, int, struct pack_plan_3d *);
                                    // which unpack function to use
  int *send_offset;                 // extraction loc for each send
  int *send_size;                   // size of each send message
  int *send_proc;                   // proc to send each message to
  int *send_bufloc;                 // if usenonblocking, offset in send buf for each isend
  struct pack_plan_3d *packplan;    // pack plan for each send message
  int *recv_offset;                 // insertion loc for each recv
  int *recv_size;                   // size of each recv message
  int *recv_proc;                   // proc to recv each message from
  int *recv_bufloc;                 // offset in scratch buf for each recv
  MPI_Request *isend_reqs;          // MPI request for each posted isend
  MPI_Request *request;             // MPI request for each posted recv
  struct pack_plan_3d *unpackplan;  // unpack plan for each recv message
  int nrecv;                        // # of recvs from other procs
  int nsend;                        // # of sends to other procs
  int self;                         // whether I send/recv with myself
  int memory;                       // user provides scratch space or not
  MPI_Comm comm;                    // group of procs performing remap
  int usecollective;                // use collective or point-to-point MPI
  int usenonblocking;               // if using point-to-point MPI, use MPI_Isend
  int usegpu_aware;                 // use GPU-Aware MPI or not
  // variables for collective MPI only
  int commringlen;                  // length of commringlist
  int *commringlist;                // ranks on communication ring of this plan
  int *sendcnts;                    // # of elements in send buffer for each rank
  int *rcvcnts;                     // # of elements in recv buffer for each rank
  int *sdispls;                     // extraction location in send buffer for each rank
  int *rdispls;                     // extraction location in recv buffer for each rank
  int selfcommringloc;              // self rank's index in commringlist
  int selfnsendloc;                 // self rank's index in send lists
  int selfnrecvloc;                 // self rank's index in recv lists
};

template<class DeviceType>
class RemapKokkos : protected Pointers {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef FFTArrayTypes<DeviceType> FFT_AT;
  RemapKokkos(class LAMMPS *);
  RemapKokkos(class LAMMPS *, MPI_Comm,int,int,int,int,int,int,
        int,int,int,int,int,int,int,int,int,int,int,int,int);
  ~RemapKokkos() override;
  void perform(typename FFT_AT::t_FFT_SCALAR_1d, typename FFT_AT::t_FFT_SCALAR_1d, typename FFT_AT::t_FFT_SCALAR_1d);

  struct remap_plan_3d_kokkos<DeviceType> *plan;

  void remap_3d_kokkos(typename FFT_AT::t_FFT_SCALAR_1d, typename FFT_AT::t_FFT_SCALAR_1d, typename FFT_AT::t_FFT_SCALAR_1d, struct remap_plan_3d_kokkos<DeviceType> *);
  struct remap_plan_3d_kokkos<DeviceType> *remap_3d_create_plan_kokkos(MPI_Comm,
                                             int, int, int, int, int, int,
                                             int, int, int, int, int, int,
                                             int, int, int, int, int, int, int);
  void remap_3d_destroy_plan_kokkos(struct remap_plan_3d_kokkos<DeviceType> *);
};

}

#endif

