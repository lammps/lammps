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

#include <mpi.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <csignal>
#include <unistd.h>
#include "kokkos.h"
#include "lammps.h"
#include "force.h"
#include "neighbor_kokkos.h"
#include "neigh_list_kokkos.h"
#include "error.h"
#include "memory_kokkos.h"

#ifdef KOKKOS_ENABLE_CUDA

// for detecting GPU-direct support:
// the function  int have_gpu_direct()
// - returns -1 if GPU-direct support is unknown
// - returns  0 if no GPU-direct support available
// - returns  1 if GPU-direct support is available

#define GPU_DIRECT_UNKNOWN static int have_gpu_direct() {return -1;}

// OpenMPI supports detecting GPU-direct as of version 2.0.0
#if OPEN_MPI

#if (OMPI_MAJOR_VERSION >= 2)
#include <mpi-ext.h>
#if defined(MPIX_CUDA_AWARE_SUPPORT)
static int have_gpu_direct() { return MPIX_Query_cuda_support(); }
#else
GPU_DIRECT_UNKNOWN
#endif

#else // old OpenMPI
GPU_DIRECT_UNKNOWN
#endif

#else // unknown MPI library
GPU_DIRECT_UNKNOWN
#endif

#endif // KOKKOS_ENABLE_CUDA

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

KokkosLMP::KokkosLMP(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  kokkos_exists = 1;
  lmp->kokkos = this;

  delete memory;
  memory = new MemoryKokkos(lmp);
  memoryKK = (MemoryKokkos*) memory;

  auto_sync = 1;

  int me = 0;
  MPI_Comm_rank(world,&me);
  if (me == 0) error->message(FLERR,"KOKKOS mode is enabled");

  // process any command-line args that invoke Kokkos settings

  ngpu = 0;
  int device = 0;
  num_threads = 1;
  numa = 1;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"d") == 0 || strcmp(arg[iarg],"device") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid Kokkos command-line args");
      device = atoi(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"g") == 0 ||
               strcmp(arg[iarg],"gpus") == 0) {
#ifndef KOKKOS_ENABLE_CUDA
      error->all(FLERR,"GPUs are requested but Kokkos has not been compiled for CUDA");
#endif
      if (iarg+2 > narg) error->all(FLERR,"Invalid Kokkos command-line args");
      ngpu = atoi(arg[iarg+1]);

      int skip_gpu = 9999;
      if (iarg+2 < narg && isdigit(arg[iarg+2][0])) {
        skip_gpu = atoi(arg[iarg+2]);
        iarg++;
      }
      iarg += 2;

      char *str;
      if ((str = getenv("SLURM_LOCALID"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpu;
        if (device >= skip_gpu) device++;
      }
      if ((str = getenv("MV2_COMM_WORLD_LOCAL_RANK"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpu;
        if (device >= skip_gpu) device++;
      }
      if ((str = getenv("OMPI_COMM_WORLD_LOCAL_RANK"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpu;
        if (device >= skip_gpu) device++;
      }

    } else if (strcmp(arg[iarg],"t") == 0 ||
               strcmp(arg[iarg],"threads") == 0) {
      num_threads = atoi(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"n") == 0 ||
               strcmp(arg[iarg],"numa") == 0) {
      numa = atoi(arg[iarg+1]);
      iarg += 2;

    } else error->all(FLERR,"Invalid Kokkos command-line args");
  }

  // initialize Kokkos

  if (me == 0) {
    if (screen) fprintf(screen,"  will use up to %d GPU(s) per node\n",ngpu);
    if (logfile) fprintf(logfile,"  will use up to %d GPU(s) per node\n",ngpu);
  }

#ifdef KOKKOS_ENABLE_CUDA
  if (ngpu <= 0)
    error->all(FLERR,"Kokkos has been compiled for CUDA but no GPUs are requested");

  // check and warn about GPU-direct availability when using multiple MPI tasks

  int nmpi = 0;
  MPI_Comm_size(world,&nmpi);
  if ((nmpi > 1) && (me == 0)) {
    if ( 1 == have_gpu_direct() ) {
      ; // all good, nothing to warn about
    } else if (-1 == have_gpu_direct() ) {
      error->warning(FLERR,"Kokkos with CUDA assumes GPU-direct is available,"
                     " but cannot determine if this is the case\n         try"
                     " '-pk kokkos gpu/direct off' when getting segmentation faults");
    } else if ( 0 == have_gpu_direct() ) {
      error->warning(FLERR,"GPU-direct is NOT available, "
                     "using '-pk kokkos gpu/direct off' by default");
    } else {
      ; // should never get here
    }
  }
#endif

#ifndef KOKKOS_ENABLE_SERIAL
  if (num_threads == 1)
    error->warning(FLERR,"When using a single thread, the Kokkos Serial backend "
                         "(i.e. Makefile.kokkos_mpi_only) gives better performance "
                         "than the OpenMP backend");
#endif

  Kokkos::InitArguments args;
  args.num_threads = num_threads;
  args.num_numa = numa;
  args.device_id = device;

  Kokkos::initialize(args);

  // default settings for package kokkos command

  binsize = 0.0;
  gpu_direct_flag = 1;
  neigh_thread = 0;
  neigh_thread_set = 0;
  neighflag_qeq_set = 0;
  if (ngpu > 0) {
    neighflag = FULL;
    neighflag_qeq = FULL;
    newtonflag = 0;
    exchange_comm_classic = forward_comm_classic = reverse_comm_classic = 0;
    exchange_comm_on_host = forward_comm_on_host = reverse_comm_on_host = 0;
  } else {
    if (num_threads > 1) {
      neighflag = HALFTHREAD;
      neighflag_qeq = HALFTHREAD;
    } else {
      neighflag = HALF;
      neighflag_qeq = HALF;
    }
    newtonflag = 1;
    exchange_comm_classic = forward_comm_classic = reverse_comm_classic = 1;
    exchange_comm_on_host = forward_comm_on_host = reverse_comm_on_host = 0;
  }

#if KOKKOS_USE_CUDA
  // only if we can safely detect, that GPU-direct is not available, change default
  if (0 == have_gpu_direct()) gpu_direct_flag = 0;
#endif

#ifdef KILL_KOKKOS_ON_SIGSEGV
  signal(SIGSEGV, my_signal_handler);
#endif
}

/* ---------------------------------------------------------------------- */

KokkosLMP::~KokkosLMP()
{
  // finalize Kokkos

  Kokkos::finalize();
}

/* ----------------------------------------------------------------------
   invoked by package kokkos command
------------------------------------------------------------------------- */

void KokkosLMP::accelerator(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"neigh") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"full") == 0) neighflag = FULL;
      else if (strcmp(arg[iarg+1],"half") == 0) {
        if (num_threads > 1 || ngpu > 0)
          neighflag = HALFTHREAD;
        else
          neighflag = HALF;
      } else if (strcmp(arg[iarg+1],"n2") == 0) neighflag = N2;
      else error->all(FLERR,"Illegal package kokkos command");
      if (!neighflag_qeq_set) neighflag_qeq = neighflag;
      iarg += 2;
    } else if (strcmp(arg[iarg],"neigh/qeq") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"full") == 0) neighflag_qeq = FULL;
      else if (strcmp(arg[iarg+1],"half") == 0) {
        if (num_threads > 1 || ngpu > 0)
          neighflag_qeq = HALFTHREAD;
        else
          neighflag_qeq = HALF;
      } else if (strcmp(arg[iarg+1],"n2") == 0) neighflag_qeq = N2;
      else error->all(FLERR,"Illegal package kokkos command");
      neighflag_qeq_set = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"binsize") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      binsize = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"newton") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"off") == 0) newtonflag = 0;
      else if (strcmp(arg[iarg+1],"on") == 0) newtonflag = 1;
      else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"comm") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"no") == 0) {
        exchange_comm_classic = forward_comm_classic = reverse_comm_classic = 1;
        exchange_comm_on_host = forward_comm_on_host = reverse_comm_on_host = 0;
      } else if (strcmp(arg[iarg+1],"host") == 0) {
        exchange_comm_classic = forward_comm_classic = reverse_comm_classic = 0;
        exchange_comm_on_host = forward_comm_on_host = reverse_comm_on_host = 1;
      } else if (strcmp(arg[iarg+1],"device") == 0) {
        exchange_comm_classic = forward_comm_classic = reverse_comm_classic = 0;
        exchange_comm_on_host = forward_comm_on_host = reverse_comm_on_host = 0;
      } else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"comm/exchange") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"no") == 0) exchange_comm_classic = 1;
      else if (strcmp(arg[iarg+1],"host") == 0) {
        exchange_comm_classic = 0;
        exchange_comm_on_host = 1;
      } else if (strcmp(arg[iarg+1],"device") == 0) {
        exchange_comm_classic = 0;
        exchange_comm_on_host = 0;
      } else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"comm/forward") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"no") == 0) forward_comm_classic = 1;
      else if (strcmp(arg[iarg+1],"host") == 0) {
        forward_comm_classic = 0;
        forward_comm_on_host = 1;
      } else if (strcmp(arg[iarg+1],"device") == 0) {
        forward_comm_classic = 0;
        forward_comm_on_host = 0;
      } else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"comm/reverse") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"no") == 0) reverse_comm_classic = 1;
      else if (strcmp(arg[iarg+1],"host") == 0) {
        reverse_comm_classic = 0;
        reverse_comm_on_host = 1;
      } else if (strcmp(arg[iarg+1],"device") == 0) {
        reverse_comm_classic = 0;
        reverse_comm_on_host = 0;
      } else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"gpu/direct") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"off") == 0) gpu_direct_flag = 0;
      else if (strcmp(arg[iarg+1],"on") == 0) gpu_direct_flag = 1;
      else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"neigh/thread") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"off") == 0) neigh_thread = 0;
      else if (strcmp(arg[iarg+1],"on") == 0) neigh_thread = 1;
      else error->all(FLERR,"Illegal package kokkos command");
      neigh_thread_set = 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal package kokkos command");
  }

  // if "gpu/direct off" and "comm device", change to "comm host"

  if (!gpu_direct_flag) {
   if (exchange_comm_classic == 0 && exchange_comm_on_host == 0)
     exchange_comm_on_host = 1;
   if (forward_comm_classic == 0 && forward_comm_on_host == 0)
     forward_comm_on_host = 1;
   if (reverse_comm_classic == 0 && reverse_comm_on_host == 0)
     reverse_comm_on_host = 1;
  }

  // set newton flags
  // set neighbor binsize, same as neigh_modify command

  force->newton = force->newton_pair = force->newton_bond = newtonflag;

  if (neigh_thread && neighflag != FULL)
    error->all(FLERR,"Must use KOKKOS package option 'neigh full' with 'neigh/thread on'");

  neighbor->binsize_user = binsize;
  if (binsize <= 0.0) neighbor->binsizeflag = 0;
  else neighbor->binsizeflag = 1;
}

/* ----------------------------------------------------------------------
   called by Finish
------------------------------------------------------------------------- */

int KokkosLMP::neigh_count(int m)
{
  int inum;
  int nneigh = 0;

  ArrayTypes<LMPHostType>::t_int_1d h_ilist;
  ArrayTypes<LMPHostType>::t_int_1d h_numneigh;

  NeighborKokkos *nk = (NeighborKokkos *) neighbor;
  if (nk->lists[m]->execution_space == Host) {
    NeighListKokkos<LMPHostType>* nlistKK = (NeighListKokkos<LMPHostType>*) nk->lists[m];
    inum = nlistKK->inum;
#ifndef KOKKOS_USE_CUDA_UVM
    h_ilist = Kokkos::create_mirror_view(nlistKK->d_ilist);
    h_numneigh = Kokkos::create_mirror_view(nlistKK->d_numneigh);
#else
    h_ilist = nlistKK->d_ilist;
    h_numneigh = nlistKK->d_numneigh;
#endif
    Kokkos::deep_copy(h_ilist,nlistKK->d_ilist);
    Kokkos::deep_copy(h_numneigh,nlistKK->d_numneigh);
  } else if (nk->lists[m]->execution_space == Device) {
    NeighListKokkos<LMPDeviceType>* nlistKK = (NeighListKokkos<LMPDeviceType>*) nk->lists[m];
    inum = nlistKK->inum;
#ifndef KOKKOS_USE_CUDA_UVM
    h_ilist = Kokkos::create_mirror_view(nlistKK->d_ilist);
    h_numneigh = Kokkos::create_mirror_view(nlistKK->d_numneigh);
#else
    h_ilist = nlistKK->d_ilist;
    h_numneigh = nlistKK->d_numneigh;
#endif
    Kokkos::deep_copy(h_ilist,nlistKK->d_ilist);
    Kokkos::deep_copy(h_numneigh,nlistKK->d_numneigh);
  }

  for (int i = 0; i < inum; i++) nneigh += h_numneigh[h_ilist[i]];

  return nneigh;
}

void KokkosLMP::my_signal_handler(int sig)
{
  if (sig == SIGSEGV) {
    kill(getpid(),SIGABRT);
  }
}
