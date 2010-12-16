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

/* ----------------------------------------------------------------------
   Contributing authors: Mike Brown (ORNL), brownw@ornl.gov
------------------------------------------------------------------------- */

#ifndef PAIR_GPU_DEVICE_H
#define PAIR_GPU_DEVICE_H

#include "pair_gpu_atom.h"
#include "pair_gpu_nbor.h"
#include "mpi.h"
#include <sstream>
#include "stdio.h"
#include <string>

template <class numtyp, class acctyp>
class PairGPUDevice {
 public:
  PairGPUDevice();
  ~PairGPUDevice(); 
 
  /// Initialize the device for use by this process
  /** Sets up a per-device MPI communicator for load balancing and initializes
    * the device (>=first_gpu and <=last_gpu) that this proc will be using **/
  bool init_device(MPI_Comm world, MPI_Comm replica, const int first_gpu, 
                   const int last_gpu, const int gpu_mode, 
                   const double particle_split, const int nthreads);

  /// Initialize the device for Atom and Neighbor storage
  /** \param rot True if quaternions need to be stored
    * \param nlocal Total number of local particles to allocate memory for
    * \param host_nlocal Initial number of host particles to allocate memory for
    * \param nall Total number of local+ghost particles
    * \param gpu_nbor True if neighboring is performed on device
    * \param gpu_host 0 if host will not perform force calculations,
    *                 1 if gpu_nbor is true, and host needs a half nbor list,
    *                 2 if gpu_nbor is true, and host needs a full nbor list
    * \param max_nbors Initial number of rows in the neighbor matrix
    * \param cell_size cutoff+skin 
    * \param pre_cut True if cutoff test will be performed in separate kernel
    *                than the force kernel **/
  bool init(const bool charge, const bool rot, const int nlocal,
            const int host_nlocal, const int nall, const int maxspecial, 
            const bool gpu_nbor, const int gpu_host, const int max_nbors,
            const double cell_size, const bool pre_cut);

  /// Output a message for pair_style acceleration with device stats
  void init_message(FILE *screen, const char *name,
                    const int first_gpu, const int last_gpu);

  /// Output a message with timing information
  void output_times(UCL_Timer &time_pair, const double avg_split, 
                    const double max_bytes, FILE *screen);

  /// Clear all memory on host and device associated with atom and nbor data
  void clear();
  
  /// Clear all memory on host and device
  void clear_device();

  /// Start timer on host
  inline void start_host_timer() { _cpu_full=MPI_Wtime(); }
  
  /// Stop timer on host
  inline void stop_host_timer() { _cpu_full=MPI_Wtime()-_cpu_full; }
  
  /// Return host time
  inline double host_time() { return _cpu_full; }

  /// Return host memory usage in bytes
  double host_memory_usage() const;

  /// Return the number of procs sharing a device (size of device commincator)
  inline int procs_per_gpu() const { return _procs_per_gpu; }
  /// Return the number of threads per proc
  inline int num_threads() const { return _nthreads; }
  /// My rank within all processes
  inline int world_me() const { return _world_me; }
  /// Total number of processes
  inline int world_size() const { return _world_size; }
  /// MPI Barrier for world
  inline void world_barrier() { MPI_Barrier(_comm_world); }
  /// Return the replica MPI communicator
  inline MPI_Comm & replica() { return _comm_replica; }
  /// My rank within replica communicator
  inline int replica_me() const { return _replica_me; }
  /// Number of procs in replica communicator
  inline int replica_size() const { return _replica_size; }
  /// Return the per-GPU MPI communicator
  inline MPI_Comm & gpu_comm() { return _comm_gpu; }
  /// Return my rank in the device communicator
  inline int gpu_rank() const { return _gpu_rank; }
  /// MPI Barrier for gpu
  inline void gpu_barrier() { MPI_Barrier(_comm_gpu); }
  /// Return the 'mode' for acceleration: GPU_FORCE or GPU_NEIGH
  inline int gpu_mode() const { return _gpu_mode; }
  /// Index of first device used by a node
  inline int first_device() const { return _first_device; }
  /// Index of last device used by a node
  inline int last_device() const { return _last_device; }
  /// Particle split defined in fix
  inline double particle_split() const { return _particle_split; }
  /// Return the initialization count for the device
  inline int init_count() const { return _init_count; }

  // -------------------------- DEVICE DATA ------------------------- 

  /// Geryon Device
  UCL_Device *gpu;

  enum{GPU_FORCE, GPU_NEIGH};

  // --------------------------- ATOM DATA -------------------------- 

  /// Atom Data
  PairGPUAtom<numtyp,acctyp> atom;

  // --------------------------- NBOR DATA ----------------------------
  
  /// Neighbor Data
  PairGPUNbor nbor;

 private:
  int _init_count;
  bool _device_init;
  MPI_Comm _comm_world, _comm_replica, _comm_gpu;
  int _procs_per_gpu, _gpu_rank, _world_me, _world_size, _replica_me, 
      _replica_size;
  int _gpu_mode, _first_device, _last_device, _nthreads;
  double _particle_split;
  double _cpu_full;

  template <class t>
  inline std::string toa(const t& in) {
    std::ostringstream o;
    o.precision(2);
    o << in;
    return o.str();
  }

};

#endif
