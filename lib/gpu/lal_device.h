/***************************************************************************
                                  device.h
                             -------------------
                            W. Michael Brown (ORNL)

  Class for management of the device where the computations are performed

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov
 ***************************************************************************/

#ifndef LAL_DEVICE_H
#define LAL_DEVICE_H

#include "lal_atom.h"
#include "lal_answer.h"
#include "lal_neighbor.h"
#include "lal_pppm.h"
#include "mpi.h"
#include <sstream>
#include "stdio.h"
#include <string>
#include <queue>

namespace LAMMPS_AL {

template <class numtyp, class acctyp, 
          class grdtyp, class grdtyp4> class PPPM;

template <class numtyp, class acctyp>
class Device {
 public:
  Device();
  ~Device(); 
 
  /// Initialize the device for use by this process
  /** Sets up a per-device MPI communicator for load balancing and initializes
    * the device (>=first_gpu and <=last_gpu) that this proc will be using 
    * Returns:
    * -  0 if successfull
    * - -2 if GPU not found
    * - -4 if GPU library not compiled for GPU
    * - -6 if GPU could not be initialized for use
    * - -7 if accelerator sharing is not currently allowed on system 
    * - -11 if vendor_string has the wrong number of parameters **/
  int init_device(MPI_Comm world, MPI_Comm replica, const int first_gpu, 
                  const int last_gpu, const int gpu_mode, 
                  const double particle_split, const int nthreads,
                  const int t_per_atom, const double cell_size, 
                  char *vendor_string);

  /// Initialize the device for Atom and Neighbor storage
  /** \param rot True if quaternions need to be stored
    * \param nlocal Total number of local particles to allocate memory for
    * \param host_nlocal Initial number of host particles to allocate memory for
    * \param nall Total number of local+ghost particles
    * \param gpu_host 0 if host will not perform force calculations,
    *                 1 if gpu_nbor is true, and host needs a half nbor list,
    *                 2 if gpu_nbor is true, and host needs a full nbor list
    * \param max_nbors Initial number of rows in the neighbor matrix
    * \param cell_size cutoff+skin 
    * \param pre_cut True if cutoff test will be performed in separate kernel
    *                than the force kernel 
    * \param threads_per_atom value to be used by the neighbor list only
    *
    * Returns:
    * -  0 if successfull
    * - -1 if fix gpu not found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card **/
  int init(Answer<numtyp,acctyp> &a, const bool charge, const bool rot,
           const int nlocal, const int host_nlocal, const int nall,
           Neighbor *nbor, const int maxspecial, const int gpu_host,
           const int max_nbors, const double cell_size, const bool pre_cut,
           const int threads_per_atom, const bool vel=false);

  /// Initialize the device for Atom storage only
  /** \param nlocal Total number of local particles to allocate memory for
    * \param nall Total number of local+ghost particles
    *
    * Returns:
    * -  0 if successfull
    * - -1 if fix gpu not found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card **/
  int init(Answer<numtyp,acctyp> &ans, const int nlocal, const int nall);

  /// Output a message for pair_style acceleration with device stats
  void init_message(FILE *screen, const char *name,
                    const int first_gpu, const int last_gpu);

  /// Perform charge assignment asynchronously for PPPM
	void set_single_precompute(PPPM<numtyp,acctyp,
	                                         float,_lgpu_float4> *pppm);

  /// Perform charge assignment asynchronously for PPPM
	void set_double_precompute(PPPM<numtyp,acctyp,
	                                         double,_lgpu_double4> *pppm);

  /// Esimate the overhead from GPU calls from multiple procs
  /** \param kernel_calls Number of kernel calls/timestep for timing estimated
    *                     overhead
    * \param gpu_overhead Estimated gpu overhead per timestep (sec)
    * \param driver_overhead Estimated overhead from driver per timestep (s) **/
  void estimate_gpu_overhead(const int kernel_calls, double &gpu_overhead,
                             double &gpu_driver_overhead);

  /// Returns true if double precision is supported on card
  inline bool double_precision() { return gpu->double_precision(); }
  
  /// Output a message with timing information
  void output_times(UCL_Timer &time_pair, Answer<numtyp,acctyp> &ans, 
                    Neighbor &nbor, const double avg_split, 
                    const double max_bytes, const double gpu_overhead,
                    const double driver_overhead, 
                    const int threads_per_atom, FILE *screen);

  /// Output a message with timing information
  void output_kspace_times(UCL_Timer &time_in, UCL_Timer &time_out,
                           UCL_Timer & time_map, UCL_Timer & time_rho,
                           UCL_Timer &time_interp, 
                           Answer<numtyp,acctyp> &ans, 
                           const double max_bytes, const double cpu_time,
                           const double cpu_idle_time, FILE *screen);

  /// Clear all memory on host and device associated with atom and nbor data
  void clear();
  
  /// Clear all memory on host and device
  void clear_device();

  /// Add an answer object for putting forces, energies, etc from GPU to LAMMPS
  inline void add_ans_object(Answer<numtyp,acctyp> *ans)
    { ans_queue.push(ans); }

  /// Add "answers" (force,energies,etc.) into LAMMPS structures
  inline double fix_gpu(double **f, double **tor, double *eatom,
                        double **vatom, double *virial, double &ecoul) {
    atom.data_unavail();
    if (ans_queue.empty()==false) {
      stop_host_timer();
      double evdw=0.0;
      while (ans_queue.empty()==false) {
        evdw+=ans_queue.front()->get_answers(f,tor,eatom,vatom,virial,ecoul);
        ans_queue.pop();
      }                                                 
      return evdw;
    }
    return 0.0;
  }

  /// Start timer on host
  inline void start_host_timer() 
    { _cpu_full=MPI_Wtime(); _host_timer_started=true; }
  
  /// Stop timer on host
  inline void stop_host_timer() { 
    if (_host_timer_started) {
      _cpu_full=MPI_Wtime()-_cpu_full; 
      _host_timer_started=false;
    }
  }
  
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
  /// Return the 'mode' for acceleration: GPU_FORCE, GPU_NEIGH or GPU_HYB_NEIGH
  inline int gpu_mode() const { return _gpu_mode; }
  /// Index of first device used by a node
  inline int first_device() const { return _first_device; }
  /// Index of last device used by a node
  inline int last_device() const { return _last_device; }
  /// Particle split defined in fix
  inline double particle_split() const { return _particle_split; }
  /// Return the initialization count for the device
  inline int init_count() const { return _init_count; }
  /// True if device is being timed
  inline bool time_device() const { return _time_device; }

  /// Return the number of threads accessing memory simulatenously
  inline int num_mem_threads() const { return _num_mem_threads; }
  /// Return the number of threads per atom for pair styles
  inline int threads_per_atom() const { return _threads_per_atom; }
  /// Return the number of threads per atom for pair styles using charge
  inline int threads_per_charge() const { return _threads_per_charge; }
  /// Return the min of the pair block size or the device max block size
  inline int pair_block_size() const { return _block_pair; }
  /// Return the maximum number of atom types that can be used with shared mem
  inline int max_shared_types() const { return _max_shared_types; }
  /// Return the maximum order for PPPM splines
  inline int pppm_max_spline() const { return _pppm_max_spline; }
  /// Return the block size for PPPM kernels
  inline int pppm_block() const { return _pppm_block; }
  /// Return the block size for neighbor binning
  inline int block_cell_2d() const { return _block_cell_2d; }
  /// Return the block size for atom mapping for neighbor builds
  inline int block_cell_id() const { return _block_cell_id; }
  /// Return the block size for neighbor build kernel
  inline int block_nbor_build() const { return _block_nbor_build; }
  /// Return the block size for "bio" pair styles
  inline int block_bio_pair() const { return _block_bio_pair; }
  /// Return the block size for "ellipse" pair styles
  inline int block_ellipse() const { return _block_ellipse; }
  /// Return the maximum number of atom types for shared mem with "bio" styles
  inline int max_bio_shared_types() const { return _max_bio_shared_types; }
  /// Architecture gpu code compiled for (returns 0 for OpenCL)
  inline double ptx_arch() const { return _ptx_arch; }
  /// Number of threads executing concurrently on same multiproc
  inline int warp_size() const { return _warp_size; }

  // -------------------- SHARED DEVICE ROUTINES -------------------- 
  // Perform asynchronous zero of integer array 
  void zero(UCL_D_Vec<int> &mem, const int numel) {
    int num_blocks=static_cast<int>(ceil(static_cast<double>(numel)/
                                    _block_pair));
    k_zero.set_size(num_blocks,_block_pair);
    k_zero.run(&mem,&numel);
  }

  // -------------------------- DEVICE DATA ------------------------- 

  /// Geryon Device
  UCL_Device *gpu;

  enum{GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH};

  // --------------------------- ATOM DATA -------------------------- 

  /// Atom Data
  Atom<numtyp,acctyp> atom;

  // --------------------------- NBOR DATA ----------------------------
  
  /// Neighbor Data
  NeighborShared _neighbor_shared;

  // ------------------------ LONG RANGE DATA -------------------------
  
  // Long Range Data
  int _long_range_precompute;
  PPPM<numtyp,acctyp,float,_lgpu_float4> *pppm_single;
  PPPM<numtyp,acctyp,double,_lgpu_double4> *pppm_double;
  /// Precomputations for long range charge assignment (asynchronously)
  inline void precompute(const int ago, const int nlocal, const int nall,
                         double **host_x, int *host_type, bool &success,
                         double *charge, double *boxlo, double *prd) {
    if (_long_range_precompute==1)
      pppm_single->precompute(ago,nlocal,nall,host_x,host_type,success,charge,
                              boxlo,prd);
    else if (_long_range_precompute==2)
      pppm_double->precompute(ago,nlocal,nall,host_x,host_type,success,charge,
                              boxlo,prd);
  }
  
  inline std::string compile_string() { return _ocl_compile_string; }

 private:
  std::queue<Answer<numtyp,acctyp> *> ans_queue;
  int _init_count;
  bool _device_init, _host_timer_started, _time_device;
  MPI_Comm _comm_world, _comm_replica, _comm_gpu;
  int _procs_per_gpu, _gpu_rank, _world_me, _world_size, _replica_me, 
      _replica_size;
  int _gpu_mode, _first_device, _last_device, _nthreads;
  double _particle_split;
  double _cpu_full;
  double _ptx_arch;
  double _cell_size; // -1 if the cutoff is used

  int _num_mem_threads, _warp_size, _threads_per_atom, _threads_per_charge;
  int _pppm_max_spline, _pppm_block;
  int _block_pair, _block_ellipse, _max_shared_types;
  int _block_cell_2d, _block_cell_id, _block_nbor_build;
  int _block_bio_pair, _max_bio_shared_types;

  UCL_Program *dev_program;
  UCL_Kernel k_zero, k_info;
  bool _compiled;
  int compile_kernels();

  int _data_in_estimate, _data_out_estimate;
  
  std::string _ocl_vendor_name, _ocl_vendor_string, _ocl_compile_string;
  int set_ocl_params(char *);
  
  template <class t>
  inline std::string toa(const t& in) {
    std::ostringstream o;
    o.precision(2);
    o << in;
    return o.str();
  }

};

}

#endif
