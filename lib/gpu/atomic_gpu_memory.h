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

#ifndef ATOMIC_GPU_MEMORY_H
#define ATOMIC_GPU_MEMORY_H

#define BLOCK_1D 64

#include "pair_gpu_device.h"
#include "pair_gpu_balance.h"
#include "mpi.h"

#ifdef USE_OPENCL
#include "geryon/ocl_texture.h"
#else
#include "geryon/nvd_texture.h"
#endif

template <class numtyp, class acctyp>
class AtomicGPUMemory {
 public:
  AtomicGPUMemory();
  virtual ~AtomicGPUMemory();

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param max_nbors initial number of rows in the neighbor matrix
    * \param cell_size cutoff + skin
    * \param gpu_split fraction of particles handled by device **/
  bool init_atomic(const int nlocal, const int nall, const int max_nbors,
                   const int maxspecial, const double cell_size, 
                   const double gpu_split, FILE *screen, 
                   const char *pair_program);

  /// Check if there is enough storage for atom arrays and realloc if not
  /** \param success set to false if insufficient memory **/
  inline void resize_atom(const int inum, const int nall, bool &success) {
    if (atom->resize(inum, nall, success))
      pos_tex.bind_float(atom->dev_x,4);
  }

  /// Check if there is enough storage for neighbors and realloc if not
  /** \param nlocal number of particles whose nbors must be stored on device
    * \param host_inum number of particles whose nbors need to copied to host
    * \param current maximum number of neighbors
    * \note olist_size=total number of local particles **/
  inline void resize_local(const int inum, const int max_nbors, bool &success) {
    nbor->resize(inum,max_nbors,success);
  }

  /// Check if there is enough storage for neighbors and realloc if not
  /** \param nlocal number of particles whose nbors must be stored on device
    * \param host_inum number of particles whose nbors need to copied to host
    * \param current maximum number of neighbors
    * \note host_inum is 0 if the host is performing neighboring
    * \note nlocal+host_inum=total number local particles
    * \note olist_size=0 **/
  inline void resize_local(const int inum, const int host_inum, 
                           const int max_nbors, bool &success) {
    nbor->resize(inum,host_inum,max_nbors,success);
  }

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear_atomic();

  /// Returns memory usage on device per atom
  int bytes_per_atom_atomic(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage_atomic() const;

  /// Accumulate timers
  inline void acc_timers() {
    if (nbor_time_avail) {
      nbor->time_nbor.add_to_total();
      nbor->time_kernel.add_to_total();
      nbor_time_avail=false;
    }
    time_pair.add_to_total();
    atom->acc_timers();
  }

  /// Zero timers
  inline void zero_timers() {
    nbor_time_avail=false;
    time_pair.zero();
    atom->zero_timers();
  }

  /// Copy neighbor list from host
  int * reset_nbors(const int nall, const int inum, int *ilist, int *numj,
                    int **firstneigh, bool &success);

  /// Build neighbor list on device
  void build_nbor_list(const int inum, const int host_inum,
                       const int nall, double **host_x, int *host_type,
                       double *boxlo, double *boxhi, int *tag, int **nspecial, 
                       int **special, bool &success);

  /// Pair loop with host neighboring
  void compute(const int timestep, const int f_ago, const int inum_full,
               const int nall, double **host_x, int *host_type,
               int *ilist, int *numj, int **firstneigh, const bool eflag,
               const bool vflag, const bool eatom, const bool vatom,
               int &host_start, const double cpu_time, bool &success);

  /// Pair loop with device neighboring
  int * compute(const int timestep, const int ago, const int inum_full,
                const int nall, double **host_x, int *host_type, double *boxlo,
                double *boxhi, int *tag, int **nspecial,
                int **special, const bool eflag, const bool vflag, 
                const bool eatom, const bool vatom, int &host_start, 
                const double cpu_time, bool &success);

  // -------------------------- DEVICE DATA ------------------------- 

  /// Device Properties and Atom and Neighbor storage
  PairGPUDevice<numtyp,acctyp> *device;

  /// Geryon device
  UCL_Device *ucl_device;

  /// Device Timers
  UCL_Timer time_pair;

  /// Host device load balancer
  PairGPUBalance<numtyp,acctyp> hd_balancer;

  /// LAMMPS pointer for screen output
  FILE *screen;

  // --------------------------- ATOM DATA --------------------------

  /// Atom Data
  PairGPUAtom<numtyp,acctyp> *atom;


  // --------------------------- NBOR DATA ----------------------------

  /// Neighbor data
  PairGPUNbor *nbor;

  /// True if we need to accumulate time for neighboring
  bool nbor_time_avail;

  // ------------------------- DEVICE KERNELS -------------------------
  UCL_Program *pair_program;
  UCL_Kernel k_pair_fast, k_pair;
  inline int block_size() { return _block_size; }

  // --------------------------- TEXTURES -----------------------------
  UCL_Texture pos_tex;

 protected:
  bool _compiled;
  int _block_size;
  double _max_bytes, _max_an_bytes;

  void compile_kernels(UCL_Device &dev, const char *pair_string);

  virtual void loop(const bool _eflag, const bool _vflag) = 0;
};

#endif


