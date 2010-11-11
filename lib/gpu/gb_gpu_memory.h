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

#ifndef GB_GPU_MEMORY_H
#define GB_GPU_MEMORY_H

#define BLOCK_1D 64

#include "pair_gpu_device.h"
#include "pair_gpu_balance.h"
#include "mpi.h"

template <class numtyp, class acctyp>
class GB_GPU_Memory {
 public:
  GB_GPU_Memory();
  ~GB_GPU_Memory(); 

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param gpu_nbor true if neighboring performed on device
    * \param max_nbors initial number of rows in the neighbor matrix
    * \param cell_size cutoff + skin
    * \param gpu_split fraction of particles handled by device 
    * \return false if there is not sufficient memory or device init prob **/
  bool init(const int ntypes, const double gamma,
            const double upsilon, const double mu, double **host_shape,
            double **host_well, double **host_cutsq, double **host_sigma, 
            double **host_epsilon, double *host_lshape, int **h_form,
            double **host_lj1, double **host_lj2, double **host_lj3, 
            double **host_lj4, double **host_offset, 
            const double *host_special_lj, const int nlocal, const int nall, 
            const int max_nbors, const double cell_size,
            const double gpu_split, FILE *screen);

  /// Check if there is enough storage for atom arrays and realloc if not
  /** \param success set to false if insufficient memory **/
  inline void resize_atom(const int inum, const int nall, bool &success) {
    atom->resize(inum, nall, success);
    if (multiple_forms) atom->dev_ans.zero();
    double bytes=atom->gpu_bytes()+nbor->gpu_bytes();
    if (bytes>_max_bytes)
      _max_bytes=bytes;
  }

  /// Check if there is enough storage for neighbors and realloc if not
  /** \param nlocal number of particles whose nbors must be stored on device
    * \param host_inum number of particles whose nbors need to copied to host
    * \param current maximum number of neighbors
    * \param olist_size size of list of particles from CPU neighboring
    * \note host_inum is 0 if the host is performing neighboring
    * \note if GPU is neighboring nlocal+host_inum=total number local particles
    * \note if CPU is neighboring olist_size=total number of local particles 
    * \note if GPU is neighboring olist_size=0 **/
  inline void resize_local(const int nlocal, const int host_inum,
                           const int max_nbors, const int olist_size,
                           bool &success) {
    if (olist_size>static_cast<int>(host_olist.numel())) {
      host_olist.clear();
      int new_size=static_cast<int>(static_cast<double>(olist_size)*1.10);
      success=success && (host_olist.alloc(new_size,*ucl_device)==UCL_SUCCESS);
    }
    nbor->resize(nlocal,host_inum,max_nbors,success);
    double bytes=atom->gpu_bytes()+nbor->gpu_bytes();
    if (bytes>_max_bytes)
      _max_bytes=bytes;
  }

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();
 
  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  /// Accumulate timers
  inline void acc_timers() {
    if (nbor_time_avail) {
      nbor->time_nbor.add_to_total();
      nbor->time_kernel.add_to_total();
      nbor_time_avail=false;
    }
    time_kernel.add_to_total();
    time_gayberne.add_to_total();
    if (multiple_forms) {
      time_kernel2.add_to_total();
      time_gayberne2.add_to_total();
      time_pair.add_to_total();
    }
    atom->acc_timers();
  }
  
  /// Accumulate timers
  inline void zero_timers() {
    nbor_time_avail=false;
    time_kernel.zero();
    time_gayberne.zero();
    if (multiple_forms) {
      time_kernel2.zero();
      time_gayberne2.zero();
      time_pair.zero();
    }
    atom->zero_timers();
  }

  // -------------------------- DEVICE DATA ------------------------- 
  
  /// Device Properties and Atom and Neighbor storage
  PairGPUDevice<numtyp,acctyp> *device;
  /// Geryon device
  UCL_Device *ucl_device;
  
  /// Device Error Flag - Set if a bad matrix inversion occurs
  UCL_D_Vec<int> dev_error;
  /// Device timers
  UCL_Timer time_kernel, time_gayberne, time_kernel2, time_gayberne2, time_pair;
  /// Host device load balancer
  PairGPUBalance<numtyp,acctyp> hd_balancer;
  /// LAMMPS pointer for screen output
  FILE *screen;
  
  // --------------------------- TYPE DATA -------------------------- 

  /// lj1.x = lj1, lj1.y = lj2, lj1.z = cutsq, lj1.w = form
  UCL_D_Vec<numtyp4> lj1;
  /// lj3.x = lj3, lj3.y = lj4, lj3.z = offset
  UCL_D_Vec<numtyp4> lj3;
  /// sigma_epsilon.x = sigma, sigma_epsilon.y = epsilon
  UCL_D_Vec<numtyp2> sigma_epsilon;
  /// cut_form.x = cutsq, cut_form.y = form
  UCL_D_Vec<numtyp2> cut_form;
  // 0 - gamma, 1-upsilon, 2-mu, 3-special_lj[0], 4-special_lj[1], ...
  UCL_D_Vec<numtyp> gamma_upsilon_mu;
  
  // True if we want to use fast GB-sphere or sphere-sphere calculations 
  bool multiple_forms;
  int **host_form;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;
  int _lj_types;
   
  // --------------------------- ATOM DATA -------------------------- 

  /// Atom Data
  PairGPUAtom<numtyp,acctyp> *atom;

  /// Aspherical Const Data for Atoms
  UCL_D_Vec<numtyp4> shape, well;
  /// Aspherical Const Data for Atoms
  UCL_D_Vec<numtyp> lshape;

  int last_ellipse, max_last_ellipse;

  // --------------------------- NBOR DATA ----------------------------

  /// Neighbor data
  PairGPUNbor *nbor;
  /// ilist with particles sorted by type
  UCL_H_Vec<int> host_olist;
  /// True if we should accumulate the neighbor timer
  bool nbor_time_avail;

  // ------------------------- DEVICE KERNELS -------------------------
  UCL_Program *pair_program, *gb_program, *gb_lj_program;
  UCL_Kernel k_gb_nbor_fast, k_gb_nbor;
  UCL_Kernel k_gayberne, k_sphere_gb, k_lj_fast, k_lj;
  inline int block_size() { return _block_size; }

 private:
  bool _allocated, _compiled;
  int _block_size;
  double _max_bytes;
  
  void compile_kernels(UCL_Device &dev);
};

#endif

