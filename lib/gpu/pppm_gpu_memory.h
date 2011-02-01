/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Charge/Molecular Massively Parallel Simulator
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

#ifndef PPPM_GPU_MEMORY_H
#define PPPM_GPU_MEMORY_H

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
class PPPMGPUMemory {
 public:
  PPPMGPUMemory();
  virtual ~PPPMGPUMemory();

  /// Clear any previous data and set up for a new LAMMPS run
  bool init(const int nlocal, const int nall, FILE *screen, const int order,
            const int nxlo_out, const int nylo_out, const int nzlo_out,
            const int nxhi_out, const int nyhi_out, const int nzhi_out,
            double **rho_coeff);

  /// Check if there is enough storage for atom arrays and realloc if not
  /** \param success set to false if insufficient memory **/
  inline void resize_atom(const int inum, const int nall, bool &success) {
    if (atom->resize(nall, success)) {
      pos_tex.bind_float(atom->dev_x,4);
      q_tex.bind_float(atom->dev_q,1);
    }
    ans->resize(inum,success);
  }

  /// Check if there is enough storage for local atoms and realloc if not
  inline void resize_local(const int inum, bool &success) {
  }

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom() const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  /// Accumulate timers
  inline void acc_timers() {
    atom->acc_timers();
    ans->acc_timers();
  }

  /// Zero timers
  inline void zero_timers() {
    atom->zero_timers();
    ans->zero_timers();
  }

  void compute(const int ago,const int nlocal,const int nall,double **host_x,
               int *host_type,bool &success,double *charge,double *boxlo,
               const double delxinv,const double delyinv,const double delzinv);

  // -------------------------- DEVICE DATA ------------------------- 

  /// Device Properties and Atom and Neighbor storage
  PairGPUDevice<numtyp,acctyp> *device;

  /// Geryon device
  UCL_Device *ucl_device;

  /// Device Timers
  UCL_Timer time_in;

  /// LAMMPS pointer for screen output
  FILE *screen;

  // --------------------------- ATOM DATA --------------------------

  /// Atom Data
  PairGPUAtom<numtyp,acctyp> *atom;


  // --------------------------- GRID DATA --------------------------

  UCL_H_Vec<numtyp> *h_brick;
  UCL_D_Vec<numtyp> *d_brick;
  
  // -------------------------- STENCIL DATA -------------------------
  UCL_D_Vec<numtyp> d_rho_coeff;
  int _order, _nxlo_out, _nylo_out, _nzlo_out, _nxhi_out, _nyhi_out, _nzhi_out;

  // ------------------------ FORCE/ENERGY DATA -----------------------

  PairGPUAns<numtyp,acctyp> *ans;

  // ------------------------- DEVICE KERNELS -------------------------
  UCL_Program *pppm_program;
  UCL_Kernel k_compute;
  inline int block_size() { return _block_size; }

  // --------------------------- TEXTURES -----------------------------
  UCL_Texture pos_tex;
  UCL_Texture q_tex;

 protected:
  bool _allocated, _compiled;
  int _block_size;
  double  _max_bytes, _max_an_bytes;

  void compile_kernels(UCL_Device &dev);
};

#endif

