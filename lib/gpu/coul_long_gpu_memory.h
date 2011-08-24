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

#ifndef CL_GPU_MEMORY_H
#define CL_GPU_MEMORY_H

#include "charge_gpu_memory.h"

template <class numtyp, class acctyp>
class CL_GPU_Memory : public ChargeGPUMemory<numtyp, acctyp> {
 public:
  CL_GPU_Memory();
  ~CL_GPU_Memory();

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param max_nbors initial number of rows in the neighbor matrix
    * \param cell_size cutoff + skin
    * \param gpu_split fraction of particles handled by device
    * 
    * Returns:
    * -  0 if successfull
    * - -1 if fix gpu not found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card **/
  int init(const int nlocal, const int nall, const int max_nbors,
           const int maxspecial, const double cell_size,
	   const double gpu_split, FILE *screen, 
	   const double host_cut_coulsq, double *host_special_coul,
	   const double qqrd2e, const double g_ewald);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// lj1 dummy
  UCL_D_Vec<numtyp4> lj1;
  /// lj3 dummy
  UCL_D_Vec<numtyp4> lj3;
  /// Special Coul values [0-3]
  UCL_D_Vec<numtyp> sp_cl;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _lj_types;

  numtyp _cut_coulsq, _qqrd2e, _g_ewald;

 private:
  bool _allocated;
  void loop(const bool _eflag, const bool _vflag);
};

#endif

