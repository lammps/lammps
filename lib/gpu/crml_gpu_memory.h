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

#ifndef CRML_GPU_MEMORY_H
#define CRML_GPU_MEMORY_H

#include "charge_gpu_memory.h"

template <class numtyp, class acctyp>
class CRML_GPU_Memory : public ChargeGPUMemory<numtyp, acctyp> {
 public:
  CRML_GPU_Memory();
  ~CRML_GPU_Memory();

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param max_nbors initial number of rows in the neighbor matrix
    * \param cell_size cutoff + skin
    * \param gpu_split fraction of particles handled by device **/
  bool init(const int ntypes, double host_cut_bothsq,
            double **host_lj1, double **host_lj2, double **host_lj3,
            double **host_lj4, double **host_offset, double *host_special_lj,
            const int nlocal, const int nall, const int max_nbors, 
            const int maxspecial, const double cell_size, 
            const double gpu_split, FILE *screen, double host_cut_ljsq,
            const double host_cut_coulsq, double *host_special_coul,
            const double qqrd2e, const double g_ewald,
            const double cut_lj_innersq, const double denom_lj, 
            double **epsilon, double **sigma, const bool mix_arithmetic);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// x = lj1, y = lj2, z = lj3, w = lj4
  UCL_D_Vec<numtyp4> lj1;
  /// x = epsilon, y = sigma
  UCL_D_Vec<numtyp2> ljd;
  /// Special LJ values [0-3] and Special Coul values [4-7]
  UCL_D_Vec<numtyp> sp_lj;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types 
  int _lj_types;

  numtyp _qqrd2e, _g_ewald, _denom_lj;

  numtyp _cut_coulsq, _cut_bothsq, _cut_ljsq, _cut_lj_innersq;

 private:
  bool _allocated;
  void loop(const bool _eflag, const bool _vflag);
};

#endif

