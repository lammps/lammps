/***************************************************************************
                                charmm_long.h
                             -------------------
                            W. Michael Brown (ORNL)

  Class for acceleration of the charmm/coul/long pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov
 ***************************************************************************/

#ifndef LAL_CHARMM_LONG_H
#define LAL_CHARMM_LONG_H

#include "lal_base_charge.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class CHARMMLong : public BaseCharge<numtyp, acctyp> {
 public:
  CHARMMLong();
  ~CHARMMLong();

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param max_nbors initial number of rows in the neighbor matrix
    * \param cell_size cutoff + skin
    * \param gpu_split fraction of particles handled by device
    *
    * Returns:
    * -  0 if successful
    * - -1 if fix gpu not found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card **/
  int init(const int ntypes, double host_cut_bothsq,
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
  int loop(const int eflag, const int vflag);
};

}

#endif
