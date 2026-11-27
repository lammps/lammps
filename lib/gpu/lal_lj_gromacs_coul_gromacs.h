/***************************************************************************
                                 lj_gromacs_coul_gromacs.h
                             -------------------

  Class for acceleration of the lj/gromacs/coul/gromacs pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#ifndef LAL_LJ_GROMACS_COUL_GROMACS_H
#define LAL_LJ_GROMACS_COUL_GROMACS_H

#include "lal_base_charge.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class LJGromacsCoulGromacs : public BaseCharge<numtyp, acctyp> {
 public:
  LJGromacsCoulGromacs();
  ~LJGromacsCoulGromacs();

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param max_nbors initial number of rows in the neighbor matrix
    * \param cell_size cutoff + skin
    *
    * Returns:
    * -  0 if successful
    * - -1 if fix gpu not found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card **/
  int init(const int ntypes, double **host_lj1, double **host_lj2, double **host_ljsw1,
           double **host_ljsw2, double **host_lj3, double **host_lj4, double **host_ljsw3,
           double **host_ljsw4, double **host_ljsw5, double *host_special_lj,
           double *host_special_coul, const int nlocal, const int nall, const int max_nbors,
           const int maxspecial, const double cell_size, FILE *screen, const double host_qqrd2e,
           const double host_cut_bothsq, const double host_cut_ljsq,
           const double host_cut_lj_innersq, const double host_cut_lj_inner,
           const double host_cut_coulsq, const double host_cut_coul_innersq,
           const double host_cut_coul_inner, const double host_coulsw1, const double host_coulsw2,
           const double host_coulsw3, const double host_coulsw4, const double host_coulsw5);

  /// Send updated coeffs from host to device (to be compatible with fix adapt)
  void reinit(const int ntypes, double **host_lj1, double **host_lj2, double **host_ljsw1,
              double **host_ljsw2, double **host_lj3, double **host_lj4, double **host_ljsw3,
              double **host_ljsw4, double **host_ljsw5);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// lj1.x = lj1, lj1.y = lj2, lj1.z = ljsw1, lj1.w = ljsw2
  UCL_D_Vec<numtyp4> lj1;
  /// lj3.x = lj3, lj3.y = lj4, lj3.z = ljsw3, lj3.w = ljsw4
  UCL_D_Vec<numtyp4> lj3;
  /// ljsw5
  UCL_D_Vec<numtyp> ljsw5;
  /// Special LJ values [0-3] and Special Coul values [4-7]
  UCL_D_Vec<numtyp> sp_lj;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _lj_types;

  numtyp _qqrd2e, _cut_bothsq, _cut_ljsq, _cut_lj_innersq, _cut_lj_inner, _cut_coulsq, _cut_coul_innersq, _cut_coul_inner, _coulsw1, _coulsw2, _coulsw3, _coulsw4, _coulsw5;

 private:
  bool _allocated;
  int loop(const int eflag, const int vflag);
};

}

#endif
