/***************************************************************************
                                 coul_cut_soft.h
                             -------------------
     Claude Opus 5.5 (Anthropic), under the direction of Axel Kohlmeyer

  Class for acceleration of the coul/cut/soft pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#ifndef LAL_COUL_CUT_SOFT_H
#define LAL_COUL_CUT_SOFT_H

#include "lal_base_charge.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class CoulCutSoft : public BaseCharge<numtyp, acctyp> {
 public:
  CoulCutSoft();
  ~CoulCutSoft();

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
  int init(const int ntypes, double **host_lam1, double **host_lam2,
           double **host_cutsq, double *host_special_coul,
           const int nlocal, const int nall, const int max_nbors,
           const int maxspecial, const double cell_size,
           FILE *screen, const double qqrd2e);

  /// Send updated coeffs from host to device (to be compatible with fix adapt)
  void reinit(const int ntypes, double **host_lam1, double **host_lam2);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// x = lambda^n, y = alpha_c (1-lambda)^2
  UCL_D_Vec<numtyp2> lam;
  /// cutsq
  UCL_D_Vec<numtyp> cutsq;
  /// Special CoulCutSoft values [0-3]
  UCL_D_Vec<numtyp> sp_cl;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _lj_types;

  numtyp _qqrd2e;

 private:
  bool _allocated;
  int loop(const int eflag, const int vflag);
};

}

#endif
