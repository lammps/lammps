/***************************************************************************
                                 coul_cut_soft_gapsys.h
                             -------------------

  Class for acceleration of the coul/cut/soft/gapsys pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#ifndef LAL_COUL_CUT_SOFT_GAPSYS_H
#define LAL_COUL_CUT_SOFT_GAPSYS_H

#include "lal_base_charge.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class CoulCutSoftGapsys : public BaseCharge<numtyp, acctyp> {
 public:
  CoulCutSoftGapsys();
  ~CoulCutSoftGapsys();

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
  int init(const int ntypes, double **host_cutsq, double **host_cut_factor, double *host_special_lj,
           double *host_special_coul, const int nlocal, const int nall, const int max_nbors,
           const int maxspecial, const double cell_size, FILE *screen, const double host_qqrd2e,
           const double host_sigmaq);

  /// Send updated coeffs from host to device (to be compatible with fix adapt)
  void reinit(const int ntypes, double **host_cutsq, double **host_cut_factor);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// cutsq.x = cutsq, cutsq.y = cut_factor = alpha_q lambda^(1/6)
  UCL_D_Vec<numtyp2> cutsq;
  /// Special LJ values [0-3] and Special Coul values [4-7]
  UCL_D_Vec<numtyp> sp_lj;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _lj_types;

  numtyp _qqrd2e, _sigmaq;

 private:
  bool _allocated;
  int loop(const int eflag, const int vflag);
};

}

#endif
