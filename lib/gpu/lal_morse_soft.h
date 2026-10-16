/***************************************************************************
                                 morse_soft.h
                             -------------------
     Claude Opus 5.5 (Anthropic), under the direction of Axel Kohlmeyer

  Class for acceleration of the morse/soft pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#ifndef LAL_MORSE_SOFT_H
#define LAL_MORSE_SOFT_H

#include "lal_base_atomic.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class MorseSoft : public BaseAtomic<numtyp, acctyp> {
 public:
  MorseSoft();
  ~MorseSoft();

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
  int init(const int ntypes, double **host_cutsq,
           double **host_morse1, double **host_r0, double **host_alpha,
           double **host_d0, double **host_offset, double **host_bs, double **host_vs, double *host_special_lj,
           const int nlocal, const int nall, const int max_nbors,
           const int maxspecial, const double cell_size,
           FILE *screen);

  /// Send updated coeffs from host to device (to be compatible with fix adapt)
  void reinit(const int ntypes, double **host_cutsq, double **host_morse1, double **host_r0,
              double **host_alpha, double **host_d0, double **host_offset, double **host_bs,
              double **host_vs);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// mor1.x = cutsq, mor1.y = morse1, mor1.z = r0, mor1.w = alpha
  UCL_D_Vec<numtyp4> mor1;
  /// mor2.x = d0, mor2.y = offset, mor2.z = bs, mor2.w = vs
  UCL_D_Vec<numtyp4> mor2;
  /// Special LJ values
  UCL_D_Vec<numtyp> sp_lj;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _types;

 private:
  bool _allocated;
  int loop(const int eflag, const int vflag);
};

}

#endif
