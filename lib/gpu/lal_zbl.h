/***************************************************************************
                                    zbl.h
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Class for acceleration of the zbl pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : ndactrung@gmail.com
 ***************************************************************************/

#ifndef LAL_ZBL_H
#define LAL_ZBL_H

#include "lal_base_atomic.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class ZBL : public BaseAtomic<numtyp, acctyp> {
 public:
  ZBL();
  ~ZBL();

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
  int init(const int ntypes, double **host_cutsq, double **host_sw1,
           double **host_sw2, double **host_sw3, double **host_sw4, double **host_sw5,
           double **host_d1a, double **host_d2a, double **host_d3a, double **host_d4a,
           double **host_zze, double cut_globalsq, double cut_innersq, double cut_inner,
           const int nlocal, const int nall, const int max_nbors,
           const int maxspecial, const double cell_size,
           const double gpu_split, FILE *screen);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// coeff1.x = sw1, coeff1.y = sw2, coeff1.z = zze, coeff1.w = cutsq
  UCL_D_Vec<numtyp4> coeff1;
  /// coeff2.x = d1a, coeff2.y = d2a, coeff2.z = d3a, coeff2.w = d4a
  UCL_D_Vec<numtyp4> coeff2;
  /// coeff3.x = sw3, coeff3.y = sw4, coeff3.z = sw5;
  UCL_D_Vec<numtyp4> coeff3;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  double _cut_globalsq;
  double _cut_innersq;
  double _cut_inner;

  /// Number of atom types
  int _lj_types;

 private:
  bool _allocated;
  void loop(const bool _eflag, const bool _vflag);
};

}

#endif
