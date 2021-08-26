/***************************************************************************
                                  amoeba.h
                             -------------------
                          Trung Dac Nguyen (Northwestern)

  Class for acceleration of the amoeba pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : trung.nguyen@northwestern.edu
 ***************************************************************************/

#ifndef LAL_AMOEBA_H
#define LAL_AMOEBA_H

#include "lal_base_amoeba.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class Amoeba : public BaseAmoeba<numtyp, acctyp> {
 public:
  Amoeba();
  ~Amoeba();

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
  int init(const int ntypes, const int max_amtype, const double *host_pdamp,
           const double *host_thole, const double *host_special_polar_wscale,
           const double *host_special_polar_piscale,
           const double *host_special_polar_pscale,
           const int nlocal, const int nall, const int max_nbors,
           const int maxspecial, const int maxspecial15, const double cell_size,
           const double gpu_split, FILE *_screen,
           const double aewald, const double felec,
           const double off2, const double polar_dscale,
           const double polar_uscale);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// pdamp = damping.x; thole = damping.y
  UCL_D_Vec<numtyp4> damping;
  /// Special polar values [0-4]: 
  ///   sp_polar.x = special_polar_wscale
  ///   sp_polar.y special_polar_pscale,
  ///   sp_polar.z = special_polar_piscale
  UCL_D_Vec<numtyp4> sp_polar;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _lj_types;

  numtyp _aewald, _felec, _off2, _polar_dscale, _polar_uscale;
  numtyp _qqrd2e;

 private:
  bool _allocated;
  int loop(const int eflag, const int vflag);
};

}

#endif
