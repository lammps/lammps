/***************************************************************************
                                   ufm.h
                             -------------------
                            Rodolfo Paula Leite (Unicamp/Brazil)
                            Maurice de Koning (Unicamp/Brazil)

   Class for acceleration of the ufm pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                :  pl.rodolfo@gmail.com
                            dekoning@ifi.unicamp.br
 ***************************************************************************/

#ifndef LAL_UFM_H
#define LAL_UFM_H

#include "lal_base_atomic.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class UFM : public BaseAtomic<numtyp, acctyp> {
 public:
  UFM();
  ~UFM();

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
  int init(const int ntypes, double **host_cutsq,
           double **host_uf1, double **host_uf2, double **host_uf3,
           double **host_uf4, double **host_offset, double *host_special_lj,
           const int nlocal, const int nall, const int max_nbors, 
           const int maxspecial, const double cell_size, 
           const double gpu_split, FILE *screen);
  
  /// Send updated coeffs from host to device (to be compatible with fix adapt)
  void reinit(const int ntypes, double **host_cutsq,
              double **host_uf1, double **host_uf2, double **host_uf3,
              double **host_uf4, double **host_offset);
  
  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// uf1.x = uf1, uf1.y = uf2, uf1.z = cutsq
  UCL_D_Vec<numtyp4> uf1;
  /// uf3.x = uf3, uf3.y = uf4, uf3.z = offset
  UCL_D_Vec<numtyp4> uf3;
  /// Special LJ values
  UCL_D_Vec<numtyp> sp_lj;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types 
  int _lj_types;

 private:
  bool _allocated;
  void loop(const bool _eflag, const bool _vflag);
};

}

#endif
