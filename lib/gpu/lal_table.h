/***************************************************************************
                                 lal_table.h
                             -------------------
                           Trung Dac Nguyen (ORNL)

  Class for acceleration of the table pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#ifndef LAL_LJ_H
#define LAL_LJ_H

#include "lal_base_atomic.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class Table : public BaseAtomic<numtyp, acctyp> {
 public:
  Table();
  ~Table(); 

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
  int init(const int ntypes, double** cutsq, double ***host_table_coeffs,
           double **host_table_data, 
           double *host_special_lj,
           const int nlocal, const int nall, const int max_nbors, 
           const int maxspecial, const double cell_size, 
           const double gpu_split, FILE *screen,
           int tabstyle, int ntables, int tablength);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;
  
  // ------------------------- DEVICE KERNELS -------------------------
  UCL_Kernel k_pair_linear, k_pair_linear_fast;
  UCL_Kernel k_pair_spline, k_pair_spline_fast;
  UCL_Kernel k_pair_bitmap, k_pair_bitmap_fast;
  
  // --------------------------- TYPE DATA --------------------------

  UCL_D_Vec<int> tabindex, nshiftbits, nmask;
  
  /// coeff2.x = innersq, coeff2.y = invdelta, coeff2.z = deltasq6, 
  UCL_D_Vec<numtyp4> coeff2;
  
  /// coeff3.x = rsq, coeff3.y = e, coeff3.z = f
  UCL_D_Vec<numtyp4> coeff3;
  
  /// coeff4.x = de, coeff4.y = df
  UCL_D_Vec<numtyp4> coeff4;
  
  UCL_D_Vec<numtyp> cutsq;
  
  /// Special LJ values
  UCL_D_Vec<numtyp> sp_lj;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types 
  int _lj_types;
  
  /// Table style, length and number of tables
  int _tabstyle,_tablength,_ntables;
  
 private:
  bool _allocated, _compiled_styles;
  
  void loop(const bool _eflag, const bool _vflag);
};

}

#endif
