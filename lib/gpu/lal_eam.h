/***************************************************************************
                              lal_eam.h
                             -------------------
                      W. Michael Brown, Trung Dac Nguyen (ORNL)

  Class for acceleration of the eam pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov nguyentd@ornl.gov
 ***************************************************************************/

#ifndef LAL_EAM_H
#define LAL_EAM_H

#include "lal_base_charge.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class EAM : public BaseCharge<numtyp, acctyp> {
 public:
  EAM();
  ~EAM();
  
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
  int init(const int ntypes, double host_cutforcesq,
          int **host_type2rhor, int **host_type2z2r,
          double ***host_rhor_spline, double ***host_z2r_spline,
          double rdr, int nrhor, int nz2r, int nr,
          const int nlocal, const int nall, const int max_nbors,
          const int maxspecial, const double cell_size,
          const double gpu_split, FILE *_screen);
  
  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;
  
  // --------------------------- TYPE DATA --------------------------
    
  UCL_D_Vec<numtyp2> type2rhor_z2r;
  
  UCL_D_Vec<numtyp> rhor_spline;
  
  UCL_D_Vec<numtyp> z2r_spline;
  
  numtyp _cutforcesq,_rdr;
  
  int _nrhor,_nz2r,_nr;
  
  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;
  
  /// Number of atom types 
  int _ntypes;
  
  // --------------------------- TEXTURES -----------------------------
  UCL_Texture fp_tex;
  
 private:
  bool _allocated;
  void loop(const bool _eflag, const bool _vflag);
};

}

#endif

