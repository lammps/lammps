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

#include "lal_base_atomic.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class EAM : public BaseAtomic<numtyp, acctyp> {
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
          int **host_type2rhor, int **host_type2z2r, int *host_type2frho,
          double ***host_rhor_spline, double ***host_z2r_spline,
          double ***host_frho_spline,
          double rdr, double rdrho, int nrhor, int nrho, 
          int nz2r, int nfrho, int nr,
          const int nlocal, const int nall, const int max_nbors,
          const int maxspecial, const double cell_size,
          const double gpu_split, FILE *_screen);
  
  // Cast fp to write buffer
  template<class cpytyp>
  inline void cast_fp_data(cpytyp *host_ptr) {
    if (_fp_avail==false) {
      double t=MPI_Wtime();
      int nall = this->atom->nall();
      if (this->ucl_device->device_type()==UCL_CPU) {
        if (sizeof(numtyp)==sizeof(double)) {
          host_fp.view((numtyp*)host_ptr,nall,*(this->ucl_device));
          dev_fp.view(host_fp);
        } else
          for (int i=0; i<nall; i++) host_fp[i]=host_ptr[i];
      } else {
        if (sizeof(numtyp)==sizeof(double))
          memcpy(host_fp.begin(),host_ptr,nall*sizeof(numtyp));
        else
          for (int i=0; i<nall; i++) host_fp[i]=host_ptr[i];
      }
    }
  }

  // Copy charges to device asynchronously
  inline void add_fp_data() {
    if (_fp_avail==false) {
      ucl_copy(dev_fp,host_fp,this->atom->nall(),true);
      _fp_avail=true;
    }
  }
  
  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;
  
  /// Pair loop with host neighboring
  void compute(const int f_ago, const int inum_full, const int nall,
               double **host_x, int *host_type, int *ilist, int *numj,
               int **firstneigh, const bool eflag, const bool vflag,
               const bool eatom, const bool vatom, int &host_start,
               const double cpu_time, bool &success, double *charge,
               const int nlocal, double *boxlo, double *prd);
               
  /// Pair loop with device neighboring
  int** compute(const int ago, const int inum_full, const int nall,
                double **host_x, int *host_type, double *sublo,
                double *subhi, int *tag, int **nspecial,
                int **special, const bool eflag, const bool vflag, 
                const bool eatom, const bool vatom, int &host_start, 
                int **ilist, int **numj, const double cpu_time, bool &success,
                double *charge, double *boxlo, double *prd, int &inum);

  /// Pair loop with host neighboring
  void compute2(int *ilist, const bool eflag, const bool vflag,
                    const bool eatom, const bool vatom, double *host_fp);

  // ------------------------- DEVICE KERNELS -------------------------
  UCL_Kernel k_energy;
  
  // --------------------------- TEXTURES -----------------------------
  UCL_Texture fp_tex;

  // --------------------------- TYPE DATA --------------------------
    
  UCL_D_Vec<numtyp2> type2rhor_z2r;
  
  UCL_D_Vec<numtyp> type2frho;
  
  UCL_D_Vec<numtyp> frho_spline, rhor_spline, z2r_spline;
  
  numtyp _cutforcesq,_rdr,_rdrho;
  
  int _nfrho,_nrhor,_nrho,_nz2r,_nr;
  
  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;
  
  /// Number of atom types 
  int _ntypes;
  
  int _max_fp_size;
  
  /// Per-atom arrays
  UCL_H_Vec<acctyp> host_fp;
  UCL_D_Vec<acctyp> dev_fp;
  
protected:
  bool _allocated, _fp_avail;
  void loop(const bool _eflag, const bool _vflag);
  void loop2(const bool _eflag, const bool _vflag);
};

}

#endif

