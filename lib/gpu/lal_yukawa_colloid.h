/***************************************************************************
                              yukawa_colloid.h
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Class for acceleration of the yukawa/colloid pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#ifndef LAL_YUKAWA_COLLOID_H
#define LAL_YUKAWA_COLLOID_H

#include "lal_base_atomic.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class YukawaColloid : public BaseAtomic<numtyp, acctyp> {
 public:
  YukawaColloid();
  ~YukawaColloid(); 

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
  int init(const int ntypes, double **host_cutsq,
           double **host_a, double **host_offset, double *host_special_lj,
           const int nlocal, const int nall, const int max_nbors, 
           const int maxspecial, const double cell_size, 
           const double gpu_split, FILE *screen, const double kappa);

  inline void cast_rad_data(double* rad) {
    int nall = this->atom->nall();
    if (_shared_view) {
      c_rad.host.view((numtyp*)rad,nall,*(this->ucl_device));
      c_rad.device.view(c_rad.host);
    } else {
      if (sizeof(numtyp)==sizeof(double))
        memcpy(c_rad.host.begin(),rad,nall*sizeof(numtyp));
      else
        for (int i=0; i<nall; i++) c_rad[i]=rad[i];
    }
  }

  // Copy rad to device asynchronously
  inline void add_rad_data() {
    c_rad.update_device(this->atom->nall(),true);
  }

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;
  
  /// Pair loop with host neighboring
  void compute(const int f_ago, const int inum_full, 
               const int nall, double **host_x, int *host_type, 
               int *ilist, int *numj, int **firstneigh, 
               const bool eflag, const bool vflag,
               const bool eatom, const bool vatom, int &host_start,
               const double cpu_time, bool &success, double *rad);
               
  /// Pair loop with device neighboring
  int** compute(const int ago, const int inum_full, const int nall,
                double **host_x, int *host_type, double *sublo,
                double *subhi, tagint *tag, int **nspecial,
                tagint **special, const bool eflag, const bool vflag, 
                const bool eatom, const bool vatom, int &host_start, 
                int **ilist, int **jnum, const double cpu_time, 
                bool &success, double *rad);

  // --------------------------- TEXTURES -----------------------------
  UCL_Texture rad_tex;

  // --------------------------- TYPE DATA --------------------------

  /// coeff.x = a, coeff.y = offset, coeff.z = cutsq
  UCL_D_Vec<numtyp4> coeff;
  /// Special LJ values
  UCL_D_Vec<numtyp> sp_lj;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types 
  int _lj_types;

  int _max_rad_size;

  numtyp _kappa;

  /// Per-atom arrays
  UCL_Vector<numtyp,numtyp> c_rad;

 private:
  bool _shared_view;
  bool _allocated;
  void loop(const bool _eflag, const bool _vflag);
};

}

#endif
