/***************************************************************************
                                    eam.h
                             -------------------
                   Trung Dac Nguyen, W. Michael Brown (ORNL)

  Class for acceleration of the eam pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov nguyentd@ornl.gov
 ***************************************************************************/

#ifndef LAL_EAM_H
#define LAL_EAM_H

#include "lal_precision.h"
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
    * -  0 if successful
    * - -1 if fix gpu not found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card **/
  int init(const int ntypes, double host_cutforcesq, int **host_type2rhor,
           int **host_type2z2r, int *host_type2frho, double ***host_rhor_spline,
           double ***host_z2r_spline, double ***host_frho_spline, double** host_cutsq,
           double rdr, double rdrho, double rhomax, int nrhor, int nrho, int nz2r,
           int nfrho, int nr, const int nlocal, const int nall,
           const int max_nbors, const int maxspecial, const double cell_size,
           const double gpu_split, FILE *_screen);

  // Copy charges to device asynchronously
  inline void add_fp_data() {
    int nghost=this->atom->nall()-_nlocal;
    if (nghost>0) {
      UCL_H_Vec<numtyp> host_view;
      UCL_D_Vec<numtyp> dev_view;
      host_view.view_offset(_nlocal,_fp.host);
      dev_view.view_offset(_nlocal,_fp.device);
      ucl_copy(dev_view,host_view,nghost,true);
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
  void compute(const int f_ago, const int inum_full, const int, const int nall,
               double **host_x, int *host_type, int *ilist, int *numj,
               int **firstneigh, const bool eflag, const bool vflag,
               const bool eatom, const bool vatom, int &host_start,
               const double cpu_time, bool &success,
               void **fp_ptr);

  /// Pair loop with device neighboring
  int** compute(const int ago, const int inum_full, const int nall,
                double **host_x, int *host_type, double *sublo,
                double *subhi, tagint *tag, int **nspecial,
                tagint **special, const bool eflag, const bool vflag,
                const bool eatom, const bool vatom, int &host_start,
                int **ilist, int **numj, const double cpu_time, bool &success,
                int &inum, void **fp_ptr);

  /// Pair loop with host neighboring
  void compute2(int *ilist, const bool eflag, const bool vflag,
                const bool eatom, const bool vatom);

  // ------------------------- DEVICE KERNELS -------------------------
  UCL_Kernel k_energy, k_energy_fast, k_energy_fast_noev, *k_energy_sel;

  // --------------------------- TEXTURES -----------------------------
  UCL_Texture fp_tex;
  UCL_Texture rhor_spline1_tex, rhor_spline2_tex;
  UCL_Texture frho_spline1_tex, frho_spline2_tex;
  UCL_Texture z2r_spline1_tex, z2r_spline2_tex;

  // --------------------------- DEVICE DATA --------------------------

  /// Device Timers
  UCL_Timer time_pair2, time_fp1, time_fp2;

  // --------------------------- TYPE DATA --------------------------

  UCL_D_Vec<int2> type2rhor_z2r;
  UCL_D_Vec<int> type2frho;

  UCL_D_Vec<numtyp4> z2r_spline1, z2r_spline2;
  UCL_D_Vec<numtyp4> frho_spline1, frho_spline2;
  UCL_D_Vec<numtyp4> rhor_spline1, rhor_spline2;

  UCL_D_Vec<numtyp> cutsq;

  numtyp _cutforcesq,_rdr,_rdrho, _rhomax;

  int _nfrho,_nrhor,_nrho,_nz2r,_nr;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _ntypes;

  int _max_fp_size;

  /// True of energy kernels are compiled
  bool _compiled_energy;

  /// Per-atom arrays
  UCL_Vector<numtyp,numtyp> _fp;

protected:
  bool _allocated;
  int _nlocal;
  int loop(const int eflag, const int vflag);
  void loop2(const bool eflag, const bool vflag);
};

}

#endif

