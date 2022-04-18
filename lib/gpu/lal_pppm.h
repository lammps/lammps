/***************************************************************************
                                   pppm.h
                             -------------------
                            W. Michael Brown (ORNL)

  Class for PPPM acceleration

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov
 ***************************************************************************/

#ifndef LAL_PPPM_H
#define LAL_PPPM_H

#include "mpi.h"
#include "lal_device.h"

#if defined(USE_OPENCL)
#include "geryon/ocl_texture.h"
#elif defined(USE_CUDART)
#include "geryon/nvc_texture.h"
#elif defined(USE_HIP)
#include "geryon/hip_texture.h"
#else
#include "geryon/nvd_texture.h"
#endif

namespace LAMMPS_AL {

template <class numtyp, class acctyp> class Device;

template <class numtyp, class acctyp, class grdtyp, class grdtyp4>
class PPPM {
 public:
  PPPM();
  virtual ~PPPM();

  /// Clear any previous data and set up for a new LAMMPS run
  /** Success will be:
    * -  0 if successful
    * - -1 if fix gpu not found
    * - -2 if GPU could not be found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card **/
  grdtyp * init(const int nlocal, const int nall, FILE *screen, const int order,
                const int nxlo_out, const int nylo_out, const int nzlo_out,
                const int nxhi_out, const int nyhi_out, const int nzhi_out,
                grdtyp **rho_coeff, grdtyp **vd_brick,
                const double slab_volfactor, const int nx_pppm,
                const int ny_pppm, const int nz_pppm, const bool split,
                int &success);

  /// Check if there is enough storage for atom arrays and realloc if not
  /** \param success set to false if insufficient memory **/
  inline void resize_atom(const int inum, const int nall, bool &success) {
    if (atom->resize(nall, success)) {
      pos_tex.bind_float(atom->x,4);
      q_tex.bind_float(atom->q,1);
    }
    ans->resize(inum,success);
  }

  /// Check if there is enough storage for local atoms and realloc if not
  inline void resize_local(const int, bool &) {
  }

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear(const double cpu_time);

  /// Returns memory usage on device per atom
  int bytes_per_atom() const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  /// Accumulate timers
  inline void acc_timers() {
    if (device->time_device()) {
      ans->acc_timers();
      time_in.add_to_total();
      time_out.add_to_total();
      time_map.add_to_total();
      time_rho.add_to_total();
      time_interp.add_to_total();
    }
  }

  /// Zero timers
  inline void zero_timers() {
    atom->zero_timers();
    ans->zero_timers();
    time_in.zero();
    time_out.zero();
    time_map.zero();
    time_rho.zero();
    time_interp.zero();
  }

  /// Precomputations for charge assignment that can be done asynchronously
  inline void precompute(const int ago, const int nlocal, const int nall,
                         double **host_x, int *host_type, bool &success,
                         double *charge, double *boxlo, double *prd) {
    double delxinv=_nx_pppm/prd[0];
    double delyinv=_ny_pppm/prd[1];
    double delzinv=_nz_pppm/(prd[2]*_slab_volfactor);
    _precompute(ago,nlocal,nall,host_x,host_type,success,charge,boxlo,delxinv,
                delyinv,delzinv);
  }

  /// Returns non-zero if out of bounds atoms
  int spread(const int ago, const int nlocal, const int nall, double **host_x,
             int *host_type, bool &success, double *charge, double *boxlo,
             const double delxinv, const double delyinv, const double delzinv);

  void interp(const grdtyp qqrd2e_scale);

  // -------------------------- DEVICE DATA -------------------------

  /// Device Properties and Atom and Neighbor storage
  Device<numtyp,acctyp> *device;

  /// Geryon device
  UCL_Device *ucl_device;

  /// Device Timers
  UCL_Timer time_in, time_out, time_map, time_rho, time_interp;

  /// LAMMPS pointer for screen output
  FILE *screen;

  // --------------------------- ATOM DATA --------------------------

  /// Atom Data
  Atom<numtyp,acctyp> *atom;


  // --------------------------- GRID DATA --------------------------

  UCL_Vector<grdtyp,grdtyp> brick;
  UCL_Vector<grdtyp,grdtyp> vd_brick;

  // Count of number of atoms assigned to each grid point
  UCL_D_Vec<int> d_brick_counts;
  // Atoms assigned to each grid point
  UCL_D_Vec<grdtyp4> d_brick_atoms;

  // Error checking for out of bounds atoms
  UCL_Vector<int,int> error_flag;

  // Number of grid points in brick (including ghost)
  int _npts_x, _npts_y, _npts_z, _npts_yx;

  // Number of local grid points in brick
  int _nlocal_x, _nlocal_y, _nlocal_z, _nlocal_yx, _atom_stride;

  // -------------------------- SPLINE DATA -------------------------
  UCL_D_Vec<grdtyp> d_rho_coeff;
  int _order, _nlower, _nupper, _order_m_1, _order2;
  int _nxlo_out, _nylo_out, _nzlo_out, _nxhi_out, _nyhi_out, _nzhi_out;

  // ------------------------ FORCE/ENERGY DATA -----------------------

  Answer<numtyp,acctyp> *ans;

  // ------------------------- DEVICE KERNELS -------------------------
  UCL_Program *pppm_program;
  UCL_Kernel k_particle_map, k_make_rho, k_interp;
  inline int block_size() { return _block_size; }

  // --------------------------- TEXTURES -----------------------------
  UCL_Texture pos_tex;
  UCL_Texture q_tex;

 protected:
  bool _allocated, _compiled, _precompute_done, _kspace_split;
  int _block_size, _block_pencils, _pencil_size, _max_brick_atoms, _max_atoms;
  double  _max_bytes, _max_an_bytes;
  double _cpu_idle_time;

  grdtyp _brick_x, _brick_y, _brick_z, _delxinv, _delyinv, _delzinv;

  double _slab_volfactor;
  int _nx_pppm, _ny_pppm, _nz_pppm;

  void compile_kernels(UCL_Device &dev);
  void _precompute(const int ago, const int nlocal, const int nall,
                   double **host_x, int *host_type, bool &success,
                   double *charge, double *boxlo, const double delxinv,
                   const double delyinv, const double delzinv);
};

}

#endif
