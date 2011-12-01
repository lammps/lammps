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

#include "lal_device.h"
#include "lal_balance.h"
#include "mpi.h"

#ifdef USE_OPENCL
#include "geryon/ocl_texture.h"
#else
#include "geryon/nvd_texture.h"
#endif

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class EAM {
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
  int init_atomic(const int nlocal, const int nall, const int max_nbors,
                  const int maxspecial, const double cell_size,
                  const double gpu_split, FILE *screen,
                  const char *pair_program);
                  
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
  
  /// Estimate the overhead for GPU context changes and CPU driver
  void estimate_gpu_overhead();
  
  /// Check if there is enough storage for atom arrays and realloc if not
  /** \param success set to false if insufficient memory **/
  inline void resize_atom(const int inum, const int nall, bool &success) {
    if (atom->resize(nall, success)) {
      pos_tex.bind_float(atom->dev_x,4);
      q_tex.bind_float(atom->dev_q,1);
    }
    
    // resize rho
    dev_fp.clear();
    host_fp.clear();
    
    bool cpuview=false;
    if (ucl_device->device_type()==UCL_CPU)
      cpuview=true;
    
    int _max_atoms=static_cast<int>(static_cast<double>(nall)*1.10);
    host_fp.alloc(_max_atoms,*ucl_device);
    if (cpuview)
      dev_fp.view(host_fp);
    else 
      dev_fp.alloc(_max_atoms,*ucl_device,UCL_WRITE_ONLY);
    
    ans->resize(inum,success);
  }
  
  /// Check if there is enough storage for neighbors and realloc if not
  /** \param nlocal number of particles whose nbors must be stored on device
    * \param host_inum number of particles whose nbors need to copied to host
    * \param current maximum number of neighbors
    * \note olist_size=total number of local particles **/
  inline void resize_local(const int inum, const int max_nbors, bool &success) {
    nbor->resize(inum,max_nbors,success);
  }

  /// Check if there is enough storage for neighbors and realloc if not
  /** \param nlocal number of particles whose nbors must be stored on device
    * \param host_inum number of particles whose nbors need to copied to host
    * \param current maximum number of neighbors
    * \note host_inum is 0 if the host is performing neighboring
    * \note nlocal+host_inum=total number local particles
    * \note olist_size=0 **/
  inline void resize_local(const int inum, const int host_inum, 
                           const int max_nbors, bool &success) {
    nbor->resize(inum,host_inum,max_nbors,success);
  }

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;
  
  /// Accumulate timers
  inline void acc_timers() {
    if (device->time_device()) {
      nbor->acc_timers();
      time_pair.add_to_total();
      atom->acc_timers();
      ans->acc_timers();
    }
  }

  /// Zero timers
  inline void zero_timers() {
    time_pair.zero();
    atom->zero_timers();
    ans->zero_timers();
  }

  /// Copy neighbor list from host
  int * reset_nbors(const int nall, const int inum, int *ilist, int *numj,
                    int **firstneigh, bool &success);

  /// Build neighbor list on device
  void build_nbor_list(const int inum, const int host_inum,
                       const int nall, double **host_x, int *host_type,
                       double *sublo, double *subhi, int *tag, int **nspecial,
                       int **special, bool &success);

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
                double *charge, double *boxlo, double *prd, int inum);
  
  /// Pair loop with host neighboring
  void compute_energy(const int f_ago, const int inum_full, const int nall,
               double **host_x, int *host_type, int *ilist, int *numj,
               int **firstneigh, const bool eflag, const bool vflag,
               const bool eatom, const bool vatom, int &host_start,
               const double cpu_time, bool &success, double *charge,
               const int nlocal, double *boxlo, double *prd, double* evdwl);
               
  /// Pair loop with device neighboring
  int** compute_energy(const int ago, const int inum_full, const int nall,
                double **host_x, int *host_type, double *sublo,
                double *subhi, int *tag, int **nspecial,
                int **special, const bool eflag, const bool vflag, 
                const bool eatom, const bool vatom, int &host_start, 
                int **ilist, int **numj, const double cpu_time, bool &success,
                double *charge, double *boxlo, double *prd, double* evdwl, int &inum);

// -------------------------- DEVICE DATA ------------------------- 

  /// Device Properties and Atom and Neighbor storage
  Device<numtyp,acctyp> *device;

  /// Geryon device
  UCL_Device *ucl_device;

  /// Device Timers
  UCL_Timer time_pair;

  /// Host device load balancer
  Balance<numtyp,acctyp> hd_balancer;

  /// LAMMPS pointer for screen output
  FILE *screen;

  // --------------------------- ATOM DATA --------------------------

  /// Atom Data
  Atom<numtyp,acctyp> *atom;


  // ------------------------ FORCE/ENERGY/DENSITY DATA --------------

  Answer<numtyp,acctyp> *ans;

  // --------------------------- NBOR DATA ----------------------------

  /// Neighbor data
  Neighbor *nbor;

  // ------------------------- DEVICE KERNELS -------------------------
  UCL_Program *pair_program;
  UCL_Kernel k_pair_fast, k_pair, k_energy;
  inline int block_size() { return _block_size; }

  // --------------------------- TEXTURES -----------------------------
  UCL_Texture pos_tex;
  UCL_Texture q_tex;

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

  /// Per-atom arrays
  UCL_H_Vec<acctyp> host_fp;
  UCL_D_Vec<acctyp> dev_fp;
  
protected:
  bool _compiled;
  int _block_size, _block_bio_size, _threads_per_atom;
  double  _max_bytes, _max_an_bytes;
  double _gpu_overhead, _driver_overhead;
  UCL_D_Vec<int> *_nbor_data;

  void compile_kernels(UCL_Device &dev, const char *pair_string);  

  bool _allocated;
  void loop(const bool _eflag, const bool _vflag);
  void energy(const bool _eflag, const bool _vflag);
};

}

#endif

