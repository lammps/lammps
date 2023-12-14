/***************************************************************************
                                  neighbor.h
                             -------------------
                            Nitin Dhamankar (Intel)
                            W. Michael Brown (ORNL)
                              Peng Wang (Nvidia)

  Class for handling neighbor lists

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov, penwang@nvidia.com
 ***************************************************************************/

#ifndef LAL_NEIGHBOR_H
#define LAL_NEIGHBOR_H

#include "lal_atom.h"
#include "lal_neighbor_shared.h"
#include <sstream>

#define IJ_SIZE 131072

#if !defined(USE_OPENCL) && !defined(USE_HIP)
#ifndef LAL_USE_OLD_NEIGHBOR
// Issue with incorrect results with CUDA >= 11.2 and pre-12.0
#if (CUDA_VERSION > 11019) && (CUDA_VERSION < 12000)
#define LAL_USE_OLD_NEIGHBOR
#endif
#endif
#endif

#if defined(USE_HIP) || defined(__APPLE__)
#define LAL_USE_OLD_NEIGHBOR
#endif

namespace LAMMPS_AL {

class Neighbor {
 public:
  Neighbor() : _allocated(false), _use_packing(false), _old_max_nbors(0), _ncells(0) {}
  ~Neighbor() { clear(); }

  /// Determine whether neighbor unpacking should be used
  /** If false, twice as much memory is reserved to allow unpacking neighbors by
    * atom for coalesced access. **/
  void packing(const bool use_packing) { _use_packing=use_packing; }

  /// Clear any old data and setup for new LAMMPS run
  /** \param inum Initial number of particles whose neighbors stored on device
    * \param host_inum Initial number of particles whose nbors copied to host
    * \param max_nbors Factor (in percentage) applied to density calculated max
    * \param gpu_nbor 0 if neighboring will be performed on host
    *        gpu_nbor 1 if neighboring will be performed on device
    *        gpu_nbor 2 if binning on host and neighboring on device
    * \param gpu_host 0 if host will not perform force calculations,
    *                 1 if gpu_nbor is true, and host needs a half nbor list,
    *                 2 if gpu_nbor is true, and host needs a full nbor list
    * \param pre_cut True if cutoff test will be performed in separate kernel
    *                than the force kernel
    * \param threads_per_atom Number of threads used per atom for force
    *                         calculation
    * \param compile_flags Flags for JIT compiling
    * \param ilist_map true if ilist mapping data structures used (3-body) **/
  bool init(NeighborShared *shared, const int inum, const int host_inum,
            const int max_nbors, const int maxspecial, UCL_Device &dev,
            const int gpu_nbor, const int gpu_host, const bool pre_cut,
            const int block_cell_2d, const int block_cell_id,
            const int block_nbor_build, const int threads_per_atom,
            const int simd_size, const bool time_device,
            const std::string &compile_flags, const bool ilist_map);

  /// Set the cutoff+skin
  inline void set_cutoff(const double cutoff) {
    _cutoff=cutoff;

    #ifndef LAL_USE_OLD_NEIGHBOR
    _cell_size=_shared->cell_size();
    _auto_cell_size=_shared->auto_cell_size();
    const int cells_in_cutoff=static_cast<int>(ceil(_cutoff/_cell_size));
    if (cells_in_cutoff > 2) _cell_size=_cutoff*0.5;
    _old_ncellx = _old_ncelly = _old_ncellz = -1;
    #else
    _cell_size=cutoff;
    _auto_cell_size=false;
    #endif
  }

  /// Get the cutoff+skin
  inline double cutoff() { return _cutoff; }

  /// Check if there is enough memory for neighbor data and realloc if not
  /** \param inum Number of particles whose nbors will be stored on device
    * \param max_nbor Current max number of neighbors for a particle
    * \param success False if insufficient memory **/
  inline void resize(const int inum, int max_nbor, bool &success) {
    if (max_nbor == 0) max_nbor = 1;
    if (inum>_max_atoms || max_nbor>_max_nbors) {
      _max_atoms=static_cast<int>(static_cast<double>(inum)*1.10);
      if (max_nbor>_max_nbors)
        _max_nbors=static_cast<int>(static_cast<double>(max_nbor)*1.10);
      alloc(success);
    }
  }

  /// Check if there is enough memory for neighbor data and realloc if not
  /** \param inum Number of particles whose nbors will be stored on device
    * \param host_inum Number of particles whose nbors will be copied to host
    * \param max_nbor Current max number of neighbors for a particle
    * \param success False if insufficient memory **/
  inline void resize(const int inum, const int host_inum, int max_nbor,
                     bool &success) {
    if (max_nbor == 0) max_nbor = 1;
    if (inum>_max_atoms || max_nbor>_max_nbors || host_inum>_max_host) {
      _max_atoms=static_cast<int>(static_cast<double>(inum)*1.10);
      _max_host=static_cast<int>(static_cast<double>(host_inum)*1.10);
      if (max_nbor>_max_nbors)
        _max_nbors=static_cast<int>(static_cast<double>(max_nbor)*1.10);
      alloc(success);
    }
  }

  inline void acc_timers(FILE *) {
    if (_nbor_time_avail) {
      if (_time_device) {
        time_nbor.add_to_total();
        if (_use_packing==false) time_kernel.add_to_total();
        if (_gpu_nbor==2) {
          time_hybrid1.add_to_total();
          time_hybrid2.add_to_total();
        }
        if (_maxspecial>0)
          time_transpose.add_to_total();
        _nbor_time_avail=false;
      }
    }
  }

  /// Free all memory on host and device
  void clear();

  /// Bytes per atom used on device
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by class
  double host_memory_usage() const;

  /// Returns the type of neighboring:
  /** - 0 if neighboring will be performed on host
    * - 1 if neighboring will be performed on device
    * - 2 if binning on host and neighboring on device **/
  inline int gpu_nbor() const { return _gpu_nbor; }

  /// Make a copy of unpacked nbor lists in the packed storage area (for gb)
  inline void copy_unpacked(const int inum, const int maxj)
    { ucl_copy(dev_packed,dev_nbor,inum*(maxj+2),true); }

  /// Copy neighbor list from host (first time or from a rebuild)
  void get_host(const int inum, int *ilist, int *numj,
                int **firstneigh, const int block_size);

  /// Copy neighbor list from host for 3-body (first time or from a rebuild)
  void get_host3(const int inum, const int nlist, int *ilist, int *numj,
                 int **firstneigh, const int block_size);

  /// Return the stride in elements for each nbor row
  inline int nbor_pitch() const { return _nbor_pitch; }

  /// Return the maximum number of atoms that can currently be stored
  inline int max_atoms() const { return _max_atoms; }

  /// Return the maximum number of nbors for a particle based on current alloc
  inline int max_nbors() const { return _max_nbors; }

  /// Return the time spent binning on the CPU for hybrid neighbor builds
  inline double bin_time() const { return _bin_time; }

  /// Loop through neighbor count array and return maximum nbors for a particle
  inline int max_nbor_loop(const int inum, int *numj, int *ilist) const {
    int mn=0;
    for (int i=0; i<inum; i++)
      mn=std::max(mn,numj[ilist[i]]);
    return mn;
  }

  /// Build nbor list on the device
  template <class numtyp, class acctyp>
  void build_nbor_list(double **x, const int inum, const int host_inum,
                       const int nall, Atom<numtyp,acctyp> &atom,
                       double *sublo, double *subhi, tagint *tag,
                       int **nspecial, tagint **special, bool &success,
                       int &max_nbors, UCL_Vector<int,int> &error_flag);

  /// Return the number of bytes used on device
  inline double gpu_bytes() {
    double res = _gpu_bytes + _c_bytes + _cell_bytes;
    if (_gpu_nbor==0)
      res += 2*IJ_SIZE*sizeof(int);

    return res;
  }

  // ------------------------------- Data -------------------------------

  /// Device neighbor matrix
  /** - 1st row is i (index into atom data)
    * - 2nd row is numj (number of neighbors)
    * - 3rd row is starting location in packed nbors
    * - Remaining rows are the neighbors arranged for coalesced access **/
  UCL_D_Vec<int> dev_nbor;
  /// Starting location in packed neighbors used only by unpack kernel
  UCL_D_Vec<int> dev_packed_begin;
  /// Packed storage for neighbor lists copied from host
  UCL_D_Vec<int> dev_packed;
  /// Host buffer for copying neighbor lists
  UCL_H_Vec<int> host_packed;
  /// Host storage for nbor counts (row 1) & accumulated neighbor counts (row2)
  UCL_H_Vec<int> host_acc;
  /// Storage for accessing atom indices from the neighbor list (3-body)
  UCL_Vector<int,int> three_ilist;

  // ----------------- Data for GPU Neighbor Calculation ---------------

  /// Host/Device storage for device calculated neighbor lists
  /** - 1st row is numj
    * - Remaining rows are by atom, columns are nbors **/
  UCL_Vector<int,int> nbor_host;
  UCL_D_Vec<int> dev_numj_host;
  UCL_H_Vec<int> host_ilist;
  UCL_H_Vec<int*> host_jlist;
  /// Device storage for special neighbor counts
  UCL_D_Vec<int> dev_nspecial;
  /// Device storage for special neighbors
  UCL_D_Vec<tagint> dev_special, dev_special_t;
  /// Host/Device storage for number of particles per cell
  UCL_Vector<int,int> cell_counts;
  #ifndef LAL_USE_OLD_NEIGHBOR
  /// Host/Device storage for number of subgroups per cell
  UCL_Vector<int,int> cell_subgroup_counts;
  /// Host/Device storage for subgroup to cell mapping
  UCL_Vector<int,int> subgroup2cell;
  #endif
  int *cell_iter;

  /// Device timers
  UCL_Timer time_nbor, time_kernel, time_hybrid1, time_hybrid2, time_transpose;

  /// Effective SIMD width of neighbor build kernel
  inline int simd_size() { return _simd_size; }

  template <class t>
    inline std::string toa(const t& in) {
    std::ostringstream o;
    o.precision(2);
    o << in;
    return o.str();
  }

  /// Helper function
  void transpose(UCL_D_Vec<tagint> &out, const UCL_D_Vec<tagint> &in,
    const int columns_in, const int rows_in);

 private:
  NeighborShared *_shared;
  UCL_Device *dev;
  bool _allocated, _use_packing, _nbor_time_avail, _time_device;
  int _gpu_nbor, _max_atoms, _max_nbors, _max_host, _nbor_pitch, _maxspecial;
  int _old_max_nbors;
  bool _gpu_host, _alloc_packed, _ilist_map, _auto_cell_size;
  double _cutoff, _bin_time, _max_neighbor_factor, _cell_size;
  enum UCL_MEMOPT _packed_permissions;

  double _gpu_bytes, _c_bytes, _cell_bytes;
  void alloc(bool &success);

  int _block_cell_2d, _block_cell_id, _max_block_nbor_build, _block_nbor_build;
  int _ncells, _threads_per_atom, _total_atoms;

  template <class numtyp, class acctyp>
  inline void resize_max_neighbors(int maxn, bool &success);

  // For viewing host arrays for data copy operations
  UCL_H_Vec<int> _host_offset;
  UCL_D_Vec<int> _nbor_offset, _acc_view, _numj_view;

  #ifndef LAL_USE_OLD_NEIGHBOR
  UCL_H_Vec<int> _host_bin_stencil;
  UCL_Const _bin_stencil;
  int _old_ncellx, _old_ncelly, _old_ncellz;
  #endif

  int _simd_size;
  #ifdef LAL_USE_OLD_NEIGHBOR
  inline void set_nbor_block_size(const int mn) {
    int desired=mn/(2*_simd_size);
    desired*=_simd_size;
    if (desired<_simd_size) desired=_simd_size;
    else if (desired>_max_block_nbor_build) desired=_max_block_nbor_build;
    _block_nbor_build=desired;
  }
  #else
  inline void set_nbor_block_size(const int) {}
  #endif
};

}

#endif
