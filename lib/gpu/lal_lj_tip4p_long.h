/**************************************************************************
                               lj_tip4p_long.h
                             -------------------
                              V. Nikolskiy (HSE)

  Class for acceleration of the lj/tip4p/long pair style

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : thevsevak@gmail.com
***************************************************************************/

#ifndef LAL_LJ_TIP4P_LONG_H
#define LAL_LJ_TIP4P_LONG_H

#include "lal_base_charge.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class LJ_TIP4PLong : public BaseCharge<numtyp, acctyp> {
public:
  LJ_TIP4PLong();
  ~LJ_TIP4PLong();

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
           double **host_lj1, double **host_lj2, double **host_lj3,
           double **host_lj4, double **host_offset, double *host_special_lj,
           const int nlocal, const int tH, const int tO,
           const double alpha, const double qdist,
           const int nall, const int max_nbors,
           const int maxspecial, const double cell_size,
           const double gpu_split, FILE *screen,
           double **host_cut_ljsq,
           const double host_cut_coulsq, const double host_cut_coulsqplus,
           double *host_special_coul, const double qqrd2e,
           const double g_ewald, int* tag,
           int *map_array, int map_size,
           int *sametag, int max_same);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  /// Copy data from LAMMPS_NS
  void copy_relations_data(int n,int* tag, int *map_array, int map_size,
      int *sametag, int max_same, int ago);

  /// Reimplement BaseCharge pair loop with host neighboring
  void compute(const int f_ago, const int inum_full, const int nall,
               double **host_x, int *host_type, int *ilist, int *numj,
               int **firstneigh, const bool eflag, const bool vflag,
               const bool eatom, const bool vatom, int &host_start,
               const double cpu_time, bool &success, double *charge,
               const int nlocal, double *boxlo, double *prd);

  /// Reimplement BaseCharge pair loop with device neighboring
  int** compute(const int ago, const int inum_full, const int nall,
                double **host_x, int *host_type, double *sublo,
                double *subhi, tagint *tag,int *map_array, int map_size, int *sametag, int max_same,
                int **nspecial,
                tagint **special, const bool eflag, const bool vflag,
                const bool eatom, const bool vatom, int &host_start,
                int **ilist, int **numj, const double cpu_time, bool &success,
                double *charge, double *boxlo, double *prd);


  // --------------------------- TYPE DATA --------------------------

  /// lj1.x = lj1, lj1.y = lj2, lj1.z = cutsq_vdw
  UCL_D_Vec<numtyp4> lj1;
  /// lj3.x = lj3, lj3.y = lj4, lj3.z = offset
  UCL_D_Vec<numtyp4> lj3;
  /// cutsq
  UCL_D_Vec<numtyp> cutsq;
  /// Special LJ values [0-3] and Special Coul values [4-7]
  UCL_D_Vec<numtyp> sp_lj;

  bool shared_types;

  /// Number of atom types
  int _lj_types;

  numtyp _qqrd2e, _g_ewald;
  /// TIP4P water parameters
  int TypeO, TypeH;
  numtyp alpha, qdist;
  numtyp cut_coulsq, cut_coulsqplus;

  UCL_D_Vec<int> hneight;
  UCL_D_Vec<numtyp4> m; // position and charge of virtual particle
  UCL_D_Vec<acctyp4> ansO; // force applied to virtual particle
  // UCL_D_Vec<acctyp4> force_comp;

  UCL_D_Vec<int> tag;
  UCL_D_Vec<int> map_array;
  UCL_D_Vec<int> atom_sametag;

  UCL_Kernel k_pair_distrib, k_pair_reneigh;

 private:
  bool _allocated;
  int t_ago;
  void loop(const bool _eflag, const bool _vflag);
};

}

#endif
