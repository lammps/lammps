/***************************************************************************
                                 mdpd.h
                             -------------------
                            Trung Dac Nguyen (U Chicago)

  Class for acceleration of the mdpd pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : December 2023
    email                : ndactrung@gmail.com
 ***************************************************************************/

#ifndef LAL_MDPD_H
#define LAL_MDPD_H

#include "lal_base_dpd.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class MDPD : public BaseDPD<numtyp, acctyp> {
 public:
  MDPD();
  ~MDPD();

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
           double **host_A_att, double **host_B_rep,
           double **host_gamma, double **host_sigma,
           double **host_cut, double **host_cut_r, double *host_special_lj,
           const int nlocal, const int nall, const int max_nbors,
           const int maxspecial, const double cell_size, const double gpu_split,
           FILE *screen);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  void get_extra_data(double *host_rho);

  // --------------------------- TYPE DATA --------------------------

  /// coeff.x = A_att, coeff.x = B_rep, coeff.z = gamma, coeff.w = sigma
  UCL_D_Vec<numtyp4> coeff;
  /// coeff2.x = cut, coeff2.y = cut_r, coeff2.z = cutsq
  UCL_D_Vec<numtyp4> coeff2;

  UCL_D_Vec<numtyp> cutsq;

  /// Special LJ values
  UCL_D_Vec<numtyp> sp_lj, sp_sqrt;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _lj_types;

  /// pointer to host data
  double *mdpd_rho;

 private:
  bool _allocated;
  int loop(const int eflag, const int vflag);
};

}

#endif
