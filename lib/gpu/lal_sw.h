/***************************************************************************
                                    sw.h
                             -------------------
                            W. Michael Brown (ORNL)

  Class for acceleration of the sw pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Tue March 26, 2013
    email                : brownw@ornl.gov
 ***************************************************************************/

#ifndef LAL_SW_H
#define LAL_SW_H

#include "lal_base_three.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class SW : public BaseThree<numtyp, acctyp> {
 public:
  SW();
  ~SW();

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
  int init(const int ntypes, const int nlocal, const int nall,
           const int max_nbors, const double cell_size,
           const double gpu_split, FILE *screen, double **ncutsq,
           double **ncut, double **sigma, double **powerp, double **powerq,
           double **sigma_gamma, double **c1, double **c2, double **c3,
           double **c4, double **c5, double **c6, double ***lambda_epsilon,
           double ***costheta, const int *map, int ***e2param);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _lj_types;

  UCL_D_Vec<numtyp> cutsq;
  /// sw_pre.x = cut, sw_pre.y = sigma, sw_pre.z = powerp, sw_pre.w = powerq
  UCL_D_Vec<numtyp4> sw_pre;
  /// c_14.x = c1, c_14.y = c2, c_14.z = c3, c_14.w = c4
  UCL_D_Vec<numtyp4> c_14;
  /// c_56.x = c5, c_56.y = c6
  UCL_D_Vec<numtyp2> c_56;
  /// cut_sigma_gamma.x = cut, cut_sigma_gamma.y = sigma_gamma
  UCL_D_Vec<numtyp2> cut_sigma_gamma;
  /// sw_pre3.x = lambda_epsilon, sw_pre3.y = costheta
  UCL_D_Vec<numtyp2> sw_pre3;

 private:
  bool _allocated;
  int loop(const int eflag, const int vflag, const int evatom, bool &success);

};

}

#endif

