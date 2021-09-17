/***************************************************************************
                                 tersoff.h
                             -------------------
                              Trung Dac Nguyen

  Class for acceleration of the tersoff pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Thu April 17, 2014
    email                : ndactrung@gmail.com
 ***************************************************************************/

#ifndef LAL_TERSOFF_H
#define LAL_TERSOFF_H

#include "lal_base_three.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class Tersoff : public BaseThree<numtyp, acctyp> {
 public:
  Tersoff();
  ~Tersoff();

  /// Clear any previous data and set up for a new LAMMPS run for generic systems
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
  int init(const int ntypes, const int nlocal, const int nall, const int max_nbors,
           const double cell_size, const double gpu_split, FILE *screen,
           int* host_map, const int nelements, int*** host_elem2param, const int nparams,
           const double* lam1, const double* lam2, const double* lam3,
           const double* powermint, const double* biga, const double* bigb,
           const double* bigr, const double* bigd, const double* c1, const double* c2,
           const double* c3, const double* c4, const double* c, const double* d,
           const double* h, const double* gamma, const double* beta,
           const double* powern, const double* cutsq);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// Number of atom types
  int _ntypes;

  /// ts1.x = lam3, ts1.y = powermint,  ts1.z = c3, ts1.w = c4
  UCL_D_Vec<numtyp4> ts1;
  /// ts2.x = biga, ts2.y = lam1,  ts2.z = bigr, ts2.w = bigd
  UCL_D_Vec<numtyp4> ts2;
  /// ts3.x = c1,   ts3.y = c2,    ts3.z = c3,   ts3.w = c4
  UCL_D_Vec<numtyp4> ts3;
  /// ts4.x = c*c,  ts4.y = d*d,   ts4.z = h,    ts4.w = gamma
  UCL_D_Vec<numtyp4> ts4;
  /// ts5.x = beta, ts5.y = powern, ts5.z = lam2, ts5.w = bigb
  UCL_D_Vec<numtyp4> ts5;

  numtyp _cutsq_max;

  UCL_D_Vec<int> elem2param;
  UCL_D_Vec<int> map;
  int _nparams,_nelements;

  /// Per-atom arrays:
  /// zetaij.x = force, zetaij.y = prefactor
  UCL_D_Vec<acctyp2>   _zetaij;
  UCL_D_Vec<acctyp> _zetaij_eng;

  UCL_Kernel k_zeta, k_zeta_noev, *k_zeta_selt;

 private:
  bool _allocated;
  int loop(const int eflag, const int vflag, const int evatom, bool &success);
};

}

#endif

