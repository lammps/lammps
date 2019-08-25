/***************************************************************************
                                vashishta.h
                             -------------------
                            Anders Hafreager (UiO9)

  Class for acceleration of the vashishta pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Mon June 12, 2017
    email                : andershaf@gmail.com
 ***************************************************************************/

#ifndef LAL_VASHISHTA_H
#define LAL_VASHISHTA_H

#include "lal_base_three.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class Vashishta : public BaseThree<numtyp, acctyp> {
 public:
  Vashishta();
  ~Vashishta();

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
  int init(const int ntypes, const int nlocal, const int nall, const int max_nbors,
           const double cell_size, const double gpu_split, FILE *screen,
           int* host_map, const int nelements, int*** host_elem2param, const int nparams,
           const double* cutsq, const double* r0, 
           const double* gamma, const double* eta,
           const double* lam1inv, const double* lam4inv,
           const double* zizj, const double* mbigd,
           const double* dvrc, const double* big6w, 
           const double* heta, const double* bigh,
           const double* bigw, const double* c0,
           const double* costheta, const double* bigb,
           const double* big2b, const double* bigc);

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

  /// param1.x = eta, param1.y = lam1inv, param1.z = lam4inv, param1.w = zizj
  UCL_D_Vec<numtyp4> param1;
  /// param2.x = mbigd, param2.y = dvrc, param2.z = big6w, param2.w = heta
  UCL_D_Vec<numtyp4> param2;
  /// param3.x = bigh, param3.y = bigw, param3.z = dvrc, param3.w = c0
  UCL_D_Vec<numtyp4> param3;
  /// param4.x = r0sq, param4.y = gamma, param4.z = cutsq, param4.w = r0
  UCL_D_Vec<numtyp4> param4;
  /// param5.x = bigc, param5.y = costheta, param5.z = bigb, param5.w = big2b
  UCL_D_Vec<numtyp4> param5;

  UCL_D_Vec<int> elem2param;
  UCL_D_Vec<int> map;
  int _nparams,_nelements;
  numtyp _cutshortsq;

  UCL_Texture param1_tex, param2_tex, param3_tex, param4_tex, param5_tex;

 private:
  bool _allocated;
  void loop(const bool _eflag, const bool _vflag, const int evatom);

};

}

#endif

