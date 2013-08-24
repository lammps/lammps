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
    * -  0 if successfull
    * - -1 if fix gpu not found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card **/
  int init(const int nlocal, const int nall, const int max_nbors, 
           const double cell_size, const double gpu_split, FILE *screen,
           const double epsilon, const double sigma,
           const double lambda, const double gamma,
           const double costheta, const double biga,
           const double bigb, const double powerp,
           const double powerq, const double cut, const double cutsq);
           
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

 private:
  bool _allocated;
  void loop(const bool _eflag, const bool _vflag, const int evatom);
  numtyp sw_epsilon, sw_sigma, sw_lambda, sw_gamma, sw_costheta;
  numtyp sw_biga, sw_bigb, sw_powerp, sw_powerq, sw_cut, sw_cutsq;
};

}

#endif

