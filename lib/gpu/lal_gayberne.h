/***************************************************************************
                                 gayberne.h
                             -------------------
                            W. Michael Brown (ORNL)

  Host code for Gay-Berne potential acceleration

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov
 ***************************************************************************/

#ifndef LAL_GAYBERNE_H
#define LAL_GAYBERNE_H

#include "lal_base_ellipsoid.h"
#include "mpi.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class GayBerne : public BaseEllipsoid<numtyp, acctyp> {
 public:
  GayBerne();
  ~GayBerne();

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param max_nbors initial number of rows in the neighbor matrix
    * \param cell_size cutoff + skin
    * \param gpu_split fraction of particles handled by device
    * \return false if there is not sufficient memory or device init prob
    *
    * Returns:
    * -  0 if successful
    * - -1 if fix gpu not found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card **/
  int init(const int ntypes, const double gamma,
           const double upsilon, const double mu, double **host_shape,
           double **host_well, double **host_cutsq, double **host_sigma,
           double **host_epsilon, double *host_lshape, int **h_form,
           double **host_lj1, double **host_lj2, double **host_lj3,
           double **host_lj4, double **host_offset,
           const double *host_special_lj, const int nlocal, const int nall,
           const int max_nbors, const int maxspecial, const double cell_size,
           const double gpu_split, FILE *screen);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  /// Device Error Flag - Set if a bad matrix inversion occurs
  UCL_D_Vec<int> dev_error;

  // --------------------------- TYPE DATA --------------------------

  /// lj1.x = lj1, lj1.y = lj2, lj1.z = cutsq, lj1.w = form
  UCL_D_Vec<numtyp4> lj1;
  /// lj3.x = lj3, lj3.y = lj4, lj3.z = offset
  UCL_D_Vec<numtyp4> lj3;
  /// sigma_epsilon.x = sigma, sigma_epsilon.y = epsilon
  UCL_D_Vec<numtyp2> sigma_epsilon;
  // 0 - gamma, 1-upsilon, 2-mu, 3-special_lj[0], 4-special_lj[1], ...
  UCL_D_Vec<numtyp> gamma_upsilon_mu;

  /// If atom type constants fit in shared memory, use fast kernels
  bool _shared_types;
  int _lj_types;

  // --------------------------- ATOM DATA --------------------------

  /// Aspherical Const Data for Atoms
  UCL_D_Vec<numtyp4> shape, well;
  /// Aspherical Const Data for Atoms
  UCL_D_Vec<numtyp> lshape;

 private:
  bool _allocated;
  int loop(const int eflag, const int vflag);
};

}

#endif
