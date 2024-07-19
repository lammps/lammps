/***************************************************************************
                                 lal_dpd_coul_slater_long.h
                             -------------------
                            Eddy BARRAUD (IFPEN/Sorbonne)

  Class for acceleration of the dpd/coul/slater/long pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : May 28, 2024
    email                : eddy.barraud@outlook.fr
 ***************************************************************************/

#ifndef LAL_DPD_CHARGED_H
#define LAL_DPD_CHARGED_H

#include "lal_base_dpd.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class DPDCoulSlaterLong : public BaseDPD<numtyp, acctyp> {
 public:
  DPDCoulSlaterLong();
  ~DPDCoulSlaterLong();

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
  int init(const int ntypes, double **host_cutsq, double **host_a0, double **host_gamma,
           double **host_sigma, double **host_cut_dpd, double **host_cut_dpdsq,
           double **host_cut_slatersq, double *host_special_lj, bool tstat_only, const int nlocal,
           const int nall, const int max_nbors, const int maxspecial, const double cell_size,
           const double gpu_split, FILE *screen, double *host_special_coul, const double qqrd2e,
           const double g_ewald, const double lamda);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  /// Update coeff if needed (tstat only)
  void update_coeff(int ntypes, double **host_a0, double **host_gamma,
                    double **host_sigma, double **host_cut_dpd );

  void get_extra_data(double *host_q);

  // --------------------------- TYPE DATA --------------------------

  /// coeff.x = a0, coeff.y = gamma, coeff.z = sigma, coeff.w = cut_dpd
  UCL_D_Vec<numtyp4> coeff;

  /// cutsq.x = cutsq, cutsq.y = cut_dpdsq, cutsq.w = cut_slatersq
  UCL_D_Vec<numtyp4> cutsq;

  /// Special LJ values
  UCL_D_Vec<numtyp> sp_lj, sp_sqrt;

  /// Special Coul values [0-3]
  UCL_D_Vec<numtyp> sp_cl;


  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _lj_types;

  /// Only used for thermostat
  int _tstat_only;

  /// Coulombic terms
  numtyp _qqrd2e, _g_ewald, _lamda;

  /// pointer to host data for atom charge
  double *q;

 private:
  bool _allocated;
  int loop(const int eflag, const int vflag);
};

}

#endif
