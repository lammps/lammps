/***************************************************************************
                              born_coul_wolf_cs.h
                             -------------------
                           Trung Dac Nguyen (Northwestern)

  Class for acceleration of the born/coul/wolf/cs pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : ndactrung@gmail.com
 ***************************************************************************/

#ifndef LAL_BORN_COUL_WOLF_CS_H
#define LAL_BORN_COUL_WOLF_CS_H

#include "lal_born_coul_wolf.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class BornCoulWolfCS : public BornCoulWolf<numtyp, acctyp> {
 public:
  BornCoulWolfCS() {}
  ~BornCoulWolfCS() {}

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
  int init(const int ntypes, double **host_cutsq, double **host_rhoinv,
           double **host_born1, double **host_born2, double **host_born3,
           double **host_a, double **host_c, double **host_d,
           double **host_sigma, double **host_offset, double *host_special_lj,
           const int nlocal, const int nall, const int max_nbors,
           const int maxspecial, const double cell_size,
           const double gpu_split, FILE *screen, double **host_cut_ljsq,
           const double host_cut_coulsq, double *host_special_coul,
           const double qqrd2e, const double alf, const double e_shift,
           const double f_shift);
};

}

#endif
