/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_FFT3D_WRAP_H
#define LMP_FFT3D_WRAP_H

#include "fft3d.h"    // IWYU pragma: export
#include "pointers.h"

namespace LAMMPS_NS {

class FFT3d : protected Pointers {
 public:
  enum { FORWARD = 1, BACKWARD = -1 };

  FFT3d(class LAMMPS *, MPI_Comm, int, int, int, int, int, int, int, int, int, int, int, int, int,
        int, int, int, int, int *, int);
  ~FFT3d() override;
  void compute(FFT_SCALAR *, FFT_SCALAR *, int);
  void timing1d(FFT_SCALAR *, int, int);

 private:
  struct fft_plan_3d *plan;
};

}    // namespace LAMMPS_NS

#endif
