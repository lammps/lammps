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

#include "fft3d.h"
#include "pointers.h"

#ifdef HEFFTE
#include "heffte.h"
#ifdef FFT_FFTW3
using heffte_backend = heffte::backend::fftw;
#else
#ifdef FFT_MKL
using heffte_backend = heffte::backend::mkl;
#endif
#endif
#endif

namespace LAMMPS_NS {

class FFT3d : protected Pointers {
 public:
  enum { FORWARD = 1, BACKWARD = -1 };

  FFT3d(class LAMMPS *, MPI_Comm, int, int, int, int, int, int, int, int, int, int, int, int, int,
        int, int, int, int, int *, int);
  ~FFT3d();
  void compute(FFT_SCALAR *, FFT_SCALAR *, int);
  void timing1d(FFT_SCALAR *, int, int);

 private:
  struct fft_plan_3d *plan;

  #ifdef HEFFTE
  std::unique_ptr<heffte::fft3d<heffte_backend>> heffte_plan;
  std::vector<std::complex<FFT_SCALAR>> heffte_workspace;
  #endif

};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Could not create 3d FFT plan

The FFT setup for the PPPM solver failed, typically due
to lack of memory.  This is an unusual error.  Check the
size of the FFT grid you are requesting.

*/
