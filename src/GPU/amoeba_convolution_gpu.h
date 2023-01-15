/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_AMOEBA_CONVOLUTION_GPU_H
#define LMP_AMOEBA_CONVOLUTION_GPU_H

#include "amoeba_convolution.h"


namespace LAMMPS_NS {

class AmoebaConvolutionGPU : public AmoebaConvolution {
 public:
  AmoebaConvolutionGPU(class LAMMPS *, class Pair *, int, int, int, int, int);

  FFT_SCALAR *pre_convolution_4d() override;
  void *post_convolution_4d() override;

};

} // namespace LAMMPS_NS
#endif
