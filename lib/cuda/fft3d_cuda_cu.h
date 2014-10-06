/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#include "cuda_shared.h"

extern "C" void initfftdata(double* in, FFT_CFLOAT* out, int nfast, int nmid, int nslow);
extern "C" void permute(FFT_DATA* in, FFT_DATA* out, int nfast, int nmid, int nslow);
extern "C" void permute_scale(FFT_DATA* in, FFT_DATA* out, int nfast, int nmid, int nslow);
extern "C" void permute_part(FFT_DATA* in, FFT_DATA* out, int nfast, int nmid, int nslow, int ihi, int ilo, int jhi, int jlo, int khi, int klo);
extern "C" void FFTsyncthreads();
