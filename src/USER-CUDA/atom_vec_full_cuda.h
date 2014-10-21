/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef ATOM_CLASS

AtomStyle(full/cuda,AtomVecFullCuda)

#else

#ifndef LMP_ATOM_VEC_FULL_CUDA_H
#define LMP_ATOM_VEC_FULL_CUDA_H

#include "atom_vec_full.h"
#include "cuda_data.h"

namespace LAMMPS_NS {

class AtomVecFullCuda : public AtomVecFull {
 public:
  AtomVecFullCuda(class LAMMPS *);
  virtual ~AtomVecFullCuda() {}
  void grow_copylist(int n);
  void grow_send(int n,double** buf_send,int flag);
  void grow_both(int n);
  int pack_comm(int, int *, double *, int, int *);
  int pack_comm_vel(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  void unpack_comm_vel(int, int, double *);
  int pack_reverse(int, int, double *);
  void unpack_reverse(int, int *, double *);
  int pack_border(int, int *, double *, int, int *);
  int pack_border_vel(int, int *, double *, int, int *);
  void unpack_border(int, int, double *);
  void unpack_border_vel(int, int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
 private:
  class Cuda *cuda;
  bool cuda_init_done;
  int* copylist;
  int* copylist2;
  cCudaData<int, int, xx >* cu_copylist;
  int max_nsend;
};

}

#endif
#endif
