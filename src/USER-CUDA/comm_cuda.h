/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_COMM_CUDA_H
#define LMP_COMM_CUDA_H

#include "pointers.h"

#include "cuda_data.h"
#include "comm.h"

namespace LAMMPS_NS {

class CommCuda : public Comm {
public:
  CommCuda(class LAMMPS *);
  ~CommCuda();

  virtual void init();
  virtual void setup();                     // setup 3d communication pattern
  virtual void forward_comm(int mode=0);              // forward communication of atom coords
  virtual void forward_comm_cuda();
  virtual void forward_comm_pack_cuda();
  virtual void forward_comm_transfer_cuda();
  virtual void forward_comm_unpack_cuda();
  virtual void forward_comm_pair(Pair *pair);
  virtual void reverse_comm();              // reverse communication of forces
  virtual void exchange();                  // move atoms to new procs
  virtual void exchange_cuda();                  // move atoms to new procs
  virtual void borders();                   // setup list of atoms to communicate
  virtual void borders_cuda();                   // setup list of atoms to communicate
  virtual void borders_cuda_overlap_forward_comm();
  virtual void forward_comm_fix(class Fix *);          // forward comm from a Fix




 protected:
  class Cuda *cuda;
  cCudaData<int, int, xy>* cu_pbc;
  cCudaData<double, X_FLOAT, x>* cu_slablo;
  cCudaData<double, X_FLOAT, x>* cu_slabhi;
  cCudaData<double, X_FLOAT, xy>* cu_multilo;
  cCudaData<double, X_FLOAT, xy>* cu_multihi;

  cCudaData<int, int, xy>* cu_sendlist;
  virtual void grow_send(int,int);          // reallocate send buffer
  virtual void grow_recv(int);              // free/allocate recv buffer
  virtual void grow_list(int, int);         // reallocate one sendlist
  virtual void grow_swap(int);              // grow swap and multi arrays
  virtual void allocate_swap(int);          // allocate swap arrays
  virtual void allocate_multi(int);         // allocate multi arrays
  virtual void free_swap();                 // free swap arrays
  virtual void free_multi();                // free multi arrays
};

}

#endif
