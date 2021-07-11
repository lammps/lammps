// clang-format off
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

#ifndef LMP_GRIDCOMM_KOKKOS_H
#define LMP_GRIDCOMM_KOKKOS_H

#include "gridcomm.h"
#include "kokkos_type.h"
#include "fftdata_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class GridCommKokkos : public GridComm {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef FFTArrayTypes<DeviceType> FFT_AT;
  GridCommKokkos(class LAMMPS *, MPI_Comm, int, int, int,
           int, int, int, int, int, int,
           int, int, int, int, int, int);
  GridCommKokkos(class LAMMPS *, MPI_Comm, int, int, int, int,
           int, int, int, int, int, int,
           int, int, int, int, int, int,
           int, int, int, int, int, int);
  virtual ~GridCommKokkos();
  void forward_comm_kspace(class KSpace *, int, int,
                           FFT_DAT::tdual_FFT_SCALAR_1d &, FFT_DAT::tdual_FFT_SCALAR_1d &, MPI_Datatype);
  void reverse_comm_kspace(class KSpace *, int, int,
                           FFT_DAT::tdual_FFT_SCALAR_1d &, FFT_DAT::tdual_FFT_SCALAR_1d &, MPI_Datatype);

 private:
  DAT::tdual_int_2d k_swap_packlist;
  DAT::tdual_int_2d k_swap_unpacklist;

  DAT::tdual_int_2d k_send_packlist;

  DAT::tdual_int_2d k_recv_unpacklist;

  DAT::tdual_int_2d k_copy_packlist;
  DAT::tdual_int_2d k_copy_unpacklist;

  // -------------------------------------------
  // internal methods
  // -------------------------------------------

  void setup_regular(int &, int &);
  void setup_tiled(int &, int &);

  void forward_comm_kspace_regular(class KSpace *, int, int,
                                   FFT_DAT::tdual_FFT_SCALAR_1d &, FFT_DAT::tdual_FFT_SCALAR_1d &, MPI_Datatype);
  void forward_comm_kspace_tiled(class KSpace *, int, int,
                                 FFT_DAT::tdual_FFT_SCALAR_1d &, FFT_DAT::tdual_FFT_SCALAR_1d &, MPI_Datatype);
  void reverse_comm_kspace_regular(class KSpace *, int, int,
                                   FFT_DAT::tdual_FFT_SCALAR_1d &, FFT_DAT::tdual_FFT_SCALAR_1d &, MPI_Datatype);
  void reverse_comm_kspace_tiled(class KSpace *, int, int,
                                 FFT_DAT::tdual_FFT_SCALAR_1d &, FFT_DAT::tdual_FFT_SCALAR_1d &, MPI_Datatype);

  void grow_swap();

  int indices(DAT::tdual_int_2d &, int, int, int, int, int, int, int);
};

}

#endif
