// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_GRID3D_KOKKOS_H
#define LMP_GRID3D_KOKKOS_H

#include "grid3d.h"
#include "kokkos_type.h"
#include "fftdata_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class Grid3dKokkos : public Grid3d {
 public:
  enum { KSPACE = 0, PAIR = 1, FIX = 2 };    // calling classes
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef FFTArrayTypes<DeviceType> FFT_AT;
  Grid3dKokkos(class LAMMPS *, MPI_Comm, int, int, int);
  Grid3dKokkos(class LAMMPS *, MPI_Comm, int, int, int,
         int, int, int, int, int, int, int, int, int, int, int, int);
  ~Grid3dKokkos() override;

  void forward_comm(int, void *, int, int, int,
                    FFT_DAT::tdual_FFT_SCALAR_1d &, FFT_DAT::tdual_FFT_SCALAR_1d &, MPI_Datatype);
  void reverse_comm(int, void *, int, int, int,
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

  void setup_comm_brick(int &, int &) override;
  void setup_comm_tiled(int &, int &) override;

  void forward_comm_kspace_brick(class KSpace *, int, int,
                                   FFT_DAT::tdual_FFT_SCALAR_1d &, FFT_DAT::tdual_FFT_SCALAR_1d &, MPI_Datatype);
  void forward_comm_kspace_tiled(class KSpace *, int, int,
                                 FFT_DAT::tdual_FFT_SCALAR_1d &, FFT_DAT::tdual_FFT_SCALAR_1d &, MPI_Datatype);
  void reverse_comm_kspace_brick(class KSpace *, int, int,
                                   FFT_DAT::tdual_FFT_SCALAR_1d &, FFT_DAT::tdual_FFT_SCALAR_1d &, MPI_Datatype);
  void reverse_comm_kspace_tiled(class KSpace *, int, int,
                                 FFT_DAT::tdual_FFT_SCALAR_1d &, FFT_DAT::tdual_FFT_SCALAR_1d &, MPI_Datatype);

  void grow_swap() override;

  int indices(DAT::tdual_int_2d &, int, int, int, int, int, int, int);
};

}    // namespace LAMMPS_NS

#endif
