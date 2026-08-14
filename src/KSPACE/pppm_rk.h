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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/rk,PPPM_RK);
// clang-format on
#else

#ifndef LMP_PPPM_RK_H
#define LMP_PPPM_RK_H

#include "pppm.h"

#include <cstring>

namespace LAMMPS_NS {

class PPPM_RK : public PPPM {
 public:
  PPPM_RK(class LAMMPS *);
  ~PPPM_RK() override;

  void init() override;

  void compute_grid_potentials(int, int) override;
  void compute(int, int) override;
  void r2k_comm(int &eflag, int &vflag) override;
  void k2r_comm(int eflag, int vflag) override;

 protected:
  virtual void compute_charge_densities(int, int);
  virtual void compute_interpolate_forces(int, int);
  FFT_SCALAR ***density_brick_buf;
  FFT_SCALAR ***vdx_brick_buf, ***vdy_brick_buf, ***vdz_brick_buf;
  FFT_SCALAR ***u_brick_buf;

  double energy_buf;
  double virial_mpi_buf[6];
  int flags_buf[2];
  double domain_boxbuf[6];

  int rproc;         // 1 if an Rspace proc, 0 if Kspace
  int me_block;      // proc ID within Rspace/Kspace block
  int ratio;         // ratio of Rspace procs to Kspace procs
  MPI_Comm block;    // communicator within one block

  MPI_Comm block_density_brick;    // communicator within one block
  MPI_Comm block_kforce_x;
  MPI_Comm block_kforce_y;
  MPI_Comm block_kforce_z;
  MPI_Comm block_kforce_u;
  MPI_Comm block_energy;
  MPI_Comm block_virial;

  MPI_Request mpi_requests_density;    //Requests for asynchronous communications
  MPI_Request mpi_requests_grid_x;
  MPI_Request mpi_requests_grid_y;
  MPI_Request mpi_requests_grid_z;
  MPI_Request mpi_requests_grid_u;
  MPI_Request mpi_requests_energy;
  MPI_Request mpi_requests_virial;
  MPI_Request mpi_requests_flags;
  MPI_Request mpi_requests_domain_box;

  // Tags for communication within MPI_Comm block
  enum MPI_Block_Tags { TAG_FLAGS = 1, TAG_BOX, N_MPI_TAGS };

  int *density_sizes, *density_disps;    // MPI gather/scatter params for block comm
  int **partitionInfoK;
  int block_size;

  void setupRKBlock();

  void prepareXYZBrickScatterBufs()
  {
    FFT_SCALAR *x_buf_loc = &vdx_brick_buf[nzlo_in][nylo_in][nxlo_in];
    FFT_SCALAR *y_buf_loc = &vdy_brick_buf[nzlo_in][nylo_in][nxlo_in];
    FFT_SCALAR *z_buf_loc = &vdz_brick_buf[nzlo_in][nylo_in][nxlo_in];
    int zlo, zhi, ylo, yhi, xlo, xhi, ncpy;
    for (int k = 1; k < block_size; k++) {
      zlo = partitionInfoK[k][1];
      zhi = partitionInfoK[k][2];
      ylo = partitionInfoK[k][3];
      yhi = partitionInfoK[k][4];
      xlo = partitionInfoK[k][5];
      xhi = partitionInfoK[k][6];

      ncpy = xhi - xlo + 1;
      for (int zz = zlo; zz <= zhi; zz++)
        for (int yy = ylo; yy <= yhi; yy++) {
          memcpy(x_buf_loc, &vdx_brick[zz][yy][xlo], ncpy * sizeof(FFT_SCALAR));
          memcpy(y_buf_loc, &vdy_brick[zz][yy][xlo], ncpy * sizeof(FFT_SCALAR));
          memcpy(z_buf_loc, &vdz_brick[zz][yy][xlo], ncpy * sizeof(FFT_SCALAR));
          x_buf_loc += ncpy;
          y_buf_loc += ncpy;
          z_buf_loc += ncpy;
        }
    }
  }
  void prepareUBrickScatterBufs()
  {
    FFT_SCALAR *u_buf_loc = &u_brick_buf[nzlo_in][nylo_in][nxlo_in];
    int zlo, zhi, ylo, yhi, xlo, xhi, ncpy;
    for (int k = 1; k < block_size; k++) {
      zlo = partitionInfoK[k][1];
      zhi = partitionInfoK[k][2];
      ylo = partitionInfoK[k][3];
      yhi = partitionInfoK[k][4];
      xlo = partitionInfoK[k][5];
      xhi = partitionInfoK[k][6];

      ncpy = xhi - xlo + 1;
      for (int zz = zlo; zz <= zhi; zz++)
        for (int yy = ylo; yy <= yhi; yy++) {
          memcpy(u_buf_loc, &u_brick[zz][yy][xlo], ncpy * sizeof(FFT_SCALAR));
          u_buf_loc += ncpy;
        }
    }
  }

  virtual void waitReceiptGridPotentialsEV(int eflag, int vflag)
  {
    if (eflag) MPI_Wait(&mpi_requests_energy, MPI_STATUS_IGNORE);
    if (vflag) MPI_Wait(&mpi_requests_virial, MPI_STATUS_IGNORE);

    if (differentiation_flag == 1) {
      MPI_Wait(&mpi_requests_grid_u, MPI_STATUS_IGNORE);
    } else {
      MPI_Wait(&mpi_requests_grid_x, MPI_STATUS_IGNORE);
      MPI_Wait(&mpi_requests_grid_y, MPI_STATUS_IGNORE);
      MPI_Wait(&mpi_requests_grid_z, MPI_STATUS_IGNORE);
    }

    // Update from buffers
    // energy and virial are updated in PPPM_RK::compute_interpolate_forces(int eflag, int vflag)
    if (differentiation_flag == 1) {
      for (int k = nzlo_in; k <= nzhi_in; k++)
        for (int j = nylo_in; j <= nyhi_in; j++)
          memcpy(&u_brick[k][j][nxlo_in], &u_brick_buf[k][j][nxlo_in],
                 (nxhi_in - nxlo_in + 1) * sizeof(FFT_SCALAR));
    } else {
      for (int k = nzlo_in; k <= nzhi_in; k++)
        for (int j = nylo_in; j <= nyhi_in; j++) {
          memcpy(&vdx_brick[k][j][nxlo_in], &vdx_brick_buf[k][j][nxlo_in],
                 (nxhi_in - nxlo_in + 1) * sizeof(FFT_SCALAR));
          memcpy(&vdy_brick[k][j][nxlo_in], &vdy_brick_buf[k][j][nxlo_in],
                 (nxhi_in - nxlo_in + 1) * sizeof(FFT_SCALAR));
          memcpy(&vdz_brick[k][j][nxlo_in], &vdz_brick_buf[k][j][nxlo_in],
                 (nxhi_in - nxlo_in + 1) * sizeof(FFT_SCALAR));
        }
    }
  }
  void post_sending_scatter_grid_potentials_ev(int eflag, int vflag)
  {
    //energy_buf and virial_mpi_buf are filled at the end of
    //PPPM_RK::compute_grid_potentials(int eflag, int vflag)
    if (eflag) { MPI_Ibcast(&energy_buf, 1, MPI_DOUBLE, 0, block_energy, &mpi_requests_energy); }
    if (vflag) { MPI_Ibcast(virial_mpi_buf, 6, MPI_DOUBLE, 0, block_virial, &mpi_requests_virial); }
    if (differentiation_flag == 1) {
      FFT_SCALAR *usrc = &u_brick_buf[nzlo_in][nylo_in][nxlo_in];
      prepareUBrickScatterBufs();
      MPI_Iscatterv(usrc, density_sizes, density_disps, MPI_FFT_SCALAR, nullptr, 0, MPI_FFT_SCALAR,
                    0, block_kforce_u, &mpi_requests_grid_u);
    } else {
      FFT_SCALAR *xsrc = &vdx_brick_buf[nzlo_in][nylo_in][nxlo_in];
      FFT_SCALAR *ysrc = &vdy_brick_buf[nzlo_in][nylo_in][nxlo_in];
      FFT_SCALAR *zsrc = &vdz_brick_buf[nzlo_in][nylo_in][nxlo_in];
      prepareXYZBrickScatterBufs();
      MPI_Iscatterv(xsrc, density_sizes, density_disps, MPI_FFT_SCALAR, nullptr, 0, MPI_FFT_SCALAR,
                    0, block_kforce_x, &mpi_requests_grid_x);
      MPI_Iscatterv(ysrc, density_sizes, density_disps, MPI_FFT_SCALAR, nullptr, 0, MPI_FFT_SCALAR,
                    0, block_kforce_y, &mpi_requests_grid_y);
      MPI_Iscatterv(zsrc, density_sizes, density_disps, MPI_FFT_SCALAR, nullptr, 0, MPI_FFT_SCALAR,
                    0, block_kforce_z, &mpi_requests_grid_z);
    }
  }
  void wait_sending_scatter_grid_potentials_ev(int eflag, int vflag)
  {
    if (eflag) MPI_Wait(&mpi_requests_energy, MPI_STATUS_IGNORE);
    if (vflag) MPI_Wait(&mpi_requests_virial, MPI_STATUS_IGNORE);

    if (differentiation_flag == 1) {
      MPI_Wait(&mpi_requests_grid_u, MPI_STATUS_IGNORE);
    } else {
      MPI_Wait(&mpi_requests_grid_x, MPI_STATUS_IGNORE);
      MPI_Wait(&mpi_requests_grid_y, MPI_STATUS_IGNORE);
      MPI_Wait(&mpi_requests_grid_z, MPI_STATUS_IGNORE);
    }
  }
};    // class PPPM_RK

}    // namespace LAMMPS_NS

#endif
#endif
