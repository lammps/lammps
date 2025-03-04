/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Sebastian Huetter (OvGU)
------------------------------------------------------------------------- */

#include "meam.h"

#include "memory.h"

using namespace LAMMPS_NS;
using namespace MEAM_NS;

/* ---------------------------------------------------------------------- */

MEAM::MEAM(Memory *mem) : memory(mem)
{
  phir = phirar = phirar1 = phirar2 = phirar3 = phirar4 = phirar5 = phirar6 = nullptr;

  nmax = 0;
  rho = rho0 = rho1 = rho2 = rho3 = frhop = nullptr;
  gamma = dgamma1 = dgamma2 = dgamma3 = arho2b = nullptr;
  arho1 = arho2 = arho3 = arho3b = t_ave = tsq_ave = nullptr;

  // msmeam arrays
  msmeamflag = 0;
  arho2mb = nullptr;
  arho1m = arho2m = arho3m = arho3mb = nullptr;

  maxneigh = 0;
  scrfcn = dscrfcn = fcpair = nullptr;
  copymode = 0;

  // clang-format on
  neltypes = 0;
  for (int i = 0; i < MAXELT; i++) {
    A_meam[i] = rho0_meam[i] = beta0_meam[i] = beta1_meam[i] = beta2_meam[i] = beta3_meam[i] =
        t0_meam[i] = t1_meam[i] = t2_meam[i] = t3_meam[i] = rho_ref_meam[i] = t1m_meam[i] =
            t2m_meam[i] = t3m_meam[i] = beta1m_meam[i] = beta2m_meam[i] = beta3m_meam[i] = 0.0;
    ibar_meam[i] = ielt_meam[i] = 0;
    for (int j = 0; j < MAXELT; j++) {
      lattce_meam[i][j] = FCC;
      Ec_meam[i][j] = re_meam[i][j] = alpha_meam[i][j] = delta_meam[i][j] = ebound_meam[i][j] =
          attrac_meam[i][j] = repuls_meam[i][j] = 0.0;
      nn2_meam[i][j] = zbl_meam[i][j] = eltind[i][j] = 0;
    }
  }
  // clang-format off

  // indices and weights for Voigt notation
  int nv2 = 0;
  int nv3 = 0;
  for (int m = 0; m < 3; m++) {
    for (int n = m; n < 3; n++) {
      vind2D[m][n] = nv2;
      vind2D[n][m] = nv2;
      nv2 = nv2 + 1;
      for (int p = n; p < 3; p++) {
        vind3D[m][n][p] = nv3;
        vind3D[m][p][n] = nv3;
        vind3D[n][m][p] = nv3;
        vind3D[n][p][m] = nv3;
        vind3D[p][m][n] = nv3;
        vind3D[p][n][m] = nv3;
        nv3 = nv3 + 1;
      }
    }
  }

  v2D[0] = 1;
  v2D[1] = 2;
  v2D[2] = 2;
  v2D[3] = 1;
  v2D[4] = 2;
  v2D[5] = 1;

  v3D[0] = 1;
  v3D[1] = 3;
  v3D[2] = 3;
  v3D[3] = 3;
  v3D[4] = 6;
  v3D[5] = 3;
  v3D[6] = 1;
  v3D[7] = 3;
  v3D[8] = 3;
  v3D[9] = 1;
}

MEAM::~MEAM()
{
  if (copymode) return;

  memory->destroy(phirar6);
  memory->destroy(phirar5);
  memory->destroy(phirar4);
  memory->destroy(phirar3);
  memory->destroy(phirar2);
  memory->destroy(phirar1);
  memory->destroy(phirar);
  memory->destroy(phir);

  memory->destroy(rho);
  memory->destroy(rho0);
  memory->destroy(rho1);
  memory->destroy(rho2);
  memory->destroy(rho3);
  memory->destroy(frhop);
  memory->destroy(gamma);
  memory->destroy(dgamma1);
  memory->destroy(dgamma2);
  memory->destroy(dgamma3);
  memory->destroy(arho2b);

  memory->destroy(arho1);
  memory->destroy(arho2);
  memory->destroy(arho3);
  memory->destroy(arho3b);
  memory->destroy(t_ave);
  memory->destroy(tsq_ave);

  memory->destroy(scrfcn);
  memory->destroy(dscrfcn);
  memory->destroy(fcpair);

  // msmeam
  if (msmeamflag) {
    memory->destroy(arho1m);
    memory->destroy(arho2m);
    memory->destroy(arho3m);
    memory->destroy(arho2mb);
    memory->destroy(arho3mb);
  }
}
