// clang-format off
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
   Contributing author: Sebastian Hütter (OvGU)
------------------------------------------------------------------------- */

#include "meam.h"

#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MEAM::MEAM(Memory* mem)
  : memory(mem)
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

  neltypes = 0;
  for (int i = 0; i < maxelt; i++) {
    A_meam[i] = rho0_meam[i] = beta0_meam[i] =
      beta1_meam[i]= beta2_meam[i] = beta3_meam[i] =
      t0_meam[i] = t1_meam[i] = t2_meam[i] = t3_meam[i] =
      rho_ref_meam[i] = ibar_meam[i] = ielt_meam[i] =
      t1m_meam[i] = t2m_meam[i] = t3m_meam[i] =
      beta1m_meam[i] = beta2m_meam[i] = beta3m_meam[i] = 0.0;
    for (int j = 0; j < maxelt; j++) {
      lattce_meam[i][j] = FCC;
      Ec_meam[i][j] = re_meam[i][j] = alpha_meam[i][j] = delta_meam[i][j] = ebound_meam[i][j] = attrac_meam[i][j] = repuls_meam[i][j] = 0.0;
      nn2_meam[i][j] = zbl_meam[i][j] = eltind[i][j] = 0;
    }
  }
}

MEAM::~MEAM()
{
  if (copymode) return;

  memory->destroy(this->phirar6);
  memory->destroy(this->phirar5);
  memory->destroy(this->phirar4);
  memory->destroy(this->phirar3);
  memory->destroy(this->phirar2);
  memory->destroy(this->phirar1);
  memory->destroy(this->phirar);
  memory->destroy(this->phir);

  memory->destroy(this->rho);
  memory->destroy(this->rho0);
  memory->destroy(this->rho1);
  memory->destroy(this->rho2);
  memory->destroy(this->rho3);
  memory->destroy(this->frhop);
  memory->destroy(this->gamma);
  memory->destroy(this->dgamma1);
  memory->destroy(this->dgamma2);
  memory->destroy(this->dgamma3);
  memory->destroy(this->arho2b);

  memory->destroy(this->arho1);
  memory->destroy(this->arho2);
  memory->destroy(this->arho3);
  memory->destroy(this->arho3b);
  memory->destroy(this->t_ave);
  memory->destroy(this->tsq_ave);

  memory->destroy(this->scrfcn);
  memory->destroy(this->dscrfcn);
  memory->destroy(this->fcpair);

  // msmeam
  if (this->msmeamflag){
    memory->destroy(this->arho1m);
    memory->destroy(this->arho2m);
    memory->destroy(this->arho3m);
    memory->destroy(this->arho2mb);
    memory->destroy(this->arho3mb);
  }
}
