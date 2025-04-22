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

#include "compute_esp_grid_local.h"

#include "atom.h"
#include "grid3d.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include "domain.h"
#include "utils.h"
#include "force.h"
#include "pair.h"

using namespace LAMMPS_NS;

ComputeESPGridLocal::ComputeESPGridLocal(LAMMPS *lmp, int narg, char **arg)
  : ComputeGridLocal(lmp, narg, arg), grid(nullptr), esp_ref(nullptr), bcut_acks2(nullptr)
{
  if (narg != 10)
    error->all(FLERR, "Illegal compute esp/grid/local command");

  nx = utils::inumeric(FLERR, arg[3], false, lmp);
  ny = utils::inumeric(FLERR, arg[4], false, lmp);
  nz = utils::inumeric(FLERR, arg[5], false, lmp);
  xlo = utils::numeric(FLERR, arg[6], false, lmp);
  ylo = utils::numeric(FLERR, arg[7], false, lmp);
  zlo = utils::numeric(FLERR, arg[8], false, lmp);
  spacing = utils::numeric(FLERR, arg[9], false, lmp);

  scalar_flag = 1;
  local_flag = 1;
}

ComputeESPGridLocal::~ComputeESPGridLocal()
{
  delete grid;
  memory->destroy(esp_ref);

  if (!reaxflag)
    memory->destroy(bcut_acks2);

}

void ComputeESPGridLocal::init()
{

  Pair *pair;

  if (pair = force->pair_match("^reaxff",0)) {
    reaxflag = 1;
    int ignore = 0;
    bcut_acks2 = (double *) pair->extract("bcut_acks2", ignore);
    if (!bcut_acks2)
      error->all(FLERR, "Compute esp/grid/local could not extract bcut_acks2 from pair reaxff");
  } else {
    reaxflag = 0;
    error->warning(FLERR,"No reaxff pair style for compute esp/grid/local, using default cutoff 5.0 angstrom.");
  }

  // memory->destroy(esp_ref);

  // memory->create3d_offset(data3d_one, nzlo_out, nzhi_out, nylo_out,
  //    nyhi_out, nxlo_out, nxhi_out, "data3d_one");

  // memory->create3d_offset(esp_ref, nx, ny, nz, ilo, ihi, jlo, jhi, klo, khi, "espgrid:esp_ref");

}


void ComputeESPGridLocal::compute_local()
{
  invoked_scalar = update->ntimestep;

  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double loss_sum = 0.0;
  double weight_sum = 0.0;

  int i, j, k;
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {

    double gx = xlo + (i + 0.5) * spacing;
    double gy = ylo + (j + 0.5) * spacing;
    double gz = zlo + (k + 0.5) * spacing;

    double r_min2 = 1e6;
    int t_nearest = -1;

    for (int a = 0; a < nlocal; ++a) {
      double dx = gx - x[a][0];
      double dy = gy - x[a][1];
      double dz = gz - x[a][2];
      double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < r_min2) {
        r_min2 = r2;
        t_nearest = type[a];
      }
    }

    if (t_nearest < 1) continue;
    double r = sqrt(r_min2);
    double rcut = bcut_acks2[t_nearest];

    if (r < 1.4 || r > rcut) continue;

    double w = 1.0 / (r * r);
    w *= (1.0 - exp(-pow((r - 1.4) / 0.3, 6)));
    w *= exp(-pow(r / rcut, 6));

    double V_model = 0.0;
    for (int a = 0; a < nlocal; ++a) {
      double dx = gx - x[a][0];
      double dy = gy - x[a][1];
      double dz = gz - x[a][2];
      double r2 = dx*dx + dy*dy + dz*dz;
      V_model += q[a] / sqrt(r2 + 1e-12);
    }

    double V_ref = esp_ref[i][j][k];
    double diff = V_model - V_ref;

    loss_sum += w * diff * diff;
    weight_sum += w;
  }

  scalar = (weight_sum > 0.0) ? loss_sum / weight_sum : 0.0;
}

void *ComputeESPGridLocal::extract_reference()
{
  return static_cast<void *>(esp_ref);
}
