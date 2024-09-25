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
   Contributing authors:
   Joel Clemmer (SNL), Thomas O'Connor (CMU)
----------------------------------------------------------------------- */

#include "pair_rheo.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_grad.h"
#include "compute_rheo_interface.h"
#include "compute_rheo_kernel.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "fix_rheo_pressure.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "utils.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace RHEO_NS;
using namespace MathExtra;

static constexpr double EPSILON = 1e-2;

/* ---------------------------------------------------------------------- */

PairRHEO::PairRHEO(LAMMPS *lmp) :
    Pair(lmp), csq(nullptr), rho0(nullptr), cs(nullptr), compute_kernel(nullptr),
    compute_grad(nullptr), compute_interface(nullptr), fix_rheo(nullptr), fix_pressure(nullptr)
{
  restartinfo = 0;
  single_enable = 0;

  artificial_visc_flag = 0;
  rho_damp_flag = 0;
  thermal_flag = 0;
  harmonic_means_flag = 0;

  comm_reverse = 3;
}

/* ---------------------------------------------------------------------- */

PairRHEO::~PairRHEO()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }

  memory->destroy(cs);
}

/* ---------------------------------------------------------------------- */

void PairRHEO::compute(int eflag, int vflag)
{
  int i, j, a, b, ii, jj, inum, jnum, itype, jtype;
  int pair_force_flag, pair_rho_flag, pair_avisc_flag;
  int fluidi, fluidj;
  double xtmp, ytmp, ztmp, wp, Ti, Tj, dT, csq_ave, cs_ave;
  double rhoi, rhoj, rho0i, rho0j, voli, volj, Pi, Pj, etai, etaj, kappai, kappaj, eta_ave,
      kappa_ave, dT_prefactor;
  double mu, q, fp_prefactor, drho_damp, fmag, psi_ij, Fij;
  double *dWij, *dWji;
  double dx[3], du[3], dv[3], fv[3], dfp[3], fsolid[3], ft[3], vi[3], vj[3];

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, rsq, r, rinv;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int dim = domain->dimension;

  ev_init(eflag, vflag);

  double **gradv = compute_grad->gradv;
  double **gradr = compute_grad->gradr;
  double **v = atom->v;
  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *drho = atom->drho;
  double *pressure = atom->pressure;
  double *viscosity = atom->viscosity;
  double *conductivity = atom->conductivity;
  double *temperature = atom->temperature;
  double *heatflow = atom->heatflow;
  int *type = atom->type;
  int *status = atom->rheo_status;

  double **fp_store, *chi;
  if (compute_interface) {
    fp_store = compute_interface->fp_store;
    chi = compute_interface->chi;

    for (i = 0; i < atom->nmax; i++) {
      fp_store[i][0] = 0.0;
      fp_store[i][1] = 0.0;
      fp_store[i][2] = 0.0;
    }
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (a = 0; a < 3; a++) {
    vi[a] = 0.0;
    vj[a] = 0.0;
    du[a] = 0.0;
    fv[a] = 0.0;
    dfp[a] = 0.0;
    fsolid[a] = 0.0;
    ft[0] = 0.0;
  }

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    if (rmass)
      imass = rmass[i];
    else
      imass = mass[itype];
    etai = viscosity[i];
    fluidi = !(status[i] & PHASECHECK);
    if (thermal_flag) {
      kappai = conductivity[i];
      Ti = temperature[i];
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = lensq3(dx);
      jtype = type[j];

      if (rsq < cutksq) {
        r = sqrt(rsq);
        rinv = 1 / r;

        if (rmass)
          jmass = rmass[i];
        else
          jmass = mass[jtype];
        etaj = viscosity[j];
        fluidj = !(status[j] & PHASECHECK);
        if (thermal_flag) {
          Tj = temperature[j];
          kappaj = conductivity[j];
        }

        cs_ave = 0.5 * (cs[itype] + cs[jtype]);
        csq_ave = cs_ave * cs_ave;

        pair_rho_flag = 0;
        pair_force_flag = 0;
        pair_avisc_flag = 0;
        if (fluidi || fluidj) {
          pair_force_flag = 1;
          if (!interface_flag) pair_avisc_flag = 1;
        }
        if (fluidi && fluidj) {
          pair_rho_flag = 1;
          pair_avisc_flag = 1;
        }

        wp = compute_kernel->calc_dw(i, j, dx[0], dx[1], dx[2], r);
        dWij = compute_kernel->dWij;
        dWji = compute_kernel->dWji;

        for (a = 0; a < dim; a++) {
          vi[a] = v[i][a];
          vj[a] = v[j][a];
        }

        // Add corrections for walls
        rhoi = rho[i];
        rhoj = rho[j];
        rho0i = rho0[itype];
        rho0j = rho0[jtype];
        Pi = pressure[i];
        Pj = pressure[j];
        fmag = 0;
        if (interface_flag) {
          if (fluidi && (!fluidj)) {
            compute_interface->correct_v(vj, vi, j, i);
            rhoj = compute_interface->correct_rho(j);
            Pj = fix_pressure->calc_pressure(rhoj, jtype);

            if ((chi[j] > 0.9) && (r < (cutk * 0.5)))
              fmag = (chi[j] - 0.9) * (cutk * 0.5 - r) * rho0j * csq_ave * cutk * rinv;

          } else if ((!fluidi) && fluidj) {
            compute_interface->correct_v(vi, vj, i, j);
            rhoi = compute_interface->correct_rho(i);
            Pi = fix_pressure->calc_pressure(rhoi, itype);

            if (chi[i] > 0.9 && r < (cutk * 0.5))
              fmag = (chi[i] - 0.9) * (cutk * 0.5 - r) * rho0i * csq_ave * cutk * rinv;

          } else if ((!fluidi) && (!fluidj)) {
            rhoi = rho0i;
            rhoj = rho0j;
          }
        }

        // Repel if close to inner solid particle
        scale3(fmag, dx, fsolid);

        // Compute volume after reconstructing
        voli = imass / rhoi;
        volj = jmass / rhoj;

        // Thermal Evolution
        if (thermal_flag) {
          if (harmonic_means_flag) {
            kappa_ave = 2.0 * kappai * kappaj / (kappai + kappaj);
          } else {
            kappa_ave = 0.5 * (kappai + kappaj);
          }
          dT_prefactor =
              2.0 * kappa_ave * (Ti - Tj) * rinv * rinv * voli * volj * 2.0 / (rhoi + rhoj);

          dT = dot3(dx, dWij);
          heatflow[i] += dT * dT_prefactor;

          if (newton_pair || j < nlocal) {
            dT = dot3(dx, dWji);
            heatflow[j] += dT * dT_prefactor;
          }
        }

        if (pair_force_flag) {

          //Hydrostatic pressure forces
          fp_prefactor = voli * volj * (Pj + Pi);
          sub3(vi, vj, dv);

          if (harmonic_means_flag) {
            eta_ave = 2.0 * etai * etaj / (etai + etaj);
          } else {
            eta_ave = 0.5 * (etai + etaj);
          }

          //Add artificial viscous pressure if required
          if (artificial_visc_flag && pair_avisc_flag) {
            //Interpolate velocities to midpoint and use this difference for artificial viscosity
            copy3(dv, du);
            for (a = 0; a < dim; a++)
              for (b = 0; b < dim; b++)
                du[a] -= 0.5 * (gradv[i][a * dim + b] + gradv[j][a * dim + b]) * dx[b];

            mu = dot3(du, dx) * cutkinv3;
            mu /= (rsq * cutkinv3 * cutkinv3 + EPSILON);
            mu = MIN(0.0, mu);
            q = av * (-2.0 * cs_ave * mu + mu * mu);
            fp_prefactor += voli * volj * q * (rhoj + rhoi);
          }

          // -Grad[P + Q]
          scale3(-fp_prefactor, dWij, dfp);

          // Now compute viscous eta*Lap[v] terms
          for (a = 0; a < dim; a++) {
            fv[a] = 0.0;
            for (b = 0; b < dim; b++) fv[a] += dv[a] * dx[b] * dWij[b];
            fv[a] *= 2.0 * eta_ave * voli * volj * rinv * rinv;
          }

          add3(fv, dfp, ft);
          add3(fsolid, ft, ft);

          f[i][0] += ft[0];
          f[i][1] += ft[1];
          f[i][2] += ft[2];

          // Note the virial's definition is hazy, e.g. viscous contributions will depend on rotation
          if (evflag)
            ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, ft[0], ft[1], ft[2], dx[0], dx[1],
                         dx[2]);

          if (newton_pair || j < nlocal) {
            for (a = 0; a < dim; a++) {
              fv[a] = 0.0;
              for (b = 0; b < dim; b++) fv[a] += (vi[a] - vj[a]) * dx[b] * dWji[b];
              fv[a] *= -2.0 * eta_ave * voli * volj * rinv * rinv;
              // flip sign here b/c -= at accummulator
            }

            scale3(fp_prefactor, dWji, dfp);
            add3(fv, dfp, ft);
            add3(fsolid, ft, ft);

            f[j][0] -= ft[0];
            f[j][1] -= ft[1];
            f[j][2] -= ft[2];

            if (evflag)
              ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, ft[0], ft[1], ft[2], -dx[0], -dx[1],
                           -dx[2]);
          }

          if (compute_interface) {
            fp_store[i][0] += dfp[0];
            fp_store[i][1] += dfp[1];
            fp_store[i][2] += dfp[2];

            if (newton_pair || j < nlocal) {
              fp_store[j][0] -= dfp[0];
              fp_store[j][1] -= dfp[1];
              fp_store[j][2] -= dfp[2];
            }
          }
        }

        // Density damping
        // conventional for low-order h
        // interpolated for RK 1 & 2  (Antuono et al., Computers & Fluids 2021)
        if (rho_damp_flag && pair_rho_flag) {
          if (laplacian_order >= 1) {
            psi_ij = rhoj - rhoi;
            Fij = -rinv * rinv * dot3(dx, dWij);
            for (a = 0; a < dim; a++) psi_ij += 0.5 * (gradr[i][a] + gradr[j][a]) * dx[a];
            drho[i] += 2 * rho_damp * psi_ij * Fij * volj;
          } else {
            drho_damp = 2 * rho_damp * ((rhoj - rho0[jtype]) - (rhoi - rho0[itype])) * rinv * wp;
            drho[i] -= drho_damp * volj;
          }

          if (newton_pair || j < nlocal) {
            if (laplacian_order >= 1) {
              Fij = rinv * rinv * dot3(dx, dWji);
              psi_ij *= -1;
              drho[j] += 2 * rho_damp * psi_ij * Fij * voli;
            } else {
              drho[j] += drho_damp * voli;
            }
          }
        }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

  if (compute_interface) {
    if (newton_pair) comm->reverse_comm(this);
    comm->forward_comm(this);
  }
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairRHEO::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairRHEO::settings(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal pair_style command");

  cutk = utils::numeric(FLERR, arg[0], false, lmp);

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "rho/damp") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "pair rheo rho/damp", error);
      rho_damp_flag = 1;
      rho_damp = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg++;
    } else if (strcmp(arg[iarg], "artificial/visc") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "pair rheo artificial/visc", error);
      artificial_visc_flag = 1;
      av = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg++;
    } else if (strcmp(arg[iarg], "harmonic/means") == 0) {
      harmonic_means_flag = 1;
    } else
      error->all(FLERR, "Illegal pair_style command, {}", arg[iarg]);
    iarg++;
  }
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairRHEO::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Incorrect number of args for pair_style rheo coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = 0; j <= atom->ntypes; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair rheo coefficients");
}

/* ----------------------------------------------------------------------
 setup specific to this pair style
 ------------------------------------------------------------------------- */

void PairRHEO::setup()
{
  auto fixes = modify->get_fix_by_style("^rheo$");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use pair rheo");
  if (fixes.size() > 1) error->all(FLERR, "Must have only one fix rheo defined");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);

  // Currently only allow one instance of fix rheo/pressure
  fixes = modify->get_fix_by_style("^rheo/pressure$");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo/pressure to use pair rheo");
  if (fixes.size() > 1) error->all(FLERR, "Must have only one fix rheo/pressure defined");
  fix_pressure = dynamic_cast<FixRHEOPressure *>(fixes[0]);

  compute_kernel = fix_rheo->compute_kernel;
  compute_grad = fix_rheo->compute_grad;
  compute_interface = fix_rheo->compute_interface;
  thermal_flag = fix_rheo->thermal_flag;
  interface_flag = fix_rheo->interface_flag;
  csq = fix_rheo->csq;
  rho0 = fix_rheo->rho0;

  if (cutk != fix_rheo->cut)
    error->all(FLERR, "Pair rheo cutoff {} does not agree with fix rheo cutoff {}", cutk,
               fix_rheo->cut);

  cutksq = cutk * cutk;
  cutkinv = 1.0 / cutk;
  cutkinv3 = cutkinv * 3.0;
  laplacian_order = -1;

  int n = atom->ntypes;
  memory->create(cs, n + 1, "rheo:cs");
  for (int i = 1; i <= n; i++) cs[i] = sqrt(csq[i]);

  if (comm->ghost_velocity == 0) error->all(FLERR, "Pair RHEO requires ghost atoms store velocity");

  if (laplacian_order == -1) {
    if (fix_rheo->kernel_style == RK2)
      laplacian_order = 2;
    else if (fix_rheo->kernel_style == RK1)
      laplacian_order = 1;
    else
      laplacian_order = 0;
  }
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairRHEO::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair rheo coeffs are not set");

  return cutk;
}

/* ---------------------------------------------------------------------- */

int PairRHEO::pack_reverse_comm(int n, int first, double *buf)
{
  double **fp_store = compute_interface->fp_store;
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    buf[m++] = fp_store[i][0];
    buf[m++] = fp_store[i][1];
    buf[m++] = fp_store[i][2];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairRHEO::unpack_reverse_comm(int n, int *list, double *buf)
{
  double **fp_store = compute_interface->fp_store;
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    fp_store[j][0] += buf[m++];
    fp_store[j][1] += buf[m++];
    fp_store[j][2] += buf[m++];
  }
}
