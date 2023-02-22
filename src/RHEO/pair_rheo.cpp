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

#include "pair_rheo.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_kernel.h"
#include "compute_rheo_grad.h"
#include "compute_rheo_interface.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"
#include "utils.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathExtra;

#define EPSILON 1e-2

/* ---------------------------------------------------------------------- */

PairRHEO::PairRHEO(LAMMPS *lmp) :
  Pair(lmp), compute_kernel(nullptr), compute_grad(nullptr),
  compute_interface(nullptr), fix_rheo(nullptr)
{
  restartinfo = 0;
  single_enable = 0;

  artificial_visc_flag = 0;
  rho_damp_flag = 0;
  thermal_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairRHEO::~PairRHEO()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairRHEO::compute(int eflag, int vflag)
{
  int i, j, a, b, ii, jj, inum, jnum, itype, jtype;
  int error_flag, pair_force_flag, pair_rho_flag, pair_avisc_flag;
  double xtmp, ytmp, ztmp;
  int fluidi, fluidj;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, rsq, r, rinv;

  double w, wp, rhoi, rhoj, voli, volj, Pi, Pj;
  double *dWij, *dWji, *d2Wij, *d2Wji, *dW1ij, *dW1ji;
  double vijeij, etai, etaj, kappai, kappaj;
  double Ti, Tj, dT;
  double drho_damp, fmag;
  double mu, q, fp_prefactor;
  double dx[3] = {0};
  double fv[3] = {0};
  double dfp[3] = {0};
  double fsolid[3] = {0};
  double du[3] = {0};
  double vi[3] = {0};
  double vj[3] = {0};
  double dv[3] = {0};
  double psi_ij = 0.0;
  double Fij = 0.0;

  ev_init(eflag, vflag);

  double **gradv = compute_grad->gradv;
  double **gradt = compute_grad->gradt;
  double **gradr = compute_grad->gradr;
  double **gradn = compute_grad->gradn;
  double **v = atom->v;
  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *drho = atom->drho;
  double *temperature = atom->temperature;
  double *heatflow = atom->heatflow;
  double *chi = compute_interface->chi;
  double **f_pressure = fix_rheo->f_pressure;
  double *viscosity = fix_rheo->viscosity;
  double *conductivity = fix_rheo->conductivity;
  double *pressure = fix_rheo->pressure;
  double *special_lj = force->special_lj;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *status = atom->status;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int dim = domain->dimension;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imass = mass[itype];

    kappai = conductivity[i];
    etai = viscosity[i];
    Ti = temperature[i];
    fluidi = status[i] & FixRHEO::STATUS_FLUID;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = lensq3(dx);
      jtype = type[j];
      jmass = mass[jtype];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        rinv = 1 / r;

        kappaj = conductivity[j];
        etaj = viscosity[j];
        Tj = temperature[j];
        fluidj = status[j] & FixRHEO::STATUS_FLUID;

        pair_rho_flag = 0;
        pair_force_flag = 0;
        pair_avisc_flag = 0;
        if (fluidi || fluidj) {
          pair_force_flag = 1;
        }
        if (fluidi && fluidj) {
          pair_avisc_flag = 1;
          pair_rho_flag = 1;
        }

        wp = compute_kernel->calc_dw(i, j, dx[0], dx[1], dx[2], r);
        dWij = compute_kernel->dWij;
        dWji = compute_kernel->dWji;

        for (a = 0; a < dim; a++) {
          vi[a] = v[i][a];
          vj[a] = v[j][a];
        }

        // TODO: search for fix pressure and add calculate P function
        // Add corrections for walls
        rhoi = rho[i];
        rhoj = rho[j];
        Pi = pressure[i];
        Pj = pressure[j];
        fmag = 0;
        if (fluidi && (!fluidj)) {
          compute_interface->correct_v(v[i], v[j], vi, i, j);
          rhoj = compute_interface->correct_rho(j, i);
          Pj = fix_pressure->calculate_P(rhoj);
          if ((chi[j] > 0.9) && (r < (h * 0.5)))
            fmag = (chi[j] - 0.9) * (h * 0.5 - r) * rho0 * csq * h * rinv;

        } else if ((!fluidi) && fluidj) {
          compute_interface->correct_v(v[j], v[i], vj, j, i);
          rhoi = compute_interface->correct_rho(i,j);
          Pi = fix_pressure->calculate_P(rhoi);
          if ((chi[i] > 0.9) && (r < (h * 0.5)))
            fmag = (chi[i] - 0.9) * (h * 0.5 - r) * rho0 * csq * h * rinv;

        } else if ((!fluidi) && (!fluidj)) {
          rhoi = 1.0;
          rhoj = 1.0;
        }

        // Repel if close to inner solid particle
        scale3(fmag, dx, fsolid);

        // Compute volume after reconstructing
        voli = imass / rhoi;
        volj = jmass / rhoj;

        //Check if Second order kernels will be used for eta * Lap(v)
        error_flag = 0;
        if (laplacian_order == 2) {
          error_flag = compute_kernel->calc_d2w(i, j, dx[0], dx[1], dx[2], r);
          d2Wij = compute_kernel->d2Wij;
          d2Wji = compute_kernel->d2Wji;
        }

        // Thermal Evolution
        if (thermal_flag) {
          dT = 0.0;
          for (a = 0; a < dim; a++) dT += dx[a] * dWij[a];
            //TODO: Assumes heat capacity and density = 1, needs to be generalized
          dT *= voli * volj * (kappai + kappaj) * (Ti - Tj) * rinv * rinv;
          heatflow[i] += dT;
        }

        // If either particle is fluid, compute hydrostatic and viscous forces
        // Compute eta*Lap(v) -  different forms depending on order of RK correction
        if (pair_force_flag) {
          //Hydrostatic pressure forces
          fp_prefactor = voli * volj * (Pj + Pi);

          //Add artificial viscous pressure if required
          if (artificial_visc_flag && pair_avisc_flag){
            //Interpolate velocities to midpoint and use this difference for artificial viscosity
            for (a = 0; a < dim; a++) {
              du[a] = vi[a] - vj[a];
              for (b = 0; b < dim; b++) {
                du[a] -= 0.5 * (gradv[i][a * dim + b] + gradv[j][a * dim + b]) * dx[b];
              }
            }
            // TODO is hinv3 supposed to be 3.0/h?
            mu = lensq3(du) * 3 * hinv;
            mu /= (rsq * 9 * hinv * hinv + EPSILON);
            mu= MIN(0.0, mu);

            // "kinematic viscous pressure"  q = Q/rho
            q = av * (-2.0 * cs * mu + 1.0 * mu * mu); //TODO is 1.0 supposed to be something?
            fp_prefactor += voli * volj * q * (rhoj + rhoi);
          }

          // -Grad[P + Q]
          scale3(-fp_prefactor, dWij, dfp);

          // Now compute viscous eta * Lap[v] terms
          for (a = 0; a < dim; a ++) {
            fv[a] = 0.0;
            for (b = 0; b < dim; b++) fv[a] += (vi[a] - vj[a]) * dx[b] * dWij[b];
            fv[a] *= voli * volj * rinv * rinv * (etai + etaj);
          }

        } else {
          zero3(fv);
          zero3(dfp);
        }

        if (pair_force_flag) {
          f[i][0] += fv[0] + dfp[0] + fsolid[0];
          f[i][1] += fv[1] + dfp[1] + fsolid[1];
          f[i][2] += fv[2] + dfp[2] + fsolid[2];
          f_pressure[i][0] += dfp[0];
          f_pressure[i][1] += dfp[1];
          f_pressure[i][2] += dfp[2];
        }

        // Density damping
        // conventional for low-order h
        // interpolated for RK 1 & 2  (Antuono et al, Computers & Fluids 2021)
        if (rho_damp_flag && pair_rho_flag) {
          if ((laplacian_order >= 1) && (error_flag == 0)) {
            psi_ij = rhoj - rhoi;
            Fij = 0.0;
            for (a = 0; a < dim; a++){
              psi_ij += 0.5 * (gradr[i][a] + gradr[j][a]) * dx[a];
              Fij -= dx[a] * dWij[a];
            }
            Fij *= rinv * rinv;
            drho[i] += 2 * rho_damp * psi_ij * Fij * volj;
          } else {
            drho_damp = 2 * rho_damp * (rhoj - rhoi) * rinv * wp;
            drho[i] -= drho_damp * volj;
          }
        }

        if (evflag) // Doesn't account for unbalanced forces
          ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fv[0]+dfp[0], fv[1]+dfp[1], fv[2]+dfp[2], dx[0], dx[1], dx[2]);

        // Newton neighbors
        if (newton_pair || j < nlocal) {

          if (thermal_flag) {
            dT = 0.0;
            for (a = 0; a < dim; a++) dT += dx[a] * dWji[a]; //TODO csq and h terms
            dT *= -voli * volj * (kappai + kappaj) * (Ti - Tj) * rinv * rinv;
            heat[j] -= dT;
          }

          for (a = 0; a < dim; a ++) {
            fv[a] = 0.0;
            for (b = 0; b < dim; b++) fv[a] += (vi[a] - vj[a]) * dx[b] * dWji[b];
            fv[a] *= -voli * volj * (etai + etaj) * rinv * rinv;
            // flip sign here b/c -= at accummulator
          }

          if (pair_force_flag) {
            for (a = 0; a < dim; a++)
              dfp[a] = fp_prefactor * dWji[a];
          }

          if (rho_damp_flag && pair_rho_flag){
            if ((laplacian_order >= 1) && (error_flag == 0)) {
              Fij = 0.0;
              for (a = 0; a < dim; a++) Fij += dx[a] * dWji[a];
              Fij *= rinv * rinv;
              psi_ij *= -1;
              drho[j] += 2 * rho_damp * psi_ij * Fij * voli;
            } else {
              drho_damp = 2 * rho_damp * (rhoj - rhoi) * rinv * wp;
              drho[j] += drho_damp * voli;
            }
          }
          if (pair_force_flag) {
            f[j][0] -= fv[0] + dfp[0] + fsolid[0];
            f[j][1] -= fv[1] + dfp[1] + fsolid[1];
            f[j][2] -= fv[2] + dfp[2] + fsolid[2];
            f_pressure[j][0] -= dfp[0];
            f_pressure[j][1] -= dfp[1];
            f_pressure[j][2] -= dfp[2];
          }
        }
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
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
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairRHEO::settings(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal pair_style command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "rho/damp") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR,"Illegal pair_style command");

      rho_damp_flag = 1;
      rho_damp = utils::numeric(FLERR,arg[iarg + 1],false,lmp);
      iarg++;
    } else if (strcmp(arg[iarg], "artificial/visc") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR,"Illegal pair_style command");

      artificial_visc_flag = 1;
      av = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg++;
    } else error->all(FLERR,"Illegal pair_style command, {}", arg[iarg]);
    iarg++;
  }
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairRHEO::coeff(int narg, char **arg)
{
  if (narg != 2)
    error->all(FLERR,"Incorrect number of args for pair_style rheo coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR,arg[0],1, atom->ntypes, ilo, ihi,error);
  utils::bounds(FLERR,arg[1],1, atom->ntypes, jlo, jhi,error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = 0; j <= atom->ntypes; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair rheo coefficients");
}

/* ----------------------------------------------------------------------
 setup specific to this pair style
 ------------------------------------------------------------------------- */

void PairRHEO::setup()
{
  auto fixes = modify->get_fix_by_style("rheo");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use fix rheo/viscosity");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);

  compute_kernel = fix_rheo->compute_kernel;
  compute_grad = fix_rheo->compute_grad;
  compute_interface = fix_rheo->compute_interface;
  thermal_flag = fix_rheo->thermal_flag;
  h = fix_rheo->h;
  csq = fix_rheo->csq;
  rho0 = fix_rheo->rho0;

  hinv = 1.0 / h;
  cs = sqrt(csq);
  laplacian_order = -1;

  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair RHEO requires ghost atoms store velocity");

  if (laplacian_order == -1) {
    if (fix_rheo->kernel_type == FixRHEO::CRK2)
      laplacian_order = 2;
    else if (fix_rheo->kernel_type == FixRHEO::CRK1)
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
  if (setflag[i][j] == 0) {
      error->all(FLERR,"All pair rheo coeffs are not set");
  }

  cut[i][j] = h;
  cut[j][i] = cut[i][j];

  return cut[i][j];
}
