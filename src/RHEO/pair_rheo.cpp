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
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"
#include "utils.h"

#include <cmath>

using namespace LAMMPS_NS;

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

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, rsq, r, ir;

  double w, wp, rhoi, rhoj, voli, volj, Pi, Pj;
  double *dWij, *dWji, *d2Wij, *d2Wji, *dW1ij, *dW1ji;
  double vijeij, etai, etaj, kappai, kappaj;
  double Ti, Tj, dT;
  double drho_damp, fmag;
  double mu, q, cs, fp_prefactor;
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
  double **v = atom->v;
  double **x = atom->x;
  double **f = atom->f;
  double **fp = atom->fp;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *drho = atom->drho;
  double *temp = atom->temp;
  double *heat = atom->heat;
  double *viscosity = atom->viscosity;
  double *conductivity = atom->conductivity;
  double *special_lj = force->special_lj;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *phase = atom->phase;

  int tmp1, tmp2;
  int index_visc = atom->find_custom("rheo_viscosity", tmp1, tmp2);
  if (index_visc == -1) error->all(FLERR, "Cannot find rheo viscosity");
  double *viscosity = atom->dvector[index_visc];

  double *conductivity;
  if (thermal_flag) {
    int index_cond = atom->find_custom("rheo_conductivity", tmp1, tmp2);
    if (index_cond == -1) error->all(FLERR, "Cannot find rheo conductivity");
    conductivity = atom->dvector[index_cond];
  }

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
    Ti = temp[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
      jtype = type[j];
      jmass = mass[jtype];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        ir = 1/r;

        kappaj = conductivity[j];
        etaj = viscosity[j];
        Tj = temp[j];

        pair_rho_flag = 0;
        pair_force_flag = 0;
        pair_avisc_flag = 0;
        if (phase[i] <= FixRHEO::FLUID_MAX || phase[j] <= FixRHEO::FLUID_MAX) {
          pair_force_flag = 1;
        }
        if (phase[i] <= FixRHEO::FLUID_MAX && phase[j] <= FixRHEO::FLUID_MAX) {
          pair_avisc_flag = 1;
          pair_rho_flag = 1;
        }

        wp = compute_kernel->calc_dw(i, j, dx[0], dx[1], dx[2], r);
        dWij = compute_kernel->dWij;
        dWji = compute_kernel->dWji;

        for (a = 0; a < dim; a ++) {
          vi[a] = v[i][a];
          vj[a] = v[j][a];
          fsolid[a] = 0.0;
        }

        // Add corrections for walls
        rhoi = rho[i];
        rhoj = rho[j];
        if (phase[i] <= FixRHEO::FLUID_MAX && phase[j] > FixRHEO::FLUID_MAX) {
          compute_sinterpolation->correct_v(v[i], v[j], vi, i, j);
          rhoj = compute_sinterpolation->correct_rho(j,i);

          // Repel if close to inner solid particle
          if (compute_sinterpolation->chi[j] > 0.9 && r < (cut[itype][jtype] * 0.5)) {
            fmag = (compute_sinterpolation->chi[j] - 0.9) * (cut[itype][jtype] * 0.5 - r);
            fmag *= rho0[itype] * csq[itype] * cut[itype][jtype] * ir;
            fsolid[0] = fmag * dx[0];
            fsolid[1] = fmag * dx[1];
            fsolid[2] = fmag * dx[2];
          }
        } else if (phase[i] > FixRHEO::FLUID_MAX && phase[j] <= FixRHEO::FLUID_MAX) {
          compute_sinterpolation->correct_v(v[j], v[i], vj, j, i);
          rhoi = compute_sinterpolation->correct_rho(i,j);

          // Repel if close to inner solid particle
          if (compute_sinterpolation->chi[i] > 0.9 && r < (cut[itype][jtype] * 0.5)) {
            fmag = (compute_sinterpolation->chi[i] - 0.9) * (cut[itype][jtype] * 0.5 - r);
            fmag *= rho0[jtype] * csq[jtype] * cut[itype][jtype] * ir;
            fsolid[0] = fmag * dx[0];
            fsolid[1] = fmag * dx[1];
            fsolid[2] = fmag * dx[2];
          }
        } else if (phase[i] > FixRHEO::FLUID_MAX && phase[j] > FixRHEO::FLUID_MAX) {
          rhoi = 1.0;
          rhoj = 1.0;
        }

        // Compute volume and pressure after reconstructing
        voli = imass / rhoi;
        volj = jmass / rhoj;
        Pj = calc_pressure(rhoj, jtype);
        Pi = calc_pressure(rhoi, itype);

        //Check if Second order kernels will be used for eta*Lap(v)
        error_flag = 0;
        if (laplacian_order == 2) {
          error_flag = compute_kernel->calc_d2w(i, j, dx[0], dx[1], dx[2], r);
          d2Wij = compute_kernel->d2Wij;
          d2Wji = compute_kernel->d2Wji;
        }

        //Thermal Evolution
        if (thermal_flag) {
          dT = 0.0;
          for (a = 0; a < dim; a++) {
            dT += (kappai + kappaj) * (Ti-Tj) * dx[a] * dWij[a] * ir * ir;
            //Assumes heat capacity and density = 1, needs to be generalized
          }
          dT *= voli * volj;
          heat[i] += dT;
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
            mu = (du[0] * dx[0] + du[1] * dx[1]+ du[2] * dx[2]) * hinv3;
            mu = mu / (rsq * hinv3 * hinv3 + EPSILON);
            mu= MIN(0.0, mu);
            cs = 0.5 * (sqrt(csq[itype]) + sqrt(csq[jtype]));
            // "kinematic viscous pressure"  q = Q/rho
            q = av*(-2.0*cs*mu + 1.0*mu*mu);
            fp_prefactor += voli*volj*q*(rhoj + rhoi);
          }

          // -Grad[P + Q]
          dfp[0] = - fp_prefactor*dWij[0];
          dfp[1] = - fp_prefactor*dWij[1];
          dfp[2] = - fp_prefactor*dWij[2];

          // Now compute viscous eta*Lap[v] terms
          for (a = 0; a < dim; a ++) {
            fv[a] = 0.0;
            for (b = 0; b < dim; b++) {
              fv[a] += (etai+etaj)*(vi[a]-vj[a])*dx[b]*dWij[b]*ir*ir;
            }
            fv[a] *= voli*volj;
          }

        } else {
          for (a = 0; a < dim; a ++) {
            fv[a] = 0;
            dfp[a] = 0;
          }
        }

        if (pair_force_flag) {
          f[i][0] += fv[0] + dfp[0] + fsolid[0];
          f[i][1] += fv[1] + dfp[1] + fsolid[1];
          f[i][2] += fv[2] + dfp[2] + fsolid[2];
          fp[i][0] += dfp[0];
          fp[i][1] += dfp[1];
          fp[i][2] += dfp[2];
        }

        // Density damping
        // conventional for low-order h
        // interpolated for RK 1 & 2  (Antuono et al, Computers & Fluids 2021)
        if (rho_damp_flag && pair_rho_flag) {
          if (laplacian_order>=1 && error_flag == 0){
            psi_ij = rhoj-rhoi;
            Fij = 0.0;
            for (a = 0; a < dim; a++){
              psi_ij += 0.5*(gradr[i][a]+gradr[j][a])*dx[a];
              Fij -= dx[a]*dWij[a];
            }
            Fij *= ir*ir;
            drho[i] += 2*rho_damp*psi_ij*Fij*volj;
          }
          else {
            drho_damp = 2*rho_damp*(rhoj-rhoi)*ir*wp;
            drho[i] -= drho_damp*volj;
          }
        }

        if (evflag) // Doesn't account for unbalanced forces
          ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fv[0]+dfp[0], fv[1]+dfp[1], fv[2]+dfp[2], dx[0], dx[1], dx[2]);

        // Newton neighbors
        if (newton_pair || j < nlocal) {

          if (thermal_flag) {
            dT = 0.0;
            for(a = 0; a < dim; a++){
              //dT += kappai*dWji[a]*gradt[i][a];
              //dT -= kappaj*dWji[a]*gradt[j][a];
              dT += 1/1*(kappai+kappaj)*(Ti-Tj)*dx[a]*dWji[a]*ir*ir; //Assumes heat capacity and density = 1, needs to be generalized
            }
            dT *= -voli*volj;
            heat[j] -= dT;
          }

          for (a = 0; a < dim; a ++) {
            fv[a] = 0.0;
            for (b = 0; b < dim; b++) {
              //fv[a] += etai*dWji[b]*(gradv[i][a*dim+b]+gradv[i][b*dim+a]);
              //fv[a] -= etaj*dWji[b]*(gradv[j][a*dim+b]+gradv[j][b*dim+a]);
              fv[a] += (etai+etaj)*(vi[a]-vj[a])*dx[b]*dWji[b]*ir*ir;
            }
            fv[a] *= -voli*volj; // flip sign here b/c -= at accummulator
          }



          if (pair_force_flag) {
            for (a = 0; a < dim; a++)
              dfp[a] = fp_prefactor*dWji[a];
          }

          if (rho_damp_flag && pair_rho_flag){
            if (laplacian_order>=1 && error_flag == 0){
              Fij = 0.0;
              for (a = 0; a < dim; a++){
                Fij += dx[a]*dWji[a];
              }
              Fij *= ir*ir;
              psi_ij *= -1;
              drho[j] += 2*rho_damp*psi_ij*Fij*voli;
            }
            else {
              drho_damp = 2*rho_damp*(rhoj-rhoi)*ir*wp;
              drho[j] += drho_damp*voli;
            }
          }
          if (pair_force_flag) {
            f[j][0] -= fv[0] + dfp[0] + fsolid[0];
            f[j][1] -= fv[1] + dfp[1] + fsolid[1];
            f[j][2] -= fv[2] + dfp[2] + fsolid[2];

            fp[j][0] -= dfp[0];
            fp[j][1] -= dfp[1];
            fp[j][2] -= dfp[2];
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
