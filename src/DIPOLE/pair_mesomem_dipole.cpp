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
   Contributing author: Pietro Sillano (TU Delft)
------------------------------------------------------------------------- */

#include "pair_mesomem_dipole.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "info.h"
#include "math_const.h"
#include "math_extra.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_PI2;

static const char cite_pair_mesomem_dipole[] =
    "pair mesomem/dipole command: doi:10.1103/4dhv-8xd7\n\n"
    "@Article{Sillano26,\n"
    " author =  {P. Sillano and S. J. Marrink and T. Idema},\n"
    " title =   {{MesoMem}: A Mesoscale Membrane Model Based on an Additive Potential},\n"
    " journal = {Phys.\\ Rev.\\ E},\n"
    " year =    2026,\n"
    " volume =  114,\n"
    " number =  3,\n"
    " pages =   {034412}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

PairMesomemDipole::PairMesomemDipole(LAMMPS *lmp) :
    Pair(lmp), cut(nullptr), sigma(nullptr), eps(nullptr), ktilt(nullptr), ksplay(nullptr),
    weight_rcut(nullptr), zeta(nullptr), c0(nullptr), gscale(nullptr), wc_inv(nullptr),
    wc_half2inv(nullptr), zpow(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_mesomem_dipole);

  writedata = 1;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairMesomemDipole::~PairMesomemDipole()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(sigma);
    memory->destroy(eps);
    memory->destroy(ktilt);
    memory->destroy(ksplay);
    memory->destroy(weight_rcut);
    memory->destroy(zeta);
    memory->destroy(c0);
    memory->destroy(gscale);
    memory->destroy(wc_inv);
    memory->destroy(wc_half2inv);
    memory->destroy(zpow);
  }
}

/* ---------------------------------------------------------------------- */

void PairMesomemDipole::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, rsq, evdwl, factor_lj, inv_mag;
  double del[3], ni[3], fi[3], ti[3], tj[3];
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double **mu = atom->mu;
  double **torque = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

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

    // particles with a zero dipole moment only have the isotropic interaction

    const double *nip = nullptr;
    if (mu[i][3] > 0.0) {
      inv_mag = 1.0 / mu[i][3];
      ni[0] = mu[i][0] * inv_mag;
      ni[1] = mu[i][1] * inv_mag;
      ni[2] = mu[i][2] * inv_mag;
      nip = ni;
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      jtype = type[j];

      del[0] = xtmp - x[j][0];
      del[1] = ytmp - x[j][1];
      del[2] = ztmp - x[j][2];
      rsq = del[0] * del[0] + del[1] * del[1] + del[2] * del[2];

      if (rsq < cutsq[itype][jtype]) {
        evdwl = mesomem_analytic(itype, jtype, rsq, del, nip, mu[j], fi, ti, tj);

        if (factor_lj != 1.0) {
          evdwl *= factor_lj;
          MathExtra::scale3(factor_lj, fi);
          MathExtra::scale3(factor_lj, ti);
          MathExtra::scale3(factor_lj, tj);
        }

        f[i][0] += fi[0];
        f[i][1] += fi[1];
        f[i][2] += fi[2];
        torque[i][0] += ti[0];
        torque[i][1] += ti[1];
        torque[i][2] += ti[2];

        if (newton_pair || j < nlocal) {
          f[j][0] -= fi[0];
          f[j][1] -= fi[1];
          f[j][2] -= fi[2];
          torque[j][0] += tj[0];
          torque[j][1] += tj[1];
          torque[j][2] += tj[2];
        }

        // forces are not central, so the virial needs the force vector
        if (evflag)
          ev_tally_xyz(i, j, nlocal, newton_pair, evdwl, 0.0, fi[0], fi[1], fi[2], del[0], del[1],
                       del[2]);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   compute energy, force on i, and torques on i and j for one pair
   rsq, del  = squared distance and distance vector x_i - x_j
   ni        = unit vector along the dipole of i or nullptr if mu_i = 0
   muj       = dipole of j with its magnitude in muj[3]
   fi        = force on i, the force on j is -fi
   ti, tj    = torques on i and j
   returns the pair energy
------------------------------------------------------------------------- */

double PairMesomemDipole::mesomem_analytic(int itype, int jtype, double rsq, const double *del,
                                           const double *ni, const double *muj, double *fi,
                                           double *ti, double *tj) const
{
  const double r = sqrt(rsq);
  const double rinv = 1.0 / r;
  const double rhat[3] = {del[0] * rinv, del[1] * rinv, del[2] * rinv};

  // isotropic part:
  // r < sigma: eps * [(sigma/r)^4 - 2 (sigma/r)^2]
  // r >= sigma: -eps * cos^(2 zeta)[pi/2 (r - sigma) / (r_c - sigma)]

  const double epsilon = eps[itype][jtype];
  const double sig = sigma[itype][jtype];
  double energy, fpair;

  if (r < sig) {
    const double sr2 = sig * sig / rsq;
    const double sr4 = sr2 * sr2;
    energy = epsilon * (sr4 - 2.0 * sr2);
    fpair = 4.0 * epsilon * rinv * (sr4 - sr2);
  } else {
    const double zt = zeta[itype][jtype];
    const double dgdr = gscale[itype][jtype];
    const double g = dgdr * (r - sig);
    const double cosg = cos(g);
    const double sing = sin(g);
    const int n = zpow[itype][jtype];
    const double cospow = (n < 0) ? pow(cosg, 2.0 * zt - 1.0) : MathSpecial::powint(cosg, n);
    energy = -epsilon * cospow * cosg;
    fpair = -2.0 * zt * epsilon * cospow * sing * dgdr;
  }

  fi[0] = fpair * rhat[0];
  fi[1] = fpair * rhat[1];
  fi[2] = fpair * rhat[2];
  ti[0] = ti[1] = ti[2] = 0.0;
  tj[0] = tj[1] = tj[2] = 0.0;

  // orientation dependent part, requires a dipole on both particles and r < w_c

  if (!ni || (muj[3] <= 0.0) || (r >= weight_rcut[itype][jtype])) return energy;

  // weight function w = exp[-r^2 / ((w_c/2)^2 (1 - (r/w_c)^4))]
  // it underflows to zero close to w_c, which also avoids a division by zero

  const double rw = r * wc_inv[itype][jtype];
  const double rw4 = rw * rw * rw * rw;
  const double dw = rw4 - 1.0;
  if (dw >= 0.0) return energy;
  const double wfac = rsq * wc_half2inv[itype][jtype];
  const double w = exp(wfac / dw);
  if (w <= 0.0) return energy;

  const double inv_mag = 1.0 / muj[3];
  const double nj[3] = {muj[0] * inv_mag, muj[1] * inv_mag, muj[2] * inv_mag};

  const double nirhat = MathExtra::dot3(ni, rhat);
  const double njrhat = MathExtra::dot3(nj, rhat);
  const double ninj = MathExtra::dot3(ni, nj);

  // tilt: 1/2 k_tilt (d_i^2 + d_j^2) with d_i = n_i.rhat + s, d_j = n_j.rhat - s
  // splay: 1/2 k_splay (n_i.n_j - 1 + 2 s^2)^2
  // with the spontaneous curvature shift s = r C_0 / 2

  const double kt = ktilt[itype][jtype];
  const double ks = ksplay[itype][jtype];
  const double curv = c0[itype][jtype];
  const double s = 0.5 * r * curv;
  const double di = nirhat + s;
  const double dj = njrhat - s;
  const double splay = ninj - 1.0 + 2.0 * s * s;
  const double uang = 0.5 * kt * (di * di + dj * dj) + 0.5 * ks * splay * splay;

  // force from the angular dependence of the tilt term, with
  // d(n.rhat)/dr = (n - (n.rhat) rhat) / r

  const double ftilt = -kt * rinv;
  for (int k = 0; k < 3; ++k)
    fi[k] += w * ftilt * (di * (ni[k] - nirhat * rhat[k]) + dj * (nj[k] - njrhat * rhat[k]));

  // radial forces from the weight function and from the r dependence of s

  const double fradial = uang * 2.0 * wfac * (rw4 + 1.0) * rinv / (dw * dw) +
      0.5 * kt * curv * (dj - di) - ks * splay * curv * curv * r;
  fi[0] += w * fradial * rhat[0];
  fi[1] += w * fradial * rhat[1];
  fi[2] += w * fradial * rhat[2];

  // torques: t = n x (-dU/dn)

  double rh_x_ni[3], rh_x_nj[3], ni_x_nj[3];
  MathExtra::cross3(rhat, ni, rh_x_ni);
  MathExtra::cross3(rhat, nj, rh_x_nj);
  MathExtra::cross3(ni, nj, ni_x_nj);

  for (int k = 0; k < 3; ++k) {
    ti[k] = w * (kt * di * rh_x_ni[k] - ks * splay * ni_x_nj[k]);
    tj[k] = w * (kt * dj * rh_x_nj[k] + ks * splay * ni_x_nj[k]);
  }

  return energy + w * uang;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMesomemDipole::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(cut, np1, np1, "pair:cut");
  memory->create(sigma, np1, np1, "pair:sigma");
  memory->create(eps, np1, np1, "pair:eps");
  memory->create(ktilt, np1, np1, "pair:ktilt");
  memory->create(ksplay, np1, np1, "pair:ksplay");
  memory->create(weight_rcut, np1, np1, "pair:weight_rcut");
  memory->create(zeta, np1, np1, "pair:zeta");
  memory->create(c0, np1, np1, "pair:c0");
  memory->create(gscale, np1, np1, "pair:gscale");
  memory->create(wc_inv, np1, np1, "pair:wc_inv");
  memory->create(wc_half2inv, np1, np1, "pair:wc_half2inv");
  memory->create(zpow, np1, np1, "pair:zpow");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMesomemDipole::settings(int narg, char **arg)
{
  // the cutoff is set for each pair of atom types with pair_coeff

  if (narg > 0)
    error->all(FLERR, 1, "Illegal pair_style mesomem/dipole command: unexpected argument {}",
               arg[0]);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   args: sigma, eps, ktilt, ksplay, cut, weight_rcut, zeta, c0
------------------------------------------------------------------------- */

void PairMesomemDipole::coeff(int narg, char **arg)
{
  if (narg != 10) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double sigma_one = utils::numeric(FLERR, arg[2], false, lmp);
  double eps_one = utils::numeric(FLERR, arg[3], false, lmp);
  double ktilt_one = utils::numeric(FLERR, arg[4], false, lmp);
  double ksplay_one = utils::numeric(FLERR, arg[5], false, lmp);
  double cut_one = utils::numeric(FLERR, arg[6], false, lmp);
  double weight_rcut_one = utils::numeric(FLERR, arg[7], false, lmp);
  double zeta_one = utils::numeric(FLERR, arg[8], false, lmp);
  double c0_one = utils::numeric(FLERR, arg[9], false, lmp);

  if (sigma_one <= 0.0)
    error->all(FLERR, 2, "Pair style mesomem/dipole requires sigma > 0.0, but sigma = {}",
               sigma_one);
  if (cut_one <= sigma_one)
    error->all(FLERR, 6, "Pair style mesomem/dipole requires r_c > sigma, but r_c = {}", cut_one);
  if ((weight_rcut_one <= 0.0) || (weight_rcut_one > cut_one))
    error->all(FLERR, 7, "Pair style mesomem/dipole requires 0.0 < w_c <= r_c, but w_c = {}",
               weight_rcut_one);
  if (zeta_one < 0.5)
    error->all(FLERR, 8, "Pair style mesomem/dipole requires zeta >= 0.5, but zeta = {}", zeta_one);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      sigma[i][j] = sigma_one;
      eps[i][j] = eps_one;
      ktilt[i][j] = ktilt_one;
      ksplay[i][j] = ksplay_one;
      cut[i][j] = cut_one;
      weight_rcut[i][j] = weight_rcut_one;
      zeta[i][j] = zeta_one;
      c0[i][j] = c0_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMesomemDipole::init_style()
{
  if (!atom->mu_flag || !atom->torque_flag)
    error->all(FLERR, Error::NOLASTLINE,
               "Pair style mesomem/dipole requires atom attributes mu and torque");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMesomemDipole::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR, Error::NOLASTLINE,
               "Pair style mesomem/dipole does not support mixing. Coefficients for all pairs of "
               "atom types must be set explicitly. Status:\n" +
                   Info::get_pair_coeff_status(lmp));

  gscale[i][j] = MY_PI2 / (cut[i][j] - sigma[i][j]);
  wc_inv[i][j] = 1.0 / weight_rcut[i][j];
  wc_half2inv[i][j] = 4.0 / (weight_rcut[i][j] * weight_rcut[i][j]);

  // use the faster integer power function when possible

  const double zexp = 2.0 * zeta[i][j] - 1.0;
  if ((zexp == floor(zexp)) && (zexp < 100.0))
    zpow[i][j] = static_cast<int>(zexp);
  else
    zpow[i][j] = -1;

  eps[j][i] = eps[i][j];
  sigma[j][i] = sigma[i][j];
  ktilt[j][i] = ktilt[i][j];
  ksplay[j][i] = ksplay[i][j];
  weight_rcut[j][i] = weight_rcut[i][j];
  zeta[j][i] = zeta[i][j];
  cut[j][i] = cut[i][j];
  c0[j][i] = c0[i][j];
  gscale[j][i] = gscale[i][j];
  wc_inv[j][i] = wc_inv[i][j];
  wc_half2inv[j][i] = wc_half2inv[i][j];
  zpow[j][i] = zpow[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMesomemDipole::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&eps[i][j], sizeof(double), 1, fp);
        fwrite(&ktilt[i][j], sizeof(double), 1, fp);
        fwrite(&ksplay[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
        fwrite(&weight_rcut[i][j], sizeof(double), 1, fp);
        fwrite(&zeta[i][j], sizeof(double), 1, fp);
        fwrite(&c0[i][j], sizeof(double), 1, fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMesomemDipole::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (comm->me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (comm->me == 0) {
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &eps[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &ktilt[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &ksplay[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &weight_rcut[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &zeta[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &c0[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&eps[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&ktilt[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&ksplay[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&weight_rcut[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&zeta[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&c0[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMesomemDipole::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMesomemDipole::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairMesomemDipole::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g %g %g %g %g %g\n", i, sigma[i][i], eps[i][i], ktilt[i][i],
            ksplay[i][i], cut[i][i], weight_rcut[i][i], zeta[i][i], c0[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairMesomemDipole::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g %g %g %g %g\n", i, j, sigma[i][j], eps[i][j], ktilt[i][j],
              ksplay[i][j], cut[i][j], weight_rcut[i][j], zeta[i][j], c0[i][j]);
}
