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

// Contributing author: Pietro Sillano (TU Delft), 2026

// 07-02-2026 branched from working pair_membrane_sillano_gemini3_radial.cpp (19/12/2025 last edit). it is the potential I used for the paper anaysis.

// in this code we will change:
// - splay to the new form
// - we will introduce curvature as distance dependent term

// - new name for the pair: membrane_sillano_v2
// - new parameter: c0
// - update tilt with new form
// - update Usplay, splay radial force and splay torques with new form

#include "pair_mesomem.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMesomem::PairMesomem(LAMMPS *lmp) : Pair(lmp), eps(nullptr), sigma(nullptr), cut(nullptr)

{
  writedata = 1;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairMesomem::~PairMesomem()
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
  }
}

/* ---------------------------------------------------------------------- */

void PairMesomem::allocate()
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
  memory->create(c0, np1, np1, "pair:c0");    // Allocate c0 array
}

/* ---------------------------------------------------------------------- */

void PairMesomem::settings(int narg, char **arg)
{
  if (narg != 1) utils::missing_cmd_args(FLERR, "pair_style mesomem", error);
  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  // Reset cutoffs that have been explicitly set
  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   Set coefficients for one or more type pairs
   Args: sigma, eps, ktilt, ksplay, cut, weight_rcut, zeta
------------------------------------------------------------------------- */

void PairMesomem::coeff(int narg, char **arg)
{
  if (narg != 10) error->all(FLERR, "Incorrect args for pair coefficients");
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
  double c0_one = utils::numeric(FLERR, arg[9], false, lmp);    // New arg

  if (weight_rcut_one > cut_one || weight_rcut_one > cut_global) {
    error->all(FLERR, "Orientation cutoff w_c > isotropic distance cutoff r_c");
  }

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
      c0[i][j] = c0_one;    // Store c0

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ---------------------------------------------------------------------- */

void PairMesomem::init_style()
{
  // Requirement: atoms must have orientation (mu) and torque
  if (!atom->q_flag || !atom->mu_flag || !atom->torque_flag)
    error->all(FLERR, "Pair mesomem requires atom attributes q, mu, torque");

  neighbor->request(this, instance_me);
}

/* ---------------------------------------------------------------------- */

double PairMesomem::init_one(int i, int j)
{
  // Strict Manual Mixing:
  // If the user did not set coefficients for this pair, we error out.
  // Automatic mixing for ktilt/ksplay/zeta is not physically defined here.

  if (setflag[i][j] == 0) {
    error->all(FLERR, "All pair coeffs must be set manually for pair_style mesomem");
  }

  eps[j][i] = eps[i][j];
  sigma[j][i] = sigma[i][j];
  ktilt[j][i] = ktilt[i][j];
  ksplay[j][i] = ksplay[i][j];
  weight_rcut[j][i] = weight_rcut[i][j];
  zeta[j][i] = zeta[i][j];
  cut[j][i] = cut[i][j];
  c0[j][i] = c0[i][j];    // Copy c0

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

void PairMesomem::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, rsq, r, inv_r;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  double Utot = 0.0;
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

  // Local vector registers
  double ni[3], nj[3], rhat[3];
  double fx, fy, fz, tx, ty, tz, tx_j, ty_j, tz_j;

  // Math vars
  double inv_mag_i, inv_mag_j;
  double factor_lj;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // OPTIMIZATION: Cache inverse magnitude of I once per neighbor list
    inv_mag_i = 1.0 / mu[i][3];
    ni[0] = mu[i][0] * inv_mag_i;
    ni[1] = mu[i][1] * inv_mag_i;
    ni[2] = mu[i][2] * inv_mag_i;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq[itype][jtype]) {
        evdwl = 0.0;
        r = sqrt(rsq);
        inv_r = 1.0 / r;

        rhat[0] = delx * inv_r;
        rhat[1] = dely * inv_r;
        rhat[2] = delz * inv_r;

        // --- 1. Isotropic (LJ/Cos) Force ---
        // Calculated for all pairs within r_cut
        double eps_val = eps[itype][jtype];
        double sigma_val = sigma
            [itype]
            [jtype];    // we are not using sigma*1.12 because we want to be transparent in the code!
        double rmin = sigma_val;
        double eps_lj = 0.0;
        double Ulj = 0.0;

        if (r < rmin) {
          double t = sigma_val * inv_r;
          double t2 = t * t;
          double t4 = t2 * t2;
          Ulj = eps_val * (t4 - 2.0 * t2);
          eps_lj = 4.0 * eps_val * inv_r * (t4 - t2);
        } else {
          double rcut = sqrt(cutsq[itype][jtype]);
          double zt = zeta[itype][jtype];

          // Precompute constant factors to avoid division in calc
          double denom = 1.0 / (rcut - rmin);
          double g = MathConst::MY_PI2 * (r - rmin) * denom;

          double cos_t = cos(g);
          double sin_t = sin(g);

          // Fast power calculation
          double cos_pow = pow(cos_t, 2.0 * zt - 1.0);
          // Alternatively, if zt is always integer, use loop for speed,
          // but pow is safer for general zeta.

          double cos_2zt = cos_pow * cos_t;

          Ulj = -eps_val * cos_2zt;

          // dU/dg * dg/dr
          double dU_dg = eps_val * (2.0 * zt) * cos_pow * sin_t;
          double dg_dr = MathConst::MY_PI2 * denom;
          eps_lj = -dU_dg * dg_dr;
        }

        // Initialize total forces/torques with just LJ part
        fx = eps_lj * rhat[0];
        fy = eps_lj * rhat[1];
        fz = eps_lj * rhat[2];
        tx = ty = tz = 0.0;
        tx_j = ty_j = tz_j = 0.0;

        // --- 2. Anisotropic (Tilt/Splay) Force ---
        double wr = weight_rcut[itype][jtype];

        if (r < wr) {
          // --- A. Weight Calculation ---
          double rga = 0.5 * wr;
          double r_wr = r / wr;

          // D = (r/wc)^4
          double r_wr_2 = r_wr * r_wr;
          double r_wr_4 = r_wr_2 * r_wr_2;
          double denom_w = r_wr_4 - 1.0;    // This is always negative for r < wr

          double w = 0.0;

          // REPLACEMENT LOGIC:
          // Only calculate if we are safely away from the singularity (denom_w < -1e-14).
          // If denom_w is closer to 0 than that, the exp() result is mathematically 0.0 anyway.
          double rga_sq = rga * rga;
          if (denom_w < -1e-14) {
            double val_exp = (r * r) / (rga_sq * denom_w);
            w = exp(val_exp);
          }
          // Else: w remains 0.0, avoiding division by tiny denom_w

          // --- B. Vector Normalization ---
          inv_mag_j = 1.0 / mu[j][3];
          nj[0] = mu[j][0] * inv_mag_j;
          nj[1] = mu[j][1] * inv_mag_j;
          nj[2] = mu[j][2] * inv_mag_j;

          // --- C. Dot Products ---
          double nirhat = ni[0] * rhat[0] + ni[1] * rhat[1] + ni[2] * rhat[2];
          double njrhat = nj[0] * rhat[0] + nj[1] * rhat[1] + nj[2] * rhat[2];
          double ninj = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];

          // --- D. Tilt Calculation (Extended with C0) ---
          double rh_x_ni[3], rh_x_nj[3];
          // Cross product rhat x ni
          rh_x_ni[0] = rhat[1] * ni[2] - rhat[2] * ni[1];
          rh_x_ni[1] = rhat[2] * ni[0] - rhat[0] * ni[2];
          rh_x_ni[2] = rhat[0] * ni[1] - rhat[1] * ni[0];
          // Cross product rhat x nj
          rh_x_nj[0] = rhat[1] * nj[2] - rhat[2] * nj[1];
          rh_x_nj[1] = rhat[2] * nj[0] - rhat[0] * nj[2];
          rh_x_nj[2] = rhat[0] * nj[1] - rhat[1] * nj[0];

          // sin_a2 = 0.5 * r * c0
          double sin_a2 = 0.5 * r * c0[itype][jtype];

          double kt = ktilt[itype][jtype];

          // Diff terms: ni.r + sin(alpha/2)
          double diff_i = nirhat + sin_a2;
          double diff_j = njrhat - sin_a2;    // check if sign needs to be the same

          double Utilt = 0.5 * kt * (diff_i * diff_i + diff_j * diff_j);

          // Vector u = n - (n.r)r
          double ui[3], uj[3];
          ui[0] = ni[0] - nirhat * rhat[0];
          ui[1] = ni[1] - nirhat * rhat[1];
          ui[2] = ni[2] - nirhat * rhat[2];

          uj[0] = nj[0] - njrhat * rhat[0];
          uj[1] = nj[1] - njrhat * rhat[1];
          uj[2] = nj[2] - njrhat * rhat[2];

          // Tilt Force term (angular part only)
          double ft_pref = -kt * inv_r;

          double tilt_fx = ft_pref * (diff_i * ui[0] + diff_j * uj[0]);
          double tilt_fy = ft_pref * (diff_i * ui[1] + diff_j * uj[1]);
          double tilt_fz = ft_pref * (diff_i * ui[2] + diff_j * uj[2]);

          // Accumulate Tilt Force scaled by weight
          fx += tilt_fx * w;
          fy += tilt_fy * w;
          fz += tilt_fz * w;

          // --- E. Splay Calculation ---
          double ks = ksplay[itype][jtype];

          double ni_x_nj[3];
          ni_x_nj[0] = ni[1] * nj[2] - ni[2] * nj[1];
          ni_x_nj[1] = ni[2] * nj[0] - ni[0] * nj[2];
          ni_x_nj[2] = ni[0] * nj[1] - ni[1] * nj[0];

          double Usplay = 0.5 * ks * (ninj - 1.0 + 2.0 * (sin_a2 * sin_a2)) *
              (ninj - 1.0 + 2.0 * (sin_a2 * sin_a2));

          // --- F. Radial Correction (Energy Conservation) ---
          // This force is purely radial (along rhat) and repulsive
          // Derived from - (Utilt + Usplay) * (dw/dr)

          double U_ang_sum = Utilt + Usplay;

          // Apply radial correction
          if (w > 0) {
            // Factor = w * [ 2 * (D+1) ] / [ rga^2 * (D-1)^2 ]

            // 8 Feb 2026 removed the minus sign in rad_numerator
            double rad_numerator = 2.0 * w * (r_wr_4 + 1.0) * r;
            double rad_denominator = rga_sq * denom_w * denom_w;

            double f_rad_mag = U_ang_sum * (rad_numerator / rad_denominator);

            fx += f_rad_mag * rhat[0];
            fy += f_rad_mag * rhat[1];
            fz += f_rad_mag * rhat[2];

            // --- Radial part coming from tilt term ---
            double f_rad_tilt = 0.5 * kt * c0[itype][jtype] * (diff_j - diff_i);
            fx += w * f_rad_tilt * rhat[0];
            fy += w * f_rad_tilt * rhat[1];
            fz += w * f_rad_tilt * rhat[2];

            // --- Radial part coming from splay term ---
            double f_rad_splay = -ks * (ninj - 1.0 + 2 * sin_a2 * sin_a2) *
                (c0[itype][jtype] * c0[itype][jtype] * r);
            fx += w * f_rad_splay * rhat[0];
            fy += w * f_rad_splay * rhat[1];
            fz += w * f_rad_splay * rhat[2];
          }

          // --- G. Torques ---
          // Torque on I
          double splay_pref = ks * (ninj - 1.0 + 2 * sin_a2 * sin_a2);

          tx += w * (kt * diff_i * rh_x_ni[0] - splay_pref * ni_x_nj[0]);
          ty += w * (kt * diff_i * rh_x_ni[1] - splay_pref * ni_x_nj[1]);
          tz += w * (kt * diff_i * rh_x_ni[2] - splay_pref * ni_x_nj[2]);

          // Torque on J
          tx_j += w * (kt * diff_j * rh_x_nj[0] + splay_pref * ni_x_nj[0]);
          ty_j += w * (kt * diff_j * rh_x_nj[1] + splay_pref * ni_x_nj[1]);
          tz_j += w * (kt * diff_j * rh_x_nj[2] + splay_pref * ni_x_nj[2]);

          // Accumulate Energy
          if (eflag) evdwl = Ulj + w * U_ang_sum;
        } else {
          if (eflag) evdwl = Ulj;
        }

        // --- FINAL APPLY ---
        if (eflag) evdwl *= factor_lj;

        // Apply Factor LJ to forces/torques
        fx *= factor_lj;
        fy *= factor_lj;
        fz *= factor_lj;
        tx *= factor_lj;
        ty *= factor_lj;
        tz *= factor_lj;
        tx_j *= factor_lj;
        ty_j *= factor_lj;
        tz_j *= factor_lj;

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        torque[i][0] += tx;
        torque[i][1] += ty;
        torque[i][2] += tz;

        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] += tx_j;
          torque[j][1] += ty_j;
          torque[j][2] += tz_j;
        }

        // Use ev_tally_xyz for non-central forces (virial correction)
        if (evflag)
          ev_tally_xyz(i, j, nlocal, newton_pair, evdwl, 0.0, fx, fy, fz, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMesomem::write_restart(FILE *fp)
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
        fwrite(&c0[i][j], sizeof(double), 1, fp);    // Write c0
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMesomem::read_restart(FILE *fp)
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
          utils::sfread(FLERR, &c0[i][j], sizeof(double), 1, fp, nullptr, error);    // Read c0
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

void PairMesomem::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMesomem::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {

    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
  }

  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairMesomem::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g %g %g %g %g %g\n", i, sigma[i][i], eps[i][i], ktilt[i][i],
            ksplay[i][i], cut[i][i], weight_rcut[i][i], zeta[i][i], c0[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairMesomem::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g %g %g %g \n", i, j, sigma[i][i], eps[i][i], ktilt[i][j],
              ksplay[i][j], cut[i][i], weight_rcut[i][j], zeta[i][i]);
}
