/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Hongyan Yuan (SUSTech)
------------------------------------------------------------------------- */

#include "pair_ylz.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_4PI;
using MathConst::MY_PI2;
using MathConst::MY_TWOBYSIXTH;

static const char cite_pair_ylz[] =
    "pair ylz command:\n\n"
    "@Article{Yuan10,\n"
    " author =  {H. Yuan, C. Huang, J. Li, G. Lykotrafitis, and S. Zhang},\n"
    " title =   {One-particle-thick, solvent-free, coarse-grained model for biological and "
    "biomimetic fluid membranes},\n"
    " journal = {Phys. Rev. E},\n"
    " year =    2010,\n"
    " volume =  82,\n"
    " pages =   {011905}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

PairYLZ::PairYLZ(LAMMPS *lmp) :
    Pair(lmp), epsilon(nullptr), sigma(nullptr), cut(nullptr), zeta(nullptr), mu(nullptr),
    beta(nullptr), avec(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_ylz);

  single_enable = 0;
  writedata = 1;
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairYLZ::~PairYLZ()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(cut);
    memory->destroy(zeta);
    memory->destroy(mu);
    memory->destroy(beta);
  }
}

/* ---------------------------------------------------------------------- */

void PairYLZ::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double evdwl, one_eng, rsq, factor_lj;
  double fforce[3], ttor[3], rtor[3], r12[3];
  double a1[3][3], a2[3][3];
  int *ilist, *jlist, *numneigh, **firstneigh;
  double *iquat, *jquat;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **x = atom->x;
  double **f = atom->f;
  double **tor = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  int flag = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    // count if pair has non-ellipsoid atom and skip over it to avoid segfault

    if (ellipsoid[i] < 0) {
      ++flag;
      continue;
    }

    iquat = bonus[ellipsoid[i]].quat;
    MathExtra::quat_to_mat_trans(iquat, a1);

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      // r12 = center to center vector

      r12[0] = x[j][0] - x[i][0];
      r12[1] = x[j][1] - x[i][1];
      r12[2] = x[j][2] - x[i][2];
      rsq = MathExtra::dot3(r12, r12);
      jtype = type[j];

      // count if pair has non-ellipsoid atom and skip over it to avoid segfault

      if (ellipsoid[j] < 0) {
        ++flag;
        continue;
      }

      // compute if less than cutoff

      if (rsq < cutsq[itype][jtype]) {

        jquat = bonus[ellipsoid[j]].quat;
        MathExtra::quat_to_mat_trans(jquat, a2);
        one_eng = ylz_analytic(i, j, a1, a2, r12, rsq, fforce, ttor, rtor);

        fforce[0] *= factor_lj;
        fforce[1] *= factor_lj;
        fforce[2] *= factor_lj;
        ttor[0] *= factor_lj;
        ttor[1] *= factor_lj;
        ttor[2] *= factor_lj;

        f[i][0] += fforce[0];
        f[i][1] += fforce[1];
        f[i][2] += fforce[2];
        tor[i][0] += ttor[0];
        tor[i][1] += ttor[1];
        tor[i][2] += ttor[2];

        if (newton_pair || j < nlocal) {
          rtor[0] *= factor_lj;
          rtor[1] *= factor_lj;
          rtor[2] *= factor_lj;
          f[j][0] -= fforce[0];
          f[j][1] -= fforce[1];
          f[j][2] -= fforce[2];
          tor[j][0] += rtor[0];
          tor[j][1] += rtor[1];
          tor[j][2] += rtor[2];
        }

        if (eflag) evdwl = factor_lj * one_eng;

        if (evflag)
          ev_tally_xyz(i, j, nlocal, newton_pair, evdwl, 0.0, fforce[0], fforce[1], fforce[2],
                       -r12[0], -r12[1], -r12[2]);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

  // error out if any pairs were skipped over because they were not both ellipsoids

  int flag_all;
  MPI_Allreduce(&flag, &flag_all, 1, MPI_INT, MPI_MAX, world);
  if (flag_all) error->all(FLERR, "All atoms for pair style ylz must be ellipsoids");
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairYLZ::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(epsilon, np1, np1, "pair:epsilon");
  memory->create(sigma, np1, np1, "pair:sigma");
  memory->create(cut, np1, np1, "pair:cut");
  memory->create(zeta, np1, np1, "pair:zeta");
  memory->create(mu, np1, np1, "pair:mu");
  memory->create(beta, np1, np1, "pair:beta");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairYLZ::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairYLZ::coeff(int narg, char **arg)
{
  if (narg < 8 || narg > 8) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);
  double zeta_one = utils::numeric(FLERR, arg[4], false, lmp);
  double mu_one = utils::numeric(FLERR, arg[5], false, lmp);
  double beta_one = utils::numeric(FLERR, arg[6], false, lmp);
  double cut_one = utils::numeric(FLERR, arg[7], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      zeta[i][j] = zeta_one;
      mu[i][j] = mu_one;
      beta[i][j] = beta_one;

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairYLZ::init_style()
{
  avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  if (!avec) error->all(FLERR, "Pair style ylz requires atom style ellipsoid");

  neighbor->request(this, instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairYLZ::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], sigma[i][i], sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i], sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  zeta[j][i] = zeta[i][j];
  mu[j][i] = mu[i][j];
  beta[j][i] = beta[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairYLZ::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
        fwrite(&zeta[i][j], sizeof(double), 1, fp);
        fwrite(&mu[i][j], sizeof(double), 1, fp);
        fwrite(&beta[i][j], sizeof(double), 1, fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairYLZ::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (comm->me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (comm->me == 0) {
          utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &zeta[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &mu[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &beta[i][j], sizeof(double), 1, fp, nullptr, error);
        }

        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&zeta[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&mu[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&beta[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairYLZ::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairYLZ::read_restart_settings(FILE *fp)
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

void PairYLZ::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g %g %g %g\n", i, epsilon[i][i], sigma[i][i], cut[i][i], zeta[i][i],
            mu[i][i], beta[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairYLZ::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g %g %g\n", i, j, epsilon[i][i], sigma[i][i], cut[i][j],
              zeta[i][j], mu[i][j], beta[i][j]);
}

/* ----------------------------------------------------------------------
   compute analytic energy, force (fforce), and torque (ttor & rtor)
   based on rotation matrices a
   if newton is off, rtor is not calculated for ghost atoms
------------------------------------------------------------------------- */

double PairYLZ::ylz_analytic(const int i, const int j, double a1[3][3], double a2[3][3],
                             double *r12, const double rsq, double *fforce, double *ttor,
                             double *rtor)

{
  int *type = atom->type;
  int newton_pair = force->newton_pair;
  int nlocal = atom->nlocal;

  double r12hat[3];
  MathExtra::normalize3(r12, r12hat);
  double r = sqrt(rsq);

  double ni1[3], nj1[3], dphi_drhat[3], dUdrhat[3], dUdni1[3], dUdnj1[3];
  double dphi_dni1[3], dphi_dnj1[3];
  double t, t1, t2, t4, cos_t, U, uR, uA, dUdr, dUdphi;
  const double energy_well = epsilon[type[i]][type[j]];
  const double rmin = MY_TWOBYSIXTH * sigma[type[i]][type[j]];
  const double rcut = cut[type[i]][type[j]];
  const double zt = zeta[type[i]][type[j]];
  const double muu = mu[type[i]][type[j]];
  const double sint = beta[type[i]][type[j]];

  ni1[0] = a1[0][0];
  ni1[1] = a1[0][1];
  ni1[2] = a1[0][2];

  nj1[0] = a2[0][0];
  nj1[1] = a2[0][1];
  nj1[2] = a2[0][2];

  const double ninj = MathExtra::dot3(ni1, nj1);
  const double ni1rhat = MathExtra::dot3(ni1, r12hat);
  const double nj1rhat = MathExtra::dot3(nj1, r12hat);

  const double a = ninj + (sint - ni1rhat) * (sint + nj1rhat) - 2.0 * sint * sint;
  const double phi = 1.0 + (a - 1.0) * muu;

  dphi_drhat[0] = muu * ((sint - ni1rhat) * nj1[0] - ni1[0] * (sint + nj1rhat));
  dphi_drhat[1] = muu * ((sint - ni1rhat) * nj1[1] - ni1[1] * (sint + nj1rhat));
  dphi_drhat[2] = muu * ((sint - ni1rhat) * nj1[2] - ni1[2] * (sint + nj1rhat));

  dphi_dni1[0] = muu * (nj1[0] - r12hat[0] * (sint + nj1rhat));
  dphi_dni1[1] = muu * (nj1[1] - r12hat[1] * (sint + nj1rhat));
  dphi_dni1[2] = muu * (nj1[2] - r12hat[2] * (sint + nj1rhat));

  dphi_dnj1[0] = muu * (ni1[0] + r12hat[0] * (sint - ni1rhat));
  dphi_dnj1[1] = muu * (ni1[1] + r12hat[1] * (sint - ni1rhat));
  dphi_dnj1[2] = muu * (ni1[2] + r12hat[2] * (sint - ni1rhat));

  if (r < rmin) {
    t = rmin / r;
    t2 = t * t;
    t4 = t2 * t2;
    uR = (t4 - 2.0 * t2) * energy_well;
    U = uR + (1.0 - phi) * energy_well;
    dUdr = 4.0 * (t2 - t4) / r * energy_well;
    dUdphi = -energy_well;
  } else {
    t = MY_PI2 * (r - rmin) / (rcut - rmin);
    cos_t = cos(t);
    t1 = cos_t;

    for (int k = 1; k <= 2 * zt - 2; k++) t1 *= cos_t;    // get cos()^(2zt-1)

    uA = -energy_well * t1 * cos_t;
    U = uA * phi;
    dUdr = MY_4PI / (rcut - rmin) * (t1) *sin(t) * phi * energy_well;
    dUdphi = uA;
  }

  dUdrhat[0] = dUdphi * dphi_drhat[0];
  dUdrhat[1] = dUdphi * dphi_drhat[1];
  dUdrhat[2] = dUdphi * dphi_drhat[2];

  double dUdrhatrhat = MathExtra::dot3(dUdrhat, r12hat);

  fforce[0] = dUdr * r12hat[0] + (dUdrhat[0] - dUdrhatrhat * r12hat[0]) / r;
  fforce[1] = dUdr * r12hat[1] + (dUdrhat[1] - dUdrhatrhat * r12hat[1]) / r;
  fforce[2] = dUdr * r12hat[2] + (dUdrhat[2] - dUdrhatrhat * r12hat[2]) / r;

  // torque i

  dUdni1[0] = dUdphi * dphi_dni1[0];
  dUdni1[1] = dUdphi * dphi_dni1[1];
  dUdni1[2] = dUdphi * dphi_dni1[2];

  MathExtra::cross3(dUdni1, ni1, ttor);    //minus sign is replace by swapping ni1 and dUdni1

  if (newton_pair || j < nlocal) {

    dUdnj1[0] = dUdphi * dphi_dnj1[0];
    dUdnj1[1] = dUdphi * dphi_dnj1[1];
    dUdnj1[2] = dUdphi * dphi_dnj1[2];

    MathExtra::cross3(dUdnj1, nj1, rtor);    //minus sign is replace by swapping ni1 and dUdni1
  }
  return U;
}
