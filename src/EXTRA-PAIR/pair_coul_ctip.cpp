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
   Contributing author: Gabriel Plummer (NASA)
------------------------------------------------------------------------- */

#include "pair_coul_ctip.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace EwaldConst;

static constexpr int DELTA = 4;

/* ---------------------------------------------------------------------- */

PairCoulCTIP::PairCoulCTIP(LAMMPS *lmp) :
    Pair(lmp), params(nullptr), shield(nullptr), shieldcu(nullptr), reffc(nullptr),
    reffcsq(nullptr), reffc4(nullptr), reffc7(nullptr), s2d_shift(nullptr), f_shift(nullptr),
    e_shift(nullptr), self_factor(nullptr), qeq_x(nullptr), qeq_j(nullptr), qeq_g(nullptr),
    qeq_z(nullptr), qeq_c(nullptr), qeq_q1(nullptr), qeq_q2(nullptr), qeq_w(nullptr)
{
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairCoulCTIP::~PairCoulCTIP()
{
  if (copymode) return;

  memory->sfree(params);
  memory->destroy(elem1param);
  memory->destroy(shield);
  memory->destroy(shieldcu);
  memory->destroy(reffc);
  memory->destroy(reffcsq);
  memory->destroy(reffc4);
  memory->destroy(reffc7);
  memory->destroy(s2d_shift);
  memory->destroy(f_shift);
  memory->destroy(e_shift);
  memory->destroy(self_factor);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(qeq_x);
    memory->destroy(qeq_j);
    memory->destroy(qeq_g);
    memory->destroy(qeq_z);
    memory->destroy(qeq_c);
    memory->destroy(qeq_q1);
    memory->destroy(qeq_q2);
    memory->destroy(qeq_w);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulCTIP::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  int iparam_i, jparam_j;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, ecoul, fpair;
  double r, rsq, reff, reffsq, reff4, forcecoul, factor_coul;
  double prefactor, erfcc, erfcd, t;
  double selfion;
  int *ilist, *jlist, *numneigh, **firstneigh;

  ecoul = 0.0;
  selfion = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  double erfcd_cut, t_cut, erfcc_cut;
  erfcd_cut = exp(-alpha * alpha * cut_coulsq);
  t_cut = 1.0 / (1.0 + EWALD_P * alpha * cut_coul);
  erfcc_cut = t_cut * (A1 + t_cut * (A2 + t_cut * (A3 + t_cut * (A4 + t_cut * A5)))) * erfcd_cut;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = map[type[i]];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    iparam_i = elem1param[itype];

    selfion = self(&params[iparam_i], qtmp);

    if (evflag) ev_tally(i, i, nlocal, 0, 0.0, selfion, 0.0, 0.0, 0.0, 0.0);

    if (eflag) {
      double e_self = self_factor[iparam_i][iparam_i] * qtmp * qtmp;
      ev_tally(i, i, nlocal, 0, 0.0, e_self, 0.0, 0.0, 0.0, 0.0);
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cut_coulsq) {
        jtype = map[type[j]];
        jparam_j = elem1param[jtype];
        r = sqrt(rsq);
        reff = cbrt(rsq * r + 1 / shieldcu[iparam_i][jparam_j]);
        reffsq = reff * reff;
        reff4 = reffsq * reffsq;
        prefactor = qqrd2e * qtmp * q[j] / r;
        erfcd = exp(-alpha * alpha * rsq);
        t = 1.0 / (1.0 + EWALD_P * alpha * r);
        erfcc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * erfcd;

        forcecoul = prefactor *
            (erfcc / r + 2.0 * alpha / MY_PIS * erfcd + rsq * r / reff4 - 1 / r -
             r * f_shift[iparam_i][jparam_j] + r * s2d_shift[iparam_i][jparam_j] * (r - cut_coul)) *
            r;
        if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;
        fpair = forcecoul / rsq;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) {
          ecoul = prefactor *
              (erfcc + r / reff - 1 - r * erfcc_cut / cut_coul - r / reffc[iparam_i][jparam_j] +
               r / cut_coul + r * f_shift[iparam_i][jparam_j] * (r - cut_coul) -
               r * s2d_shift[iparam_i][jparam_j] * (r - cut_coul) * (r - cut_coul) * 0.5);
          if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
        } else
          ecoul = 0.0;

        if (evflag) ev_tally(i, j, nlocal, newton_pair, 0.0, ecoul, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

double PairCoulCTIP::self(Param *param, double qi)
{
  double s1 = param->chi, s2 = param->eta, s3 = param->qmin, s4 = param->qmax, s5 = param->omega;

  if (qi < s3) {
    return qi * ((s1 - 2 * s3 * s5) + qi * (0.50 * s2 + s5)) + s3 * s3 * s5;
  } else if (qi < s4) {
    return qi * (s1 + qi * (0.50 * s2));
  } else {
    return qi * ((s1 - 2 * s4 * s5) + qi * (0.50 * s2 + s5)) + s4 * s4 * s5;
  }

  return 0.0;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairCoulCTIP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(qeq_x, n + 1, "pair:qeq_x");
  memory->create(qeq_j, n + 1, "pair:qeq_j");
  memory->create(qeq_g, n + 1, "pair:qeq_g");
  memory->create(qeq_z, n + 1, "pair:qeq_z");
  memory->create(qeq_c, n + 1, "pair:qeq_c");
  memory->create(qeq_q1, n + 1, "pair:qeq_q1");
  memory->create(qeq_q2, n + 1, "pair:qeq_q2");
  memory->create(qeq_w, n + 1, "pair:qeq_w");

  map = new int[n + 1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulCTIP::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Illegal pair_style command");

  alpha = utils::numeric(FLERR, arg[0], false, lmp);
  cut_coul = utils::numeric(FLERR, arg[1], false, lmp);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCoulCTIP::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  map_element2type(narg - 3, arg + 3);

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulCTIP::init_style()
{
  if (!atom->q_flag) error->all(FLERR, "Pair style coul/ctip requires atom attribute q");

  neighbor->add_request(this);

  cut_coulsq = cut_coul * cut_coul;

  memory->destroy(shield);
  memory->destroy(shieldcu);
  memory->destroy(reffc);
  memory->destroy(reffcsq);
  memory->destroy(reffc4);
  memory->destroy(reffc7);
  memory->destroy(s2d_shift);
  memory->destroy(f_shift);
  memory->destroy(e_shift);
  memory->destroy(self_factor);

  memory->create(shield, nelements, nelements, "pair:shield");
  memory->create(shieldcu, nelements, nelements, "pair:shieldcu");
  memory->create(reffc, nelements, nelements, "pair:reffc");
  memory->create(reffcsq, nelements, nelements, "pair:reffcsq");
  memory->create(reffc4, nelements, nelements, "pair:reffc4");
  memory->create(reffc7, nelements, nelements, "pair:reffc7");
  memory->create(s2d_shift, nelements, nelements, "pair:s2d_shift");
  memory->create(f_shift, nelements, nelements, "pair:f_shift");
  memory->create(e_shift, nelements, nelements, "pair:e_shift");
  memory->create(self_factor, nelements, nelements, "pair:self_factor");

  double qqrd2e = force->qqrd2e;
  double cut_coulcu, cut_coul4, alphacu, erfcd_cut, t_cut, erfcc_cut;
  cut_coulcu = cut_coulsq * cut_coul;
  cut_coul4 = cut_coulsq * cut_coulsq;
  alphacu = alpha * alpha * alpha;
  erfcd_cut = exp(-alpha * alpha * cut_coulsq);
  t_cut = 1.0 / (1.0 + EWALD_P * alpha * cut_coul);
  erfcc_cut = t_cut * (A1 + t_cut * (A2 + t_cut * (A3 + t_cut * (A4 + t_cut * A5)))) * erfcd_cut;

  for (int elt1 = 0; elt1 < nelements; elt1++) {
    for (int elt2 = 0; elt2 < nelements; elt2++) {
      shield[elt1][elt2] = sqrt(params[elt1].gamma * params[elt2].gamma);
      shieldcu[elt1][elt2] = shield[elt1][elt2] * shield[elt1][elt2] * shield[elt1][elt2];
      reffc[elt1][elt2] = std::cbrt(cut_coulcu + 1.0 / shieldcu[elt1][elt2]);
      reffcsq[elt1][elt2] = reffc[elt1][elt2] * reffc[elt1][elt2];
      reffc4[elt1][elt2] = reffcsq[elt1][elt2] * reffcsq[elt1][elt2];
      reffc7[elt1][elt2] = reffc4[elt1][elt2] * reffcsq[elt1][elt2] * reffc[elt1][elt2];
      s2d_shift[elt1][elt2] = 2.0 * erfcc_cut / cut_coulcu +
          4.0 * alpha / MY_PIS * erfcd_cut / cut_coulsq + 4.0 * alphacu / MY_PIS * erfcd_cut -
          2.0 / cut_coulcu + 4.0 * cut_coul4 / reffc7[elt1][elt2] -
          2.0 * cut_coul / reffc4[elt1][elt2];
      f_shift[elt1][elt2] = erfcc_cut / cut_coulsq + 2.0 * alpha / MY_PIS * erfcd_cut / cut_coul -
          1.0 / cut_coulsq + cut_coulsq / reffc4[elt1][elt2];
      e_shift[elt1][elt2] = 2.0 * erfcc_cut / cut_coul + 2.0 * alpha / MY_PIS * erfcd_cut -
          2.0 / cut_coul + 1.0 / reffc[elt1][elt2] + cut_coulcu / reffc4[elt1][elt2] +
          s2d_shift[elt1][elt2] * cut_coulsq * 0.5;
      self_factor[elt1][elt2] = -(e_shift[elt1][elt2] * 0.5 + alpha / MY_PIS) * qqrd2e;
    }
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulCTIP::init_one(int /*i*/, int /*j*/)
{
  return cut_coul;
}

/* ---------------------------------------------------------------------- */

void PairCoulCTIP::read_file(char *file)
{
  memory->sfree(params);
  params = nullptr;
  nparams = 0;
  maxparam = 0;

  // open file on proc 0
  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "coul/ctip");
    char *line;

    while ((line = reader.next_line(NPARAMS_PER_LINE))) {
      try {
        ValueTokenizer values(line);

        std::string iname = values.next_string();

        // ielement = 1st args
        int ielement;

        for (ielement = 0; ielement < nelements; ielement++)
          if (iname == elements[ielement]) break;

        // load up parameter settings and error check their values

        if (nparams == maxparam) {
          maxparam += DELTA;
          params = (Param *) memory->srealloc(params, maxparam * sizeof(Param), "pair:params");

          // make certain all addional allocated storage is initialized
          // to avoid false positives when checking with valgrind

          memset(params + nparams, 0, DELTA * sizeof(Param));
        }

        params[nparams].ielement = ielement;
        params[nparams].chi = values.next_double();
        params[nparams].eta = values.next_double();
        params[nparams].gamma = values.next_double();
        params[nparams].zeta = values.next_double();
        params[nparams].zcore = values.next_double();
        params[nparams].qmin = values.next_double();
        params[nparams].qmax = values.next_double();
        params[nparams].omega = values.next_double();

      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }

      // parameter sanity check

      if (params[nparams].eta < 0.0 || params[nparams].zeta < 0.0 || params[nparams].zcore < 0.0 ||
          params[nparams].gamma < 0.0 || params[nparams].qmin > params[nparams].qmax ||
          params[nparams].omega < 0.0)
        error->one(FLERR, "Illegal coul/ctip parameter");

      nparams++;
    }
  }

  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxparam, 1, MPI_INT, 0, world);

  if (comm->me != 0) {
    params = (Param *) memory->srealloc(params, maxparam * sizeof(Param), "pair:params");
  }

  MPI_Bcast(params, maxparam * sizeof(Param), MPI_BYTE, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairCoulCTIP::setup_params()
{
  int i, m, n;

  // set elem1param

  memory->destroy(elem1param);
  memory->create(elem1param, nelements, "pair:elem1param");

  for (i = 0; i < nelements; i++) {
    n = -1;
    for (m = 0; m < nparams; m++) {
      if (i == params[m].ielement) {
        if (n >= 0) error->all(FLERR, "Potential file has duplicate entry");
        n = m;
      }
    }
    if (n < 0) error->all(FLERR, "Potential file is missing an entry");
    elem1param[i] = n;
  }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulCTIP::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) { fwrite(&setflag[i][j], sizeof(int), 1, fp); }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulCTIP::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulCTIP::write_restart_settings(FILE *fp)
{
  fwrite(&alpha, sizeof(double), 1, fp);
  fwrite(&cut_coul, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulCTIP::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &alpha, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &cut_coul, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_coul, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

void *PairCoulCTIP::extract(const char *str, int &dim)
{
  if (strcmp(str, "cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_coul;
  }
  if (strcmp(str, "chi") == 0 && qeq_x) {
    dim = 1;
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0)
        qeq_x[i] = params[map[i]].chi;
      else
        qeq_x[i] = 0.0;
    return (void *) qeq_x;
  }
  if (strcmp(str, "eta") == 0 && qeq_j) {
    dim = 1;
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0)
        qeq_j[i] = params[map[i]].eta;
      else
        qeq_j[i] = 0.0;
    return (void *) qeq_j;
  }
  if (strcmp(str, "gamma") == 0 && qeq_g) {
    dim = 1;
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0)
        qeq_g[i] = params[map[i]].gamma;
      else
        qeq_g[i] = 0.0;
    return (void *) qeq_g;
  }
  if (strcmp(str, "zeta") == 0 && qeq_z) {
    dim = 1;
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0)
        qeq_z[i] = params[map[i]].zeta;
      else
        qeq_z[i] = 0.0;
    return (void *) qeq_z;
  }
  if (strcmp(str, "zcore") == 0 && qeq_c) {
    dim = 1;
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0)
        qeq_c[i] = params[map[i]].zcore;
      else
        qeq_c[i] = 0.0;
    return (void *) qeq_c;
  }
  if (strcmp(str, "qmin") == 0 && qeq_q1) {
    dim = 1;
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0)
        qeq_q1[i] = params[map[i]].qmin;
      else
        qeq_q1[i] = 0.0;
    return (void *) qeq_q1;
  }
  if (strcmp(str, "qmax") == 0 && qeq_q2) {
    dim = 1;
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0)
        qeq_q2[i] = params[map[i]].qmax;
      else
        qeq_q2[i] = 0.0;
    return (void *) qeq_q2;
  }
  if (strcmp(str, "omega") == 0 && qeq_w) {
    dim = 1;
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0)
        qeq_w[i] = params[map[i]].omega;
      else
        qeq_w[i] = 0.0;
    return (void *) qeq_w;
  }
  return nullptr;
}
