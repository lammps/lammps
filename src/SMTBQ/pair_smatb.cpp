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
   Contributing author: Daniele Rapetti (iximiel@gmail.com)
------------------------------------------------------------------------- */

#include "pair_smatb.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSMATB::PairSMATB(LAMMPS *_lmp) :
    Pair(_lmp), nmax(0), on_eb(nullptr), r0(nullptr), p(nullptr), A(nullptr), q(nullptr),
    QSI(nullptr), cutOffStart(nullptr), cutOffEnd(nullptr), cutOffEnd2(nullptr), a3(nullptr),
    a4(nullptr), a5(nullptr), x3(nullptr), x4(nullptr), x5(nullptr)
{
  single_enable = 0;              // 1 if single() routine exists
  restartinfo = 1;                // 1 if pair style writes restart info
  respa_enable = 0;               // 1 if inner/middle/outer rRESPA routines
  one_coeff = 0;                  // 1 if allows only one coeff * * call
  manybody_flag = 1;              // 1 if a manybody potential
  no_virial_fdotr_compute = 0;    // 1 if does not invoke virial_fdotr_compute()
  writedata = 1;                  // 1 if writes coeffs to data file
  ghostneigh = 0;                 // 1 if pair style needs neighbors of ghosts

  // set comm size needed by this Pair
  comm_forward = 1;
  comm_reverse = 1;
}

/* ---------------------------------------------------------------------- */

PairSMATB::~PairSMATB()
{
  if (copymode) { return; }
  memory->destroy(on_eb);
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(r0);
    memory->destroy(p);
    memory->destroy(A);
    memory->destroy(q);
    memory->destroy(QSI);
    memory->destroy(cutOffStart);
    memory->destroy(cutOffEnd);
    memory->destroy(cutOffEnd2);
    memory->destroy(a3);
    memory->destroy(a4);
    memory->destroy(a5);
    memory->destroy(x5);
    memory->destroy(x4);
    memory->destroy(x3);
  }
}

/* ---------------------------------------------------------------------- */

void PairSMATB::compute(int eflag, int vflag)
{
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, del[3], fpair;
  double dijsq, dij;
  double espo, aexpp, qsiexpq, eb_i, Fb, Fr;
  double polyval, polyval2, polyval3, polyval4, polyval5;

  if (eflag || vflag) {
    ev_setup(eflag, vflag);
    eng_vdwl = 0;
  } else {
    evflag = vflag_fdotr = eflag_global = eflag_atom = 0;
  }

  // grow on_eb array if necessary

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->grow(on_eb, nmax, "pair_smatb:on_eb");
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  int newton_pair = force->newton_pair;

  // zero out on_eb
  memset(on_eb, 0, nall * sizeof(double));

  int inum = list->inum;
  int *ilist = list->ilist;
  int *jlist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // FIRST LOOP: CALCULATES the squared bonding energy and accumulate it in on_eb for each atom
  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; ++jj) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      del[0] = xtmp - x[j][0];
      del[1] = ytmp - x[j][1];
      del[2] = ztmp - x[j][2];
      dijsq = del[0] * del[0] + del[1] * del[1] + del[2] * del[2];

      if (dijsq < cutOffEnd2[itype][jtype]) {
        dij = sqrt(dijsq);
        if (dij < cutOffStart[itype][jtype]) {
          qsiexpq = (QSI[itype][jtype] * QSI[itype][jtype]) *
              exp(2.0 * q[itype][jtype] * (1.0 - dij / r0[itype][jtype]));
        } else {
          polyval = dij - cutOffEnd[itype][jtype];
          polyval3 = polyval * polyval * polyval;
          polyval4 = polyval3 * polyval;
          polyval5 = polyval4 * polyval;
          qsiexpq = x5[itype][jtype] * polyval5 + x4[itype][jtype] * polyval4 +
              x3[itype][jtype] * polyval3;
          qsiexpq = qsiexpq * qsiexpq;
        }
        on_eb[i] += qsiexpq;
        if (newton_pair) on_eb[j] += qsiexpq;
      }
    }
  }

  // communicate the squared bonding energy between the various bins

  if (newton_pair) comm->reverse_comm(this);

  // Support Loop: take the square root of the bonding energy and
  // accumulate it in the energy accumulator if needed the store the
  // reciprocal in on_eb in order to not do it in the SECOND LOOP

  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    if (i < nlocal) {
      eb_i = sqrt(on_eb[i]);
      if (eb_i != 0.0) {
        on_eb[i] = 1.0 / eb_i;
      } else {
        on_eb[i] = 0.0;
      }
      // if needed the bonding energy is accumulated:
      if (eflag_either) {
        if (eflag_atom) { eatom[i] -= eb_i; }
        if (eflag_global) { eng_vdwl -= eb_i; }
      }
    }
  }
  // this communication stores the denominators in the ghosts atoms,
  // this is needed because of how forces are calculated
  comm->forward_comm(this);

  // SECOND LOOP: given on_eb[i] calculates forces and energies
  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      del[0] = xtmp - x[j][0];
      del[1] = ytmp - x[j][1];
      del[2] = ztmp - x[j][2];

      dijsq = del[0] * del[0] + del[1] * del[1] + del[2] * del[2];
      if (dijsq < cutOffEnd2[itype][jtype]) {
        dij = sqrt(dijsq);
        if (dij < cutOffStart[itype][jtype]) {
          espo = 1.0 - dij / r0[itype][jtype];
          aexpp = exp(p[itype][jtype] * espo) * A[itype][jtype];
          Fr = (2.0 * aexpp) * (p[itype][jtype] / r0[itype][jtype]);
          qsiexpq = (QSI[itype][jtype] * QSI[itype][jtype]) * exp(2.0 * q[itype][jtype] * espo);
          Fb = -qsiexpq * q[itype][jtype] / r0[itype][jtype];
        } else {
          polyval = dij - cutOffEnd[itype][jtype];
          polyval2 = polyval * polyval;
          polyval3 = polyval2 * polyval;
          polyval4 = polyval3 * polyval;
          polyval5 = polyval4 * polyval;
          aexpp = a5[itype][jtype] * polyval5 + a4[itype][jtype] * polyval4 +
              a3[itype][jtype] * polyval3;
          Fr = -2.0 *
              (5.0 * a5[itype][jtype] * polyval4 + 4.0 * a4[itype][jtype] * polyval3 +
               3.0 * a3[itype][jtype] * polyval2);
          qsiexpq = x5[itype][jtype] * polyval5 + x4[itype][jtype] * polyval4 +
              x3[itype][jtype] * polyval3;
          Fb = ((5.0 * x5[itype][jtype] * polyval4 + 4.0 * x4[itype][jtype] * polyval3 +
                 3.0 * x3[itype][jtype] * polyval2)) *
              qsiexpq;
        }

        // calculates the module of the pair energy between i and j
        fpair = (Fb * (on_eb[i] + on_eb[j]) + Fr) / dij;

        f[i][0] += del[0] * fpair;
        f[i][1] += del[1] * fpair;
        f[i][2] += del[2] * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= del[0] * fpair;
          f[j][1] -= del[1] * fpair;
          f[j][2] -= del[2] * fpair;
        }
        if (evflag) {
          ev_tally(i, j, nlocal, newton_pair, 2.0 * aexpp, 0.0, fpair, del[0], del[1], del[2]);
        }
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairSMATB::settings(int narg, char **)
{
  if (narg > 0) error->all(FLERR, "Illegal pair_style command: smatb accepts no options");
}

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairSMATB::allocate()
{
  const int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair_smatb:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair_smatb:cutsq");
  memory->create(r0, np1, np1, "pair_smatb:r0");
  memory->create(p, np1, np1, "pair_smatb:p");
  memory->create(A, np1, np1, "pair_smatb:A");
  memory->create(q, np1, np1, "pair_smatb:q");
  memory->create(QSI, np1, np1, "pair_smatb:QSI");
  memory->create(cutOffStart, np1, np1, "pair_smatb:cutOffStart");
  memory->create(cutOffEnd, np1, np1, "pair_smatb:cutOffEnd");
  memory->create(cutOffEnd2, np1, np1, "pair_smatb:cutOffEnd2");
  memory->create(a3, np1, np1, "pair_smatb:a1");
  memory->create(a4, np1, np1, "pair_smatb:a2");
  memory->create(a5, np1, np1, "pair_smatb:a5");
  memory->create(x3, np1, np1, "pair_smatb:x1");
  memory->create(x4, np1, np1, "pair_smatb:x2");
  memory->create(x5, np1, np1, "pair_smatb:x3");

  allocated = 1;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairSMATB::coeff(int narg, char **arg)
{
  if (!allocated) { allocate(); }
  if (narg != 9) utils::missing_cmd_args(FLERR, "pair_style smatb", error);

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double myr0 = utils::numeric(FLERR, arg[2], false, lmp);
  double myp = utils::numeric(FLERR, arg[3], false, lmp);
  double myq = utils::numeric(FLERR, arg[4], false, lmp);
  double myA = utils::numeric(FLERR, arg[5], false, lmp);
  double myQSI = utils::numeric(FLERR, arg[6], false, lmp);
  double mycutOffStart = utils::numeric(FLERR, arg[7], false, lmp);
  double mycutOffEnd = utils::numeric(FLERR, arg[8], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      r0[i][j] = myr0;
      p[i][j] = myp;
      A[i][j] = myA;
      q[i][j] = myq;
      QSI[i][j] = myQSI;
      cutOffStart[i][j] = mycutOffStart;
      cutOffEnd[i][j] = mycutOffEnd;

      setflag[i][j] = 1;

      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ------------------------------------------------------------------------ */

void PairSMATB::init_style()
{
  if (force->newton_pair == 0) error->all(FLERR, "Pair style smatb requires newton pair on");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */

double PairSMATB::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    ///@todo implement smatb mixing rules
    cutOffStart[i][j] = MIN(cutOffStart[i][i], cutOffStart[j][j]);
    cutOffEnd[i][j] = MAX(cutOffEnd[i][i], cutOffEnd[j][j]);

    error->all(FLERR, "All pair coeffs are not set");
  }

  double es = cutOffEnd[i][j] - cutOffStart[i][j];
  double es2 = es * es;
  double es3 = es2 * es;

  // variables for poly for p and A
  double expp = A[i][j] * exp(p[i][j] * (1. - cutOffStart[i][j] / r0[i][j]));
  double ap = -1. / es3;
  double bp = p[i][j] / (r0[i][j] * es2);
  double cp = -(p[i][j] * p[i][j]) / (es * r0[i][j] * r0[i][j]);

  a5[i][j] = expp * (12. * ap + 6. * bp + cp) / (2. * es2);
  a4[i][j] = expp * (15. * ap + 7. * bp + cp) / es;
  a3[i][j] = expp * (20. * ap + 8. * bp + cp) / 2.;

  // variables for poly for q and qsi
  double expq = QSI[i][j] * exp(q[i][j] * (1. - cutOffStart[i][j] / r0[i][j]));
  double aq = -1 / es3;
  double bq = q[i][j] / (es2 * r0[i][j]);
  double cq = -(q[i][j] * q[i][j]) / (es * r0[i][j] * r0[i][j]);

  x5[i][j] = expq * (12. * aq + 6. * bq + cq) / (2. * es2);
  x4[i][j] = expq * (15. * aq + 7. * bq + cq) / es;
  x3[i][j] = expq * (20. * aq + 8. * bq + cq) / 2.;

  cutOffEnd2[i][j] = cutOffEnd[i][j] * cutOffEnd[i][j];
  if (i != j) {
    setflag[j][i] = 1;
    cutOffEnd2[j][i] = cutOffEnd2[i][j];

    r0[j][i] = r0[i][j];
    p[j][i] = p[i][j];
    q[j][i] = q[i][j];
    A[j][i] = A[i][j];
    QSI[j][i] = QSI[i][j];
    cutOffStart[j][i] = cutOffStart[i][j];
    cutOffEnd[j][i] = cutOffEnd[i][j];

    a3[j][i] = a3[i][j];
    a4[j][i] = a4[i][j];
    a5[j][i] = a5[i][j];
    x3[j][i] = x3[i][j];
    x4[j][i] = x4[i][j];
    x5[j][i] = x5[i][j];
  }

  return cutOffEnd[i][j];
}

/* ---------------------------------------------------------------------- */

int PairSMATB::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; ++i) {
    j = list[i];
    buf[m++] = on_eb[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSMATB::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; ++i) { on_eb[i] = buf[m++]; }
}

/* ---------------------------------------------------------------------- */

int PairSMATB::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; ++i) { buf[m++] = on_eb[i]; }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSMATB::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    on_eb[j] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void PairSMATB::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  fwrite(&tail_flag, sizeof(int), 1, fp);
}

/* ---------------------------------------------------------------------- */

void PairSMATB::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairSMATB::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&r0[i][j], sizeof(double), 1, fp);
        fwrite(&p[i][j], sizeof(double), 1, fp);
        fwrite(&q[i][j], sizeof(double), 1, fp);
        fwrite(&A[i][j], sizeof(double), 1, fp);
        fwrite(&QSI[i][j], sizeof(double), 1, fp);
        fwrite(&cutOffStart[i][j], sizeof(double), 1, fp);
        fwrite(&cutOffEnd[i][j], sizeof(double), 1, fp);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairSMATB::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &r0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &p[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &q[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &A[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &QSI[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cutOffStart[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cutOffEnd[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&r0[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&p[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&q[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&A[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&QSI[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cutOffStart[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cutOffEnd[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ---------------------------------------------------------------------- */

void PairSMATB::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) {
    fprintf(fp, "%d %g %g %g %g %g %g %g\n", i, r0[i][i], p[i][i], q[i][i], A[i][i], QSI[i][i],
            cutOffStart[i][i], cutOffEnd[i][i]);
  }
}

/* ---------------------------------------------------------------------- */

void PairSMATB::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      fprintf(fp, "%d %d %g %g %g %g %g %g %g\n", i, j, r0[i][j], p[i][j], q[i][j], A[i][j],
              QSI[i][j], cutOffStart[i][j], cutOffEnd[i][j]);
    }
  }
}
